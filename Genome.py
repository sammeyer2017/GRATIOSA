#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
import operator
from useful_functions_genome import *
from itertools import groupby
from globvar import *
from Gene import Gene
#from TSS import TSS
#from TTS import TTS
#from TU import TU
from TSS_TTS_TU import TSS, TTS, TU
from datetime import datetime
#from btssfinder import *
from scipy import stats
import Transcriptome
import HiC
import Chipseq
from pathlib import Path

#=============================================================================#

class Genome:

    def __init__(self, name):
        """ 
        Called when a HiC instance is created,
        initializes the attribute name.
        Example:
            >>> g = Genome.Genome("dickeya")
        """
        self.name = name


    def load_seq(self, filename="sequence.fasta"):
        """
        Load_seq loads DNA sequence from a .fasta file present in the 
        main directory of the organism using useful_functions_genome.load_seq 
        function. Adds this sequence, its complement, and its length to a
        Genome instance.

        Args:
            self (Genome instance)
            filename (Optional [str.}): name of the file containing the DNA
                                        sequence in FASTA format. 
                                        Default: "sequence.fasta"

        Outputs: creates 3 new attributes to the Genome instance
            seq (str.): genomic sequence compose of A,T,G and C
            seqcompl (str.): complement sequence to seq
            length (int.): length of the genomic sequence

        Example:
            >>> g = Genome.Genome("dickeya")
            >>> g.load_seq()
            >>> g.seq[100:151]
            'AATGTCGATCTTCAACATATCGCCGATCCGACGGGCACCCAGATCCTGCAG'
            >>> g.seqcompl[100:151]
            'TTACAGCTAGAAGTTGTATAGCGGCTAGGCTGCCCGTGGGTCTAGGACGTC'
        """
        seq_dir = f"{basedir}data/{self.name}/sequence.fasta"
        self.seq = load_seq(seq_dir).upper()
        self.seqcompl = ''
        self.length = len(self.seq)
        l = self.seq
        l = l.replace('A','t').replace('T','a')
        l = l.replace('C','g').replace('G','c')
        l = l.replace('a','A').replace('t','T')
        l = l.replace('c','C').replace('g','G')
        self.seqcompl = l


    def load_annotation(self, annot_file="sequence.gff3",*args, **kwargs):
        """ 
        load_annotation loads a gene annotation (coordinates, length, 
        name...) from a file present in the /annotation/ directory. 
        If the file is a .gff3 or .gff information will be loaded using 
        useful_functions_genome.load_gff. 
        Else, the information importation requires an annotation.info file, 
        containing column indices of each information in the data file and some 
        additional information, in the following order: 
        [0] Filename [1] Separator [2] Locus_tag column [3] Name column  
        [4] ID column [5] Strand column [6] Left coordinate column
        [7] Right coordinate column [8] File start line

        Args:
            self (Genome instance)
            filename (Optional [str.]): name of the file containing the genomic
                                        annotation. Default: "sequence.gff3" 
        Output: 
            self.genes (dict.): new attribute of the Genome instance. 
                                self.genes is a dictionary of shape 
                                {locus_tags: Gene object}. Each Gene object is 
                                initialized with the following attributes: 
                                locus_tag, ID, name,strand, left, right, start,
                                end, middle and length
                                See __init__ in the Gene class for more details
                                about each attribute.
        Example: 
            >>> g = Genome.Genome("dickeya")
            >>> g.load_annotation()
            >>> g.genes["Dda3937_00005"].name
            'guaA'
        """
        path2dir = f"{basedir}data/{self.name}/annotation/"
        path2file = f"{path2dir}{annot_file}"
        file_type = Path(path2file).suffix

        print(f"Trying to load annotation from: {path2file}")
        if file_type in [".gff", ".gff3"]:
            self.genes = load_gff(path2file)
            load = True
        else:
            load = False
            with open(f"{path2dir}annotation.info", "r") as f:
                skiphead = next(f)  
                for line in f:
                    line = line.strip().split('\t')
                    if line[0] == annot_file:
                        self.genes = load_annot_general(
                                                    f"{path2dir}{line[0]}",
                                                    line[1], int(line[2]), 
                                                    int(line[3]), int(line[4]), 
                                                    int(line[5]), int(line[6]),
                                                    int(line[7]), int(line[8]))
                        load = True
                f.close()
            if not load:
                print(f"Unable to load annotation from: {path2file}."
                      "Please add information in annotation.info or use gff3 or gff format.")
        if load == True :
            print("Annotation loaded")


    def load_genes_per_pos(self, window=0):
        """
        load_genes_per_pos associates each position with a list of genes 
        overlapping this and its surrounding positions delimited 
        by the window size given as an argument. 

        Args: 
            self (Genome instance)
            window (Optional [int.]): window size in b. load_genes_per_pos finds 
                                      genes between pos-window/2 and pos+window/2 
                                      (inclusive).
                                      Default: window = 0 

        Output:
            self.genes_per_pos (dict.): new attribute of the Genome instance. 
                                        self.genes is a dictionary of shape 
                                        {position: list of genes}. It contains,
                                        for each position p, the list of 
                                        genes overlapping any position between
                                        p-window/2 and p+window/2 (inclusive)
        
        N.B.: This method needs a genomic annotation. If no annotation is 
        loaded, the load_annotation method with the default "sequence.gff3" 
        file is computed. To use another annotation, please load an 
        annotation before using the load_genes_per_pos method.
        
        Example:
            >>> g = Genome.Genome("dickeya")
            >>> g.load_annotation(annot_file="sequence.gff3")
            >>> g.load_genes_per_pos(window=1000)
            >>> g.genes_per_pos[120]
            ['Dda3937_00155', 'Dda3937_00156', 'Dda3937_00154']
        """
        if not hasattr(self, "genes"):
            self.load_annotation()
        if not hasattr(self, "seq"):
            self.load_seq()

        self.genes_per_pos = {}

        for locus in self.genes.keys():
            gene = self.genes[locus]
            for pos in np.arange(gene.left - int(window / 2),
                                 gene.right + 1 + int(window / 2)):
                if pos < 0 : pos += self.length 
                if pos > self.length : pos -= self.length
                if pos in self.genes_per_pos.keys():
                    self.genes_per_pos[pos].append(locus)
                else:
                    self.genes_per_pos[pos] = [locus]

        for pos in np.arange(self.length):
            if pos not in self.genes_per_pos.keys():
                self.genes_per_pos[pos] = [None]

        
    def load_neighbor_all(self):
        """
        For each gene and positions, load_neighbor_all finds nearest 
        neighbors (left and right) on genome, whatever their strand.

        Arg: 
            self (Genome instance)

        Output: 
            2 new attributes of Gene instances related to the Genome 
            instance given as argument:
                self.genes[locus].left_neighbor: locus of the nearest 
                                                 left-side neighbor gene
                self.genes[locus].right_neighbor: locus of the nearest  
                                                  right-side neighbor gene
            4 new attributes of Genome instance : 
                self.genomic_situation (dict.): dict of shape {position:
                                                situation} with situation
                                                either "intergenic" or 
                                                "intragenic"
                self.left_neighbor (dict.): dict of shape {position:
                                            locus of the nearest 
                                            left-side neighbor gene}
                self.right_neighbor (dict.): dict of shape {position:
                                             locus of the nearest 
                                             right-side neighbor gene}
                self.gene (dict.): dict of shape {position: gene} with
                                   gene = "NA" if the position is 
                                   intergenic.

        N.B.: This method needs a genomic annotation. If no annotation is 
        loaded, the load_annotation method with the default "sequence.gff3" 
        file is computed. To use another annotation, please load an 
        annotation before using the load_neighbor_all method.

        Example:
            >>> g = Genome.Genome("dickeya")
            >>> g.load_annotation(annot_file="sequence.gff3")
            >>> g.load_neighbor_all()
            >>> g.genes["Dda3937_00005"].left_neighbor
            'Dda3937_00006'
            >>> g.genes["Dda3937_00005"].right_neighbor
            'Dda3937_00004'
            >>> g.g.genomic_situation[569] 
            'intragenic'
            >>> g.gene[569]
            'Dda3937_00156'
            >>> g.left_neighbor[569]
            'Dda3937_00155'
            >>> g.right_neighbor[569]
            'Dda3937_00157'
        """
        if not hasattr(self, "genes"):
            self.load_annotation()
        if not hasattr(self, "length"):
            self.load_seq()

        # Sorts genes according to their position.
        list_start = []
        for i in self.genes:
            list_start.append(int(self.genes[i].start))
        genes_order = [x for _, x in sorted(zip(list_start, self.genes.keys()))]

        self.genomic_situation = ["intergenic"]*self.length
        self.left_neighbor = {}
        self.right_neighbor = {}
        self.gene = ["NA"]*self.length

        for i in np.arange(len(genes_order)):
            # Adds right and left neighbors to each Gene object
            if i != 0:
                self.genes[genes_order[i]].add_left_neighbor(genes_order[i-1])
            else : 
                self.genes[genes_order[i]].add_left_neighbor(genes_order[-1])
            if i != len(genes_order)-1:
                self.genes[genes_order[i]].add_right_neighbor(genes_order[i+1])
            else : 
                self.genes[genes_order[i]].add_right_neighbor(genes_order[0])
            
            # Adds right and left neighbors to each genomic position
            for p in np.arange(self.genes[genes_order[i]].left, self.genes[genes_order[i]].right) :
                self.genomic_situation[p] = "intragenic"
                self.gene[p] = genes_order[i]
                if i != 0 :
                    self.left_neighbor[p] = genes_order[i-1]
                else : 
                    self.left_neighbor[p] = genes_order[-1]
                if i != len(genes_order)-1 :
                    self.right_neighbor[p] = genes_order[i+1]
                else : 
                    self.right_neighbor[p] = genes_order[0]
            if i != len(genes_order) - 1:
                for p in np.arange(self.genes[genes_order[i]].right, self.genes[genes_order[i+1]].left) :
                    self.left_neighbor[p] = genes_order[i]
                    self.right_neighbor[p] = genes_order[i+1]
            else : 
                for p in np.arange(self.genes[genes_order[i]].right, self.length) :
                    self.left_neighbor[p] = genes_order[i]
                    self.right_neighbor[p] = genes_order[0]
                for p in np.arange(0,self.genes[genes_order[0]].left) :
                    self.left_neighbor[p] = genes_order[i]
                    self.right_neighbor[p] = genes_order[0]


    def load_gene_orientation(self, couple=3, max_dist=5000):
        """
        Compute gene orientation with the following criteria : 
        If couple = 3, gene is considered :
            divergent if left neighbor on - strand 
                     and right neighbor on + strand,
            convergent if left neighbor on + strand 
                      and right neighbor on - strand, 
            tandem if left and right neighbors on same strand 
                   (whatever the strand of the given gene is)
            isolated if the distance between neighbors is higher than the  
                     maximal distance given as argument
        If couple = 2, gene is considered 
            tandem if predecessor (left neighbor for gene on + strand, 
                   right neighbor for gene on - strand) is on same strand, 
            divergent if the predecessor is on opposite strand.
        
        Arg: 
            self (Genome instance)
            couple (Optional [int.]): number of genes to consider in a "couple"
                           if couple = 2 : computes the orientation of a 
                                           gene relative to its predecessor
                           if couple = 3 : computes the orientation of a 
                                           gene relative to its two
                                           neighbors
                           Default: 3
            max_dist (Optional [int.]): maximal distance between 2 genes start 
                                  positions for seeking neighbor (Default: 5kb)

        Outputs:
            self.orientation (dict.): new attribute of the Genome instance.
                                      self.orientation is a dictionary of 
                                      shape {orientation: list of genes}. 
                                      It contains the list of genes for each 
                                      orientation (tandem, divergent, 
                                      convergent, and isolated if couple=3,
                                      tandem and divergent if couple=2)                          
            self.genes[locus].orientation (str.): new attribute of Gene 
                                                instances related to the Genome 
                                                instance given as argument

        N.B.: This method needs a genomic annotation. If no annotation is 
        loaded, the load_annotation method with the default "sequence.gff3" 
        file is computed. To use another annotation, please load an 
        annotation before using this method.

        Example:
            >>> g = Genome.Genome("dickeya")
            >>> g.load_annotation(annot_file="sequence.gff3")
            >>> g.load_gene_orientation()
            >>> g.genes["Dda3937_00005"].orientation
            'tandem'
            >>> g.orientation["isolated"]
            ['Dda3937_02360', 'Dda3937_01107', 'Dda3937_03898', 'Dda3937_04216',
            'Dda3937_00244', 'Dda3937_04441', 'Dda3937_01704', 'Dda3937_01705',
            ...
            'Dda3937_02126', 'Dda3937_01530', 'Dda3937_04419', 'Dda3937_02081']
        """
        self.load_neighbor_all()
        res = {"tandem": [], "divergent": [], "convergent": [], "isolated": []}
        for gene in self.genes:
            orient = ""
            try:
                g = self.genes[gene]
                lg = self.genes[g.left_neighbor]
                rg = self.genes[g.right_neighbor]

                #circular DNA
                if lg.start > g.start :
                    rstart = rg.start + self.length
                    gstart = g.start + self.length
                elif g.start > rg.start :
                    rstart = rg.start + self.length
                    gstart = g.start
                else : 
                    rstart = rg.start
                    gstart = g.start


                if couple == 3:
                    if ((gstart - lg.start) < max_dist 
                        and (rstart - gstart) < max_dist):
                        if not lg.strand and not rg.strand:
                            orient = "tandem"
                        elif lg.strand and rg.strand:
                            orient = "tandem"
                        elif lg.strand and not rg.strand:
                            orient = "convergent"
                        elif not lg.strand and rg.strand:
                            orient = "divergent"
                    else:
                        orient = "isolated"

                elif couple == 2:  
                    if g.strand:
                        if (gstart - lg.start) < max_dist:
                            if lg.strand:
                                orient = "tandem"
                            elif not lg.strand:
                                orient = "divergent"
                    if not g.strand:
                        if (rstart - gstart) < max_dist:
                            if rg.strand:
                                orient = "divergent"
                            elif not rg.strand:
                                orient = "tandem"  

                self.genes[gene].add_orientation(orient)
                res[orient].append(gene)                   
            except BaseException as e:
                print(f"Warning with locus {gene}: {e}")
        self.orientation = res

    def load_pos_orientation(self, max_dist=5000):
        """
        Compute gene orientation with the following criteria : 
            divergent if left neighbor on - strand and 
                         right neighbor on + strand,
            convergent if left neighbor on + strand and 
                          right neighbor on - strand, 
            tandem if left and right neighbors on same strand 
            isolated if the distance between neighbors is higher than the  
                          maximal distance given as argument

        Arg: 
            self (Genome instance)
            max_dist (Optional [int.]): maximal distance between 2 genes start 
                                  positions for seeking neighbor (Default: 5kb)

        Outputs:
            self.pos_orientation (dict of dic): new attribute of the Genome instance.
                                      self.pos_orientation is a dictionary of 
                                      containing 2 subdictionaries. One dictionary 
                                      for "intergenic" positions and one for 
                                      "intragenic" positions. 
                                      Each subdictionary contains the list of 
                                      position for each orientation.  
                                      {"intergenic":{orientation: list of positions}},
                                       "intragenic":{orientation: list of positions}}
                                       with orientation in ["divergent","convergent",
                                       "tandem","isolated"]
            self.orientation_per_pos (dict.): new attribute of the Genome instance.
                                      self.orientation_per_pos is a dictionary of 
                                      shape {position: orientation}.              

        N.B.: This method needs a genomic annotation. If no annotation is 
        loaded, the load_annotation method with the default "sequence.gff3" 
        file is computed. To use another annotation, please load an 
        annotation before using this method.

        Example:
            >>> g = Genome.Genome("dickeya")
            >>> g.load_annotation(annot_file="sequence.gff3")
            >>> g.load_pos_orientation()
            >>> g.orientation_per_pos[30]
            'tandem'
            >>> g.pos_orientation["intragenic"]["isolated"]
            [11534,11535,11536,11537,11538,11539,11540,11541,11542,11543,...]
        """
        self.load_neighbor_all()
        res_inter = {"tandem": [], "divergent": [], "convergent": [], "isolated": []}
        res_intra = {"tandem": [], "divergent": [], "convergent": [], "isolated": []}
        self.orientation_per_pos = {}
        for pos in np.arange(self.length):
            orient = ""
            try:
                lg = self.genes[self.left_neighbor[pos]]
                rg = self.genes[self.right_neighbor[pos]]
                
                #circular DNA
                if lg.start > pos :
                    pstart = pos + self.length
                    rstart = rg.start + self.length
                elif pos > rg.start :
                    pstart = pos
                    rstart = rg.start + self.length
                else : 
                    pstart = pos
                    rstart = rg.start
   
                if ((pstart - lg.start) < max_dist and (rstart - pstart) < max_dist):
                    if not lg.strand and not rg.strand:
                        orient = "tandem"
                    elif lg.strand and rg.strand:
                        orient = "tandem"
                    elif lg.strand and not rg.strand:
                        orient = "convergent"
                    elif not lg.strand and rg.strand:
                        orient = "divergent"
                else:
                    orient = "isolated"  

                self.orientation_per_pos[pos] = orient
                if self.genomic_situation[pos] == "intragenic":
                    res_intra[orient].append(pos)
                else : 
                    res_inter[orient].append(pos)
            except BaseException as e:
                print(f"Warning with position {pos}: {e}")
        self.pos_orientation = {"intragenic":res_intra,"intergenic":res_inter}

    def load_TSS(self):
        """ 
        load_TSS loads a TSS annotation from a file present in the /TSS/ 
        directory. The information importation requires a TSS.info file, 
        containing column indices of each information in the data file and some 
        additional information, in the following order: 
        [0] Condition [1] Filename [2] Locus tag column [3] TSS_column
        [4] File start line [5] Separator [6] Strand column
        [7] Sigma column [8] Sites column
        
        Arg:
            self (Genome instance)

        Output: 
            self.TSSs (dict. of dict.): new attribute of the Genome instance. 
                                        self.TSSs is a dictionary of shape 
                                        {Condition: {TSS position: TSS object}}. 
                                        One subdictionary is created for each 
                                        condition listed in TSS.info file.
                                        Each TSS object is initialized with 
                                        the following attributes: 
                                        pos, genes, promoter, score, strand
                                        The promoter attribute is a dictionary
                                        containing, for each sigma factor (keys)
                                        another dictionary (value). The first 
                                        created key of this second dictionary is 
                                        "sites" and the associated value is 
                                        a tuple containing the positions of 
                                        promoter elements.
                                        See __init__ and add_promoter in the 
                                        TSS class for more details about each 
                                        attribute.
            self.TSSs['all_TSSs'] (subdictionary of self.TSSs): 
                                        additional subdictionary of self.TSSs, 
                                        of shape : 
                                        self.TSSs['all_TSSs']={TSSpos: [TSScond]}
                                        with [TSScond] the list of TSS conditions 
                                        where this TSS was found.
        
        N.B.: This method needs a genomic annotation. If no annotation is 
        loaded, the load_annotation method with the default "sequence.gff3" 
        file is computed. To use another annotation, please load an 
        annotation before using this method.

        Example: 
            >>> g = Genome.Genome("dickeya")
            >>> g.load_TSS()
            >>> g.TSSs["dickeya-btss"][4707030].promoter
            {'sigma70': {'sites': (4707039, 4707044, 4707060, 4707065)},
             'sigma32': {'sites': (4707037, 4707046, 4707062, 4707068)}}
            >>> g.TSSs["dickeya-btss"][4707030].strand
            False
        """
        self.TSSs = {} 
        self.TSSs['all_TSS'] = {}     
        if not hasattr(self, "genes"):
            self.load_annotation()

        path2dir = f"{basedir}data/{self.name}/TSS/"
        if Path(f"{path2dir}TSS.info").exists():
            with open(f"{path2dir}TSS.info", "r") as f:
                skiphead = next(f)  
                for line in f:
                    line = line.strip('\n').split('\t')
                    filename = line[1]
                    startline = int(line[4])
                    sep = line[5]
                    strand = int(line[6])
                    TSScol = int(line[3])
                    genescol = int(line[2]) if line[2] != "" else None
                    sigcol = int(line[7]) if line[7] != "" else None
                    sitescol = int(line[8]) if line[8] != "" else None
                    scorecol = int(line[9]) if line[9] != "" else None
                    try:  
                        self.TSSs[line[0]] = load_TSS_cond(self.genes,
                                                           path2dir+filename,
                                                           TSScol, startline,
                                                           sep, strand,
                                                           genescol, sigcol,
                                                           sitescol, scorecol)
                        # appends all entries to the "all_TSS" subdictionary
                        for TSSpos in self.TSSs[line[0]].keys():
                            if TSSpos in self.TSSs['all_TSS'].keys():  
                                self.TSSs['all_TSS'][TSSpos].append(line[0])
                            else :
                                self.TSSs['all_TSS'][TSSpos] = [line[0]]
                    except Exception as e:
                        print("Error loading", line[0], e)
        else:
            print("No TSS.info, unable to load TSS")


    def load_prom_elements(self,shift=0,prom_region=[0,0]):
        """
        load_prom_elements extracts sequences of the different promoter elements 
        (spacer, -10, -35, discriminator, region around TSS) based on -35 and 
        -10 coordinates loaded with load_TSS method, for all TSS conditions.
        
        Args: 
            self (Genome instance)
            shift (Optional [int.]): number of nucleotides to include beyond each 
                                     region on either side (Default: 0nt)
            prom_region (Optional [int.,int.]): region upstream and downstream TSSs  
                                                to extract. Argument with the shape: 
                                                [length before TSS, length after TSS] 
                                                Default: [0,0] ie no sequence will 
                                                be extracted around TSS

        Outputs:
            For each sigma factor associated to a TSS annotation, creates a 
            subdictionary in self.TSSs[condTSS][TSS].promoter with the shape 
            promoter[sigma] = {element: sequence} with element in 
            ["spacer", "minus10", "minus35", "discriminator", "region"]
            See load_TSS description to understand the structure of the 
            subdictionary self.TSSs[condTSS][TSS].promoter

        N.B.: This method needs a genomic annotation. If no annotation is 
        loaded, the load_annotation method with the default "sequence.gff3" 
        file is computed. To use another annotation, please load an 
        annotation before using this method.

        Example:
            >>> g = Genome.Genome("dickeya")
            >>> g.load_prom_elements()
            >>> g.TSSs["dickeya-btss"][4707030].promoter
            {'sigma70': {'sites': (4707039, 4707044, 4707060, 4707065),
                         'spacer': 'TCGCCCACCCTCAAT',
                         'minus10': 'CATCAT',
                         'minus35': 'TACCCC',
                         'discriminator': 'GAATAACC'},
             'sigma32': {'sites': (4707037, 4707046, 4707062, 4707068),
                         'spacer': 'CCTCGCCCACCCTCA',
                         'minus10': 'ATCATCATGA',
                         'minus35': 'CCGTACC',
                         'discriminator': 'ATAACC'}}
        """
        if not hasattr(self, 'TSSs'): 
            self.load_TSS()
        if not hasattr(self, 'seq'): 
            self.load_seq()
        for cond_TSS in self.TSSs.keys():
            try:
                if cond_TSS != 'all_TSS':
                    for TSS in self.TSSs[cond_TSS].keys():
                        self.TSSs[cond_TSS][TSS].add_prom_elements(
                            self.seq, self.seqcompl, shift=shift, prom=prom_region)
            
            except BaseException:
                print('Unable to load prom elements for:', cond_TSS)


    def load_TU(self):
        """ 
        load_TU loads a TU annotation from a file present in the /TU/ 
        directory. The information importation requires a TU.info file, 
        containing column indices of each information in the data file and some 
        additional information, in the following order: 
        [0] Condition [1] Filename [2] Start column [3] Stop column
        [4] Strand column [5] Gene column [6] File start line [7] Separator 
        
        Arg:
            self (Genome instance)

        Output: 
            self.TUs (dict. of dict.): new attribute of the Genome instance. 
                                       self.TUs is a dictionary of shape 
                                       {Condition: {TU start pos: TU object}} 
                                       One subdictionary is created for each 
                                       condition listed in TU.info file.
                                       Each TU object is initialized with 
                                       the following attributes: 
                                       start, stop, orientation, genes,
                                       left,right.                                   
                                       See __init__ in the TU class for more 
                                       details about each attribute.

        Example: 
            >>> g = Genome.Genome("dickeya")
            >>> g.load_TU()
            >>> g.TUs["TU_Forquet"][2336702].right
            2339864
            >>> g.TUs["TU_Forquet"][2336702].strand
            True
        """
        self.TUs = {}
        path2dir = f"{basedir}data/{self.name}/TU/"
        if Path(f"{path2dir}TU.info").exists():
            with open(f"{path2dir}TU.info", "r") as f:
                skiphead = next(f) 
                for line in f:
                    line = line.strip().split('\t')
                    try:
                        self.TUs[line[0]] = load_TU_cond(path2dir+line[1],
                                                int(line[2]), int(line[3]), 
                                                int(line[4]), int(line[5]), 
                                                int(line[6]), line[7])
                    except BaseException:
                        print("Error loading cond", line[0])
            f.close()
        else:
            print("No TU.info file, please create one")


    def load_TTS(self):
        """ 
        load_TTS loads a TTS annotation from a file present in the /TTS/ 
        directory. The information importation requires a TTS.info file, 
        containing column indices of each information in the data file and  
        some additional information, in the following order: 
        [0] Condition [1] Filename [2] Left coordinate column
        [3] Right coordinate column [4] Strand column [5] Startline  
        [6] Separator [7] Sequence column 
        and optionally :
        [8] Score column [9] Genes column [10] Rho dependency column

        Arg:
            self (Genome instance)

        Output: 
            self.TTSs (dict. of dict.): new attribute of the Genome instance. 
                                        self.TTSs is a dictionary of shape 
                                        {Condition: {TTS position: TTS object}}. 
                                        One subdictionary is created for each 
                                        condition listed in TTS.info file.
                                        Each TTS object is initialized with 
                                        the following attributes: 
                                        left, right,  start, end, strand, 
                                        rho_dpdt, genes, seq, score. 
                                        If the data do not contain information 
                                        about associated genes, sequence, or rho 
                                        dependency, the corresponding attributes
                                        will be initialized as "None".
                                        See __init__ in the 
                                        TTS class for more details about each 
                                        attribute.
            self.TTSs['all_TTSs'] (subdictionary of self.TTSs): 
                                        additional subdictionary of self.TTSs, 
                                        of shape : 
                                        self.TTSs['all_TTSs']={TTSpos: [TTScond]}
                                        with [TTScond] the list of TTS conditions 
                                        where this TTS was found.
        Example: 
            >>> g = Genome.Genome("dickeya")
            >>> g.load_TTS()
            >>> g.TTSs["RhoTerm"][2791688].left
            2791688
        """
        self.TTSs = {}

        path2dir = f"{basedir}data/{self.name}/TTS/"
        if Path(f"{path2dir}TTS.info").exists():
            with open(f"{path2dir}TTS.info", "r") as f:
                skiphead = next(f) 
                for line in f:
                    line = line.strip().split('\t')
                    seqcol = int(line[7]) if line[7] != "" else None
                    scorecol = int(line[8]) if line[8] != "" else None
                    genescol = int(line[9]) if line[9] != "" else None
                    try:
                        self.TTSs[line[0]] = load_TTS_cond(
                                                path2dir + line[1],
                                                line[6], int(line[5]), 
                                                int(line[2]), int(line[3]), 
                                                int(line[4]), int(line[10]), 
                                                seqcol, scorecol, genescol)
                    except Exception as e:
                        print("Error loading cond", line[0])
            f.close()
            self.TTSs["all"] = {}
            for x in self.TTSs.keys():
                self.TTSs["all"].update(self.TTSs[x])
        else:
            print("No TTS.info file, please create one")


    def load_GO(self):
        """
        load_GO loads file specified in GO.info to assign GO terms to genes.
        Other annotation systems such as COG, or domain assignment can also 
        be used. 
        The information importation requires a GO.info file, 
        containing column indices of each information in the data file and  
        some additional information, in the following order: 
        [0] Annotation system [1] Filename [2] Locus tag column 
        [3] GOterm column [4] Separator

        Arg:
            self (Genome instance)

        Output: 
            self.GO (dict. of dict.): new attribute of the Genome instance. 
                                      self.GO is a dictionary of shape 
                                      {annot_syst: {GOterm: list of genes}} 
                                      i.e. one subdictionary is created for 
                                      each annotation system (such as GOc,
                                      COG or domain) listed in GO.info file.
            self.genes[locus].GO (list): new attribute of Gene instances 
                                         related to the Genome instance given 
                                         as argument. List of terms (such as 
                                         GO terms) associated to the gene.
        
        N.B.: This method needs a genomic annotation. If no annotation is 
        loaded, the load_annotation method with the default "sequence.gff3" 
        file is computed. To use another annotation, please load an 
        annotation before using this method.

        Example: 
            >>> g = Genome.Genome("dickeya")
            >>> g.load_GO()
            >>> g.GO['GO']['GO:0000100']
            ['Dda3937_01944', 'Dda3937_03618']
            >>> g.genes["Dda3937_00004"].GO
            {'GO': ['GO:0005829','GO:0055114','GO:0046872','GO:0042802',
                    'GO:0000166','GO:0006177','GO:0003938','GO:0005829',
                    'GO:0055114','GO:0046872','GO:0042802','GO:0000166',
                    'GO:0006177','GO:0003938']}
        """
        if not hasattr(self, "genes"): 
            self.load_annotation()
        
        warn_locus = [] 
        self.GO = {}
        path2dir = f"{basedir}data/{self.name}/GO_analysis/"
        
        # Opens .info file to obtain the information required to correctly 
        # load the data for each annotation system. 
        if Path(f"{path2dir}GO.info").exists():
            with open(f"{path2dir}GO.info", "r") as info_file:
                skiphead = next(info_file) 
                for info_line in info_file:
                    info_line = info_line.strip().split('\t')
                    system = info_line[0]  
                    filename = info_line[1]
                    tagtype = info_line[2]
                    tagcol = int(info_line[3])
                    GOcol = int(info_line[4])
                    separator = info_line[5]
                    startline = int(info_line[6])
                    
                    # Loads the data: associates GO terms to each gene 
                    print(f"Loading {system}...")
                    GOterms = []
                    with open(path2dir + filename, 'r') as GO_file:
                        dict_GO = {}
                        i=0
                        while i < startline:
                            header = next(GO_file)
                            i+=1
                        for line in GO_file:
                            line = line.strip('\n')
                            if separator == '\\t':
                                line = line.split('\t')
                            else:
                                line = line.split(separator)
                            GOterms.append(line[GOcol])
                            if line[tagcol] in dict_GO.keys() : 
                                dict_GO[line[tagcol]].append(line[GOcol])
                            else :
                                dict_GO[line[tagcol]] = [line[GOcol]]

                    self.GO[system] = {}
                    GOterms = list(set(GOterms))
                    for t in GOterms:
                        self.GO[system][t] = []

                    for locus in self.genes.keys():
                        g = self.genes[locus]
                        try :
                            if tagtype == "ASAP" :
                                name = g.ASAP
                            elif tagtype == "locus_tag" :  
                                name = locus
                            else:
                                sys.exit("Unknown \"tagtype\" in GO.info file:"
                                         "accepted tagtypes are ASAP and locus_tag")
                            self.genes[locus].add_GO(system,dict_GO[name])
                            for t in dict_GO[name]:
                                if locus not in self.GO[system][t] :
                                    self.GO[system][t].append(locus)
                        except Exception :
                            warn_locus.append(locus)
                                
                    GO_file.close()
                    # Locuses associated with no GO term are listed in warn_locus and are
                    # printed as a warning for each annotation system
                    success = len(self.genes.keys())-len(warn_locus)
                    print(f"\t{success} genes were successfully associated with some GO terms")
                    if len(warn_locus) > 20:
                        print(f"\t{len(warn_locus)} genes were not associated with any GO term")
                    else:
                        print(f"\t{warn_locus} are associated with no GO term")
            info_file.close()
        else:
            print("No GO.info file, please create one")


    def load_genomic_RNASeq_cov(self, cond="all"):
        """
        Load RNASeq coverage to a Genome instance. The importation requires 
        a cov.info or cov_txt.info file. See load_cov_rnaseq method 
        (Transcriptome class) for the details.

        Args: 
            self (Genome instance)
            cond (Optional [list of str.]): selection of one or several conditions 
                                    (1 condition corresponds to 1 data file).
                                    By default : cond ='all' ie all available 
                                    coverages are loaded.
        Outputs:
            self.cov_rnaseq_pos: dictionary of shape {condition: cov+}
                                 with cov+ an array containing one signal data 
                                 per genomic position for the + strand
            self.cov_rnaseq_neg: idem for the - strand

        Example: 
            >>> g = Genome.Genome("ecoli)
            >>> g.load_genomic_RNASeq_cov()
            >>> g.GO["GOc"]["GO:0033177"]
            ['Dda3937_00149', 'Dda3937_00150', 'Dda3937_00151']
            >>>  g.cov_rnaseq_neg["WT"]
            array([0., 0., 0., ..., 0., 0., 0.])
        """
        tr = Transcriptome.Transcriptome(self.name)
        tr.load_cov_rnaseq(cond)
        self.cov_rnaseq_pos = tr.cov_rnaseq_pos
        self.cov_rnaseq_neg = tr.cov_rnaseq_neg


    def load_genomic_signal(self, cond='all',data_treatment=None,
                            replicates=False,signal_per_genes=False,
                            *args, **kwargs):
        """ 
        Load_genomic_signals loads a 1D distribution along the chromosome 
        (typically a CHIPSeq distribution) from a data file (such as a 
        .bedGraph file obtained with bamCompare) containing bin starting and 
        ending positions and signal in each bin. See Chipseq class for more
        details.

        Args: 
            self (Genome instance)
            cond (Optional [list of str.]): selection of one or several 
                                            conditions (1 condition corresponds
                                            to 1 data file).
                                            By default cond ='all' ie all 
                                            available signals are loaded.
            data_treatment (Optional [str.]): "binning", "smoothing" or None 
                                               Default: None
            replicates (Optional [Bool.]): if True, the average between the 
                                           conditions will be computed.
                                           Default: False
            signal_per_genes(Optional [Bool.]): if True, compute the mean
                                                signal for each Gene

        Keyword Args:
            average_name (str.)
            window ([int.] only required if data_treatment == "smoothing"): 
                window size to compute the moving average             
            binsize ([int.] only required if data_treatment == "binning"):
                binsize to compute the binning
        
        Outputs: 
            self.all_signals (dict.): dictionary of shape {condition : 
                array containing one value per genomic position}.
            A second dictionary of shape {condition : array containing one 
                value per genomic position}. Depending on the data_treatment,
                the output is either: 
                - self.signals with the selected condition(s) as key(s)
                - self.signals_average with the average_name given in input
                    as key
                - self.binned_signal with the condition(s) merged with 
                  the binsize as key(s) (example: WT_bin200b)
                - self.smoothed_signal with the condition(s) merged with 
                  the smoothing window as key(s) (example: WT_smooth200b)
            self.genes[locus].signal (float.): new attribute of Gene instances 
                                               related to the Genome instance 
                                               given as argument.Contains the 
                                               Gene mean signal. 
            self.signals_gene (dict. of dict.): Dictionary containing one 
                                                subdictionary per condition 
                                                given in input. Each subdictionary
                                                contains the signal (value) for 
                                                each gene (locus_tag as key).
        Example:
            >>> g = Genome.Genome(name="ecoli")
            >>> g.load_genomic_signal(cond = ["Bates_WT_R1","Bates_WT_R2","Bates_WT_R3"],
                                      average_name = "Bates_WT_2kb",data_treatment ="binning"
                                      replicates = True, signal_per_genes = True,
                                      binsize = 2000)
            {'Bates_WT_2kb': 0.03255951846079419}
            >>> g.signals_average['Bates_WT_2kb_test']
            array([-0.06270156, -0.06270156, -0.06270156, ...,  0.25471515,
            0.25471515,  0.25471515])
            """

        Chip = Chipseq.Chipseq(self.name)

        if signal_per_genes:
            if not hasattr(self, "signals_gene"):
                self.signals_gene = {}
            if not hasattr(self, "genes"):
                self.load_annotation()

        if not hasattr(self, "all_signals"):
                self.all_signals = {}

        if replicates == True:
            if cond == "all":
                sys.exit("Please select conditions")
            binsize = kwargs.get('binsize', 1)
            window = kwargs.get('window', 1)
            average_name = kwargs.get('average_name', str(
                datetime.now())[:-10].replace(" ", "_"))
            Chip.load_signals_average(
                list_cond=cond,
                average_name=average_name,
                data_treatment=data_treatment,
                binsize=binsize,
                window=window)
            if not hasattr(self, "signals_average"):
                self.signals_average = {}
            self.signals_average[average_name] = Chip.signals_average[average_name]
            self.all_signals[average_name] = self.signals_average[average_name]
            if hasattr(self,"length"):
                if self.length != len(self.signals_average[average_name]) :
                    print(f"WARNING : {average_name} signal length isn't equal to genomic sequence length")
            if signal_per_genes == True:
                if not hasattr(self, "signals_gene"):
                    self.signals_gene = {average_name: {}}
                else:
                    self.signals_gene[average_name] = {}
                for locus in self.genes.keys():
                    g = self.genes[locus]
                    s = np.mean(
                       self.signals_average[average_name][g.left:g.right + 1])
                    self.genes[locus].add_signal(average_name, s)
                    self.signals_gene[average_name][locus] = s

        elif data_treatment == "binning":
            binsize = kwargs.get('binsize')
            Chip.load_binned_signal(binsize, cond)
            if not hasattr(self, "binned_signal"):
                self.binned_signal = {}
            for cbin in Chip.binned_signal.keys():
                self.binned_signal[cbin] = Chip.binned_signal[cbin]
                self.all_signals[cbin] = self.binned_signal[cbin] 
                if hasattr(self,"length"):
                    if self.length != len(self.binned_signal[cbin]) :
                        print(f"WARNING : {csmoo} signal length isn't equal to genomic sequence length")
                if signal_per_genes == True:
                    if not hasattr(self, "signals_gene"):
                        self.signals_gene = {cbin: {}}
                    else:
                        self.signals_gene[cbin] = {}
                    for locus in self.genes.keys():
                        g = self.genes[locus]
                        s = np.mean(
                            self.binned_signal[cbin][g.left:g.right + 1])
                        self.genes[locus].add_signal(c_bin, s)
                        self.signals_gene[cbin][locus] = s

        elif data_treatment == "smoothing":
            window = kwargs.get('window')
            Chip.load_smoothed_signal(window, cond)
            if not hasattr(self, "smoothed_signal"):
                self.smoothed_signal = {}
            for csmoo in Chip.smoothed_signal.keys():
                self.smoothed_signal[csmoo] = Chip.smoothed_signal[csmoo]
                self.all_signals[csmoo] = self.smoothed_signal[csmoo]
                if hasattr(self,"length"):
                    if self.length != len(self.smoothed_signal[csmoo]) :
                        print(f"WARNING : {csmoo} signal length isn't equal to genomic sequence length")
                if signal_per_genes == True:
                    if not hasattr(self, "signals_gene"):
                        self.signals_gene = {csmoo: {}}
                    else:
                        self.signals_gene[csmoo] = {}
                    for locus in self.genes.keys():
                        g = self.genes[locus]
                        s = np.mean(
                            self.smoothed_signal[csmoo][g.left:g.right + 1])
                        self.genes[locus].add_signal(csmoo, s)
                        self.signals_gene[csmoo][locus] = s

        else:
            print("No data treatment selected")
            Chip.load_signal(cond)
            if not hasattr(self, "signal"):
                self.signal = {}
            for c in Chip.signal.keys():
                self.signal[c] = Chip.signal[c]
                self.all_signals[c] = self.signal[c]
                if hasattr(self,"length"):
                    if self.length != len(self.signal[c]) :
                        print(f"WARNING : {c} signal length isn't equal to genomic sequence length")
                if signal_per_genes == True:
                    if not hasattr(self, "signals_gene"):
                        self.signals_gene = {c: {}}
                    else:
                        self.signals_gene[c] = {}
                    for locus in self.genes.keys():
                        g = self.genes[locus]
                        s = np.mean(self.signal[c][g.left:g.right + 1])
                        self.genes[locus].add_signal(c, s)
                        self.signals_gene[c][locus] = s

    def load_state_from_FC(self, thresh_pval=0.05,thresh_fc=0):
        """
        Loads Fold-changes and p-values data and computes genes state from these data.
        If no p-values is available in the data file, the default value of 0 is assigned 
        to each gene pvalue. 
        A gene is considered :
        - activated if its FC is above the FC threshold given as argument 
                   and its pvalue is below the pvalue threshold given as argument
                   ie. FC > thresh_FC and pval < thresh_pval
        - repressed if its FC is below the opposite of the FC threshold 
                    and its pvalue is below the pvalue threshold 
                    ie. FC < -thresh_FC and pval < thresh_pval
        - not affected either if its pvalue is above the threshold, 
                           or if its FC is between the - thresh_FC and + thresh_FC
        The FC and pvalues data importation requires an FC.info file, in the fold_changes 
        directory, containing column indices of each information in the data file, in 
        the following order: 
        [0] Condition [1] Filename [2] Locus_tag column [3] FC columns 
        [4] Separator [5] File start line [6] P-value column

        Args: 
            self (Genome instance)
            thresh_pval (Float): pvalue threshold used for the genes classification
            thresh_fc (Float): fold-change threshold used for the genes classification

        
        Outputs:
            self.genes[locus].state (dict.): new attribute of Gene instances 
                                            related to the Genome instance given 
                                            as argument. Dictionary of shape 
                                            {condition: state} with state either
                                            - 'act' if the gene is activated
                                            - 'rep' if the gene is repressed
                                            - 'non' if the gene is not affected
                                            - 'null' if the gene is not present 
                                                in the data
            self.statesFC (dict. of dict.): new attribute of the Genome instance.
                                            Dictionary containing one subdictionary 
                                            per condition listed in fc.info. Each 
                                            subdictionary contains the list of genes
                                            corresponding to each state ('act', 'rep',
                                             'non' or 'null').
            self.genes[locus].fc_pval (dict.): new attribute of Gene instances 
                                               related to the Genome instance given 
                                               as argument. Dictionary of shape 
                                               {condition: (FC,pvalue)}. 

        N.B.: This method needs a genomic annotation. If no annotation is 
        loaded, the load_annotation method with the default "sequence.gff3" 
        file is computed. To use another annotation, please load an 
        annotation before using this method.

        Example: 
            >>> g = Genome.Genome("ecoli)
            >>> g.load_state_from_FC()
            >>> g.genes['b0001'].state
            {'osmotic': 'rep', 'acidic_1mn': 'act'}
            >>> g.genes['b0001'].fc_pval
            {'osmotic': (-1.73717046009437, 0.0), 'acidic_1mn': (1.73, 0.0)}
            >>> g.statesFC['osmotic']['rep']
            ['b0001', 'b0002', 'b0003', 'b0004',...]
        """
        tr = Transcriptome.Transcriptome(self.name)
        if not hasattr(self, 'genes'):
            self.load_annotation()

        tr.compute_state_from_fc(thresh_pval=thresh_pval, thresh_fc=thresh_fc)
        for g in self.genes.keys():
            self.genes[g].state = tr.genes[g].state
            try : 
                self.genes[g].fc_pval = tr.genes[g].fc_pval
            except :
                pass

        self.statesFC = tr.statesFC

    def load_genomic_expression(self):
        """
        load_genomic_expression loads expression data from files present 
        in the expression directory. The information importation requires an
        expression.info file, containing column indices of each information in 
        the data file and some additional information, in the following order: 
        [0] Filename [1] Locus_tag column [2] Expression column
        [3] is Log ? [4] Strand column    [5] Separator 

        Arg:
            self (Genome instance)

        Output: 
            self.genes[locus].expression (dict.): new attribute of Gene instances 
                                                  related to the Genome instance 
                                                  given as argument.Dictionary of 
                                                  shape {condition : expression
                                                  level (float.)}

        Example: 
            >>> g = Genome.Genome("dickeya")
            >>> g.load_genomic_expression()
            >>> gg.genes["Dda3937_00004"].expression
            {'WT_nov_stat(E10)_rpkm': 66.6413789332929,
            'WT_stat(E3)': 6.980245227157392,
            'WT_PGA_stat(E13)_rpkm': 13.9428053948966}
        """
        tr = Transcriptome.Transcriptome(self.name)

        tr.load_expression()

        for g in self.genes.keys():
            try:
                self.genes[g].expression = tr.genes[g].expression
            except BaseException:
                print(f"{g} has no \"expression\" attribute ")

    def load_loops(self, cond='all', per_genes=True, window=0):
        """
        Load loops genomic positions determined with HiC

        """
        if not hasattr(self, "seq"):
            self.load_seq()

        if per_genes:
            if not hasattr(self, "genes"):
                self.load_annotation()
            if not hasattr(self, "loops_genes"):
                self.loops_genes = {}

        if not hasattr(self, "loops_pos"):
            self.loops_pos = {}

        HC = HiC.HiC(self.name)
        HC.load_hic_loops()

        if isinstance(cond, str):
            if cond == 'all':
                cond = list(HC.loops.keys())
            else:
                cond = [cond]

        for c in cond:
            print(f"Loading loops in: {c}")

            self.loops_pos[c] = {"loops": [], "no_loops": []}
            loops_list = list(HC.loops[c].keys())
            binsize = HC.loops[c][loops_list[0]]["binsize"]

            is_loop = [False] * self.length

            for loop in loops_list:
                for pos in list(np.arange(
                        loop[0], loop[0] + binsize + 1)) + list(np.arange(loop[1], loop[1] + binsize + 1)):
                    is_loop[pos] = True

            for x in np.arange(self.length):
                if is_loop[x]:
                    self.loops_pos[c]["loops"].append(x)
                else:
                    self.loops_pos[c]["no_loops"].append(x)

            if per_genes:
                self.loops_genes[f"{c}_w{window}b"] = {
                    "loops": [], "no_loops": []}
                for locus in self.genes.keys():
                    self.genes[locus].add_is_loop(c, False)
                    gene = self.genes[locus]
                    for pos in np.arange(gene.left - window,
                                         gene.right + 1 + window):
                        if gene.right + 1 + window >= self.length:
                            pos = pos - self.length
                        if is_loop[pos]:
                            self.genes[locus].add_is_loop(c, True)
                    if self.genes[locus].is_loop[c]:
                        self.loops_genes[f"{c}_w{window}b"]["loops"].append(
                            locus)
                    else:
                        self.loops_genes[f"{c}_w{window}b"]["no_loops"].append(
                            locus)

    def load_borders(self, cond="all", per_genes=True, window=0):
        """
        Load borders genomic positions determined with HiC
        """

        if not hasattr(self, "seq"):
            self.load_seq()

        if per_genes:
            if not hasattr(self, "genes"):
                self.load_annotation()
            if not hasattr(self, "borders_genes"):
                self.borders_genes = {}

        if not hasattr(self, "borders_pos"):
            self.borders_pos = {}

        HC = HiC.HiC(self.name)
        HC.load_hic_borders()
        if isinstance(cond, str):
            if cond == 'all':
                cond = list(HC.borders.keys())
            else:
                cond = [cond]
        for c in cond:
            print(f"Loading borders in: {c}")         
            borders_list = list(HC.borders[c].keys())
            binsize = HC.borders[c][borders_list[0]]["binsize"]
            pos_borders = []
            is_border = [False] * self.length
            for border in borders_list:
                ps = np.arange(border, border + binsize + 1)
                pos_borders.extend(ps)
                for p in ps:
                    is_border[p] = True
            pos_no_borders = set(np.arange(self.length)) - set(pos_borders)
            self.borders_pos[c] = {"borders":pos_borders, "no_borders": pos_no_borders}

            if per_genes:
                self.borders_genes[f"{c}_w{window}b"] = {
                    "borders": [], "no_borders": []}
                for locus in self.genes.keys():
                    self.genes[locus].add_is_border(c, False)
                    gene = self.genes[locus]
                    for pos in np.arange(gene.left - window,
                                         gene.right + 1 + window):
                        if gene.right + 1 + window >= self.length:
                            pos = pos - self.length
                        if is_border[pos]:
                            self.genes[locus].add_is_border(c, True)
                    if self.genes[locus].is_border[c]:
                        self.borders_genes[f"{c}_w{window}b"]["borders"].append(
                            locus)
                    else:
                        self.borders_genes[f"{c}_w{window}b"]["no_borders"].append(
                            locus)


    def load_peaks_genome(self, cond="all", per_genes=True, window=0):
        """
        
        """
        if not hasattr(self, "seq"):
            self.load_seq()

        if per_genes:
            if not hasattr(self, "genes"):
                self.load_annotation()
            if not hasattr(self, "peaks_genes"):
                self.peaks_genes = {}

        if not hasattr(self, "peaks_pos"):
            self.peaks_pos = {}

        Chip = Chipseq.Chipseq(self.name)
        Chip.load_peaks()

        if isinstance(cond, str):
            if cond == 'all':
                cond = list(Chip.peaks.keys())
            else:
                cond = [cond]

        for c in cond:
            print(f"Loading peaks in: {c}")
            pos_peaks = []
            is_peak = [False]*self.length
            for start,end in Chip.peaks[c]:
                pos_peaks.extend(np.arange(start,end+1))
                for p in np.arange(start,end+1) :
                    is_peak[p] = True
            pos_no_peaks = set(np.arange(self.length)) - set(pos_peaks)
            self.peaks_pos[c] = {"peaks":pos_peaks, "no_peaks": pos_no_peaks}

            if per_genes:
                self.peaks_genes[c] = {"peaks": [], "no_peaks": []}
                for locus in self.genes.keys():
                    self.genes[locus].add_is_peak(c, False)
                    gene = self.genes[locus]
                    for pos in np.arange(gene.left - window,
                                         gene.right + 1 + window):
                        if gene.right + 1 + window >= self.length:
                            pos = pos - self.length
                        if is_peak[pos]:
                            self.genes[locus].add_is_peak(c, True)
                    if self.genes[locus].is_peak[c]:
                        self.peaks_genes[c]["peaks"].append(locus)
                    else:
                        self.peaks_genes[c]["no_peaks"].append(locus)









###################### A ENLVER POUR PUBLI #######################

    """
    def predict_promoter_from_TSS(self,list_TSS,*args,**kwargs): #running bTSSfinder
        freedom = kwargs.get('free',0)
        nameOut = kwargs.get('out',list_TSS+'-btss-'+str(freedom))
        '''
        kwags possible : out & free
        bTSSfinder need to be install on this computer
        for more informations you can read btssfinder.py
        obj is Genome object
        list_TSS is the TSS list eg. biocyc
        out is  the name of fasta file (whiout extension .fasta) or if don't exist will become file's name of all out
        free is number of additionnal base at the TSS region
        NEXT
        convert out-bTSSfinder.gff on basedir/data/[nom_liste_TSS]
        AND FINALY
        write the localization of new csv in file TSS.info to the next load.TSS()
        '''
        try:
            test = self.TSSs[list_TSS]
            test = self.seq[0:4]
        except:
            self.load_TSS()
            self.load_seq()

        run_btssfinder(self,list_TSS,nameOut,freedom)

        gff2csv(self,list_TSS,nameOut,freedom)
        TSSinfo = basedir+"data/"+self.name+"/TSS/TSS.info"
        if Path(TSSinfo).exists():
            exist = False
            f = open(TSSinfo,"r")
            for i in f.readlines():
                line = i.split('\t')
                if line[0] == nameOut:
                    exist = True
            f.close()
            if not exist:
                f = open(TSSinfo,"a")
                f.write(nameOut+'\t'+nameOut+'.csv'+'\t'+"2"+'\t'+"0"+'\t'+"2"+'\t'+"\\t"+'\t'+"1"+'\t'+"3"+'\t'+"4"+'\n')
                f.close()
        else:
            print("TSS info not found")
        print("Finished…"+'\n'+"Now, you can visualise file "+TSSinfo+" or you can just reload TSS list.")
            

    def load_SIST(self, start, end,*args, **kwargs):
        if not hasattr(self, 'SIST_profile'):
            self.SIST_profile={}
            self.load_seq()
        option = kwargs.get('option')
        if option:
            if option == 'A':
                self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq, option=option)
            elif option == 'Z':
                self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq, option=option)
            elif option == 'C':
                self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq, option=option)
            else:
                print("This option doesn't exist")
        else:
            self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq)



    def load_sites(self, *args, **kwargs):
        '''
        #Load sites from sites.info: left and right binding coordinates on genome, sequence of site, score, strand
        '''
        self.sites = {}
        path2dir = f"{basedir}data/{self.name}/sites/"
        if Path(f"{path2dir}sites.info").exists():
            with open(f"{path2dir}sites.info", "r") as f:
                skiphead = next(f)  # skip head
                for header in f:
                    header = header.strip().split('\t')
                    self.sites[header[0]] = load_sites_cond(path2dir+ header[1], header[0], header[2], int(
                        header[3]), int(header[4]), int(header[5]), int(header[6]), int(header[7]), int(header[8]))
            f.close()

    def load_rpkm(self):
        Load a RPKM file information where indice 0 = Condition
        1 = filename type, 2 = RPKM  column, 3 = Start line,
        4 = type of separator, 5=locus_column 
        if not (self.genes):
            self.load_annotation()

        if os.path.exists(basedir+"data/"+self.name+"/rpkm/rpkm.info"):
            with open(basedir+"data/"+self.name+"/rpkm/rpkm.info","r") as f:
                for line in f:
                    line = line.strip('\n')
                    line = line.split('\t')
                    self.genes=add_single_rpkm_to_genes(self.genes, basedir+"data/"+self.name+"/rpkm/"+line[1],line[0],int(line[2]),int(line[3]),line[4],int(line[5]))
        else:
            print(" no rpkm file in this folder ")


    def load_reads(self):
        '''
        Load paired end reads from .npz files that have been generated using process_bam_paired_end in useful_functions
        and which are described in reads.info
        New attribute reads : reads_pos & reads_neg, of shape {[condition] : .npy}, e.g. self.reads_pos[cond1]
        '''
        self.reads_pos = {} # reads on + strand
        self.reads_neg = {} # reads on - strand
        if not os.path.exists(basedir+"data/"+self.name+'/rnaseq_reads/reads.info'):
            print('Unable to locate reads.info in /rnaseq_reads/')
        else:
            # open info file
            with open(basedir+"data/"+self.name+'/rnaseq_reads/reads.info',"r") as f:
                header = next(f) # first line = header
                # one line / condition
                for line in f:
                    line=line.strip()
                    line=line.split('\t')
                    print('Loading condition',line[0])
                    self.reads_pos[line[0]] = np.load(basedir+"data/"+self.name+'/rnaseq_reads/'+line[1])["Rpos"]
                    self.reads_neg[line[0]] = np.load(basedir+"data/"+self.name+'/rnaseq_reads/'+line[1])["Rneg"]
            print('Done')

    """