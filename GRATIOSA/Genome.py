#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import operator
import numpy as np
from datetime import datetime
from pathlib import Path
from GRATIOSA.useful_functions_genome import *
from GRATIOSA.globvar import *
from GRATIOSA.Gene import Gene
from GRATIOSA.TSS_TTS_TU import TSS, TTS, TU


class Genome:

    '''
    The Genome class is one of the main classes of this package. It gathers all
    the attributes of a genome such as the sequence, the set of genes (with their
    functional annotations and their orientation) but also all the annotations of
    TSS, TU and TTS.

    Each Genome instance has to be initialized with an organism name
    
        >>> g = Genome.Genome("ecoli")
    '''

    def __init__(self, name):
        self.name = name

    def load_seq(self, filename="sequence.fasta"):
        """
        Load_seq loads DNA sequence from a .fasta file present in the
        main directory of the organism using useful_functions_genome.load_seq
        function. Adds this sequence, its complement, and its length to a
        Genome instance.
        
        Creates 3 new attributes to the Genome instance
            * seq (str.): genomic sequence compose of A,T,G and C
            * seqcompl (str.): complement sequence to seq
            * length (int.): length of the genomic sequence

        Args:
            filename (Optional [str.}): name of the file containing the DNA
                    sequence in FASTA format. Default: "sequence.fasta"

        Example:
            >>> g = Genome.Genome("dickeya")
            >>> g.load_seq()
            >>> g.seq[100:151]
            'AATGTCGATCTTCAACATATCGCCGATCCGACGGGCACCCAGATCCTGCAG'
            >>> g.seqcompl[100:151]
            'TTACAGCTAGAAGTTGTATAGCGGCTAGGCTGCCCGTGGGTCTAGGACGTC'
        """
        seq_dir = f"{basedir}data/{self.name}/sequence.fasta"
        self.seq = read_seq(seq_dir).upper()
        self.seqcompl = ''
        self.length = len(self.seq)
        l = self.seq
        l = l.replace('A', 't').replace('T', 'a')
        l = l.replace('C', 'g').replace('G', 'c')
        l = l.replace('a', 'A').replace('t', 'T')
        l = l.replace('c', 'C').replace('g', 'G')
        self.seqcompl = l


    def load_annotation(self, annot_file="sequence.gff3"):
        """
        load_annotation loads a gene annotation (coordinates, length,
        name...) from a file present in the /annotation/ directory.

        Creates a new attribute of the Genome instance: 
            * self.genes (dict.) 
                self.genes is a dictionary of shape 
                {locus_tags: Gene object}. Each Gene object is 
                initialized with the following attributes:
                locus_tag, ID, name,strand, left, right, start, end, 
                middle, length and ASAP. 

        Args:
            filename (Optional [str.]): name of the file containing the 
                    genomic annotation. Default: "sequence.gff3"

        Note: 
            If the file is a .gff3 or .gff information will be loaded using
            useful_functions_genome.load_gff.
            Else, the information importation requires an annotation.info file,
            containing column indices of each information in the data file and 
            some additional information, in the following order:
            [0] Filename [1] Separator [2] Locus_tag column [3] Name column
            [4] ID column [5] Strand column [6] Left coordinate column
            [7] Right coordinate column [8] File start line

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
        if load == True:
            print("Annotation loaded")


    def load_genes_per_pos(self, window=0):
        """
        load_genes_per_pos associates each position with a list of genes
        overlapping this and its surrounding positions delimited
        by the window size given as an argument. 

        Creates a new attribute of the Genome instance: 
            * self.genes_per_pos (dict.)
                    self.genes is a dictionary of shape {position: list of 
                    genes}. It contains, for each position p, the list of
                    genes overlapping any position between p-window/2 and 
                    p+window/2 (inclusive)

        Args:
            window (Optional [int.]): window size in b. load_genes_per_pos 
                    finds genes between pos-window/2 and pos+window/2
                    (inclusive). Default: window = 0

        Warning:
            This method needs a genomic annotation. If no annotation is
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
                if pos < 0:
                    pos += self.length
                if pos > self.length:
                    pos -= self.length
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
        
        Creates 2 new attributes of Gene instances:
            * self.genes[locus].left_neighbor 
                locus of the nearest left-side neighbor gene
            * self.genes[locus].right_neighbor 
                locus of the nearest right-side neighbor gene

        and 4 new attributes of Genome instance:
            * self.genomic_situation (dict.) 
                dict of shape {position: situation} with situation either 
                "intergenic" or "intragenic".
            * self.left_neighbor (dict.) 
                dict of shape {position: locus of the nearest left-side 
                neighbor gene}
            * self.right_neighbor (dict.) 
                dict of shape {position: locus of the nearest right-side 
                neighbor gene}
            * self.gene (dict.) 
                dict of shape {position: gene} with gene = "NA" if the 
                position is intergenic.
            
        Warning:
            This method needs a genomic annotation. If no annotation is 
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
        genes_order = [x for _, x in sorted(
            zip(list_start, self.genes.keys()))]

        self.genomic_situation = ["intergenic"] * self.length
        self.left_neighbor = {}
        self.right_neighbor = {}
        self.gene = ["NA"] * self.length

        for i in np.arange(len(genes_order)):
            # Adds right and left neighbors to each Gene object
            if i != 0:
                self.genes[genes_order[i]].add_left_neighbor(
                    genes_order[i - 1])
            else:
                self.genes[genes_order[i]].add_left_neighbor(genes_order[-1])
            if i != len(genes_order) - 1:
                self.genes[genes_order[i]].add_right_neighbor(
                    genes_order[i + 1])
            else:
                self.genes[genes_order[i]].add_right_neighbor(genes_order[0])

            # Adds right and left neighbors to each genomic position
            for p in np.arange(self.genes[genes_order[i]].left, self.genes[genes_order[i]].right):
                self.genomic_situation[p] = "intragenic"
                self.gene[p] = genes_order[i]
                if i != 0:
                    self.left_neighbor[p] = genes_order[i - 1]
                else:
                    self.left_neighbor[p] = genes_order[-1]
                if i != len(genes_order) - 1:
                    self.right_neighbor[p] = genes_order[i + 1]
                else:
                    self.right_neighbor[p] = genes_order[0]
            if i != len(genes_order) - 1:
                for p in np.arange(self.genes[genes_order[i]].right, self.genes[genes_order[i + 1]].left):
                    self.left_neighbor[p] = genes_order[i]
                    self.right_neighbor[p] = genes_order[i + 1]
            else:
                for p in np.arange(self.genes[genes_order[i]].right, self.length):
                    self.left_neighbor[p] = genes_order[i]
                    self.right_neighbor[p] = genes_order[0]
                for p in np.arange(0, self.genes[genes_order[0]].left):
                    self.left_neighbor[p] = genes_order[i]
                    self.right_neighbor[p] = genes_order[0]


    def load_gene_orientation(self, couple=3, max_dist=5000):
        """
        Compute gene orientation with the following criteria: 

        * If couple = 3, gene is considered:

            * `divergent` if left neighbor on - strand and right neighbor on + strand,
            * `convergent` if left neighbor on + strand and right neighbor on - strand,
            * `tandem` if left and right neighbors on same strand (whatever the strand of the given gene is),
            * `isolated` if the distance between neighbors is higher than the maximal distance given as argument.
        * If couple = 2, gene is considered:

            * `tandem` if predecessor (left neighbor for gene on + strand, right neighbor for gene on - strand) is on same strand,
            * `divergent` if the predecessor is on opposite strand.

        Creates new attributes: 
            * self.orientation (dict.) 
                    new attribute of the Genome instance.
                    self.orientation is a dictionary of shape {orientation: 
                    list of genes}. It contains the list of genes for each
                    orientation (tandem, divergent, convergent, and isolated 
                    if couple=3, tandem and divergent if couple=2)
            * self.genes[locus].orientation (str.) 
                    new attribute of Gene instances related to the Genome instance 
                    given as argument

        Args:
            couple (int.): number of genes to consider in a "couple". 

                * If couple = 2: computes the orientation of a gene 
                  relative to its predecessor
                * If couple = 3: computes the orientation of a gene 
                  relative to its two neighbors

                Default: 3

            max_dist (Optional [int.]): maximal distance between 2 genes 
                    start positions for seeking neighbor (Default: 5kb)            

        Warning:
            This method needs a genomic annotation. If no annotation is 
            loaded, the load_annotation method with the default "sequence.gff3"
            file is computed. To use another annotation, please load an
            annotation before using the load_neighbor_all method.

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

                # circular DNA
                if lg.start > g.start:
                    rstart = rg.start + self.length
                    gstart = g.start + self.length
                elif g.start > rg.start:
                    rstart = rg.start + self.length
                    gstart = g.start
                else:
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
        Computes gene orientation with the following criteria:

            * `divergent` if left neighbor on - strand and right neighbor on + strand,
            * `convergent` if left neighbor on + strand and right neighbor on - strand,
            * `tandem` if left and right neighbors on same strand
            * `isolated` if the distance between neighbors is higher than the maximal distance given as argument
        
        Creates 2 new attributes of the Genome instance:
            * self.pos_orientation (dict of dict)
                    dictionary containing 2 subdictionaries. 
                    One subdictionary for "intergenic" positions and one for 
                    "intragenic" positions.
                    Each subdictionary contains the list of position for 
                    each orientation. {"intergenic":{orientation: list of 
                    positions}}, "intragenic":{orientation: list of positions}}
                    with orientation in ["divergent","convergent","tandem",
                    "isolated"]
            * self.orientation_per_pos (dict.)
                    dictionary of shape {position: orientation}.

        Args:
            max_dist (Optional [int.]): maximal distance between 2 genes 
                    start positions for seeking neighbor (Default: 5kb)            

        Warning:
            This method needs a genomic annotation. If no annotation is
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
        res_inter = {
            "tandem": [],
            "divergent": [],
            "convergent": [],
            "isolated": []}
        res_intra = {
            "tandem": [],
            "divergent": [],
            "convergent": [],
            "isolated": []}
        self.orientation_per_pos = {}
        for pos in np.arange(self.length):
            orient = ""
            try:
                lg = self.genes[self.left_neighbor[pos]]
                rg = self.genes[self.right_neighbor[pos]]

                # circular DNA
                if lg.start > pos:
                    pstart = pos + self.length
                    rstart = rg.start + self.length
                elif pos > rg.start:
                    pstart = pos
                    rstart = rg.start + self.length
                else:
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
                else:
                    res_inter[orient].append(pos)
            except BaseException as e:
                print(f"Warning with position {pos}: {e}")
        self.pos_orientation = {
            "intragenic": res_intra,
            "intergenic": res_inter}


    def load_TSS(self):
        """
        load_TSS loads a TSS annotation from a file present in the /TSS/
        directory.

        Creates:
            * self.TSSs (dict. of dict.)
                    New attribute of the Genome instance
                    self.TSSs is a dictionary of shape 
                    {Condition: {TSS position: TSS object}}. 
                    One subdictionary is created for each condition listed in 
                    TSS.info file. Each TSS object is initialized with the 
                    following attributes: pos, genes, promoter, score, strand
                    The promoter attribute is a dictionary containing, for 
                    each sigma factor (keys) a subdictionary (value). The 
                    first created key of this subdictionary is "sites" and 
                    the associated value is a tuple containing the positions 
                    of promoter elements. See __init__ and add_promoter in 
                    the TSS class for more details about each attribute.

            * self.TSSs['all_TSSs'] (subdictionary of self.TSSs)
                    additional subdictionary of self.TSSs, of shape:
                    self.TSSs['all_TSSs']={TSSpos: [TSScond]} with [TSScond] 
                    the list of TSS conditions where this TSS was found.        

        Note:
            The information importation requires a TSS.info file,
            containing column indices of each information in the data file and 
            some additional information, in the following order:
            [0] Condition [1] Filename [2] Locus tag column [3] TSS position
            [4] File start line [5] Separator [6] Strand column
            [7] Sigma factor column [8] Sites column [9] Score column

        Note:
            If some of the data types are missing (locus tags, sigma factors,
            scores or sites), an empty space can be left in the .info file.

        Warning:
            This method needs a genomic annotation. If no annotation is
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
                                                           path2dir + filename,
                                                           TSScol, startline,
                                                           sep, strand,
                                                           genescol, sigcol,
                                                           sitescol, scorecol)
                        # appends all entries to the "all_TSS" subdictionary
                        for TSSpos in self.TSSs[line[0]].keys():
                            if TSSpos in self.TSSs['all_TSS'].keys():
                                self.TSSs['all_TSS'][TSSpos].append(line[0])
                            else:
                                self.TSSs['all_TSS'][TSSpos] = [line[0]]
                    except Exception as e:
                        print("Error loading", line[0], e)
        else:
            print("No TSS.info, unable to load TSS")


    def load_prom_elements(self, shift=0, prom_region=[0, 0]):
        """
        load_prom_elements extracts sequences of the different promoter 
        elements (spacer, -10, -35, discriminator, region around TSS) based 
        on -35 and -10 coordinates loaded with load_TSS method, for all TSS 
        conditions. For each sigma factor associated to a TSS annotation, creates 
        a subdictionary in self.TSSs[condTSS][TSS].promoter with the shape
        promoter[sigma] = {element: sequence of the element} with element 
        in ["spacer", "minus10", "minus35", "discriminator", "region"].

        Args:
            shift (Optional [int.]): number of nucleotides to include beyond 
                    each region on either side (Default: 0nt)
            prom_region (Optional [int.,int.]): region upstream and 
                    downstream TSSs to extract. Argument with the shape:
                    [length before TSS, length after TSS]. 
                    Default: [0,0] ie no sequence will be extracted around 
                    TSS.
        Note:
            See load_TSS description to understand the structure of the
            subdictionary self.TSSs[condTSS][TSS].promoter

        Warning:
            This method needs a genomic annotation. If no annotation is
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
        directory.

        Creates a new attribute of the Genome instance:
            * self.TUs (dict. of dict.)
                    self.TUs is a dictionary of shape
                    {Condition: {TU start pos: TU object}}
                    One subdictionary is created for each condition listed 
                    in TU.info file. Each TU object is initialized with the 
                    following attributes:  start, stop, orientation, genes,
                    left,right, expression. 
                    See __init__ in the TU class for more details 
                    about each attribute.

        Note:
            The information importation requires a TU.info file,
            containing column indices of each information in the data file and 
            some additional information, in the following order:
            [0] Condition [1] Filename [2] TU ID [3] Start column [4] Stop column
            [5] Strand column [6] Gene column [7] File start line [8] Separator
            [9] Expression (optionnal)
            See the load_TU_cond function in useful_functions_genome for more
            details.

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
                        exprcol = None
                        if len(line) > 9 :
                            try:
                                exprcol = int(line[9])
                            except:
                                pass 
                        self.TUs[line[0]] = load_TU_cond(path2dir + line[1],
                                                         int(line[2]), 
                                                         int(line[3]),
                                                         int(line[4]), 
                                                         int(line[5]),
                                                         int(line[6]), 
                                                         int(line[7]),
                                                         line[8],
                                                         exprcol)
                    except BaseException:
                        print("Error loading cond", line[0])
            f.close()
        else:
            print("No TU.info file, please create one")


    def load_TTS(self):
        """
        load_TTS loads a TTS annotation from a file present in the /TTS/
        directory.

        Creates:
            * self.TTSs (dict. of dict.)
                New attribute of the Genome instance
                self.TTSs is a dictionary of shape
                {Condition: {TTS position: TTS object}}.
                One subdictionary is created for each condition listed 
                in TTS.info file. Each TTS object is initialized with the 
                following attributes: left, right,  start, end, strand,                    
                rho_dpdt, genes, seq, score.
                If the data do not contain information about associated 
                genes, sequence, or rho dependency, the corresponding 
                attributes will be initialized as "None".
                See __init__ in the TTS class for more details about 
                each attribute.
            * self.TTSs['all_TTSs'] (subdictionary of self.TTSs)
                additional subdictionary of self.TTSs, of shape:
                self.TTSs['all_TTSs']={TTSpos: [TTScond]}
                with [TTScond] the list of TTS conditions where this TTS 
                was found.

        Note: 
            The information importation requires a TTS.info file,
            containing column indices of each information in the data file and
            some additional information, in the following order:
            [0] Condition [1] Filename [2] Left coordinate column
            [3] Right coordinate column [4] Strand column [5] Startline
            [6] Separator [7] Sequence column
            and optionally:
            [8] Score column [9] Genes column [10] Rho dependency column

        Note: 
            See the load_TTS_cond function in useful_functions_genome for more
            details.

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

        Creates:      
            * self.GO (dict. of dict)
                new attribute of the Genome instance. self.GO is a dictionary of 
                shape {annot_syst: {GOterm: list of genes}} i.e. one subdictionary 
                is created for each annotation system (such as GOc, COG or domain) 
                listed in GO.info file.
            * self.genes[locus].GO (list)
                new attribute of Gene instances related to the Genome instance given
                as argument. List of terms (such as GO terms) associated to the gene.            
        
        Note:           
            The information importation requires a GO.info file,
            containing column indices of each information in the data file and
            some additional information, in the following order:
            [0] Annotation system [1] Filename [2] Locus tag column
            [3] GOterm column [4] Separator
            GO.info and the data files have to be in the /GO/ directory

        Note:
            Other annotation systems such as COG, or domain assignment can also
            be used.

        Warning:
            This method needs a genomic annotation. If no annotation is
            loaded, the load_annotation method with the default "sequence.gff3"
            file is computed. To use another annotation, please load an
            annotation before using this method.

        Example
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
        path2dir = f"{basedir}data/{self.name}/GO/"

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
                        i = 0
                        while i < startline:
                            header = next(GO_file)
                            i += 1
                        for line in GO_file:
                            line = line.strip('\n')
                            if separator == '\\t':
                                line = line.split('\t')
                            else:
                                line = line.split(separator)
                            GOterms.append(line[GOcol])
                            if line[tagcol] in dict_GO.keys():
                                dict_GO[line[tagcol]].append(line[GOcol])
                            else:
                                dict_GO[line[tagcol]] = [line[GOcol]]

                    self.GO[system] = {}
                    GOterms = list(set(GOterms))
                    for t in GOterms:
                        self.GO[system][t] = []

                    for locus in self.genes.keys():
                        g = self.genes[locus]
                        try:
                            if tagtype == "ASAP":
                                name = g.ASAP
                            elif tagtype == "locus_tag":
                                name = locus
                            else:
                                sys.exit("Unknown \"tagtype\" in GO.info file:"
                                         "accepted tagtypes are ASAP and locus_tag")
                            self.genes[locus].add_GO(system, dict_GO[name])
                            for t in dict_GO[name]:
                                if locus not in self.GO[system][t]:
                                    self.GO[system][t].append(locus)
                        except Exception:
                            warn_locus.append(locus)

                    GO_file.close()
                    # Locuses associated with no GO term are listed in warn_locus and are
                    # printed as a warning for each annotation system
                    success = len(self.genes.keys()) - len(warn_locus)
                    print(f"\t{success} genes were successfully associated with some GO terms")
                    if len(warn_locus) > 20:
                        print(f"\t{len(warn_locus)} genes were not associated with any GO term")
                    else:
                        print(f"\t{warn_locus} are associated with no GO term")
            info_file.close()
        else:
            print("No GO.info file, please create one")
