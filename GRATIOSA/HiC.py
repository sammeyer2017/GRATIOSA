#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path
from GRATIOSA.globvar import *
from GRATIOSA.useful_functions_HiC import *
from GRATIOSA import Genome

class HiC:
    '''
    The HiC class is used to gather spatial information obtained from HiC data 
    analysis with tools such as Chromosight. The two patterns analyzed here are 
    borders (between two topological domains) and loops (more specifically the 
    positions of the two regions that come into contact).

    Each HiC instance has to be initialized with an organism name
    
        >>> hc = HiC.HiC("dickeya")
    '''

    def __init__(self, name):
        self.name = name

    def load_hic_borders(self, cond="all"):
        """
        load_hic_borders imports a list of borders from a data file (typically
        a borders.tsv data file obtained with Chromosight) containing at least,
        for each border, its position (bin number) and optionally its identification 
        score (pearson correlation coefficient between the border kernel and the 
        detected pattern), pvalue and qvalue.
        Creates 2 new attributes of the HiC instance:
            * self.borders (dict. of dict.) 
                    dictionary containing one subdictionary per condition 
                    with the shape self.borders[cond]={bin_nb: {"score": 
                    score, "pval": pvalue, "qval": qvalue,"binsize":binsize}}.
                    N.B.: If data are binned at 2kb, the bin with number 10 
                    corresponds to data between 20000 and 22000.    
            * self.borders_pos (dict. of dict)
                    dictionary containing one subdictionary per 
                    condition, with the shape self.borders_pos[cond] =
                    {'borders':list of positions that are in a border
                    'no_borders':list of positions that are not in a border}

        Args:
            cond (Optional [list of str.]): selection of one or several
                    conditions (1 condition corresponds to 1 data file and 
                    each condition has to be listed in the border.info file).
                    By default: cond ='all' ie all available conditions are 
                    loaded.
        Note:
            The data importation requires a borders.info file that contains the
            column indices of each information in the data file and some
            additional information, in the following order:
                * (required) [0] Condition, [1] Filename, [2] Startline,
                  [3] Separator, [4] Bin, [5] Binsize (in b)
                * (optional) [6] Score, [7] Pvalue,[8] Qvalue
            borders.info and data file have to be in the /HiC/Borders/ directory

        Note: 
            See Chromosight documentation for more details about the input format.

        Warning:
            Make sure that the version of the genome is the same as
            the version used in the HiC analysis workflow

        Example:
            >>> hc = HiC.HiC("dickeya")
            >>> hc.load_hic_borders("WT")
            >>> hc.borders_pos['WT']['no_borders'][0:5]
            [2001,2002,2003,2004,2005]
        """

        # gets the path to data and .info files
        path2dir = f"{basedir}data/{self.name}/HiC/Borders/"

        if isinstance(cond, str):
            cond = [cond]

        # tries to open borders.info file
        if Path(f"{path2dir}borders.info").exists():
            with open(f"{path2dir}borders.info", "r") as f:
                skiphead = next(f)
                # 1 line in borders.info corresponds to 1 condition (1 data
                # file)
                for line in f:
                    line = line.strip('\n').split('\t')
                    # loads data file information for this condition
                    if cond == ["all"] or line[0] in cond:

                        path2file = f"{basedir}data/{self.name}/HiC/Borders/{line[1]}"
                        startline, sep = int(line[2]), line[3]
                        bin_col = line[4]

                        # loads binsize
                        binsize = int(line[5])

                        # loads column indices of each optional information
                        print(line[0])
                        score_col, pval_col, qval_col = None, None, None
                        if len(line) > 6:
                            try : 
                                score_col = int(line[6])
                            except :
                                pass
                        if len(line) > 7:
                            try : 
                                pval_col = int(line[7])
                            except :
                                pass
                        if len(line) > 8:
                            try : 
                                qval_col = int(line[8])
                            except :
                                pass

                        # creates the new attribute using the load_HiC_site_cond function
                        # from useful_functions_HiC
                        if not hasattr(self, "borders"):
                            self.borders = {}  # creates the new attribute
                        self.borders[line[0]] = load_HiC_site_cond(
                            'borders', path2file, startline, sep, bin_col, binsize,
                            score_col=score_col, pval_col=pval_col, qval_col=qval_col)
        else:
            print("No borders.info, unable to load borders")

        if cond == ['all']:
            cond = list(self.borders.keys())

        if not hasattr(self, "length"):
            g = Genome.Genome(self.name)
            g.load_seq()
            self.length = g.length

        if not hasattr(self, "borders_pos"):
            self.borders_pos = {}

        for name in cond:
            print(f"Loading borders in: {name}")
            borders_list = list(self.borders[name].keys())
            binsize = self.borders[name][borders_list[0]]["binsize"]
            pos_borders = []
            is_border = [False] * self.length
            for border in borders_list:
                ps = np.arange(border, border + binsize + 1)
                pos_borders.extend(ps)
                for p in ps:
                    is_border[p] = True
            pos_no_borders = set(np.arange(self.length)) - set(pos_borders)
            self.borders_pos[name] = {
                "borders": list(pos_borders),
                "no_borders": list(pos_no_borders)}

    def load_hic_loops(self, cond="all"):
        """
        load_hic_loops imports a list of loops from a data file (typically
        a loops.tsv data file obtained with Chromosight) containing, for each 
        loop, at least the loop coordinates (2 bin numbers) and optionally the 
        loop score (pearson correlation coefficient between the loop kernel and 
        the detected pattern), the p-value and the q-value.

        Creates 2 new attributes of the HiC instance:
            * self.loops (dict. of dict.)
                dictionary containing one subdictionary per condition with 
                the shape self.loops[cond] =
                {(bin1_nb,bin2_nb): {"score": score, "pval": pvalue,
                "qval": qvalue, "binsize":binsize}}
                N.B.: If data are binned at 2kb, the bin with number 10 
                corresponds to data between 20000 and 22000.
            * self.loops_pos (dict. of dict)
                dictionary containing one subdictionary per condition, with 
                the shape self.loops_pos[cond] =
                {'loops':list of positions that are in a loop,
                'no_loops':list of positions that are not in a loop}

        Args:
            self (HiC instance)
            cond (Optional [list of str.]): selection of one or several 
                    conditions (1 condition corresponds to 1 data file and 
                    has to be in the loops.info file). By default: 
                    cond ='all' ie all available data are loaded.

        Note:
            The data importation requires a loops.info file that contains the
            column indices of each information in the data file and some 
            additional information, in the following
            order:
                * (required) [0] Condition, [1] Filename, [2] Startline,
                  [3] Separator, [4] Bin1, [5] Bin2, [6] Binsize (in b),
                * (optional) [7] Score, [8] Pvalue,[9] Qvalue
            loops.info and data file have to be in the /HiC/Loops/ directory
        
        Note: 
            See Chromosight documentation for more details about the input format.

        Warning:
            Make sure that the version of the genome is the same as
            the version used in the HiC analysis workflow

        Example:
            >>> hc = HiC.HiC("dickeya")
            >>> hc.load_hic_loops()
            >>> hc.loops_pos['WT']['no_borders'][0:5]
            [2001,2002,2003,2004,2005]
            >>> hc.loops['WT1']
            {(42000, 80000): {'binsize': 2000},
             (238000, 252000): {'binsize': 2000}
             ...}
        """
        # gets the path to data and .info files
        path2dir = f"{basedir}data/{self.name}/HiC/Loops/"

        if isinstance(cond, str):
            cond = [cond]

        # tries to open loops.info file
        if Path(f"{path2dir}loops.info").exists():
            with open(f"{path2dir}loops.info", "r") as f:
                skiphead = next(f)  # skip head

                # 1 line in loops.info corresponds to 1 condition (1 file)
                for line in f:
                    line = line.strip('\n')
                    line = line.split('\t')

                    if cond == ["all"] or line[0] in cond:

                        # loads data file information for this condition
                        path2file = f"{path2dir}{line[1]}"
                        startline, sep = int(line[2]), line[3]
                        bin1_col, bin2_col = line[4], line[5]

                        # loads binsize
                        binsize = int(line[6])

                        # loads column indices of each optional information
                        score_col, pval_col, qval_col = None, None, None
                        if len(line) > 7:
                            try : 
                                score_col = int(line[7])
                            except :
                                pass
                        if len(line) > 8:
                            try : 
                                pval_col = int(line[8])
                            except :
                                pass
                        if len(line) > 9:
                            try : 
                                qval_col = int(line[9])
                            except :
                                pass
                        # creates the new attribute using the load_HiC_site_cond 
                        # function from useful_functions_HiC

                        if not hasattr(self, "loops"):
                            self.loops = {}  # creates the new attribute
                        self.loops[line[0]] = load_HiC_site_cond(
                            "loops", path2file, startline, sep, bin1_col, binsize,
                            bin2_col=bin2_col, score_col=score_col,
                            pval_col=pval_col, qval_col=qval_col)

        else:
            print("No loops.info, unable to load loops")
        if cond == ['all']:
            cond = list(self.loops.keys())

        if not hasattr(self, "loops_pos"):
            self.loops_pos = {}

        if not hasattr(self, "length"):
            g = Genome.Genome(self.name)
            g.load_seq()
            self.length = g.length

        # register loops genomic positions in the loops_pos attribute
        for c in cond:
            print(f"Loading loops in: {c}")

            self.loops_pos[c] = {"loops": [], "no_loops": []}
            loops_list = list(self.loops[c].keys())
            binsize = self.loops[c][loops_list[0]]["binsize"]

            is_loop = [False] * self.length

            for loop in loops_list:
                for pos in list(np.arange(loop[0], loop[0] + binsize + 1)) + list(np.arange(loop[1], loop[1] + binsize + 1)):
                    is_loop[pos] = True

            for x in np.arange(self.length):
                if is_loop[x]:
                    self.loops_pos[c]["loops"].append(x)
                else:
                    self.loops_pos[c]["no_loops"].append(x)

    def load_loops_genes(self, cond='all', window=0):
        """
        load_loops_genes imports a list of loops determined with HiC from a 
        data file (typically a loops.tsv data file obtained with Chromosight)
        containing at least the loops coordinates (2 bin numbers) and 
        determines the list of genes overlapping one of the loops positions.

        Creates:
            self.loops_pos (dict. of dict.) 
                    new attribute of the Genome 
                    instance. Dictionary containing one dictionary per 
                    condition listed in loops.info. Each subdictionary has 
                    2 keys: "loops" and "no_loops" and contains the 
                    corresponding lists of genomic positions.
            self.loops_genes (dict. of dict.) 
                    New attribute of the Genome instance. 
                    Dictionary containing one dictionary per condition listed 
                    in loops.info. The key of this dictionary contains the 
                    condition name and the chosen window size (for example 
                    "cond_w100b" for a 100b window).
                    Each subdictionary has 2 keys: "loops" and "no_loops" and
                    contains the corresponding lists of genes.
                    For example, self.loops_genes[cond_w0b]["loops"] returns
                    the list of genes that overlap the loops positions.
            self.genes[locus].is_loop (dict.) 
                    New attribute of Gene instances related to the Genome
                    instance given as argument.Dictionary of shape
                    {condition:boolean}. The boolean is True if the gene
                    overlaps the loop positions (window included).

        Args:
            cond (Optional [list of str.]): selection of one or several 
                    conditions (1 condition corresponds to 1 data file and 
                    has to be in the loops.info file). By default: 
                    cond ='all' ie all available data are loaded.
            per_genes(Optional [Bool.]): if True, determines the list of 
                    genes overlapping the loops positions and returns the 
                    outputs self.loops_genes and self.genes[locus].is_loop 
                    described below.
            window (Optional [int.]): window around the loop positions (ie
                    around the 2 bins) for the seeking of overlapping genes
                    (Default: 0). All genes overlaping any position between:
                        * loops_start_bin1 - window and loops_end_bin1 + window or 
                        * loops_start_bin2 - window and loops_end_bin2 + window
                    are considered "loops" genes.

        Note:
            The data importation requires a loops.info file that contains 
            the column indices of each information in the data file and 
            some additional information, in the following order:
                * (required) [0] Condition, [1] Filename, [2] Startline, 
                  [3] Separator, [4] Bin1, [5] Bin2, [6] Binsize (in b),
                * (optional) [7] Score, [8] Pvalue,[9] Qvalue
            loops.info and data file have to be in the /HiC/Loops/ directory
        
        Note:    
            The position (pos) is the bin number. If data are binned at 2kb,
            the bin with number 10 corresponds to data between 20000 and 22000.

        Warning:
            This method needs a genomic sequence and a genomic annotation.
            If no annotation is loaded, the load_annotation method with the 
            default "sequence.gff3" file is computed. If no sequence is loaded, 
            the load_seq method with de defualt "sequence.fasta" is computed.
            To use another annotation or sequence, please load them to your
            Transcriptome instance with the following commands before using this 
            method:
                >>> from GRATIOSA import Genome, HiC
                >>> HC = HiC.HiC("ecoli")
                >>> g = Genome.Genome(HC.name)
                >>> g.load_annotation(annot_file=chosen_file)
                >>> g.load_seq(filename=chosen_file2
                >>> HC.genes = g.genes
                >>> HC.length = g.length

        Warning:
            Make sure that the version of the genome is the same as
            the version used in the HiC analysis workflow

        Example:
            >>> g = Genome.Genome("dickeya")
            >>> g.load_loops_genes()
            >>> g.loops_pos['WT2kb_bin2000b']['loops']
            [254000,254001,254002,254003,...]
            >>> g.loops_genes['WT2kb_bin2000b_w0b']['loops']
            ['Dda3937_00221','Dda3937_04438','Dda3937_03673',...]
            >>> g.genes['Dda3937_00221'].is_loop
            {'acid2kb_bin2000b': False, 'WT2kb_bin2000b': True}
        """
        if not hasattr(self, "seq"):
            gen = Genome.Genome(self.name)
            gen.load_seq()
            self.length = gen.length

        if not hasattr(self, "genes"):
            gen = Genome.Genome(self.name)
            gen.load_annotation()
            self.genes = gen.genes

        if not hasattr(self, "loops_genes"):
            self.loops_genes = {}

        if not hasattr(self, "loops_pos"):
            self.loops_pos = {}

        if isinstance(cond, str):
            cond = [cond]

        if not hasattr(self, "loops_pos") or cond == ["all"]:
            self.load_hic_loops(cond)
        else:
            unloaded_cond = set(cond) - set(self.loops_pos.keys())
            if unloaded_cond:
                self.load_hic_loops(list(unloaded_cond))

        if cond == ['all']:
            cond = self.loops_pos.keys()
        for c in cond:
            print(f"Loading loops genes in {c}_w{window}b")
            loops_pos = self.loops_pos[c]["loops"]
            is_loop = [False] * self.length
            for pos in loops_pos:
                is_loop[pos] = True

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

    def load_borders_genes(self, cond="all", window=0):
        """
        load_hic_borders_genes imports a list of borders from a data file 
        (typically a borders.tsv data file obtained with Chromosight) 
        containing at least the border position (bin number) and determines 
        the list of genes overlapping one of the borders.

        Creates:
            * self.borders_genes (dict. of dict.) 
                    created only if per_genes = True. New attribute of the 
                    Genome instance. Dictionary containing one dictionary 
                    per condition listed in borders.info. The key of this 
                    dictionary contains the condition name and the chosen
                    window size (for example "cond_w100b" for a 100b window).
                    Each subdictionary has 2 keys: "borders" and "no_borders" 
                    and contains the corresponding lists of genes. For 
                    example, self.borders_genes[cond_w0b]["borders"] returns
                    the list of genes that overlap the borders positions.
            * self.borders_pos (dict. of dict) 
                    new attribute of the HiC instance, dictionary containing 
                    one subdictionary per condition, with the shape 
                    self.borders_pos[cond] =
                    {'borders':list of positions that are in a border,
                    'no_borders':list of positions that are not in a border}
            * self.genes[locus].is_border (dict.) 
                    created only if per_genes = True.
                    New attribute of Gene instances related to the Genome
                    instance given as argument.Dictionary of shape
                    {condition:boolean}. The boolean is True if the gene
                    overlaps the border positions (window included).


        Args:
            cond (Optional [list of str.]): selection of one or several 
                    conditions (1 condition corresponds to 1 data file and 
                    has to be in the borders.info file). By default: 
                    cond ='all' ie all available data are loaded.
            per_genes(Optional [Bool.]): if True, determines the list of 
                    genes overlapping the borders positions and returns the 
                    outputs self.borders_genes and self.genes[locus].is_border 
                    described below.
            window (Optional [int.]): window around the border positions ie 
                    around the bin) for the seeking of overlapping genes
                    (Default: 0). All genes overlaping any position between 
                    border_start_bin - window and border_end_bin + window 
                    are considered "borders" genes.

        Note:
            The data importation requires a borders.info file that contains 
            the column indices of each information in the data file and some
            additional information, in the following order:
                * (required) [0] Condition, [1] Filename, [2] Startline,
                  [3] Separator, [4] Bin, [5] Binsize (in b)
                * (optional) [6] Score, [7] Pvalue,[8] Qvalue
            borders.info and data file have to be in the /HiC/Borders/ directory

        Note:
            The position (pos) is the bin number. If data are binned at 2kb,
            the bin with number 10 corresponds to data between 20000 and 22000.

        Warning:
            This method needs a genomic sequence and a genomic annotation.
            If no annotation is loaded, the load_annotation method with the 
            default "sequence.gff3" file is computed. If no sequence is loaded, 
            the load_seq method with de defualt "sequence.fasta" is computed.
            To use another annotation or sequence, please load them to your 
            Transcriptome instance with the following commands before using this 
            method:
                >>> from GRATIOSA import Genome, HiC
                >>> HC = HiC.HiC("ecoli")
                >>> g = Genome.Genome(HC.name)
                >>> g.load_annotation(annot_file=chosen_file)
                >>> g.load_seq(filename=chosen_file2
                >>> HC.genes = g.genes
                >>> HC.length = g.length


        Warning:
            make sure that the version of the genome is the same as
            the version used in the HiC analysis workflow

        Example:
            >>> g = Genome.Genome("dickeya")
            >>> g.load_borders_genes()
            >>> g.borders_pos["WT2_2kb_bin2000b"]['borders']
            [20000,20001,20002,20003,20004,...]
            >>> g.borders_genes['acid1kb_bin1000b_w0b']['borders']
            ['Dda3937_00158','Dda3937_01107','Dda3937_01108',...]
            >>> g.genes['Dda3937_00221'].is_border
            {'acid1kb_bin1000b': False,'WT2_2kb_bin2000b': False}
        """

        if not hasattr(self, "seq"):
            gen = Genome.Genome(self.name)
            gen.load_seq()
            self.length = gen.length

        if not hasattr(self, "genes"):
            gen = Genome.Genome(self.name)
            gen.load_annotation()
            self.genes = gen.genes

        if not hasattr(self, "borders_genes"):
            self.borders_genes = {}

        if isinstance(cond, str):
            cond = [cond]

        if not hasattr(self, "borders_pos") or cond == ["all"]:
            self.load_hic_borders(cond)
        else:
            unloaded_cond = set(cond) - set(self.borders_pos.keys())
            if unloaded_cond:
                self.load_hic_borders(list(unloaded_cond))

        if cond == ["all"]:
            cond = self.borders_pos.keys()

        for c in cond:
            print(f"Loading borders genes in {c}_w{window}b")
            borders_pos = self.borders_pos[c]["borders"]
            is_border = [False] * self.length
            for pos in borders_pos:
                is_border[pos] = True

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
