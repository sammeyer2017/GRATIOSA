#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
from useful_functions_HiC import *
from pathlib import Path
from globvar import *

#==============================================================================#

class HiC:

    def __init__(self, name):
        """ Called when an HiC instance is created, 
        initializes the attribute name. 
        """
        self.name = name

    def load_hic_borders(self):  
        """ 
        load_hic_borders imports a list of borders from a data file (typically 
        a borders.tsv data file obtained with Chromosight) containing :
        + (required) border coordinate (bin number),
        + (optional) border score (pearson correlation coefficient between the 
        border kernel and the detected pattern), border pvalue and border qvalue
        See Chromosight documentation for more details.

        Requirements:
        + The data importation requires a borders.info file that contains the 
        column indices of each information in the data file and some additional 
        information, in the following order :
        (required) [0] Condition, [1] Filename, [2] Startline, [3] Separator, 
                   [4] Bin, [5] Binsize (in b)
        (optional) [6] Score, [7] Pvalue,[8] Qvalue
        + borders.info and data file have to be in the /HiC/Borders/ directory 
        
        Output: 
        a new attribute borders, of shape self.borders = 
        {cond:{pos: {"score": score, "pval": pvalue, "qval": qvalue}}}

        N.B.: The position (pos) is the bin number. If data are binned at 2kb,
        the bin with number 10 corresponds to data between 20000 and 22000.

        !!!Caution: make sure that the version of the genome is the same as 
        the version used in the HiC analysis workflow!!!
        """

        # gets the path to data and .info files
        path2dir = f"{basedir}data/{self.name}/HiC/Borders/"

        # tries to open borders.info file
        if Path(f"{path2dir}borders.info").exists():
            with open(f"{path2dir}borders.info", "r") as f:
                skiphead = next(f)  
                # 1 line in borders.info corresponds to 1 condition (1 data file)
                for line in f:
                    line = line.strip('\n').split('\t')
                    # loads data file information for this condition
                    cond = line[0]
                    path2file = f"{basedir}data/{self.name}/HiC/Borders/{line[1]}"
                    startline, sep = int(line[2]), line[3]
                    bin_col = line[4]

                    # loads binsize
                    if len(line) > 7 :
                        binsize = int(line[5]) 
                    else : 
                        binsize = None

                    # loads column indices of each optional information             
                    score_col = line[6] if len(line) > 6 else None
                    pval_col = line[7] if len(line) > 7 else None
                    qval_col = line[8] if len(line) > 8 else None

                    # creates the new attribute using the load_HiC_site_cond function
                    # from useful_functions_HiC 
                    if not hasattr(self, "borders"):
                        self.borders = {}  # creates the new attribute 
                    self.borders[f"{cond}_bin{str(binsize)}b"] = load_HiC_site_cond(
                        'borders', path2file, startline, sep, bin_col,binsize, 
                        score_col=score_col, pval_col=pval_col, qval_col=qval_col)
        else:
            print("No borders.info, unable to load borders")


    def load_hic_loops(self):
        """ 
        load_hic_loops imports a list of loops from a data file (typically 
        a loops.tsv data file obtained with Chromosight) containing :
        + (required) loops coordinates (bin numbers),
        + (optional) loops score (pearson correlation coefficient between the 
        loop kernel and the detected pattern), loop pvalue and loop qvalue
        See Chromosight documentation for more details.

        Requirements:
        + The data importation requires a loops.info file that contains the 
        column indices of each information in the data file and some additional 
        information, in the following 
        order :
        (required) [0] Condition, [1] Filename, [2] Startline, [3] Separator, 
                   [4] Bin1, [5] Bin2, [6] Binsize (in b),
        (optional) [7] Score, [8] Pvalue,[9] Qvalue
        + loops.info and data file have to be in the /HiC/Loops/ directory 
        
        Output: 
        a new attribute borders, of shape self.borders = 
        {cond:{pos: {"score": score, "pval": pvalue, "qval": qvalue}}}

        N.B.: The position (pos) is the bin number. If data are binned at 2kb,
        the bin with number 10 corresponds to data between 20000 and 22000.

        !!!Caution: make sure that the version of the genome is the same as 
        the version used in the HiC analysis workflow!!!
        """
        # gets the path to data and .info files
        path2dir = f"{basedir}data/{self.name}/HiC/Loops/"

        # tries to open loops.info file
        if Path(f"{path2dir}loops.info").exists():
            with open(f"{path2dir}loops.info", "r") as f:
                skiphead = next(f)  # skip head

                # 1 line in loops.info corresponds to 1 condition (1 data file)
                for line in f:
                    line = line.strip('\n')
                    line = line.split('\t')

                    # loads data file information for this condition
                    cond = line[0]
                    path2file = f"{path2dir}{line[1]}"
                    startline, sep = int(line[2]), line[3]
                    bin1_col, bin2_col = line[4], line[5]

                    # loads binsize
                    if len(line) > 7 :
                        binsize = int(line[6]) 
                    else : 
                        binsize = None

                    # loads column indices of each optional information             
                    score_col = line[7] if len(line) > 7 else None
                    pval_col = line[8] if len(line) > 8 else None
                    qval_col = line[9] if len(line) > 9 else None

                    # creates the new attribute using the load_HiC_site_cond function
                    # from useful_functions_HiC 

                    if not hasattr(self, "loops"):
                        self.loops = {}  # creates the new attribute 
                    self.loops[f"{cond}_bin{str(binsize)}b"] = load_HiC_site_cond(
                        "loops",path2file, startline, sep, bin1_col,binsize,
                        bin2_col=bin2_col,score_col=score_col, 
                        pval_col=pval_col, qval_col=qval_col)
                    
        else:
            print("No loops.info, unable to load loops")