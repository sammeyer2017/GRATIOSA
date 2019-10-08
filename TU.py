#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd


#==============================================================================#

class TU:

    def __init__(self, *args, **kwargs):
        """ Possible kwargs arguments: start, stop, strand
        """
        self.start=kwargs.get('start')
        self.stop=kwargs.get('stop')
        orient = kwargs.get('orientation')
        if type(orient) == str:
            if orient.lower() in ["true","1","+"]:         
                self.orientation = True
            elif orient.lower() in ["false","0","-"]:         
                self.orientation = False
            else:
                self.orientation = None
        else:
            self.orientation = orient

        self.genes = kwargs.get('genes')
        if self.orientation:
            self.left = self.start
            self.right = self.stop
        else:
            self.left = self.stop
            self.right = self.start

    def add_genes(self, genes):
        """ Adds a list of gene IDs to the TU.
        """
        self.genes=genes

    def add_TSS_cond(self, condition, TSS):
        """ Adds a list of TSSs to the TU corresponding to given condition.
        """
        if not hasattr(self, 'TSSs_cond'):
            self.TSS_cond = {}
        self.TSS_cond[condition]=TSS

    def add_TTS_cond(self, condition, TTS):
        """ Adds a list of TTSs to the TU corresponding to a given condition.
        """
        if not hasattr(self, 'TTSs'):
            self.TTS_cond = {}
        self.TTS_cond[condition]=TTS

    def add_correlation(self, correlations):
        """ Adds a list of expression correlation values among genes of TU
            Shape [(gene1,gene2,correlation among conditions),...]
        """
        self.correlation = correlations
        self.mean_correlation = np.mean([x[2] for x in correlations])

    def add_intergenic_cov(self, cov):
        """ Adds a list of coverage values from RNAseq data between intergenic regions among successive genes of TU
            Shape [(gene1,gene2,mean coverage among conditions),...]
        """
        self.intergenic_cov = cov
        self.mean_intergenic_cov = np.mean([np.mean(x[2]) for x in cov])

    def add_expression_ratio(self, expr_ratio):
        """ Adds a list of expression ratio values from RNAseq data (log2RPKM) among genes of TU
            Shape [(gene1,gene2,mean expression ratio among conditions),...]
        """
        self.expression_ratio = expr_ratio
        self.mean_expression_ratio = np.mean([x[2] for x in expr_ratio])

    def add_idx_corr_ratio(self, idx):
        """ 
        """
        self.idx_corr_ratio = idx
        self.mean_idx_corr_ratio = np.mean([x[2] for x in idx])

    def add_TSS(self, TSS):
        """
        Attributes potential TSS i.e. those located 100 bp upstream TU first start 
        """
        self.TSS = TSS

    def add_TTS(self, TTS):
        """
        Attributes potential TTS i.e. those located 100 bp downstream TU last stop 
        """
        self.TTS = TTS
