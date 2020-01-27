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



    def add_TSS_primary(self, TSS):
        """
        Attributes potential primary TSS = TSS of TU
        """
        self.TSS_prim = TSS

    def add_TTS_primary(self, TTS):
        """
        Attributes potential primary TTS = TTS of TU
        """
        self.TTS_prim = TTS
   
    def add_TSS_internal(self, x, TSS):
        """
        Attributes potential internal TSS
        """
        if not hasattr(self, 'TSS_intern'):
            self.TSS_intern = {}
        self.TSS_intern[x] = TSS

    def add_TTS_internal(self, x, TTS):
        """
        Attributes potential internal TTS
        """
        if not hasattr(self, 'TTS_intern'):
            self.TTS_intern = {}
        self.TTS_intern[x] = TTS



    def add_TSS_cov(self, x, TSS):
        """
        Attributes potential TSS detected by cov
        x = idx of gene in TU, 1 = first gene of TU = start of TU
        """
        if not hasattr(self, 'TSS_cov'):
            self.TSS_cov = {}
        
        self.TSS_cov[x] = TSS

    def add_TTS_cov(self, x, TTS):
        """
        Attributes potential TTS detected by cov
        x = idx of gene in TU, 1 = first gene of TU = start of TU
        """
        if not hasattr(self, 'TTS_cov'):
            self.TTS_cov = {}
        
        self.TTS_cov[x] = TTS

    def add_TSS(self, x, TSS):
        """
        Attributes potential TSS from list
        x = idx of gene in TU, 1 = first gene of TU = start of TU
        TSS = tuple (position, proportion of total starts)
        """
        if not hasattr(self, 'TSS'):
            self.TSS = {}
        
        self.TSS[x] = TSS

    def add_TTS(self, x, TTS):
        """
        Attributes potential TTS from list
        x = idx of gene in TU, 1 = first gene of TU = start of TU
        TTS = tuple (position, type of terminator)
        """
        if not hasattr(self, 'TTS'):
            self.TTS = {}
        
        self.TTS[x] = TTS
