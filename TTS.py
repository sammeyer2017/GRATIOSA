#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd


#==============================================================================#

class TTS:

    def __init__(self, *args, **kwargs):
        """ Possible kwargs arguments: start, stop, strand
        """
        self.left= kwargs.get('left',None)
        self.right= kwargs.get('right',None)
        self.strand = kwargs.get('strand',None)
        self.rho_dpdt = kwargs.get('rho_dpdt',None)

        self.genes= kwargs.get('genes',[])
        self.seq= kwargs.get('seq',"")
        self.score= kwargs.get('score',None)

        if self.strand:
            self.start = self.left
            self.end = self.right
        elif not self.strand:
            self.start = self.right
            self.end = self.left


    def add_genes(self, genes):
        """ Adds a list of gene IDs to the TTS.
        """
        self.genes=genes
