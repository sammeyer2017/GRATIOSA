#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np



#==============================================================================#

class TU:

    def __init__(self, *args, **kwargs):
        """ Possible kwargs arguments: start, stop, strand
        """
        self.start=kwargs.get('start')
        self.stop=kwargs.get('stop')
        self.orientation = kwargs.get('orientation')
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
