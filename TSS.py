#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
# -------------------
# useful function


#==============================================================================#

class TSS:
    def __init__(self, *args, **kwargs):
        self.id = kwargs.get('id')
        self.pos = kwargs.get('pos')
        self.condition=[]
        self.genes={}

    def add_condition(self,condition):
        if condition not in self.condition:
            self.condition.append(condition)
            self.genes[condition]=[]

    def add_sig(self,sig):
        self.sig = sig

    def add_strand(self,strand):
        self.orientation = strand
        self.strand = strand
        if self.strand == '+':
            self.strand= True
        else:
            self.strand=False

    def add_genes(self,genes,condition):
        self.genes[condition].append(genes)
