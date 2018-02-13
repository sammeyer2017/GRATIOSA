#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np

#==============================================================================#

class Operon:
    def __init__(self, left, right, strand):
        self.left = left
        self.right = right
        self.strand = strand
        if strand == '-':
            self.start=self.right
            self.end=self.left
        else:
            self.start=self.left
            self.end=self.right
        self.genes=[]

    def add_genes(self,gene):
        self.genes.append(gene)

    def add_ID(self,ID):
        self.ID=ID
