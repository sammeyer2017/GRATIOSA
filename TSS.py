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
        self.pos = kwargs.get('pos')
        self.strand = None
        self.genes=[]
        self.promoter={}

    def add_strand(self,strand):
        self.strand = strand
        if self.strand == '+':
            self.strand= True
        elif self.strand == '-':
            self.strand=False
        else:
            self.strand = None

    def add_genes(self,tags,genes_dict):
        # tags must have the shape : [gene1,gene2]
        for tag in tags.strip().replace(' ','').split(','):
            if tag != '':
                if tag in genes_dict.keys():
                    self.genes.append(tag)
                else:
                    print(tag + " not in annotations")

    def add_promoter(self, sig, *arg, **kwargs):
        # shape TSS.promoter = {[sig]:[sites]}
        # sites must have the shape : [-10l,-10r,-35l,-35r] with l = left, r = right
        self.promoter[sig] = []
        sites = kwargs.get('sites')
        if sites:
            for site in sites.strip().replace(' ','').split(','):
                self.promoter[sig].append(site)

