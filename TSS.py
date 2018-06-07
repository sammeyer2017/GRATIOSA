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
        # shape TSS.promoter = {[sig]:(sites)}
        # sites must have the shape : [-10l,-10r,-35l,-35r] with l = left, r = right
        sites = kwargs.get('sites')
        self.promoter[sig] = {}
        if sites:
            self.promoter[sig]['sites'] = tuple(map(int, sites.strip().replace(' ','').split(',')))

    def compute_magic_prom(self,gen_seq,gen_seqcompl,*arg,**kwargs):
        '''
        Compute spacer, -10 and - 35 sequences starting from boxes (sites) coordinates and genome sequence
        gen_seq : + strand, gen_seqcompl : - strand (complementary, /!\ 3' -> 5'). shift : nb of nucleotides to
        include on either side.
        '''
        shift = kwargs.get('shift',0) # number of nt to include beyond each region on either side, e.g. to compute angle
        try:
            for sig in self.promoter.keys(): # for each sigma factor
                try:                    
                    if self.strand == True:
        # if promoter on + strand, spacer length = -10L - -35R -1
                        self.promoter[sig]['spacer'] = gen_seq[self.promoter[sig]['sites'][3]-shift:self.promoter[sig]['sites'][0]-1+shift]
                        self.promoter[sig]['minus10'] = gen_seq[self.promoter[sig]['sites'][0]-1-shift:self.promoter[sig]['sites'][1]+shift]
                        self.promoter[sig]['minus35'] = gen_seq[self.promoter[sig]['sites'][2]-1-shift:self.promoter[sig]['sites'][3]+shift]
                        self.promoter[sig]['discriminator'] = gen_seq[self.promoter[sig]['sites'][1]-shift:self.pos-1+shift]

                    elif self.strand == False:
        # if promoter on - strand, spacer length = -35L - -10R -1
                        self.promoter[sig]['spacer'] = gen_seqcompl[self.promoter[sig]['sites'][1]-shift:self.promoter[sig]['sites'][2]-1+shift][::-1]
                        self.promoter[sig]['minus10'] = gen_seqcompl[self.promoter[sig]['sites'][0]-1-shift:self.promoter[sig]['sites'][1]+shift][::-1]
                        self.promoter[sig]['minus35'] = gen_seqcompl[self.promoter[sig]['sites'][2]-1-shift:self.promoter[sig]['sites'][3]+shift][::-1]
                        self.promoter[sig]['discriminator'] = gen_seqcompl[self.pos-shift:self.promoter[sig]['sites'][0]-1+shift][::-1]

                except: # sigma factor without site coordinates or invalid
                    pass

        except Exception as e:
            print 'Error for TSS',self.pos,':',e