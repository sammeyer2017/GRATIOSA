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
        prom_region = kwargs.get('prom',False)
        if prom_region != False:
            if self.strand == True:
                self.promoter['region'] = gen_seq[self.pos-1-prom_region[0]:self.pos+prom_region[1]+shift]

            elif self.strand == False:
                self.promoter['region'] = gen_seqcompl[self.pos-prom_region[1]-1-shift:self.pos+prom_region[0]+shift][::-1]

        if self.strand == True:
            self.promoter['discr_model'] = gen_seq[self.pos-12-1:self.pos+8]

        elif self.strand == False:
            self.promoter['discr_model'] = gen_seqcompl[self.pos-8-1:self.pos+12][::-1]


        try:
            for sig in self.promoter.keys(): # for each sigma factor
                try:                    
                    if self.strand == True:
        # if promoter on + strand, spacer length = -10L - -35R -1
                        self.promoter[sig]['spacer'] = gen_seq[self.promoter[sig]['sites'][3]-shift:self.promoter[sig]['sites'][0]-1+shift]
                        self.promoter[sig]['minus10'] = gen_seq[self.promoter[sig]['sites'][0]-1-shift:self.promoter[sig]['sites'][1]+shift]
                        self.promoter[sig]['minus35'] = gen_seq[self.promoter[sig]['sites'][2]-1-shift:self.promoter[sig]['sites'][3]+shift]
                        self.promoter[sig]['discriminator'] = gen_seq[self.promoter[sig]['sites'][1]-shift:self.pos-1+shift]

                        lenbef = self.pos - self.promoter[sig]['sites'][0] + 1 ; rest = 21 - lenbef
                        self.promoter[sig]['discr_model'] = gen_seq[self.promoter[sig]['sites'][0]-1:self.promoter[sig]['sites'][0]+20]

                    elif self.strand == False:
        # if promoter on - strand, spacer length = -35L - -10R -1
                        self.promoter[sig]['spacer'] = gen_seqcompl[self.promoter[sig]['sites'][1]-shift:self.promoter[sig]['sites'][2]-1+shift][::-1]
                        self.promoter[sig]['minus10'] = gen_seqcompl[self.promoter[sig]['sites'][0]-1-shift:self.promoter[sig]['sites'][1]+shift][::-1]
                        self.promoter[sig]['minus35'] = gen_seqcompl[self.promoter[sig]['sites'][2]-1-shift:self.promoter[sig]['sites'][3]+shift][::-1]
                        self.promoter[sig]['discriminator'] = gen_seqcompl[self.pos-shift:self.promoter[sig]['sites'][0]-1+shift][::-1]

                        lenbef = self.promoter[sig]['sites'][1] - self.pos + 1 ; rest = 21 - lenbef
                        #self.promoter[sig]['discr_model'] = gen_seqcompl[self.pos-rest-1:self.promoter[sig]['sites'][1]][::-1]
                        self.promoter[sig]['discr_model'] = gen_seqcompl[self.promoter[sig]['sites'][1]-20-1:self.promoter[sig]['sites'][1]][::-1]

                except Exception as e: # sigma factor without site coordinates or invalid
                    pass

        except Exception as e:
            print 'Error for TSS',self.pos,':',e
