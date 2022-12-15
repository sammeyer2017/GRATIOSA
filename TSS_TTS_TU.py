#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd

#==============================================================================#
class TSS:
    def __init__(self, *args, **kwargs):
        self.pos = kwargs.get('pos')
        self.strand = kwargs.get('strand',None)
        self.genes=kwargs.get('genes',[])
        self.promoter={}
        self.score= kwargs.get('score',None)
        
    def add_prom_elements(self,gen_seq,gen_seqcompl,*arg,**kwargs):
        '''
        Compute spacer, -10 and - 35 sequences starting from boxes (sites) coordinates and genome sequence
        gen_seq : + strand, gen_seqcompl : - strand (complementary, /!\ 3' -> 5'). shift : nb of nucleotides to
        include on either side.
        '''
        shift = kwargs.get('shift',0) # number of nt to include beyond each region on either side, e.g. to compute angle
        prom_region = kwargs.get('prom',[0,0]) # specify the length of the region upstream TSSs that need to be extracted if required
        # shape [length before TSS, length after TSS] e.g. [250,100]
        if prom_region != [0,0]:
            if self.strand == True:
                self.promoter['region'] = gen_seq[self.pos-1-prom_region[0]:self.pos+prom_region[1]+shift]

            elif self.strand == False:
                self.promoter['region'] = gen_seqcompl[self.pos-prom_region[1]-1-shift:self.pos+prom_region[0]+shift][::-1]

        try:

            for sig in self.promoter.keys(): # for each sigma factor
                try:    
                    sites_pos = self.promoter[sig]['sites']                
                    if self.strand :
        # if promoter on + strand, spacer length = -10L - -35R -1
                        self.promoter[sig]['spacer'] = gen_seq[sites_pos[3]-shift:sites_pos[0]-1+shift]
                        self.promoter[sig]['minus10'] = gen_seq[sites_pos[0]-1-shift:sites_pos[1]+shift]
                        self.promoter[sig]['minus35'] = gen_seq[sites_pos[2]-1-shift:sites_pos[3]+shift]
                        self.promoter[sig]['discriminator'] = gen_seq[sites_pos[1]-shift:self.pos-1+shift]
                    elif not self.strand:
        # if promoter on - strand, spacer length = -35L - -10R -1
                        self.promoter[sig]['spacer'] = gen_seqcompl[sites_pos[1]-shift:sites_pos[2]-1+shift][::-1]
                        self.promoter[sig]['minus10'] = gen_seqcompl[sites_pos[0]-1-shift:sites_pos[1]+shift][::-1]
                        self.promoter[sig]['minus35'] = gen_seqcompl[sites_pos[2]-1-shift:sites_pos[3]+shift][::-1]
                        self.promoter[sig]['discriminator'] = gen_seqcompl[self.pos-shift:sites_pos[0]-1+shift][::-1]
                except Exception as e: # sigma factor without site coordinates or invalid
                    pass

        except Exception as e:
            print('Error for TSS',self.pos,':',e)

""" A ENLEVER POUR PUBLI ?
    def add_strand(self,strand):
        self.strand = strand

    def add_score(self,score):
        self.score = score

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
"""

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

        self.start = kwargs.get('start',None)
        self.end = kwargs.get('end',None)
       
        if self.strand:
            self.start = self.left
            self.end = self.right
        elif not self.strand:
            self.start = self.right
            self.end = self.left

"""A ENLEVER POUR PUBLI ? 
    def add_genes(self, genes):
        #Adds a list of gene IDs to the TTS.
        self.genes=genes
"""

class TU:

    def __init__(self, *args, **kwargs):
        """ Possible kwargs arguments: start, stop, strand
        """
        self.start=kwargs.get('start')
        self.stop=kwargs.get('stop')
        strand = kwargs.get('strand')
        if type(strand) == str:
            if strand.lower() in ["true","1","+"]:         
                self.strand = True
            elif strand.lower() in ["false","0","-"]:         
                self.strand = False
            else:
                self.strand = None
        else:
            print("Unknown 'strand' type")

        self.genes = kwargs.get('genes')
        if self.strand:
            self.left = self.start
            self.right = self.stop
        else:
            self.left = self.stop
            self.right = self.start

'''A ENLEVER POUR PUBLI ??
    def add_genes(self, genes):
        """ #Adds a list of gene IDs to the TU.
        """
        self.genes=genes

    def add_correlation(self, correlations):
        """ #Adds a list of expression correlation values among genes of TU
            #Shape [(gene1,gene2,correlation among conditions),...]
        """
        self.correlation = correlations
        self.mean_correlation = np.mean([x[2] for x in correlations])

    def add_intergenic_cov(self, cov):
        """ #Adds a list of coverage values from RNAseq data between intergenic regions among successive genes of TU
            #Shape [(gene1,gene2,mean coverage among conditions),...]
        """
        self.intergenic_cov = cov
        self.mean_intergenic_cov = np.mean([np.mean(x[2]) for x in cov])

    def add_expression_ratio(self, expr_ratio):
        """ #Adds a list of expression ratio values from RNAseq data (log2RPKM) among genes of TU
            #Shape [(gene1,gene2,mean expression ratio among conditions),...]
        """
        self.expression_ratio = expr_ratio
        self.mean_expression_ratio = np.mean([x[2] for x in expr_ratio])

    def add_TSS(self, x, TSS):
        """
        #Attributes potential TSS from list
        #x = idx of gene in TU, 1 = first gene of TU = start of TU
        #TSS = tuple (position, proportion of total starts)
        """
        if not hasattr(self, 'TSS'):
            self.TSS = {}

        if x not in self.TSS.keys():
            self.TSS[x] = []

        self.TSS[x] += TSS

    def add_TTS(self, x, TTS):
        """
        #Attributes potential TTS from list
        #x = idx of gene in TU, 1 = first gene of TU = start of TU
        #TTS = tuple (position, type of terminator)
        """
        if not hasattr(self, 'TTS'):
            self.TTS = {}
        
        if x not in self.TTS.keys():
            self.TTS[x] = []

        self.TTS[x] += TTS

    def add_TSS_treated(self, x, TSS):
        """
        #Attributes treated TSS from list
        #x = idx of gene in TU, 1 = first gene of TU = start of TU
        """
        if not hasattr(self, 'TSS_treated'):
            self.TSS_treated = {}

        self.TSS_treated[x] = TSS

    def add_TTS_treated(self, x, TTS):
        """
        #Attributes treated TTS from list
        #x = idx of gene in TU, 1 = first gene of TU = start of TU
        """
        if not hasattr(self, 'TTS_treated'):
            self.TTS_treated = {}
        
        self.TTS_treated[x] = TTS


    def add_TSS_strong(self, x, TSS):
        """
        #Attributes strong enough TSS from list
        #x = idx of gene in TU, 1 = first gene of TU = start of TU
        """
        if not hasattr(self, 'TSS_strong'):
            self.TSS_strong = {}

        self.TSS_strong[x] = TSS

    def add_TTS_strong(self, x, TTS):
        """
        #Attributes strong enough TTS from list
        #x = idx of gene in TU, 1 = first gene of TU = start of TU
        """
        if not hasattr(self, 'TTS_strong'):
            self.TTS_strong = {}
        
        self.TTS_strong[x] = TTS
'''