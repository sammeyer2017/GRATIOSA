#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd

'''
The TTS (Transcription Termination Site), TSS (Transcription Start Site) and
TU (Transcription Unit) classes allow to gather information on the positions, 
reliability and sequences of these elements.
'''

class TSS:
    """ 
    Each TTS instance has to be initialized with the following attributes:
        * pos (int.): position of the TSS
        * strand (bool.): DNA strand (True for forward strand, False for 
          complementary strand)
        * genes (list of str.): list of locus tags of genes associated with 
          the TSS
        * score (float.): TSS score
    """

    def __init__(self, pos=None, strand=None, genes=[], score=None):
        self.pos = pos 
        self.strand = strand
        self.genes = genes  
        self.promoter = {}
        self.score = score


    def add_genes(self, tags, genes_dict):
        """
        add_genes completes the genes attribute of the TSS instance
        with the genes that are given as tags and that are in the
        genes_dict.

        Args:
            tags (str.): locus tags separated by commas
            genes_dict: dictionary of shape {locus tag : Gene object}

        Example:
            >>> from GRATIOSA import Genome
            >>> g = Genome.Genome("dickeya")
            >>> g.load_annotation()
            >>> from GRATIOSA.TSS_TTS_TU import TSS
            >>> tss = TSS()
            >>> tss.add_genes("Dda3937_00001,Dda3937_00002",g.genes)
        """
        for tag in tags.strip().replace(' ', '').split(','):
            if tag != '':
                if tag in genes_dict.keys():
                    if tag not in self.genes:
                        self.genes.append(tag)
                else:
                    print(tag + " not in annotations")


    def add_promoter(self, sig, *arg, **kwargs):
        """
        add_promoter completes the promoter attribute of the TSS instance.
        Creates or complete the 'promoter' of the TSS instance with the shape {[sig]:(sites)}

        Args:
            sig (str.): sigma factor
            sites (Optional [str.]): the left and right coordinates of the -10 and -35 
                    elements with the shape "-10l,-10r,-35r,-35l"
            genes_dict (Optional [dict.]): dictionary of shape {locus tag : Gene object}

        Example:
            >>> from GRATIOSA.TSS_TTS_TU import TSS
            >>> tss = TSS(pos=4707030,strand=False)
            >>> tss.add_promoter(sig = "sigma70", sites="4707037,4707046,4707062,4707068")
            >>> tss.promoter
            {'sigma70': {'sites': (4707037, 4707046, 4707062, 4707068)}}
        """

        # sites must have the shape : [-10l,-10r,-35l,-35r] with l = left, r =
        # right
        sites = kwargs.get('sites')
        self.promoter[sig] = {}
        if sites:
            self.promoter[sig]['sites'] = tuple(
                map(int, sites.strip().replace(' ', '').split(',')))


    def add_prom_elements(self, gen_seq, gen_seqcompl, shift=0, prom=[0, 0]):
        '''
        add_prom_elements extracts sequences of the different promoter 
        elements (spacer, -10, -35, discriminator, region around TSS) based 
        on -35 and -10 coordinates that are included in the promoter attribute.

        Completes the 'promoter' attribute of the TSS instance:
                Before using this method, this attribute was a dictionary of 
                shape {sigma factor: subdictionary} containing, for each 
                sigma factor, a subdictionary of shape "sites":(-10l,-10r,-35r,
                -35l)} with -10l, -10r, -35r and -35l the left and right 
                coordinates of the -10 and -35 elements. Now, this 
                subdictionary contains 4 new keys : 
                "spacer", "minus10","minus35" and "discriminator". 
                Each associated value is the sequence of the element.

        Args:
            gen_seq (str.): sequence of the forward strand
            gen_seqcompl (str.): sequence of the complementary strand
                    (WARNING: 3' -> 5')
            shift (Optional [int.]): number of nucleotides to include beyond 
                    each region on either side (Default: 0nt)
            prom (Optional [int.,int.]): region upstream and downstream TSS
                    to extract. Argument with the shape:
                    [length before TSS, length after TSS]
                    Default: [0,0] ie no sequence will be extracted around 
                    TSS.
                    
        Note: 
            If prom !=[0,0], self.promoter dictionary has also a new 
            "region" key. self.promoter['region'] returns the sequence 
            of the region chosen with the prom argument.
                

        Warning: 
            the promoter attribute has to be loaded prior using 
            this method. Use add_promoter method.
            

        Example:
            >>> from GRATIOSA import Genome
            >>> g = Genome.Genome("dickeya")
            >>> g.load_seq()
            >>> from GRATIOSA.TSS_TTS_TU import TSS
            >>> tss = TSS(pos=4707030,strand=False)
            >>> tss.add_promoter(sig = "sigma70", sites="4707037,4707046,4707062,4707068")
            >>> tss.add_prom_elements(g.seq,
            ...                       g.seqcompl,
            ...                       sig = "sigma70",
            ...                       sites="4707039,4707044,4707060,4707065",
            ...                       prom=[20,20])
            >>> tss.promoter
            {'sigma70': {'sites': (4707037, 4707046, 4707062, 4707068),
             'spacer': 'CCTCGCCCACCCTCA',
             'minus10': 'ATCATCATGA',
             'minus35': 'CCGTACC',
             'discriminator': 'ATAACC'},
             'region': 'CTCAATCATCATGAATAACCCCCCTCCTTGTGTCTTTCTTA'}
        '''
        if prom != [0, 0]:
            if self.strand:
                self.promoter['region'] = gen_seq[self.pos -
                                                  1 - prom[0]:self.pos + prom[1] + shift]

            elif self.strand == False:
                self.promoter['region'] = gen_seqcompl[self.pos -
                                                       prom[1] - 1 - shift:self.pos + prom[0] + shift][::-1]

        try:

            for sig in self.promoter.keys():  # for each sigma factor
                try:
                    sites_pos = self.promoter[sig]['sites']
                    if self.strand:
                        self.promoter[sig]['spacer'] = gen_seq[sites_pos[3] - shift:sites_pos[0] - 1 + shift]
                        self.promoter[sig]['minus10'] = gen_seq[sites_pos[0] - 1 - shift:sites_pos[1] + shift]
                        self.promoter[sig]['minus35'] = gen_seq[sites_pos[2] - 1 - shift:sites_pos[3] + shift]
                        self.promoter[sig]['discriminator'] = gen_seq[sites_pos[1] - shift:self.pos - 1 + shift]
                    elif not self.strand:
                        self.promoter[sig]['spacer'] = gen_seqcompl[sites_pos[1] - shift:sites_pos[2] - 1 + shift][::-1]
                        self.promoter[sig]['minus10'] = gen_seqcompl[sites_pos[0] - 1 - shift:sites_pos[1] + shift][::-1]
                        self.promoter[sig]['minus35'] = gen_seqcompl[sites_pos[2] - 1 - shift:sites_pos[3] + shift][::-1]
                        self.promoter[sig]['discriminator'] = gen_seqcompl[self.pos - shift:sites_pos[0] - 1 + shift][::-1]
                except Exception as e:  # sigma factor without site coordinates or invalid
                    pass

        except Exception as e:
            print('Error for TSS', self.pos, ':', e)


class TTS:
    """
    Each TTS instance has to be initialized with the following attributes:
        * left (int.) and right (int.): TTS coordinates (does not take into
          account the strand, ie right > left)
        * start (int.) and end (int.): positions of the beginning and the 
          end of the TTS
        * strand (bool.): DNA strand (True for forward strand, False for 
          complementary strand)
        * rho_dpdt (bool.): rho dependency of the TTS
        * genes (list of str.): list of locus tags of genes associated with 
          the TTS
        * score (float.): TSS score
    """

    def __init__(self, start=None, end=None, left=None, right=None,
                 strand=None, rho_dpdt=None, genes=[], seq="", score=None):
        self.start = start
        self.end = end
        self.left = left
        self.right = right
        self.strand = strand
        self.rho_dpdt = rho_dpdt
        self.genes = genes
        self.seq = seq
        self.score = score

        if strand is not None:
            if self.left is not None and self.right is not None and (
                    self.start is None or self.stop is None):
                if self.strand:
                    self.start = self.left
                    self.end = self.right
                elif not self.strand:
                    self.start = self.right
                    self.end = self.left
            elif self.start is not None and self.end is not None and (self.left is None or self.right is None):
                if self.strand:
                    self.left = self.start
                    self.right = self.end
                elif not self.strand:
                    self.right = self.start
                    self.left = self.end


    def add_genes(self, tags, genes_dict):
        """
        add_genes completes the genes attribute of the TTS instance
        with the genes that are given as tags and that are in the genes_dict.

        Args:
            tags (str.): locus tags separated by commas
            genes_dict: dictionary of shape {locus tag : Gene object}

        Example:
            >>> from GRATIOSA import Genome
            >>> g = Genome.Genome("dickeya")
            >>> g.load_annotation()
            >>> from GRATIOSA.TSS_TTS_TU import TTS
            >>> tts = TTS()
            >>> tts.add_genes("Dda3937_00001,Dda3937_00002",g.genes)
            >>> tts.genes
            ["Dda3937_00001,Dda3937_00002"]
        """
        for tag in tags.strip().replace(' ', '').split(','):
            if tag != '':
                if tag in genes_dict.keys():
                    if tag not in self.genes:
                        self.genes.append(tag)
                else:
                    print(tag + " not in annotations")


class TU:
    """
    Each TU instance has to be initialized with the following attributes:
        * left (int.) and right (int.): TU coordinates (does not take into
          account the strand, ie right > left)
        * start (int.) and end (int.): positions of the beginning and the end of the TU
        * strand (bool.): DNA strand (True for forward strand, False for 
          complementary strand)
        * genes (list of str.): list of locus tags of genes associated with the TU

    """
    def __init__(self, start=None, end=None, left=None,
                 right=None, strand=None, genes=[],expression=None):
        self.start = start
        self.end = end
        self.left = left
        self.right = right
        self.strand = strand

        if strand is not None:
            if self.left is not None and self.right is not None and (
                    self.start is None or self.stop is None):
                if self.strand:
                    self.start = self.left
                    self.end = self.right
                elif not self.strand:
                    self.start = self.right
                    self.end = self.left
            elif self.start is not None and self.end is not None and (self.left is None or self.right is None):
                if self.strand:
                    self.left = self.start
                    self.right = self.end
                elif not self.strand:
                    self.right = self.start
                    self.left = self.end
        self.genes = genes
        self.expression = expression


    def add_genes(self, tags, genes_dict):
        """
        add_genes completes the genes attribute of the TU instance
        with the genes that are given as tags and that are in the genes_dict.

        Args:
            tags (str.): locus tags separated by commas
            genes_dict: dictionary of shape {locus tag : Gene object}

        Example:
            >>> from GRATIOSA import Genome
            >>> g = Genome.Genome("dickeya")
            >>> g.load_annotation()
            >>> from GRATIOSA.TSS_TTS_TU import TU
            >>> tu = TU()
            >>> tu.add_genes("Dda3937_00001,Dda3937_00002",g.genes)
            >>> tu.genes
            ["Dda3937_00001,Dda3937_00002"]
        """
        for tag in tags.strip().replace(' ', '').split(','):
            if tag != '':
                if tag in genes_dict.keys():
                    if tag not in self.genes:
                        self.genes.append(tag)
                else:
                    print(tag + " not in annotations")


    def add_TU_expression(self, condition, TU_expr):
        """
        add_TU_expression adds an expression value associated to a new 
        condition to the expression attribute of the TU instance. If a list 
        of expressions is given in input (with the TU_expr argument), the 
        mean is computed and the mean value is added to the attribute. Elif 
        only one value is given in input, this value is directly added to the 
        attribute.

        Creates or completes the "expression" attribute of the TU instance. 
        This attribute is a dictionary of shape {condition:expression value}

        Args:
            condition (str.): name of the condition
            TU_expr (float or list of float): expression value(s)
                             
        Example:
            >>> from GRATIOSA.TSS_TTS_TU import TU
            >>> tu = TU()
            >>> tu.add_TU_expression("Control",0.24)
            >>> tu.add_TU_expression("test",[0.2,0.4])
            >>> tu.expression
            {'Control': 0.24, 'test': 0.30}
        """
        if not hasattr(self, 'expression'):
            self.expression = {}

        if isinstance(TU_expr, float):
            self.expression[condition] = TU_expr
        else:
            self.expression[condition] = np.mean(TU_expr)


    def add_TSS(self, TSS):
        """
        add_TSS attributes potential TSS to the TU object. 
        Creates or completes the 'TSS' attributes (list of tuples) of the 
        TU instance. 

        Args:
            TSS: tuple of shape (TSS position (int.),
                                 propotion of total starts (float))

        Example:
            >>> from GRATIOSA.TSS_TTS_TU import TU
            >>> tu = TU()
            >>> tu.add_TU_expression("Control",0.24)
            >>> tu.add_TU_expression("test",[0.2,0.4])
            >>> tu.expression
            {'Control': 0.24, 'test': 0.30}
        """
        if not hasattr(self, 'TSS'):
            self.TSS = []

        if TSS not in self.TSS:
            self.TSS.append(TSS)


    def add_TTS(self, TTS):
        """
        add_TTS attributes potential TTS to the TU object.
        Creates or completes the 'TTS' attributes (list of tuples) of the 
        TU instance. 

        Args:
            TTS: tuple of shape 
                    (TTS position (int.), propotion of total starts (float))
        """
        if not hasattr(self, 'TTS'):
            self.TTS = []

        if TTS not in self.TTS:
            self.TTS.append(TTS)