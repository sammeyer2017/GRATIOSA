#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd

'''
The TTS (Transcription Termination Site), TSS (Transcription Start Site) and
TU (Transcription Unit) classes allow to gather information on the positions, 
reliability
and sequences of these elements.
'''


class TSS:
    def __init__(self, pos=None, strand=None, genes=[], score=None):
        """
        Called when a TSS instance is created,initializes the following 
        attributes: pos, strand, genes, promoter and score.
        Args:
            self: TSS instance
            pos (int.): position of the TSS
            strand (bool.): DNA strand (True for forward strand, False for 
                    complementary strand)
            genes (list of str.): list of locus tags of genes associated with 
                    the TSS
            score (float.): TSS score

        Example:
            >>> from TSS_TTS_TU import TSS
            >>> tss = TSS(pos=4707030,strand=False)
            >>> tss.pos
            4707030
        """
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
            self: TSS instance
            tags (str.): locus tags separated by commas
            genes_dict: dictionary of shape {locus tag : Gene object}

        Example:
            >>> import Genome
            >>> g = Genome.Genome("dickeya")
            >>> g.load_annotation()
            >>> from TSS_TTS_TU import TSS
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

        Required args:
            self: TSS instance
            sig (str.): sigma factor

        Optional arg:
            sites (str.): the left and right coordinates of the -10 and -35 
                    elements with the shape "-10l,-10r,-35r,-35l"
            genes_dict: dictionary of shape {locus tag : Gene object}

        Output:
            self.promoter (dict.): completed attribute of the TSS instance.
                                   with the shape {[sig]:(sites)}

        Example:
            >>> from TSS_TTS_TU import TSS
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
        on -35 and -10 coordinates that are included in the promoter 
        attribute (which must first be loaded using the add_promoter method).

        Args:
            self (TSS instance with a promoter attribute)
            gen_seq (str.): sequence of the forward strand
            gen_seqcompl (str.): sequence of the complementary strand
                    (WARNING: 3' -> 5')
            shift (Optional [int.]): number of nucleotides to include beyond 
                    each region on either side (Default: 0nt)
            prom (Optional [int.,int.]): region upstream and downstream TSS
                    to extract. Argument with the shape:
                    [length before TSS, length after TSS]
                    Default: [0,0] ie no sequence will be extracted around 
                    TSS

        Outputs:
            self.promoter (dict.): completed attribute of the TSS instance.
                    Before this function, this dictionary of shape 
                    {sigma factor: subdictionary} contained, for each sigma 
                    factor, a subdictionary of shape "sites":(-10l,-10r,-35r,
                    -35l)} with -10l, -10r, -35r and -35l the left and right 
                    coordinates of the -10 and -35 elements. Now, this 
                    subdictionary contains 4 new keys : 
                    "spacer", "minus10","minus35" and "discriminator". 
                    Each associated value is the sequence of the element.

                if prom !=[0,0], self.promoter dictionary has also a new 
                "region" key. self.promoter['region'] returns the sequence 
                of the region chosen with the prom argument.

            New TSS attributes promoter
            Dictionnary of shape {element: sequence} with element in
            ["spacer", "minus10", "minus35", "discriminator", "region"]
            See load_TSS description to understand the structure of the
            subdictionary self.TSSs[condTSS][TSS].promoter


        Example:
            >>> import Genome
            >>> g = Genome.Genome("dickeya")
            >>> g.load_seq()
            >>> from TSS_TTS_TU import TSS
            >>> tss = TSS(pos=4707030,strand=False)
            >>> tss.add_promoter(sig = "sigma70", sites="4707037,4707046,4707062,4707068")
            >>> tss.add_prom_elements(g.seq,
                                      g.seqcompl,
                                      sig = "sigma70",
                                      sites="4707039,4707044,4707060,4707065",
                                      prom=[20,20])
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

    def __init__(self, start=None, end=None, left=None, right=None,
                 strand=None, rho_dpdt=None, genes=[], seq="", score=None):
        """
        Called when a TTS instance is created,initializes the following 
        attributes: left, right, start, end, strand, rho_dpdt, genes, score 
        and seq.
        Args:
            self: TTS instance
            left (int.) and right (int.): TTS coordinates (does not take into
                    account the strand, ie right > left)
            start (int.) and end (int.): positions of the beginning and the 
                    end of the TTS
            strand (bool.): DNA strand (True for forward strand, False for 
                    complementary strand)
            rho_dpdt (bool.): rho dependency of the TTS
            genes (list of str.): list of locus tags of genes associated with 
                    the TTS
            score (float.): TSS score

        Example:
            >>> from TSS_TTS_TU import TSS
            >>> tts = TTS(left=1,right=10,strand=True)
            >>> tss.left
            1
        """
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
            self: TTS instance
            tags (str.): locus tags separated by commas
            genes_dict: dictionary of shape {locus tag : Gene object}

        Example:
            >>> import Genome
            >>> g = Genome.Genome("dickeya")
            >>> g.load_annotation()
            >>> from TSS_TTS_TU import TTS
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

    def __init__(self, start=None, end=None, left=None,
                 right=None, strand=None, genes=[]):
        """
        Called when a TU instance is created,initializes the following 
        attributes: left, right, start, end, strand, and genes.

        Args:
            self: TU instance
            left (int.) and right (int.): TU coordinates (does not take into
                    account the strand, ie right > left)
            start (int.) and end (int.): positions of the beginning and the 
                    end of the TU
            strand (bool.): DNA strand (True for forward strand, False for 
                    complementary strand)
            genes (list of str.): list of locus tags of genes associated with 
                    the TU

        Example:
            >>> from TSS_TTS_TU import TU
            >>> tu = TU(start=100,end=10,strand=False)
            >>> tu.left
            10
        """
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


    def add_genes(self, tags, genes_dict):
        """
        add_genes completes the genes attribute of the TU instance
        with the genes that are given as tags and that are in the genes_dict.

        Args:
            self: TU instance
            tags (str.): locus tags separated by commas
            genes_dict: dictionary of shape {locus tag : Gene object}

        Example:
            >>> import Genome
            >>> g = Genome.Genome("dickeya")
            >>> g.load_annotation()
            >>> from TSS_TTS_TU import TU
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

        Args:
            self: TU instance
            condition (str.): name of the condition
            TU_expr (float or list of float): expression value(s)

        Output:
            self.expression (dict.): New (or completed) attribute
                                     of the TU instance with the shape
                                     {condition:expression value}
        Example:
            >>> from TSS_TTS_TU import TU
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

        Args:
            self: TU instance
            TSS: tuple of shape (TSS position (int.),
                                 propotion of total starts (float))

        Output:
            self.TSS (list of tuples): New (or completed) attribute of the TU
                    instance containing the TSS tuple given as argument.

        Example:
            >>> from TSS_TTS_TU import TU
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

        Args:
            self: TU instance
            TTS: tuple of shape 
                    (TTS position (int.), propotion of total starts (float))

        Output:
            self.TTS (list of tuples): New (or completed) attribute of the TU
                    instance containing the TTS tuple given as argument.
        """
        if not hasattr(self, 'TTS'):
            self.TTS = []

        if TTS not in self.TTS:
            self.TTS.append(TTS)