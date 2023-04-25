#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

class Gene:
    '''
    The Gene class allows to gather information at the level of a single gene:
    its identifiers, its coordinates, its expression levels under different
    conditions, its functional annotations...
    '''
    def __init__(self, locus_tag, name, ID, left, right, strand, ASAP_name):
        self.locus_tag = locus_tag
        self.ID = ID
        self.name = name
        self.left = int(left)
        self.right = int(right)
        self.strand = strand
        self.middle = float(left + right) / 2
        self.length = int(right) - int(left) + 1
        if self.strand:
            self.start = int(left)
            self.end = int(right)
        else:
            self.start = int(right)
            self.end = int(left)
        self.ASAP = ASAP_name

    def add_single_rpkm(self, condition, expression_value, factor):
        if not hasattr(self, 'rpkm'):
            self.rpkm = {}

    def add_id_TSS(self, id_TSS):
        if not hasattr(self, 'id_TSS'):
            self.id_TSS = []
        if not (id_TSS in self.id_TSS):
            self.id_TSS.append(id_TSS)

    def add_fc_pval_cond(self, fc_value, condition, p_value):
        if not hasattr(self, 'fc_pval'):
            self.fc_pval = {}
        self.fc_pval[condition] = (fc_value, p_value)

    def add_state_cond(self, cond, state):
        if not hasattr(self, 'state'):
            self.state = {}
        self.state[cond] = state

    def add_left_neighbor(self, neighbor):
        self.left_neighbor = neighbor

    def add_right_neighbor(self, neighbor):
        self.right_neighbor = neighbor

    def add_orientation(self, orient):
        self.orientation = orient

    def add_expression(self, condition, expression_value):
        if not hasattr(self, 'expression'):
            self.expression = {}
        self.expression[condition] = expression_value

    def set_mean_expression(self, expression_value=None):
        if expression_value:
            self.mean_expression = expression_value
        else:
            self.mean_expression = np.mean(self.expression.values())

    def add_signal(self, cond, signal):
        if not hasattr(self, "signal"):
            self.signal = {}
        self.signal[cond] = signal

    def add_is_loop(self, cond, is_loop):
        if not hasattr(self, "is_loop"):
            self.is_loop = {}
        self.is_loop[cond] = is_loop

    def add_is_border(self, cond, is_border):
        if not hasattr(self, "is_border"):
            self.is_border = {}
        self.is_border[cond] = is_border

    def add_is_peak(self, cond, is_peak):
        if not hasattr(self, "is_peak"):
            self.is_peak = {}
        self.is_peak[cond] = is_peak

    def add_GO(self, annot, GO):
        if not hasattr(self, "GO"):
            self.GO = {}
        self.GO[annot] = GO