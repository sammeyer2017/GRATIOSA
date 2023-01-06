#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np

#==============================================================================#

class Gene:

    def __init__(self, locus_tag,name,ID,left,right,strand,ASAP_name):
        """
        DECRIRE CHAQUE ATTRIBUT !!
        """
        self.locus_tag = locus_tag
        self.ID = ID 
        self.name = name
        self.left = int(left)
        self.right = int(right)
        self.strand = strand
        self.middle = float(left+right)/2
        self.length = int(right)-int(left)+1
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
        self.rpkm[condition] = (expression_value*np.power(10,9)) / (factor)

    def set_model_parameters(self, k_on, p_off, condition):
        """ Sets the gene's k_on and p_off parameters. p_off is the probability
        that the polymerase that transcribed the previous gene does not
        transcribe this one as well (readthrough). k_on is the quantity of new
        polymerases binding to and transcribing this gene.
        """
        if not hasattr(self, 'model_parameters'):
            self.model_parameters = {}
        self.model_parameters[condition] = (k_on, p_off)

    def __eq__(self, other):
        return (self.name == other.name) and \
               (self.left == other.left) and \
               (self.right == other.right) and \
               (self.orientation == other.orientation)

    def add_id_TSS(self, id_TSS):
        if not hasattr(self,'id_TSS'):
            self.id_TSS=[]
        if not (id_TSS in self.id_TSS):
            self.id_TSS.append(id_TSS)

    def add_fc_pval_cond(self,fc_value, condition, p_value):
        if not hasattr(self,'fc_pval'):
            self.fc_pval={}
        self.fc_pval[condition]=(fc_value,p_value)

    def add_state_cond(self,cond,state):
        if not hasattr(self,'state'):
            self.state={}
        self.state[cond]=state

    def add_left_neighbor(self, neighbor):
        self.left_neighbor = neighbor

    def add_right_neighbor(self, neighbor):
        self.right_neighbor = neighbor

    def add_orientation(self,orient):
        self.orientation = orient

    def add_id_operon(self, operon):
        if not hasattr(self,'id_operon'):
            self.id_operon=[]
        if not(operon in self.id_operon):
            self.id_operon.append(operon)

    def add_expression_data(self, conditions, expression_values):
        """ Adds expression in the form of a dictionnary where the keys
        correspond to the different conditions. The conditions and expression
        values are passed as lists.
        """
        if not hasattr(self,'expression'):
            self.expression = {}
        
        self.expression.update(dict(zip(conditions, expression_values)))

        #self.mean_expression = np.mean(self.expression.values())


    def add_single_expression(self, condition, expression_value):
        if not hasattr(self, 'expression'):
            self.expression = {}
        self.expression[condition] = expression_value

    def set_mean_expression(self, expression_value=None):
        if expression_value:
            self.mean_expression = expression_value
        else:
            self.mean_expression = np.mean(self.expression.values())

    def add_list_expression(self):
        self.list_expression=[]
        if hasattr(self,'expression'):
            for i in self.expression:
                self.list_expression.append(self.expression[i])

    def add_signal(self,cond,signal) :
        if not hasattr(self, "signal"):
            self.signal = {}
        self.signal[cond] = signal

    def add_is_loop(self,cond,is_loop) :
        if not hasattr(self, "is_loop"):
            self.is_loop = {}
        self.is_loop[cond] = is_loop

    def add_is_border(self,cond,is_border) :
        if not hasattr(self, "is_border"):
            self.is_border= {}
        self.is_border[cond] = is_border
    
    def add_is_peak(self,cond,is_peak) :
        if not hasattr(self, "is_peak"):
            self.is_peak= {}
        self.is_peak[cond] = is_peak


    def add_GO(self,annot,GO) :
        if not hasattr(self, "GO"):
            self.GO= {}
        self.GO[annot] = GO