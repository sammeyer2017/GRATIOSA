#! /usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import os
import numpy as np
from Gene import Gene

#==============================================================================#

def add_list_condition(list_genes):
    condition=[]
    for i in list_genes:
        if hasattr(i,'expression'):
            for j in i.expression:
                condition.append(j)
            return condition


class Domain:

    def __init__(self, *args, **kwargs):
        self.start=kwargs.get("start")
        self.end=kwargs.get("end")
        self.index=kwargs.get("index")
        self.genes=kwargs.get("genes")
        if not self.genes:
            self.genes=[]
            self.condition=[]
            self.length=self.end-self.start
            if self.start > self.end:
                self.start_deb= True
            else:
                self.start_deb= False

    def add_genes(self,dict_genes):
        for i in dict_genes:
            if self.start_deb == False:
                if dict_genes[i].start > self.start:
                    if dict_genes[i].end < self.end:
                        self.genes.append(dict_genes[i])
            else:
                if dict_genes[i].start > self.start:
                    self.genes.append(dict_genes[i])
                if dict_genes[i].start < self.end:
                    if dict_genes[i].end > self.start:
                        self.genes.append(dict_genes[i])
                    if dict_genes[i].end < self.end:
                        self.genes.append(dict_genes[i])

    def add_domain_expression(self):
        if self.genes==[]:
            print('load genes in domain please')
        else:
            self.expression={}
            self.condition=add_list_condition(self.genes)
            for i in self.condition:
                expression=0.0
                n=0.0
                for j in self.genes:
                    if hasattr(j,'expression'):
                        expression+=float(j.expression[i])
                        n+=1.0
                self.expression[i]= expression/n

    def add_list_expression(self):
        self.list_expression=[]
        for i in self.condition:
            self.list_expression.append(self.expression[i])
        for i in self.genes:
            i.add_list_expression()
