#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Vincent CABELI
"""

import sys
import os
import numpy as np
from Gene import Gene

#==============================================================================#

class Domain:

    def __init__(self, genes_list):
        self.genes_list = genes_list
        for gene in genes_list:
            if not isinstance(gene, Gene):
                raise TypeError("The genes list must contain Gene objects.")
        self.number_of_genes = len(genes_list)
        self.genes_list.sort(key = lambda x: x.start)
        self.genes_orientations = np.asarray([gene.orientation for gene in genes_list])
        self.genes_mean_expressions = [gene.mean_expression for gene in genes_list]
        #self.genes_numbers = [gene.number for gene in genes_list]
        self.genes_middle = [(gene.left+gene.right)/2 for gene in genes_list]
        self.left = np.min([gene.left for gene in genes_list])
        self.right = np.max([gene.right for gene in genes_list])

    def reverse_genes_list(self):
        self.genes_list = self.genes_list[::-1]
        #!!
        self.genes_orientations = np.invert(self.genes_orientations)[::-1]
        #self.genes_numbers = [gene.number for gene in self.genes_list]
        self.genes_middle = [(gene.left+gene.right)/2 for gene in self.genes_list]

    def find_neighbours(self, all_genes_list, n=1):
        """ Find nearest genes to the left and right from a list, where n is 
        the number of genes to add at each extremity.
        """
        sorted_genes = all_genes_list
        sorted_genes.sort(key = lambda x: x.left)
        sorted_pos = np.asarray([gene.left for gene in sorted_genes])
        self.genes_list.append(sorted_genes[np.where(sorted_pos < self.left)[0][-n:]])
        self.genes_list.append(sorted_genes[np.where(sorted_pos > self.right)[0][:n]])
        self.genes_list.sort(key = lambda x: x.start)
        self.number_of_genes = len(self.genes_list)
        self.genes_orientations = np.asarray([gene.orientation for gene in self.genes_list])
        self.genes_mean_expressions = [gene.mean_expression for gene in self.genes_list]
        #self.genes_numbers = [gene.number for gene in self.genes_list]
        self.genes_middle = [(gene.left+gene.right)/2 for gene in self.genes_list]
        self.left = np.min([gene.left for gene in self.genes_list])
        self.right = np.max([gene.right for gene in self.genes_list])

    def includes(self, other):
        """ Checks whether self include a lesser domain of genes, other.
        """
        if self.number_of_genes < other.number_of_genes:
            return False
        for gene in other.genes_list:
            if gene not in self.genes_list:
                return False
        return True

    def __eq__(self, other):
        if self.number_of_genes != other.number_of_genes:
            return False
        for gene in self.genes_list:
            if gene not in other.genes_list:
                return False
        return True
