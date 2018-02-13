#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np

#==============================================================================#

class Terminator:
    def __init__(self, left, right, strand, ID, dict_terminator):
        self.start = left
        self.end = right
        self.strand = strand
        self.ID=ID
        if strand == '-':
            self.left=self.end
            self.right=self.start
        else:
            self.left=self.start
            self.right=self.end
        for x in dict_terminator:
            if 'Hp' in x:
                self.hp=dict_terminator[x]
            if 'Tail' in x:
                self.tail=dict_terminator[x]
            if 'Regionconf' in x:
                self.regionconf=dict_terminator[x]
            if 'opp_overlap' in dict_terminator[x]:
                line = dict_terminator[x]
                self.opp_overlap = line.split()[1::]
