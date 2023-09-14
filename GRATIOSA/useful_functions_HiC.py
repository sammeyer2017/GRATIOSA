#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions called by HiC methods
"""


def load_HiC_site_cond(site_type, 
                       path2file, 
                       startline, 
                       sep, 
                       bin1_col,
                       binsize, 
                       bin2_col=None, 
                       score_col=None, 
                       pval_col=None, 
                       qval_col=None):
    '''
    Called by load_HiC_borders and load_HiC_loops, allows borders and loops 
    data to be loaded by specifying files, and where each information is (one 
    column for each type of information).

    Args:
        type (str.): HiC detected motif type ("borders" or "loops")
        path2file (str.): path to the file containing the data
        startline (int.): file start line
        sep (str.): file separator
        bin1_col (int.): index of the column containing the position of the 
                first bin
        binsize (int.): binsize (in b)
        bin2_col (Optional [int.]): index of the column containing the position of the
                second bin
        score_col (Optional [int.]): index of the column containing the position of the score
        pval_col (Optional [int.]): index of the column containing the position of the pvalue
        qval_col (Optional [int.]): index of the column containing the position of the qvalue

    Note: 
        Column numbering starts at 0.
    '''
    if sep == '\\t':
        sep = '\t'

    if site_type == "borders":
        Borders = {}
    elif site_type == "loops":
        Loops = {}

    with open(path2file, 'r') as f:
        i = 0
        while i < startline:
            header = next(f)
            i += 1
        for line in f:
            line = line.strip('\n').split(sep)
            try:
                if site_type.lower() == "borders":
                    b = int(line[int(bin1_col)]) * binsize
                    Borders[b] = {"binsize": binsize}
                    if score_col:
                        Borders[b]['score'] = line[int(score_col)]
                    if pval_col:
                        Borders[b]['pval'] = line[int(pval_col)]
                    if qval_col:
                        Borders[b]['qval'] = line[int(qval_col)]
                elif site_type.lower() == "loops":
                    pos_tuple = (
                        int(line[int(bin1_col)]) * binsize, int(line[int(bin2_col)]) * binsize)
                    Loops[pos_tuple] = {"binsize": binsize}
                    if score_col:
                        Loops[pos_tuple]['score'] = line[int(score_col)]
                    if pval_col:
                        Loops[pos_tuple]['pval'] = line[int(pval_col)]
                    if qval_col:
                        Loops[pos_tuple]['qval'] = line[int(qval_col)]
            except Exception as e:
                print(e)
    f.close()

    if site_type == "borders":
        return Borders

    elif site_type == "loops":
        return Loops
