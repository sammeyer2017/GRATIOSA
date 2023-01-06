#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import math
from globvar import *

#==============================================================================#

##### functions called by transcriptome methods #####
def add_single_rpkm_to_genes(genes_dict, expression_filename, condition, TSS_column, start_line, separator,tag_column):
    """ Adds rpkm data to Gene objects by parsing a file with two columns:
    gene name and value
    """
    if separator == '\\t' : separator = '\t' 
    with open(expression_filename, 'r') as f:
        i=1
        while i != start_line:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n').split(separator)
            try:
                genes_dict[line[tag_column]].add_single_rpkm(
                    #log2(x+1)
                    condition, float(line[TSS_column]))
            except:
                if line[tag_column] not in list(genes_dict.keys()):
                    # the rpkm value corresponds to an un-annoted gene
                    print("rpkm : gene " + line[tag_column] + " not in annotation")
                else:
                    # the rpkm value cannot be converted
                    genes_dict[line[tag_column]].add_single_rpkm(condition, float("NaN"))
    # look if some genes are not in rpkm file: add nan
    for g in list(genes_dict.keys()):
        if isinstance(genes_dict[g],Gene):
            if not hasattr(genes_dict[g], 'rpkm'):
                genes_dict[g].add_single_rpkm(condition, float("NaN"))
    return genes_dict


def add_expression_to_genes(genes_dict, filename, tag_col, first_expression_col, is_log, separator):
    """ Adds expression data to Gene objects by parsing a file with as many
    columns as there are different conditions in the experiment, plus one for
    the gene names (first column).
    """        
    genes_valid = {} 
    except_locus = []
    if separator == '\\t' : separator = '\t' 
    with open(filename, 'r') as f:
        header=next(f)
        header=header.strip().split(separator)        
        header=header[first_expression_col:]
        genes_valid["conditions"] = header
        for line in f:
            line=line.strip().split(separator)
            try:
                if is_log == 'no':
                    genes_dict[line[tag_col]].add_expression_data(header,[math.log(float(i),2) for i in line[first_expression_col:]])
                else:
                    genes_dict[line[tag_col]].add_expression_data(header,[float(i) for i in line[first_expression_col:]])
                genes_valid[line[tag_col]] = genes_dict[line[tag_col]]
            except KeyError:
                if line[tag_col] == 'none':
                    print("expressions without locus tag")
                else:
                    except_locus.append(line[tag_col])
    if except_locus :
        if len(except_locus) > 20 :
            print(f"{len(except_locus)} locus are not in annotation")
        else : 
            print(f"{except_locus} are not in annotation" )
    return genes_valid


def load_fc_pval_cond(genes_dict, filename, condition, tag_col, fc_col, separator, start_line, *args, **kwargs):
    ''' Called by load_fc_pval, allows expression data to be loaded by specifying files, and where each
    information is (tag, fc, pval...). If no p-value column, assigns pval = 0 to each gene
    '''
    genes_valid = [] # list containing all genes having valid FC / pval
    p_val_col= kwargs.get('p_value')
    if separator == '\\t' : separator = '\t' 
    except_locus = []
    except_fc = 0 
    with open(filename, 'r') as f:
        i=0
        while i < start_line:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n').split(separator)
            try:
                if p_val_col:
                    genes_dict[line[tag_col]].add_fc_pval_cond(float(line[fc_col]),condition, float(line[p_val_col]))
                else:
                    genes_dict[line[tag_col]].add_fc_pval_cond(float(line[fc_col]),condition, float(0))
                genes_valid.append(line[tag_col])
            except:
                if line[tag_col] not in genes_dict.keys():
                    if line[tag_col] != '':
                        except_locus.append(line[tag_col])
                    else:
                        except_fc += 1
    f.close()
    if except_locus :
        if len(except_locus) > 20 :
            print(f"\t{len(except_locus)} locus are not in annotation")
        else : 
            print(f"\t{except_locus} are not in annotation" )
    if except_fc :
        print(f"\t{except_fc} fc without locus")
    return genes_valid

def load_expr_cond(genes_dict, filename, condition, tag_col, nb_replicates, expr_col, separator, start_line):
    ''' Called by load_expr_cond, allows expression data to be loaded by specifying files, and where each
    information is (tag, expr...).
    '''
    if separator == '\\t' : separator = '\t'
    except_locus = []
    with open(filename, 'r') as f:
        i=1
        while i < start_line:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n').split(separator)
            st = int(expr_col) ; nb = int(nb_replicates)
            for i in range(st,st+nb):
                try:
                    genes_dict[line[tag_col]].add_expr_cond(float(line[i]),condition)
                except:
                    if line[tag_col] not in genes_dict.keys():
                        if line[tag_col] != '':
                            except_locus.append(line[tag_col])
                        else:
                            print("fc without locus")
    f.close()
    if except_locus :
        if len(except_locus) > 20 :
            print(f"{len(except_locus)} locus are not in annotation")
        else : 
            print(f"{except_locus} are not in annotation" )

def set_mean_expression(genes_dict, expression_filename):
    """ 
    For each gene of genome object, set mean expression value based on loaded expression data in various conditions
    """
    with open(expression_filename, 'r') as f:
        header = next(f)
        header = header.strip('\n').split(',')
        header = header[1:]
        for line in f:
            line = line.strip('\n')
            line = line.split(',')
            try:
                genes_dict[line[0]].set_mean_expression(line[1])
            except KeyError:
                print("Expressions : Could not find gene " + line[0])
    return genes_dict