#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd
import inspect as inspect
import math
import matplotlib.pyplot as plt
import operator
from useful_function import *
from itertools import groupby
from globvar import *
from Gene import Gene
from TSS import TSS
from TU import TU
from Operon import Operon
from Terminator import Terminator
from domain_phuc import Domain
from datetime import datetime
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from matplotlib import patches
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import sqrt


#==============================================================================#

# -------------------
# useful function


def annotations_parser(annotations_filename, tag_column=1):
    """ Returns a dictionary of Gene objects which has the shape
    {gene_tag_name:gene_object}, where the tag name column can be specified.
    """
    genes_dict = {}
    with open(annotations_filename, 'r') as f:
        header = next(f)
        for line in f:
            line = line.strip('\n')
            line = line.split(',')
            try:
                if (line[tag_column] in genes_dict):
                    print("Warning! Overwriting value for gene " + \
                         line[tag_column])
                # if key already exists this will overwrite the value
                genes_dict[line[tag_column]] = Gene(annotations_list=line)
            except (IndexError, ValueError):
                print("Annotations : could not read line " + ','.join(line))
    return genes_dict

def annotations_parser_general_old(annotations_filename,separator,tag_column,strand_column,left_column,         right_column,start_line):
    genes_dict = {}
    with open(annotations_filename, 'r') as f:
        i=1
        j=0
        while i < start_line:
            header=next(f)
            i+=1
        for line in f:
            list=[]
            line=line.strip()
            if      separator == '\\t':
                line=line.split('\t')
            else:
                line=line.split(separator)
            list.append(line[tag_column])
            if line[strand_column]=="complement":
                list.append('-')
            elif line[strand_column]=="forward":
                list.append('+')
            elif line[strand_column]== "1":
                list.append('+')
            elif line[strand_column]== "-1":
                list.append('-')
            else:
                list.append(line[strand_column])
            list.append(line[left_column])
            list.append(line[right_column])
            if line[tag_column] in genes_dict:
                print("Warning! Overwriting value for gene ")
                print(line[tag_column])
            genes_dict[line[tag_column]]=Gene(annotations_general=list)

    return genes_dict


def annotations_parser_general(annotations_filename,separator,tag_column,strand_column,left_column,right_column,start_line):
    genes_dict = {}
    with open(annotations_filename, 'r') as f:
        i=1
        j=0
        head=next(f)
        headl=head.strip()
        if separator == '\\t':
            headl=headl.split('\t')
        else:
            headl=headl.split(separator)
        i=2
        while i < start_line:
            header=next(f)
            i+=1
        for line in f:
            list=[]
            line=line.strip()
            if      separator == '\\t':
                line=line.split('\t')
            else:
                line=line.split(separator)
            list.append(line[tag_column])
            if line[strand_column]=="complement":
                list.append('-')
            elif line[strand_column]=="forward":
                list.append('+')
            elif line[strand_column]== "1":
                list.append('+')
            elif line[strand_column]== "-1":
                list.append('-')
            else:
                list.append(line[strand_column])
            list.append(line[left_column])
            list.append(line[right_column])
            list.append(line)
            if line[tag_column] in genes_dict:
                print("Warning! Overwriting value for gene ")
                print(line[tag_column])
            genes_dict[line[tag_column]]=Gene(annotations_general=list, head=headl)

    return genes_dict



def annotations_parser_gff(annotations_filename):
    genes_dict = {}
    under_dict={}
    my_file = open(annotations_filename, "r")
    for line in my_file.readlines():
        if line[0]!= '#':
            if line != '\n':
                line=line.split('\t')
                underline=line[8]
                underline=underline.split(';')
                for x in underline:
                    x=x.strip()
                    x=x.split('=')
                    under_dict[x[0]]=x[1]
                line[8]=under_dict
                under_dict={}
                try:
                    if('locus_tag' in line[8]):
                        if(line[8]['locus_tag'] in genes_dict):
                            print("Warning! Overwriting value for gene ")
                            print(line[8]['locus_tag'])
                            if('gene_biotype' in line[8]):
                                print(line[8]['gene_biotype'])
                        genes_dict[line[8]['locus_tag']]= Gene(annotations_list_gff=line)
                except (IndexError, ValueError):
                    print("Annotations : could not read line ")
    return genes_dict

def add_operon(dict_genes,file):
    dict_operon={}
    under_dict={}
    i = 1
    my_file = open(file, "r")
    for line in my_file.readlines():
        if line[0]!= '#':
            line=line.strip()
            line=line.split('\t')
            dict_operon[i]=Operon(int(line[3]),int(line[4]),line[6])
            underline=line[8]
            underline=underline.split(';')
            for x in underline:
                x=x.split('=')
                under_dict[x[0]]=x[1]
            for x in under_dict:
                if 'Genes' in x:
                    dict_operon[i].add_genes(under_dict[x])
                if 'ID' in x:
                    dict_operon[i].add_ID(under_dict[x])
            for x in dict_operon[i].genes:
                if x in list(dict_genes.keys()):
                    dict_genes[x].add_id_operon(i)
                else:
                    print(x+" not in annotation")
            i+=1
            under_dict={}
    return dict_operon

def add_terminator(file):
    dict_terminator={}
    under_dict={}
    my_file = open(file, "r")
    for line in my_file.readlines():
        if line[0]!= '#':
            line=line.strip()
            line=line.split('\t')
            underline=line[8]
            underline=underline.split(';')
            for x in underline:
                x=x.split('=')
                under_dict[x[0]]=x[1]
            for x in under_dict:
                if 'ID' in x:
                    dict_terminator[int(under_dict[x])]=Terminator(int(line[3]),int(line[4]),line[6], int(under_dict[x]), under_dict)
    return dict_terminator


def feature_table_parser(table_filename):
    """ Returns a dictionary of Gene objects which has the shape
    {gene_tag_name:gene_object}, where the tag name column can be specified.
    """
    genes_dict = {}
    with open(table_filename, 'r') as f:
        next(f)
        for line in f:
            line = line.strip('\n')
            line = line.split('\t')
            if len(line) != 3:
                continue
            elif line[2] == 'gene':
                start = int(line[0].strip('<').strip('>'))
                end = int(line[1].strip('<').strip('>'))
                orientation = (start < end)
                while 'locus_tag' not in line :
                    line = next(f)
                    line = line.strip('\n')
                    line = line.split('\t')
                if 'gene' in line:
                    continue
                if line[-1] in genes_dict:
                    print("Warning! Overwriting value for gene " + \
                         line[-1])
                # if key already exists this will overwrite the value
                genes_dict[line[-1]] = Gene(name = line[-1],
                                           left = min(start,end),
                                           right = max(start, end),
                                           orientation = orientation)
    return genes_dict


def add_expression_to_genes(genes_dict, expression_filename, tag_col, first_expression_col, is_log):
    """ Adds expression data to Gene objects by parsing a file with as many
    columns as there are different conditions in the experiment, plus one for
    the gene names (first column).
    """
    with open(expression_filename, 'r') as f:
        header=next(f)
        header=header.strip()
        header=header.split('\t')
        header=header[first_expression_col:]
        for line in f:
            line=line.strip()
            line = line.split('\t')
            try:
                if is_log == 'no':
                    genes_dict[line[tag_col]].add_expression_data(header,[math.log(float(i),2) for i in line[first_expression_col:]])
                else:
                    genes_dict[line[tag_col]].add_expression_data(header,[float(i) for i in line[first_expression_col:]])
            except KeyError:
                if line[tag_col] == 'none':
                    print("expressions without locus tag")
                else:
                    print(line[tag_col] + " this locus not in annotation")
        return genes_dict


def add_single_expression_to_genes(genes_dict, expression_filename):
    """ Adds expression data to Gene objects by parsing a file with as many
    columns as there are different conditions in the experiment, plus one for
    the gene names (first column).
    """
    condition = expression_filename[expression_filename.rfind('/')+1:
                                    expression_filename.rfind('/')+3]
    with open(expression_filename, 'r') as f:
        header = next(f)
        header = header.strip('\n').split(',')
        header = header[1:]
        for line in f:
            line = line.strip('\n')
            line = line.split(',')
            try:
                genes_dict[line[0]].add_single_expression(
                    #log2(x+1)
                    condition, np.log2(float(line[1])+1))
            except KeyError:
                print("Expressions : Could not find gene " + line[0])
                #genes_dict[line[0]] = Gene(name=line[0], left=int(line[2]),
                #                          right=int(line[3]),
                #                          orientation=(line[4]=='+'))
                #genes_dict[line[0]].add_single_expression(
                #    expression_filename[i:i+3], float(line[1]))
    return genes_dict



def add_single_rpkm_to_genes(genes_dict, expression_filename, condition, TSS_column, start_line, separator,tag_column):
    """ Adds rpkm data to Gene objects by parsing a file with two columns:
    gene name and value
    """
    with open(expression_filename, 'r') as f:
        i=1
        while i != start_line:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n')
            if separator == '\\t':
                line = line.split('\t')
            else:
                line=line.split(separator)
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


def set_mean_expression(genes_dict, expression_filename):
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



def filter_TSS_old(xxx_todo_changeme,filt,win):
    (plus,minus) = xxx_todo_changeme
    isreal=np.any(plus[:,1:]>=filt,axis=1)
    plus=plus[isreal]
    ordre=np.argsort(plus[:,0])
    plus=plus[ordre]
    print(len(plus))
    # --- group together close ones
    finplus=[]
    if isinstance(win,int):
        for ip,p in enumerate(plus[:-1]):
            if (plus[ip+1,0]-p[0])<win:
                if np.mean(p[1:])<np.mean(plus[ip+1,1:]):
                    # replace position of TSS to group them
                    ind=plus[ip+1,0]
                else:
                    ind=plus[ip,0]
                finplus.append(ind)
        if plus[-1,0]!=ind:
            finplus.append(ind)
                #line=p+plus[ip+1]
                #line[0]=ind
                #finplus.append(line)
            #plus=np.array(finplus)
    else:
        finplus=plus[:,0]
    # ---------------------------------
    # -------- - strand
    isreal=np.any(minus[:,1:]>=filt,axis=1)
    minus=minus[isreal]
    ordre=np.argsort(minus[:,0])
    minus=minus[ordre]
    # ---- group together close ones
    finminus=[]
    if isinstance(win,int):
        for ip,p in enumerate(minus[:-1]):
            if (minus[ip+1,0]-p[0])<win:
                #print p
                if np.mean(p[1:])<np.mean(minus[ip+1,1:]):
                    # replace position of TSS to group them
                    ind=minus[ip+1,0]
                else:
                    ind=minus[ip,0]
                finminus.append(ind)
        if minus[-1,0]!=ind:
            finminus.append(ind)
                #line=p+plus[ip+1]
                #line[0]=ind
                #finplus.append(line)
                #plus=np.array(finplus)
    else:
        finminus=minus[:,0]
    return np.array(finplus),np.array(finminus)

def filter_TSS(xxx_todo_changeme1,filt,win):
    (plus,minus) = xxx_todo_changeme1
    isreal=np.any(plus[:,1:]>=filt,axis=1)
    plus=plus[isreal]
    ordre=np.argsort(plus[:,0])
    plus=plus[ordre]
    # --- group together close ones
    dftrue=plus[1:,0]-plus[:-1,0]<win
    a=[(k,len(list(g))) for k,g in groupby(dftrue)]
    # ind is the first index
    ind=0
    if a[0][0]:
        pluslist=[]
    else:
        pluslist=[plus[0,0]]
    for k,l in a:
        if k:
            # l separations are small
            # between index ind and index ind+l+1
            pluslist.append(plus[ind+np.argmax(np.mean(plus[ind:(ind+l+1),1:],axis=1)),0])
            ind+=l
        else:
            if l>1:
                pluslist+=plus[(ind+1):(ind+l),0].tolist()
            ind+=l
    # -------- minus strand$
    isreal=np.any(minus[:,1:]>=filt,axis=1)
    plus=minus[isreal]
    ordre=np.argsort(plus[:,0])
    plus=plus[ordre]
    # --- group together close ones
    dftrue=abs(plus[1:,0]-plus[:-1,0])<win
    a=[(k,len(list(g))) for k,g in groupby(dftrue)]
    # ind is the first index
    ind=0
    if a[0][0]:
        minuslist=[]
    else:
        minuslist=[plus[0,0]]
    for k,l in a:
        if k:
            # l separations are small
            # between index ind and index ind+l+1
            minuslist.append(plus[ind+np.argmax(np.mean(plus[ind:(ind+l+1),1:],axis=1)),0])
            ind+=l
        else:
            if l>1:
                minuslist+=plus[(ind+1):(ind+l),0].tolist()
            ind+=l
    return np.array(pluslist),np.array(minuslist)



def add_fc_to_genes(genes_dict, fc_filename):
    """ Adds expression data to Gene objects by parsing a file with as many
    columns as there are different conditions in the experiment, plus one for
    the gene names (first column).
    """
    with open(fc_filename, 'r') as f:
        header = next(f)
        header = header.strip('\n').split(',')
        header = header[1:]
        for line in f:
            line = line.strip('\n')
            line = line.split(',')
            try:
                genes_dict[line[0]].add_fc_data(header,
                    [float(i) for i in line[1:]])
            except KeyError:
                print("Fold change : Could not find gene " + line[0])
    return genes_dict

def add_single_fc_to_genes(genes_dict, filename, condition, tag_col, fc_col, separator, start_line, n, *args, **kwargs):
    p_val_col= kwargs.get('p_value')
    list_valid_genes=[]
    with open(filename, 'r') as f:
        i=1
        while i < start_line:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n')
            if separator == '\\t':
                line = line.split('\t')
            else:
                line=line.split(separator)
            try:
                if p_val_col:
                    genes_dict[line[tag_col]].add_fc(float(line[fc_col]),condition, p_value=float(line[p_val_col]))
                else:
                    genes_dict[line[tag_col]].add_fc(float(line[fc_col]),condition)
                list_valid_genes.append(line[tag_col])
            except:
                if line[tag_col] not in genes_dict.keys():
                    if line[tag_col] != '':
                        print(line[tag_col] + " not in annotation ")
                    else:
                        print("fc without locus")
    f.close()
    return list_valid_genes



def domain_parser(genes_dict, domain_filename):
    """ Creates a list of domains from a dictionnary of Gene objects and a text
    file.
    """
    domains = []
    with open(domain_filename, 'r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.split(',')
            try:
                genes_list = [genes_dict[gene_name] for gene_name in line]
                domains.append(Domain(genes_list))
            except KeyError:
                print("Domains : Could not find gene " + line)
    return domains


def operon_domains(genes_dict):
    """ Creates a list of domains from the operon attribute of Gene objects.
    """
    domains = []
    genes = list(genes_dict.values())
    operons_list = [x.operon for x in list(genes_dict.values())]
    operons_list = np.unique(operons_list)
    for operon in operons_list:
        operon_genes = [genes[i] for i in np.where([x.operon==operon
            for x in genes])[0]]
        domains.append(Domain(operon_genes))
    return domains


def load_seq(filename):
    seq=str()
    my_file = open(filename, "r")
    for i in my_file.readlines():
        line=i.strip() #Removes \n
        if line != '':#Inspect if empty line
            if line[0]!=">":
                seq+=line
    my_file.close
    return seq


def load_single_TSS(genes_dict,TSS_dict, filename, TSS_column, start_line , separator, condition, strand, dic_plus, dic_minus, *args, **kwargs):
    sig = kwargs.get('sig')
    tag_column = kwargs.get('tag_column')
    if TSS_dict == {}:
        indice = 1
    else:
        indice = list(TSS_dict.keys())
        indice=len(indice) +1
    with open(filename, 'r') as f:
        i=1
        while i < start_line:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n')
            if separator == '\\t':
                line = line.split('\t')
            else:
                line=line.split(separator)
            try:
                if line[strand] == '+':
                    if not (line[TSS_column] in list(dic_plus.keys())):
                        TSS_dict[indice] = TSS(pos=int(line[TSS_column]),id=indice)
                        TSS_dict[indice].add_condition(condition)
                        TSS_dict[indice].add_strand(line[strand])
                        if sig:
                            if 'Sig' in line[sig]:
                                TSS_dict[indice].add_sig(line[int(sig)])
                        dic_plus[line[TSS_column]]=indice
                        indice +=1
                    else:
                        if not (condition in TSS_dict[dic_plus[line[TSS_column]]].condition):
                            TSS_dict[dic_plus[line[TSS_column]]].add_condition(condition)
                else:
                    if not (line[TSS_column] in list(dic_minus.keys())):
                        TSS_dict[indice] = TSS(pos=int(line[TSS_column]),id=indice)
                        TSS_dict[indice].add_condition(condition)
                        TSS_dict[indice].add_strand(line[strand])
                        if sig:
                            if 'Sig' in line[sig]:
                                TSS_dict[indice].add_sig(line[int(sig)])
                        dic_minus[line[TSS_column]]=indice
                        indice +=1
                    else:
                        if not (condition in TSS_dict[dic_minus[line[TSS_column]]].condition):
                            TSS_dict[dic_minus[line[TSS_column]]].add_condition(condition)
            except:
                pass
            try:
                if len(line[tag_column].split(',')):
                    for tag in line[tag_column].strip().replace(' ','').split(','):
                        try:
                            if line[strand] == '+':
                                if tag !='':
                                    TSS_dict[dic_plus[line[TSS_column]]].add_genes(tag,condition)
                                genes_dict[tag].add_id_TSS(dic_plus[line[TSS_column]])
                            else:
                                if tag !='':
                                    TSS_dict[dic_minus[line[TSS_column]]].add_genes(tag,condition)
                                genes_dict[tag].add_id_TSS(dic_minus[line[TSS_column]])

                        except:
                            if tag not in list(genes_dict.keys()):
                                if tag != '':
                                    print(tag + " not in annotations")
                else:
                    try:
                        if line[strand] == '+':
                            TSS_dict[dic_plus[line[TSS_column]]].add_genes(line[tag_column],condition)
                            genes_dict[tag_column].add_id_TSS(dic_plus[line[TSS_column]])
                        else:
                            TSS_dict[dic_plus[line[TSS_column]]].add_genes(line[tag_column],condition)
                            genes_dict[tag_column].add_id_TSS(dic_minus[line[TSS_column]])
                    except:
                        if line[tag_column] not in list(genes_dict.keys()):
                            if line[tag_column] != '':
                                print(line[tag_column] + " not in annotations")
            except:
                pass

    return genes_dict


def proportion_of(dict, list_genes, composition, seq_plus, seq_minus):
    dict_proportion={}
    dict_proportion['plus']=0.0
    dict_proportion['minus']=0.0
    activated=0
    repressed=0
    for i in list_genes:
        if dict[i].relax_log_fc > 0.0:
            if dict[i].strand:
                dict_proportion['plus']+= float(seq_plus[int(dict[i].left):int(dict[i].right)].count(composition))/float(dict[i].length)
            else:
                dict_proportion['plus']+= float(seq_minus[int(dict[i].left):int(dict[i].right)].count(composition))/float(dict[i].length)
            activated +=1
        if dict[i].relax_log_fc < 0.0:
            if dict[i].strand:
                dict_proportion['minus']+= float(seq_plus[int(dict[i].left):int(dict[i].right)].count(composition))/float(dict[i].length)
            else:
                dict_proportion['minus']+= float(seq_minus[int(dict[i].left):int(dict[i].right)].count(composition))/float(dict[i].length)
            repressed+=1
    if activated !=0:
        dict_proportion['plus']=dict_proportion['plus']/float(activated)
    else:
        dict_proportion['plus']=0.0
    if repressed !=0:
        dict_proportion['minus']=dict_proportion['minus']/float(repressed)
    else:
        dict_proportion['minus']=0.0
    return dict_proportion

def add_neighbour(dict_genes,list):
    for i in range(len(list)):
        if i != 0:
            dict_genes[list[i][1]].add_left_neighbour(list[i-1][1])
        if i != len(list)-1:
            dict_genes[list[i][1]].add_right_neighbour(list[i+1][1])
    return dict_genes

# ----------------------




class Genome:

    def __init__(self, *args, **kwargs):
        """ Possible kwargs arguments: name, seq, length,
        """
        self.name = kwargs.get('name')
        self.length = kwargs.get('length')
        self.genes=kwargs.get('genes')
        self.TSS_complete = {}
        self.TSS_plus={}
        self.TSS_minus={}

    def load_seq(self):
        self.seq=load_seq(basedir+"data/"+self.name+"/sequence.fasta")
        self.seq_reverse=''
        if(self.length):
            if(self.length != len(self.seq)):
                print("Warning not the same length, the new sequence will be removed")
                self.seq=''
        self.length=len(self.seq)
        l=self.seq[::-1]
        l=l.replace('A','t')
        l=l.replace('T','a')
        l=l.replace('C','g')
        l=l.replace('G','c')
        l=l.replace('a','A')
        l=l.replace('t','T')
        l=l.replace('c','C')
        l=l.replace('g','G')
        self.seq_reverse=l

    def give_proportion_of(self,*args, **kwargs):
        self.load_fc()
        self.load_seq()
        composition = kwargs.get('composition')
        if composition:
            values=[0.0,0.0]
            for i in composition:
                values[0]+=proportion_of(self.genes, self.list_genes_fc, i, self.seq, self.seq_reverse)['plus']
                values[1]+=proportion_of(self.genes, self.list_genes_fc, i, self.seq, self.seq_reverse)['minus']
            return values
        else:
            val_plus=[0.0,0.0,0.0,0.0]
            val_minus=[0.0,0.0,0.0,0.0]
            explode = (0,0,0,0)
            name = ['A','T','G','C']
            val_plus[0]+=proportion_of(self.genes, self.list_genes_fc, 'A', self.seq, self.seq_reverse)['plus']
            val_plus[1]+=proportion_of(self.genes, self.list_genes_fc, 'T', self.seq, self.seq_reverse)['plus']
            val_plus[2]+=proportion_of(self.genes, self.list_genes_fc, 'G', self.seq, self.seq_reverse)['plus']
            val_plus[3]+=proportion_of(self.genes, self.list_genes_fc, 'C', self.seq, self.seq_reverse)['plus']
            val_minus[0]+=proportion_of(self.genes, self.list_genes_fc, 'A', self.seq, self.seq_reverse)['minus']
            val_minus[1]+=proportion_of(self.genes, self.list_genes_fc, 'T', self.seq, self.seq_reverse)['minus']
            val_minus[2]+=proportion_of(self.genes, self.list_genes_fc, 'G', self.seq, self.seq_reverse)['minus']
            val_minus[3]+=proportion_of(self.genes, self.list_genes_fc, 'C', self.seq, self.seq_reverse)['minus']
            plt.figure(1)
            plt.subplot(211)
            plt.pie(val_plus, explode=explode, labels=name, autopct='%1.1f%%', startangle=90, shadow=True)
            plt.axis('equal')
            plt.subplot(212)
            plt.pie(val_minus, explode=explode, labels=name, autopct='%1.1f%%', startangle=90, shadow=True)
            plt.axis('equal')
            plt.show()



    def load_annotation(self):
        """ Laod a annotations file information where indice
                0 = file, 1 = separator ,2 = Locus column,3 = Strand column,
                4,5 Left Rigth column, 6 start line """
        if os.path.exists(basedir+"data/"+self.name+"/annotation/annotation.info"):
            with open(basedir+"data/"+self.name+"/annotation/annotation.info","r") as f:
                for line in f:
                    line=line.strip()
                    line=line.split('\t')
                    self.genes=annotations_parser_general(basedir+"data/"+self.name+'/annotation/'+line[0],line[1],int(line[2]),int(line[3]),int(line[4]),int(line[5]),int(line[6]))

                f.close()
            return True
        else:
            print(" Please load annotation of GFF file ")
            return False


    def load_annotation_gff(self):
        self.genes=annotations_parser_gff(basedir+"data/"+self.name+"/annotation/sequence.gff3")


    def load_neighbour(self):
        if not self.genes:
            if os.path.exists(basedir+"data/"+self.name+"/annotation/annotation.info"):
                self.load_annotation()
            else:
                self.load_annotation_gff()
        dict_plus={}
        dict_minus={}
        for i in self.genes:
            if self.genes[i].strand == True:
                dict_plus[int(self.genes[i].left)]=i
            else:
                dict_minus[int(self.genes[i].left)]=i
        l_plus=sorted(list(dict_plus.items()), key=operator.itemgetter(0))
        l_minus=sorted(list(dict_minus.items()), key=operator.itemgetter(0))
        self.genes=add_neighbour(self.genes,l_plus)
        self.genes=add_neighbour(self.genes,l_minus)


    def load_expression_level(self):
        """ Add expression level for all genes in dictionary """
        if not self.genes:
            if os.path.exists(basedir+"data/"+self.name+"/annotation/annotation.info"):
                self.load_annotation()
            else:
                self.load_annotation_gff()
        if os.path.exists(basedir+"data/"+self.name+"/expression/expression.info"):
            with open(basedir+"data/"+self.name+"/expression/expression.info","r") as f:
                for line in f:
                    line=line.strip()
                    line=line.split('\t')
                    self.genes=add_expression_to_genes(self.genes,basedir+"data/"+self.name+"/expression/"+line[0], int(line[1]), int(line[2]), line[3])
        else:
            print(" not found expression file information")

    def load_genes_positions(self):
        if not hasattr(self, 'genepos'):
            self.genepos = {}
        l=pd.read_csv(basedir+"data/"+self.name+"/genes_annotations.csv",header=0)
        gp=l[l.Strand=='forward']
        gp=gp[["Left.End.ASAP","Right.End.ASAP"]]
        gpq=np.array(gp)
        self.genepos["+"]=gpq[np.argsort(gpq[:,0])]
        gm=l[l.Strand=='complement']
        gm=gp[["Right.End.ASAP","Left.End.ASAP"]]
        gmq=np.array(gm)
        self.genepos["-"]=gmq[np.argsort(gmq[:,0])]


    def load_complete_TSS_list(self, *args, **kwargs):
        """ Load TSS list from file.
        Two filters can be applied:
        - filt: a minimum number of starts (default 20)
        - buffer: if two TSS are closer than this buffer size, we keep only the strongest one
        """
        filt=kwargs.get("filt")
        if not hasattr(self, 'TSS_complete'):
            self.TSS_complete = {}
        if not filt:
            filt=20
        self.TSS_complete["filter"]=filt
        plus=np.loadtxt(basedir+"data/"+self.name+"/TSS/TSS+.dat",skiprows=1)
        minus=np.loadtxt(basedir+"data/"+self.name+"/TSS/TSS-.dat",skiprows=1)
        print(plus)
        win=kwargs.get("buffer")
        self.TSS_complete["buffer"]=win
        pl,mi=filter_TSS((plus,minus),filt,win)
        self.TSS_complete["+"]=np.array(pl)
        self.TSS_complete["-"]=np.array(mi)

    def add_TSS_data(self,TSS_list,TSS_condition):
        """ Adds TSS from list generated from database as numpy array.
        TSS_condition describes the list, in the case where filters were applied etc."""
        if not hasattr(self, 'expression'):
            self.expression = {}
        self.TSS[TSS_condition]=np.array(TSS_list)

    def load_TSS(self, *args, **kwargs):
        """ Laod a TSS file information where indice 0 = condition, 1 = filename,
        2 = locus_tag, 3 = TSS_column, 4 = start_line, 5 = separator, 6 = strand column, 7 = Sig column
        if much other condition give it in the seconde line of file and change TSS column """
        if not (self.genes):
            if os.path.exists(basedir+"data/"+self.name+"/annotation/annotation.info"):
                self.load_annotation()
            else:
                self.load_annotation_gff()
        if os.path.exists(basedir+"data/"+self.name+"/TSS/TSS.info"):
            with open(basedir+"data/"+self.name+"/TSS/TSS.info","r") as f:
                for line in f:
                    line = line.strip('\n')
                    line = line.split('\t')
                    try:
                        test=int(line[2])# If no locus tag column in file information put None, test if locus tag column is number
                        try:
                            test=int(line[7])# Same test if no sig column
                            self.genes=load_single_TSS(self.genes, self.TSS_complete, basedir+"data/"+self.name+"/TSS/"+line[1], int(line[3]), int(line[4]), line[5], line[0], int(line[6]), self.TSS_plus, self.TSS_minus, tag_column = int(line[2]), sig=int(line[7]))
                        except:
                            self.genes=load_single_TSS(self.genes, self.TSS_complete, basedir+"data/"+self.name+"/TSS/"+line[1], int(line[3]), int(line[4]), line[5], line[0], int(line[6]), self.TSS_plus, self.TSS_minus,tag_column = int(line[2]))
                    except:
                        try:
                            test=int(line[7])
                            self.genes=load_single_TSS(self.genes, self.TSS_complete, basedir+"data/"+self.name+"/TSS/"+line[1], int(line[3]), int(line[4]), line[5], line[0], int(line[6]), self.TSS_plus, self.TSS_minus, sig=int(line[7]))
                        except:
                            self.genes=load_single_TSS(self.genes, self.TSS_complete, basedir+"data/"+self.name+"/TSS/"+line[1], int(line[3]), int(line[4]), line[5], line[0], int(line[6]), self.TSS_plus, self.TSS_minus)
        else:
            print(" no TSS.info maybe try use add_TSS_data()")


    def load_rpkm(self):
        """ Laod a RPKM file information where indice 0 = Condition
        1 = filename type, 2 = RPKM  column, 3 = Start line,
        4 = type of separator, 5=locus_column """
        if not (self.genes):
            if os.path.exists(basedir+"data/"+self.name+"/annotation/annotation.info"):
                self.load_annotation_gff()
            else:
                self.load_annotation()
        if os.path.exists(basedir+"data/"+self.name+"/rpkm/rpkm.info"):
            with open(basedir+"data/"+self.name+"/rpkm/rpkm.info","r") as f:
                for line in f:
                    line = line.strip('\n')
                    line = line.split('\t')
                    self.genes=add_single_rpkm_to_genes(self.genes, basedir+"data/"+self.name+"/rpkm/"+line[1],line[0],int(line[2]),int(line[3]),line[4],int(line[5]))
        else:
            print(" no rpkm file in this folder ")


    def load_genes_in_TU(self, TU):
        """ adds genes to TU accoding to annotation, ordered along orientation of the TU.
        Condition= gene entirely in TU
        """
        if not hasattr(self, 'genes'):
            print("loading annotation")
            load_annotation(self)
        glist=[]
        for g in list(self.genes.values()):
            if g.orientation==TU.orientation and g.left>=TU.left and g.right<=TU.right:
                glist.append(g.name)
        # sort list
        if len(glist)==0:
            TU.add_genes([])
        else:
            TU.add_genes(glist)
            lefts=[self.genes[x].left for x in glist]
            leftsort=sorted(lefts)
            if TU.orientation:
                TU.add_genes([glist[lefts.index(x)] for x in leftsort])
            else:
                TU.add_genes([glist[lefts.index(x)] for x in leftsort[::-1]])

    def add_single_TU(self, TU, index):
        """ adds TU to genome
        """
        if not hasattr(self, 'TU'):
            self.TU={}
        self.TU[index]=TU


    def load_genes_in_TUs(self):
        """ adds genes to all existing TUs.
        """
        if not hasattr(self, "TU"):
            print("no TU in genome")
        else:
            for TU in list(self.TU.values()):
                load_genes_in_TU(self, TU)



    def load_operon(self):
        if not self.genes:
            if os.path.exists(basedir+"data/"+self.name+"/annotation/annotation.info"):
                self.load_annotation()
            else:
                self.load_annotation_gff()
        self.operon_complete=add_operon(self.genes,basedir+"data/"+self.name+"/operon")


    def load_terminator(self):
        if not self.genes:
            if os.path.exists(basedir+"data/"+self.name+"/annotation/annotation.info"):
                self.load_annotation()
            else:
                self.load_annotation_gff()
        self.terminator_complete=add_terminator(basedir+"data/"+self.name+"/terminator")


    def get_cov_from_accession_files(self):
        if os.path.exists(basedir+"data/"+self.name+"/rnaseq_cov/description.info"):
            with open(basedir+"data/"+self.name+"/rnaseq_cov/description.info","r") as f:
                for line in f:
                    line = line.strip('\n')
                    line=line.split('\t')
                    if line[1] == '2':
                        list_cov=download_pair(basedir+"data/"+self.name+"/rnaseq_cov/"+line[0],self.name)
                    else:
                        list_cov=download_single(basedir+"data/"+self.name+"/rnaseq_cov/"+line[0],self.name)
                    create_cov_info(list_cov,self.name)
            f.close()
        else:
            print("Warning no description.info")


    def get_profile_from_matrix(self,factor_name):
        if not hasattr(self, 'profile'):
            self.profile={}
            self.load_seq()
        a = IUPAC.unambiguous_dna
        self.profile[factor_name]={}
        matrix=create_matrix_from_file(basedir+"data/pwm.txt",factor_name)
        #matrix=create_matrix_from_file_2(basedir+"data/pwm.txt",factor_name)
        seq=Seq(self.seq,a)
        for o in matrix.search_pwm(seq):
            self.profile[factor_name][o[0]]=o[1]

    def load_SIST(self, start, end,*args, **kwargs):
        if not hasattr(self, 'SIST_profile'):
            self.SIST_profile={}
            self.load_seq()
        option = kwargs.get('option')
        if option:
            if option == 'A':
                self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq, option=option)
            elif option == 'Z':
                self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq, option=option)
            elif option == 'C':
                self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq, option=option)
            else:
                print("This option doesn't exist")
        else:
            self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq)

    def load_dom(self):
        i=1
        if not hasattr(self, 'dom'):
            self.dom={}
            if not self.load_annotation():
                self.load_annotation_gff()
            self.load_expression_level()
        with open(basedir+"data/"+self.name+"/domains.txt","r") as f:
            header=next(f)
            for line in f:
                line=line.strip()
                line=line.split('\t')
                self.dom[i]=Domain(start=int(line[1]),end=int(line[2]),index=i)
                self.dom[i].add_genes(self.genes)
                self.dom[i].add_domain_expression()
                self.dom[i].add_list_expression()
                i+=1

###################### RAPH #############################

    def load_reads(self): # new attribute reads : reads_pos & reads_neg, of shape {[condition] : .npy}, e.g. self.reads_pos[cond1]
        self.reads_pos = {} # reads on + strand
        self.reads_neg = {} # reads on - strand
        if not os.path.exists(basedir+"data/"+self.name+'/rnaseq_reads/reads.info'):
            print 'Unable to locate reads.info in /rnaseq_reads/'
        else:
            # open info file
            with open(basedir+"data/"+self.name+'/rnaseq_reads/reads.info',"r") as f:
                header = next(f) # first line = header
                # one line / condition
                for line in f:
                    line=line.strip()
                    line=line.split('\t')
                    print 'Loading condition',line[0]
                    self.reads_pos[line[0]] = np.load(basedir+"data/"+self.name+'/rnaseq_reads/'+line[1])["Rpos"]
                    self.reads_neg[line[0]] = np.load(basedir+"data/"+self.name+'/rnaseq_reads/'+line[1])["Rneg"]
            print 'Done'


    def load_cov(self): # new attribute cov : cov_pos & cov_neg, of shape {[condition] : .npy}, e.g. self.cov_pos[cond1]
        self.cov_pos = {} # cov on + strand
        self.cov_neg = {} # cov on - strand
        # function tries first to deal with cov_info and .npy files directly, if cov_info not available then
        # tries to open cov_txt.info, convert .txt files into .npy, create cov.info and load them
        if os.path.exists(basedir+"data/"+self.name+'/rnaseq_cov/cov.info'): # cov.info available, cov.info opening instead of cov_txt.info
            with open(basedir+"data/"+self.name+"/rnaseq_cov/cov.info","r") as f:
                header = next(f)       
                for line in f: # for each condition
                    line=line.strip()
                    line=line.split('\t')
                    print 'Loading condition',line[0]
                    # load attributes
                    self.cov_neg[line[0]]= np.load(basedir+"data/"+self.name+'/rnaseq_cov/'+line[1])["cov_neg"]
                    self.cov_pos[line[0]]= np.load(basedir+"data/"+self.name+'/rnaseq_cov/'+line[1])["cov_pos"]
        
        if not os.path.exists(basedir+"data/"+self.name+'/rnaseq_cov/cov.info') and os.path.exists(basedir+"data/"+self.name+'/rnaseq_cov/cov_txt.info'):
            print 'Unable to locate cov.info in /rnaseq_cov/'
            print 'Working with .txt file (cov_txt.info)'
            file = open(basedir+"data/"+self.name+'/rnaseq_cov/cov.info','w') 
            file.write('Condition\tCov file\tDate\tReads file')
            file.close() 
            with open(basedir+"data/"+self.name+"/rnaseq_cov/cov_txt.info","r") as f: # load cov_txt.info
                for line in f: # for each condition (each .txt file in cov_txt.info)
                    line = line.strip('\n') 
                    line = line.split('\t') 
                    print 'Loading condition:',line[0]
                    # create .npy from .txt
                    cov_neg=np.loadtxt(basedir+"data/"+self.name+"/rnaseq_cov/"+line[1],usecols=[int(line[2])]) 
                    cov_pos=np.loadtxt(basedir+"data/"+self.name+"/rnaseq_cov/"+line[3],usecols=[int(line[4])]) 
                    # load attributes
                    self.cov_neg[line[0]]= cov_neg
                    self.cov_pos[line[0]]= cov_pos
                    # save .npy into .npz
                    np.savez(basedir+"data/"+self.name+"/rnaseq_cov/"+line[0]+'_cov.npz', cov_pos=cov_pos, cov_neg=cov_neg)
                    # update cov.info
                    file=open(basedir+"data/"+self.name+"/rnaseq_cov/cov.info","a")                   
                    file.write('\n'+line[0]+'\t'+line[0]+'_cov.npz\t'+str(datetime.now())+'\tUnknown')
                    file.close()
            f.close()
        if not os.path.exists(basedir+"data/"+self.name+'/rnaseq_cov/cov.info') and not os.path.exists(basedir+"data/"+self.name+'/rnaseq_cov/cov_txt.info'):
            print 'cov.info not available nor cov_txt.info, please check /rnaseq_cov/ folder'
    
        print 'Done'
        

    def compute_rpkm_from_cov(self, before=100):
        """Adds rpkm values from coverage: along whole genes Before= number of bps to add before = to take into account 
        DNA region upstream of the coding sequence of the gene """
        if not self.genes: # if no genes loaded
            # try to load them
            if os.path.exists(basedir+"data/"+self.name+"/annotation/annotation.info"):
                self.load_annotation()
            else:
                self.load_annotation_gff()
        try:
            for g in list(self.genes.keys()): # for each gene
                if hasattr(self.genes[g],'orientation'): # which has a known orientation
                    if self.genes[g].orientation==1: 
            # gene in + strand
                        for cond in list(self.cov_pos.keys()): # for each condition of cov
                            self.genes[g].add_single_rpkm(cond, np.mean(self.cov_pos[cond][(self.genes[g].left-100):self.genes[g].right]))
                    else:
            # gene in - strand
                        for cond in list(self.cov_neg.keys()):
                            self.genes[g].add_single_rpkm(cond, np.mean(self.cov_neg[cond][self.genes[g].left:(self.genes[g].right+100)]))
        except:
            print("You need to load coverage pls")


    def load_fc(self,*args, **kwargs):
        """ Fc info file, 0=condition 1=file_name, 2=tag_column, 3=fc_column
                4=type of separator, 5 = start line, 6 = p_value( if no write nothing),
                if other source give the line of the source in fc information file """
        self.list_genes_fc = {} # for each condition, list of genes having a FC
        self.list_genes_fc_pval = {} # for each condition, list of genes having a FC and a pvalue
        if not self.genes: # if no genes loaded
            # try to load them
            if os.path.exists(basedir+"data/"+self.name+"/annotation/annotation.info"):
                self.load_annotation()
            else:
                self.load_annotation_gff()
        n=0 
        if os.path.exists(basedir+"data/"+self.name+"/fold_changes/fc.info"): 
            with open(basedir+"data/"+self.name+"/fold_changes/fc.info","r") as f:
                skiphead = next(f) # skip head
                for header in f:
                    header=header.strip()
                    header=header.split('\t')
                    try: # if p-value in file, works
                        self.list_genes_fc_pval[header[0]]=add_single_fc_to_genes(self.genes,basedir+"data/"+self.name+"/fold_changes/"+header[1],header[0],int(header[2]),int(header[3]),header[4],int(header[5]),n,p_value=int(header[6]))
                        self.list_genes_fc[header[0]] = self.list_genes_fc_pval[header[0]]
                    except: # without p-value otherwise
                        self.list_genes_fc[header[0]]=add_single_fc_to_genes(self.genes,basedir+"data/"+self.name+"/fold_changes/"+header[1],header[0], int(header[2]),int(header[3]),header[4],int(header[5]),n)

                    n+=1
            f.close()
        else:
            print("No fc.info file, please create one")


    def compute_melting_energy(self,windows=500000, increment=4000):
        # compute melting energy on genome windows with a specific increment
        self.melting_energy = []
        if not hasattr(self, 'seq'): # if no seq loaded
            try:
                print 'Trying to load seq...'
                self.load_seq()
                print 'seq loaded'
            except:
                print'Unable to load seq'
                sys.exit()

        bins = [] # bins = windows of the genome : [start coordinate,end coordinate,melting energy]
        bins_overlap = [] # bins where end coordinate < start coordinate (overlap circular chromosome)

        for i in range(1,self.length,increment): # create bins depending on windows size and increment value
            if (i+windows) <= self.length: # enough length to create a bin
                bins.append([i,i+windows,0])
            else: # i + windows > genome size, overlap with the beginning (circular chromosome)
                bins_overlap.append([i, windows - (self.length-i),0])
        
        bins = np.array(bins) # convert to .npy
        bins_overlap = np.array(bins_overlap)

        # compute melting energy on bins
        for start,end,melting_energy in bins:
            seq = Seq(self.seq[start-1:end])
            melting_energy = mt.Tm_NN(seq)
        
        for start,end,melting_energy in bins_overlap:
            seq = Seq(self.seq[start-1:] + self.seq[0:end])
            melting_energy = mt.Tm_NN(seq)
        
        bins = np.concatenate((bins,bins_overlap))

        self.melting_energy = list(bins[:,2])

    def draw_melting_energy_circle(self, *args, **kwargs):
        # draw melting energy circle from melting energy
        # opt arg : colormap, vmin, vmax
        if not hasattr(self, 'melting_energy'): # if melting energy not computed
            try:
                print 'Trying to compute melting energy with default values...'   
                self.compute_melting_energy()
                print 'Melting energy computed'
            except:
                print'Unable to compute melting_energy'
                sys.exit()

        colormap= kwargs.get('colormap','seismic') # default value seismic
        try:
            cScale_fc = plt.get_cmap(colormap)
        except:
            print 'Incorrect colormap, please check https://matplotlib.org/users/colormaps.html'
            print 'Loading default seismic'
            cScale_fc = plt.get_cmap('seismic')

        vmin = kwargs.get('vmin', min(self.melting_energy))
        vmax = kwargs.get('vmax', max(self.melting_energy))          
        # normalisation of colors
        cNorm_fc  = colors.Normalize(vmin=vmin, vmax=vmax) 
        # map which assigns a colour depending on value between vmin and vmax
        cMap_fc = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_fc) 
        # config, see globvar for more
        # init plots
        fig, ax = plt.subplots(1,1) ; fig.set_size_inches(fig_width, fig_height)
        plt.axis([0, fig_width, 0, fig_height]) ; ax.set_axis_off()

        angle = 360.0/len(self.melting_energy) # angle between two fragments
        i=0
        for value in self.melting_energy:
            # edgecolor = assign a colour depending on value using cMap
            # draw arc on circle
            arc = patches.Arc((center_x,center_y), radius, radius, angle=0,theta1=i, theta2=i+angle, edgecolor=cMap_fc.to_rgba(value), lw=10)
            ax.add_patch(arc)
            i+= angle

        cMap_fc._A = [] # fake array to print map
        plt.colorbar(cMap_fc).set_label("Melting energy")
        fig.suptitle('Melting energy '+self.name, fontweight='bold') #, fontsize=14, fontweight='bold')
        fig.savefig(basedir+"data/"+self.name+"/annotation/melting_energy.png", format='png', dpi=400, transparent=False) # png (72,300,400 dpi) or svg       

    def draw_expression_circles(self, *arg, **kwargs):
        # generate density circles based on FC and pvalues
        # opt arguments : colormap, vmin vmax (color normalisation), windows, increment
        windows= kwargs.get('windows', 500000)
        increment= kwargs.get('increment', 4000)

        path = basedir+"data/"+self.name+"/fold_changes/circles-"+str(datetime.now())
        os.makedirs(path)

        if not hasattr(self, 'list_genes_fc_pval'): # if no fc loaded 
            try:
                print 'Trying to load FC...'
                self.load_fc()
                print 'FC loaded'
            except:
                print 'Unable to load fc'
                sys.exit()

        if self.list_genes_fc_pval.keys() == []:
            print 'No condition where genes have valids FC and p-value, unable to continue'
            sys.exit()
        else:
            for cond in self.list_genes_fc_pval.keys():
                print 'Computing condition',cond
                gen_states = self.compute_table_genes(cond)
                bins = self.count_genes_in_windows(cond, gen_states, windows, increment)
                print 'Windows on genome :\n',bins
                tot_act = len(gen_states[np.where(gen_states[:,1] == 1)])
                print 'Total activated genes on genome :',tot_act
                tot_repr = len(gen_states[np.where(gen_states[:,1] == -1)])
                print 'Total repressed genes on genome :',tot_repr
                zscores = self.compute_zscores(tot_act,tot_repr,bins)
                # Colormap for fold change
                colormap= kwargs.get('colormap','seismic') # default value seismic
                try:
                    cScale_fc = plt.get_cmap(colormap)
                except:
                    print 'Incorrect colormap, please check https://matplotlib.org/users/colormaps.html'
                    print 'Loading default seismic'
                    cScale_fc = plt.get_cmap('seismic')
                # default values = smallest zscore, highest
                vmin = kwargs.get('vmin', min(zscores))
                vmax = kwargs.get('vmax', max(zscores))
                # normalisation of colours
                cNorm_fc  = colors.Normalize(vmin=vmin, vmax=vmax) 
                # map which assigns a colour depending on value between vmin and vmax
                cMap_fc = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_fc) 
                # init plots
                fig, ax = plt.subplots(1,1) ; fig.set_size_inches(fig_width, fig_height)
                plt.axis([0, fig_width, 0, fig_height]) ; ax.set_axis_off()
                
                angle = 360.0/len(zscores) # angle between two fragments
                i=0
                for value in zscores:
                    # edgecolor = assign a colour depending on value using cMap
                    # draw arc on circle
                    arc = patches.Arc((center_x,center_y), radius, radius, angle=0,theta1=i, theta2=i+angle, edgecolor=cMap_fc.to_rgba(value), lw=5)
                    ax.add_patch(arc)
                    i+= angle

                cMap_fc._A = [] # fake array to print map
                plt.colorbar(cMap_fc).set_label("Z-score")
                fig.suptitle(cond, fontweight='bold') #, fontsize=14, fontweight='bold')
                fig.savefig(path+"/circle-"+cond+".png", format='png', dpi=400, transparent=False) # png (72,300,400 dpi) or svg

    def compute_table_genes(self, cond): 
        # returns a npy where a row is a gene caracterised by a start pos and a gene state
        # gene is considered activated above a given fc, repressed below a given fc
        gen_states = []
        for gene in self.list_genes_fc_pval[cond]:
            # if activated
            if self.genes[gene].all_fc[cond] >= fc_treshold_pos and self.genes[gene].all_pval[cond] <= pval_treshold:
                gen_states.append([self.genes[gene].start,1])
            # if repressed
            elif self.genes[gene].all_fc[cond] <= fc_treshold_neg and self.genes[gene].all_pval[cond] <= pval_treshold:
                gen_states.append([self.genes[gene].start,-1])
            # if not affected
            else:
                gen_states.append([self.genes[gene].start,0])
        
        gen_states = np.array(gen_states)
        return gen_states

    def count_genes_in_windows(self, cond, gen_states, windows, increment):
        # compute bins on the genome depending on windows size and increment, and
        # calculate the nb of activated / repressed / non affected genes in each bin for further zscore
        if not hasattr(self, 'seq'): # if no seq loaded
            try:
                print 'Trying to load seq...'
                self.load_seq()
                print 'seq loaded'
            except:
                print'Unable to load seq'
                sys.exit()

        bins = [] # bins = windows of the genome : [start coordinate,end coordinate,nb of activated genes,nb of repressed genes,nb of genes not affected]
        bins_overlap = [] # bins where end coordinate < start coordinate (overlap circular chromosome)
        for i in range(1,self.length,increment): # create bins depending on windows size and increment value
            if (i+windows) <= self.length: # enough length to create a bin
                bins.append([i,i+windows,0,0,0])
            else: # i + windows > genome size, overlap with the beginning (circular chromosome)
                bins_overlap.append([i, windows - (self.length-i),0,0,0])
        bins_overlap = np.array(bins_overlap)
        bins = np.array(bins) # convert to .npy
        
        for start,state in gen_states: # reminder gene state : a row = beginning of gene, state (activated, repressed or not affected)
            if state == 1: # activated
            # test to which bins the gene belongs to, and add one to the nb of activated genes of these bins
                bins[np.where((bins[:,0] < start) & (bins[:,1] > start)),2] += 1
                bins_overlap[np.where((bins_overlap[:,0] < start) | (bins_overlap[:,1] > start)),2] += 1
            elif state == -1: # repressed gene
                bins[np.where((bins[:,0] < start) & (bins[:,1] > start)),3] += 1
                bins_overlap[np.where((bins_overlap[:,0] < start) | (bins_overlap[:,1] > start)),3] += 1
            elif state == 0: # not affected gene
                bins[np.where((bins[:,0] < start) & (bins[:,1] > start)),4] += 1
                bins_overlap[np.where((bins_overlap[:,0] < start) | (bins_overlap[:,1] > start)),4] += 1
        
        bins = np.concatenate((bins,bins_overlap))
        return bins

    def compute_zscores(self, tot_act, tot_repr, bins):
        # p_exp = nb of total activated genes / nb of total activated + repressed genes on the whole genome
        zscores = []
        p_exp = float(tot_act) / float((tot_act + tot_repr))
        print 'g+ / (g+ + g-) on genome :',p_exp
        # compute zscore for each bin
        for start,end,nb_act,nb_repr,nb_null in bins:
            nb_act = float(nb_act) ; nb_repr = float(nb_repr)
            try:
                zscore =(nb_act - (nb_act+nb_repr)*p_exp) / (sqrt((nb_act+nb_repr)*p_exp*(1-p_exp)))
            except: # division by zero if no genes activated nor repressed
                zscore = 0
            zscores.append(zscore)
        return zscores
