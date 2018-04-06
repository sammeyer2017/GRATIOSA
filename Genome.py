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
# useful functions

def annotations_parser_general(annotations_filename,separator,tag_column,strand_column,left_column,right_column,start_line):
    ''' Called by load annotation, allows genes to be loaded from info file
    '''
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
            gene=[]
            line=line.strip()
            if separator == '\\t':
                line=line.split('\t')
            else:
                line=line.split(separator)

            gene.append(line[tag_column])
            if line[strand_column]=="complement":
                gene.append('-')
            elif line[strand_column]=="forward":
                gene.append('+')
            elif line[strand_column]== "1":
                gene.append('+')
            elif line[strand_column]== "-1":
                gene.append('-')
            else:
                gene.append(line[strand_column])

            gene.append(line[left_column])
            gene.append(line[right_column])
            gene.append(line)
            if line[tag_column] in genes_dict:
                print("Warning! Overwriting value for gene ")
                print(line[tag_column])
            genes_dict[line[tag_column]]=Gene(annotations_general=gene, head=headl)
    return genes_dict

def annotations_parser_gff(annotations_filename):
    ''' Called by load annotation, allows genes to be loaded from gff
    '''
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

def load_fc_pval_cond(genes_dict, filename, condition, tag_col, fc_col, separator, start_line, *args, **kwargs):
    ''' Called by load_fc_pval, allows expression data to be loaded by specifying files, and where each
    information is (tag, fc, pval...). If no p-value column, assigns pval = 0 to each gene
    '''
    genes_valid = [] # list containing all genes having valid FC / pval
    p_val_col= kwargs.get('p_value')

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
                    genes_dict[line[tag_col]].add_fc_pval_cond(float(line[fc_col]),condition, float(line[p_val_col]))
                else:
                    genes_dict[line[tag_col]].add_fc_pval_cond(float(line[fc_col]),condition, float(0))
                genes_valid.append(line[tag_col])
            except:
                if line[tag_col] not in genes_dict.keys():
                    if line[tag_col] != '':
                        print(line[tag_col] + " not in annotation ")
                    else:
                        print("fc without locus")
    f.close()
    return genes_valid



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

def load_TSS_cond(genes_dict,filename, TSS_column, start_line , separator, strand, *args, **kwargs):
    TSS_dict = {} # dict of TSS objects
    sig = kwargs.get('sig')
    tag_column = kwargs.get('tag_column')
    sites = kwargs.get('sites')
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
                pos = int(line[TSS_column])
                # if TSS need to be init
                if pos not in TSS_dict.keys():
                    # init object tss
                    TSS_dict[pos] = TSS(pos = pos)
                    TSS_dict[pos].add_strand(line[strand])
                    if tag_column: # if tag column
                        TSS_dict[pos].add_genes(line[tag_column],genes_dict)
                        for gene in TSS_dict[pos].genes: # add TSS to gene attributes
                            genes_dict[gene].add_id_TSS(pos)
                # Add sigma factor and binding sites to the promoter dict
                if sig: # if sigma column
                    if line[sig] != '':
                        if sites:
                            if line[sites] != '':
                                TSS_dict[pos].add_promoter(line[sig], sites = line[sites])
                        else:
                            TSS_dict[pos].add_promoter(line[sig])
            
            except Exception as e:
                print 'Error in line, wrong information type :',e

    return TSS_dict


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
        self.seqcompl=''
        if(self.length):
            if(self.length != len(self.seq)):
                print("Warning not the same length, the new sequence will be removed")
                self.seq=''
        self.length=len(self.seq)
        l=self.seq
        l=l.replace('A','t')
        l=l.replace('T','a')
        l=l.replace('C','g')
        l=l.replace('G','c')
        l=l.replace('a','A')
        l=l.replace('t','T')
        l=l.replace('c','C')
        l=l.replace('g','G')
        self.seqcompl=l

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
        """ Load annotation. Two options : if gff file in directory -> load annotation from gff
        if no gff file in directory -> tries to load annotation.info (0 = file, 1 = separator ,2 = 
        Locus column,3 = Strand column, 4,5 Left Rigth column, 6 start line)
        """
        if os.path.exists(basedir+"data/"+self.name+"/annotation/sequence.gff3"):
            self.genes=annotations_parser_gff(basedir+"data/"+self.name+"/annotation/sequence.gff3")
        elif os.path.exists(basedir+"data/"+self.name+"/annotation/annotation.info"):
            with open(basedir+"data/"+self.name+"/annotation/annotation.info","r") as f:
                for line in f:
                    line=line.strip()
                    line=line.split('\t')
                    self.genes=annotations_parser_general(basedir+"data/"+self.name+'/annotation/'+line[0],line[1],int(line[2]),int(line[3]),int(line[4]),int(line[5]),int(line[6]))
                f.close()
        else:
            print('No GFF file nor annotation.info, unable to load annotation')


    def load_neighbour(self):
        if not self.genes:
            self.load_annotation()

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


    def load_TSS(self, *args, **kwargs):
        """ Load a TSS file info where indice 0 = condition, 1 = filename,
        2 = locus_tag, 3 = TSS_column, 4 = start_line, 5 = separator, 6 = strand column, 7 = Sig column
        8 = Sites column if much other condition give it in the seconde line of file and change TSS column """
        self.TSSs = {} # shape (dict of dict) : TSSs = {TSScond : {TSS:attr}}
        self.TSSs['all_TSS'] = {} # dict containing all TSS and where they appear (shape all_TSS = {pos:[conditions]})
        if not (self.genes):
            self.load_annotation()

        if os.path.exists(basedir+"data/"+self.name+"/TSS/TSS.info"):
            with open(basedir+"data/"+self.name+"/TSS/TSS.info","r") as f:
                skiphead = next(f) # skip head
                for line in f:
                    line = line.strip('\n')
                    line = line.split('\t')
                    try: # successively try to load :
                        try: # sites + tag non tested because sites without sig not possible
                            try: # tag + sig + sites
                                self.TSSs[line[0]]=load_TSS_cond(self.genes, basedir+"data/"+self.name+"/TSS/"+line[1], int(line[3]), int(line[4]), line[5], int(line[6]), tag_column = int(line[2]), sig=int(line[7]), sites =int(line[8]))
                            except: # tag + sig
                                self.TSSs[line[0]]=load_TSS_cond(self.genes, basedir+"data/"+self.name+"/TSS/"+line[1], int(line[3]), int(line[4]), line[5], int(line[6]), tag_column = int(line[2]), sig=int(line[7]))
                        except:
                            try: # sig + sites
                                self.TSSs[line[0]]=load_TSS_cond(self.genes, basedir+"data/"+self.name+"/TSS/"+line[1], int(line[3]), int(line[4]), line[5], int(line[6]), sig=int(line[7]), sites =int(line[8]))
                            except: # tag
                                self.TSSs[line[0]]=load_TSS_cond(self.genes, basedir+"data/"+self.name+"/TSS/"+line[1], int(line[3]), int(line[4]), line[5], int(line[6]), tag_column = int(line[2]))
                    except: # if no tag, no sig, no sites
                            self.TSSs[line[0]]=load_TSS_cond(self.genes, basedir+"data/"+self.name+"/TSS/"+line[1], int(line[3]), int(line[4]), line[5], int(line[6]))
                    # append all entries to all TSS dict
                    for entry in self.TSSs[line[0]].keys():
                        try: # works if entry already in dict
                            self.TSSs['all_TSS'][entry].append(line[0])
                        except: # init list of conditions for entry
                            self.TSSs['all_TSS'][entry] = []
                            self.TSSs['all_TSS'][entry].append(line[0])
        else:       
            print("No TSS.info, unable to load TSS")


    def load_rpkm(self):
        """ Laod a RPKM file information where indice 0 = Condition
        1 = filename type, 2 = RPKM  column, 3 = Start line,
        4 = type of separator, 5=locus_column """
        if not (self.genes):
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
            self.load_annotation()

        self.operon_complete=add_operon(self.genes,basedir+"data/"+self.name+"/operon")


    def load_terminator(self):
        if not self.genes:
            self.load_annotation()

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

    def load_reads(self): 
        '''
        new attribute reads : reads_pos & reads_neg, of shape {[condition] : .npy}, e.g. self.reads_pos[cond1]
        '''
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

    def load_cov(self): 
        '''new attribute cov : cov_pos & cov_neg, of shape {[condition] : .npy}, e.g. self.cov_pos[cond1]
        '''
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
        '''
        Adds rpkm values from coverage: along whole genes Before= number of bps to add before = to take into account 
        DNA region upstream of the coding sequence of the gene
        '''
        if not self.genes: # if no genes loaded
            # try to load them
            self.load_annotation()
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


    def load_fc_pval(self,*args, **kwargs):
        ''' Load fc and pval specified in fc.info. Fc info file columns (tab-separated) : condition, filename, tag column
        , fc column, file separator, startline, pvalcolumn (if not available leave empty)
        '''
        self.genes_valid = {} # for each condition, list of genes having a FC /pval value
        if not self.genes: # if no genes loaded
            # try to load them
            self.load_annotation()

        if os.path.exists(basedir+"data/"+self.name+"/fold_changes/fc.info"): 
            with open(basedir+"data/"+self.name+"/fold_changes/fc.info","r") as f:
                skiphead = next(f) # skip head
                for header in f:
                    header=header.strip()
                    header=header.split('\t')
                    try: # if p-value column specified works
                        self.genes_valid[header[0]]=load_fc_pval_cond(self.genes,basedir+"data/"+self.name+"/fold_changes/"+header[1],header[0],int(header[2]),int(header[3]),header[4],int(header[5]), p_value=int(header[6]))
                    except: # otherwise, set pvalue = 0
                        self.genes_valid[header[0]]=load_fc_pval_cond(self.genes,basedir+"data/"+self.name+"/fold_changes/"+header[1],header[0],int(header[2]),int(header[3]),header[4],int(header[5]))
            f.close()
        else:
            print("No fc.info file, please create one")

    def compute_melting_energy(self,windows=500000, increment=4000):
        ''' 
        Compute melting energy on genome windows with a specific increment
        '''
        self.melting_energy = []
        if not hasattr(self, 'seq'): # if no seq loaded
            try:
                print 'Trying to load seq...'
                self.load_seq()
                print 'seq loaded'
            except:
                print'Unable to load seq'
                sys.exit()

        bins = [] # bins = windows of the genome : [start coordinate,end coordinate]
        bins_overlap = [] # bins where end coordinate < start coordinate (overlap circular chromosome)

        for i in range(1,self.length,increment): # create bins depending on windows size and increment value
            if (i+windows) <= self.length: # enough length to create a bin
                bins.append([i,i+windows])
            else: # i + windows > genome size, overlap with the beginning (circular chromosome)
                bins_overlap.append([i, windows - (self.length-i)])
        
        bins = np.array(bins) # convert to .npy
        bins_overlap = np.array(bins_overlap)

        # compute melting energy on bins
        for start,end in bins:
            seq = Seq(self.seq[start-1:end])
            self.melting_energy.append(mt.Tm_GC(seq, strict=False))
            # if len(seq) < 1000:
            #     melting_energy = mt.Tm_NN(seq,strict=False)
            # else:
            #     melting_energy = mt.Tm_GC(seq,strict=False)
    
        for start,end in bins_overlap:
            seq = Seq(self.seq[start-1:] + self.seq[0:end])
            self.melting_energy.append(mt.Tm_GC(seq, strict=False))
        
        self.melting_energy_wind = windows
        self.melting_energy_incr = increment

    def draw_melting_energy_circle(self, *args, **kwargs):
        '''
        Draw melting energy circle from melting energy
        opt arg : colormap, vmin, vmax
        '''
        if not hasattr(self, 'melting_energy'): # if melting energy not computed
            try:
                print 'Trying to compute melting energy with default values...'   
                self.compute_melting_energy()
                print 'Melting energy computed'
            except:
                print'Unable to compute melting_energy'
                sys.exit()

        colormap= kwargs.get('colormap','jet') # default value
        try:
            cScale_fc = plt.get_cmap(colormap)
        except:
            print 'Incorrect colormap, please check https://matplotlib.org/users/colormaps.html'
            print 'Loading default'
            cScale_fc = plt.get_cmap('jet')

        vmin = kwargs.get('vmin', min(self.melting_energy))
        vmax = kwargs.get('vmax', max(self.melting_energy))          
        # normalisation of colors
        cNorm_fc  = colors.Normalize(vmin=vmin, vmax=vmax) 
        # map which assigns a colour depending on value between vmin and vmax
        cMap_fc = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_fc) 
        # config, see globvar for more
        # init plots
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.subplot(1,1,1)
        plt.axis([0, fig_width, 0, fig_height]) ; ax.set_axis_off() 

        angle = 360.0/len(self.melting_energy) # angle between two fragments
        # display colour in the middle of the windows
        start_angle = angle * (self.melting_energy_wind/2 - self.melting_energy_incr/2)  / self.melting_energy_incr

        i=0
        for value in self.melting_energy:
            # edgecolor = assign a colour depending on value using cMap
            # draw arc on circle
            arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_fc.to_rgba(value), lw=7)
            ax.add_patch(arc)
            i+= angle

        cMap_fc._A = [] # fake array to print map
        cbar = plt.colorbar(cMap_fc,shrink=0.3)
        cbar.set_label("Melting temp °C",size=9)
        cbar.ax.tick_params(labelsize=7)
        plt.axis("equal")
        plt.annotate('OriC', xy=(center_x,radius), xytext=(center_x,radius+0.4),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=8)
        plt.annotate('Ter', xy=(center_x,0), xytext=(center_x,0),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=8)
        plt.annotate('Melting energy '+self.name, xy=(center_x,radius), xytext=(center_x,center_y),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=7)
        plt.savefig(basedir+"data/"+self.name+"/annotation/melting_energy.pdf", transparent=False)        
        plt.close()

    def draw_expression_circles(self, *arg, **kwargs):
        '''
        generate density circles based on FC and pvalues
        opt arguments : colormap = see matplotlib doc, vmin vmax (color normalisation), windows, increment
        meth (= act/repr/actrepr/overlay), format (= png/pdf/svg)
        '''
        windows= kwargs.get('windows', 500000)
        increment= kwargs.get('increment', 4000)

        path = basedir+"data/"+self.name+"/fold_changes/circles-"+str(datetime.now())
        os.makedirs(path)

        if not hasattr(self, 'genes_valid'): # if no fc loaded 
            try:
                print 'Trying to load FC...'
                self.load_fc_pval()
                print 'FC loaded'
            except:
                print 'Unable to load fc'
                sys.exit()

        for cond in self.genes_valid.keys():
            print 'Computing condition',cond
            gen_states = self.compute_state_genes(cond)
            bins = self.count_genes_in_windows(cond, gen_states, windows, increment)
            print 'Windows on genome :\n',bins
            try:
                tot_act = len(gen_states[np.where(gen_states[:,1] == 1)])
                print 'Total activated genes on genome :',tot_act
                tot_repr = len(gen_states[np.where(gen_states[:,1] == -1)])
                print 'Total repressed genes on genome :',tot_repr
                tot_non = len(gen_states[np.where(gen_states[:,1] == 0)])
                print 'Total non affected genes on genome :',tot_non

                meth = kwargs.get('meth', 'actrepr')

                if meth == "act":
                    zscores = self.compute_zscores_act(tot_act,tot_repr,tot_non,bins)
                elif meth == "repr":
                    zscores = self.compute_zscores_repr(tot_act,tot_repr,tot_non,bins)
                elif meth == "actrepr":
                    zscores = self.compute_zscores_actrepr(tot_act,tot_repr,tot_non,bins)
                elif meth == "overlay":
                    zscores_act = self.compute_zscores_act(tot_act,tot_repr,tot_non,bins)
                    zscores_repr = self.compute_zscores_repr(tot_act,tot_repr,tot_non,bins)
                else:
                    print "Unknown method, computing default..."
                    zscores = self.compute_zscores_actrepr(tot_act,tot_repr,tot_non,bins)

            except:
                print 'Invalid data (e.g. no genes repressed nor activated)'
                sys.exit()

            # init plots
            fig = plt.figure(figsize=(fig_width,fig_height))
            ax = plt.subplot(1,1,1)
            plt.axis([0, fig_width, 0, fig_height]) ; ax.set_axis_off() 
           
            vmin = kwargs.get('vmin', -4)
            vmax = kwargs.get('vmax', +4)

            if meth == "overlay":
                # normalisation of colours
                try:
                    cNorm_fc = colors.Normalize(vmin=vmin, vmax=vmax)
                except:
                    print 'Incorrect values for vmin / vmax'
                    print 'Loading default'
                    cNorm_fc = colors.Normalize(vmin=-4, vmax=4)

                cScale_act = plt.get_cmap("Reds")
                cScale_repr = plt.get_cmap("Blues")
                cMap_act = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_act) 
                cMap_repr = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_repr) 
                angle = 360.0/len(zscores_act) # angle between two fragments
                # display colour in the middle of the windows
                start_angle = angle * (windows/2 - increment/2)  / increment 
                i=0
                for value in zscores_act:
                    # edgecolor = assign a colour depending on value using cMap
                    # draw arc on circle
                    arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_act.to_rgba(value), lw=7, alpha = 0.7)
                    ax.add_patch(arc)
                    i+= angle
                angle = 360.0/len(zscores_repr) 
                start_angle = angle * (windows/2 - increment/2)  / increment
                i=0
                for value in zscores_repr:
                    # edgecolor = assign a colour depending on value using cMap
                    # draw arc on circle
                    arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_repr.to_rgba(value), lw=7, alpha = 0.7)
                    ax.add_patch(arc)
                    i+= angle

                cMap_act._A = [] # fake array to print map
                cMap_repr._A = [] # fake array to print map
                cbar = plt.colorbar(cMap_act,shrink=0.3)
                cbar.set_label("Z-score act",size=9)
                cbar.ax.tick_params(labelsize=7)
                
                cbar = plt.colorbar(cMap_repr,shrink=0.3)
                cbar.set_label("Z-score repr",size=9)
                cbar.ax.tick_params(labelsize=7)

            else:
                # Colormap for Zscore
                colormap= kwargs.get('colormap','jet') # default value
                try:
                    cScale_fc = plt.get_cmap(colormap)
                except:
                    print 'Incorrect colormap, please check https://matplotlib.org/users/colormaps.html'
                    print 'Loading default'
                    cScale_fc = plt.get_cmap('jet')
                # normalisation of colours
                try:
                    cNorm_fc = colors.Normalize(vmin=vmin, vmax=vmax)
                except:
                    print 'Incorrect values for vmin / vmax'
                    print 'Loading default'
                    cNorm_fc = colors.Normalize(vmin=-4, vmax=4)

                # map which assigns a colour depending on value between vmin and vmax
                cMap_fc = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_fc) 
                angle = 360.0/len(zscores) # angle between two fragments
                # display colour in the middle of the windows
                start_angle = angle * (windows/2 - increment/2)  / increment
                i=0
                for value in zscores:
                    # edgecolor = assign a colour depending on value using cMap
                    # draw arc on circle
                    arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_fc.to_rgba(value), lw=7)
                    ax.add_patch(arc)
                    i+= angle

                cMap_fc._A = [] # fake array to print map
                cbar = plt.colorbar(cMap_fc,shrink=0.3)
                cbar.set_label("Z-score",size=9)
                cbar.ax.tick_params(labelsize=7)

            plt.axis("equal")

            plt.annotate('OriC', xy=(center_x,radius), xytext=(center_x,radius+0.4),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=8)
            plt.annotate('Ter', xy=(center_x,0), xytext=(center_x,0),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=8)
            plt.annotate(cond, xy=(center_x,radius), xytext=(center_x,center_y),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=7)
            form = kwargs.get('format','pdf')
            if form == "png":
                plt.savefig(path+"/circle-"+cond+".png", dpi=400, transparent=False)
            elif form == "svg":
                plt.savefig(path+"/circle-"+cond+".svg", transparent=False) 
            elif form == "pdf":
                plt.savefig(path+"/circle-"+cond+".pdf", transparent=False)
            else:
                print 'Unknown format, computing default...'
                plt.savefig(path+"/circle-"+cond+".pdf",  transparent=False)
            plt.close()

    def compute_state_genes(self, cond): 
        '''
        returns a npy where a row is a gene caracterised by a start pos and a gene state
        gene is considered activated above a given fc, repressed below a given fc
        '''
        gen_states = []
        for gene in self.genes_valid[cond]:
            # if activated
            if self.genes[gene].fc_pval[cond][0] >= fc_treshold_pos and self.genes[gene].fc_pval[cond][1] <= pval_treshold:
                gen_states.append([self.genes[gene].start,1])
            # if repressed
            elif self.genes[gene].fc_pval[cond][0] <= fc_treshold_neg and self.genes[gene].fc_pval[cond][1] <= pval_treshold:
                gen_states.append([self.genes[gene].start,-1])
            # if not affected
            else:
                gen_states.append([self.genes[gene].start,0])
        
        gen_states = np.array(gen_states)
        return gen_states

    def count_genes_in_windows(self, cond, gen_states, windows, increment):
        '''
        compute bins on the genome depending on windows size and increment, and
        calculate the nb of activated / repressed / non affected genes in each bin for further zscore
        '''
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
                bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),2] += 1
                bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),2] += 1
            elif state == -1: # repressed gene
                bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),3] += 1
                bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),3] += 1
            elif state == 0: # not affected gene
                bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),4] += 1
                bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),4] += 1
        
        bins = np.concatenate((bins,bins_overlap))
        return bins

    def compute_zscores_actrepr(self, tot_act, tot_repr, tot_non, bins):
        '''
        p_exp = nb of total activated genes / nb of total activated + repressed genes on the whole genome
        '''
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
   
    def compute_zscores_act(self, tot_act, tot_repr, tot_non, bins):
        '''
        p_exp = nb of total activated genes / nb of total genes on the whole genome
        '''
        zscores = []
        p_exp = float(tot_act) / float((tot_act + tot_repr + tot_non))
        print 'g+ / gtot on genome :',p_exp
        # compute zscore for each bin
        for start,end,nb_act,nb_repr,nb_null in bins:
            nb_act = float(nb_act) ; nb_tot = float((nb_act+nb_repr+nb_null))
            try:
                zscore =(nb_act - (nb_tot*p_exp)) / sqrt(nb_tot*p_exp*(1-p_exp))
            except: # division by zero if no genes activated nor repressed
                zscore = 0
            zscores.append(zscore)
        return zscores

    def compute_zscores_repr(self, tot_act, tot_repr, tot_non, bins):
        '''
        p_exp = nb of total repressed genes / nb of total genes on the whole genome
        '''
        zscores = []
        p_exp = float(tot_repr) / float((tot_act + tot_repr + tot_non))
        print 'g- / gtot on genome :',p_exp
        # compute zscore for each bin
        for start,end,nb_act,nb_repr,nb_null in bins:
            nb_repr = float(nb_repr) ; nb_tot = float((nb_act+nb_repr+nb_null))
            try:
                zscore =(nb_repr - (nb_tot*p_exp)) / sqrt(nb_tot*p_exp*(1-p_exp))
            except: # division by zero if no genes activated nor repressed
                zscore = 0
            zscores.append(zscore)
        return zscores

    def compute_spacer_response(self,cond_fc,cond_tss,*arg,**kwargs):
        '''For all TSS in TSSs[cond_tss], compute spacer length for each sigma factor from promoter site coordinates.
        Then, associate expression data available in cond_fc for genes associated to TSS to spacer lengths, which allows
        the analysis of the link between FC and spacer length. kwargs = thresh_pval = only consider genes with pval below, default 0.05. Returns a dict oh shape {[sigma_factor]:{[sp_length]:{[genes]:[g1,g2],[expr]:[1,2]}}}
        '''
        if not hasattr(self, 'genes_valid'): # if no fc loaded 
            self.load_fc_pval()
        if not hasattr(self, 'TSSs'): # if no TSS loaded
            self.load_TSS() 
        if not hasattr(self, 'seq'): # if no seq loaded
            self.load_seq() 

        thresh_pval = kwargs.get('thresh_pval', 0.05)
        try:
            thresh_pval = float(thresh_pval)
        except:
            print 'Invalid thresh_pval, setting default 0.05...'
            thresh_pval = 0.05

        spacer_sigma = {} # dict of shape {[sigma_factor]:{[sp_length]:{[genes]:[g1,g2],[expr]:[1,2]}}}

        try :
            for TSSpos in self.TSSs[cond_tss].keys(): # for all TSS in cond_tss
                TSSu = self.TSSs[cond_tss][TSSpos] # single TS
                for sig in TSSu.promoter.keys(): # for all sigma factors binding to this TSS
                    try:
                        if sig not in spacer_sigma.keys(): # add sigma factor to sigma dict
                            spacer_sigma[sig] = {'genes_valid':[],'genes_tot':[],'expr':[]}
                        try: # if site coordinates
                            TSSu.compute_magic_prom(self.seq,self.seqcompl)
                            spacer = len(TSSu.promoter[sig]['spacer'])
                        except: # if spacer length instead of sites coordinates
                            spacer = TSSu.promoter[sig]['sites'][0]

                        # valid gene : gene with pval < thresh_pval
                        expr = [] # expression values for valid genes in that TSS 
                        valid_genes = [] # valid genes in that TSS
                        tot_genes = [] # genes that appear in FC + TSS data but with pval > thresh_pval

                        for gene in TSSu.genes:
                            try:
                                if gene in self.genes_valid[cond_fc]:
                                    tot_genes.append(gene) # gene in FC data + TSS data

                                    if self.genes[gene].fc_pval[cond_fc][1] <= thresh_pval:
                                    # gene in FC data + TSS data and valid (pval < thresh)
                                        valid_genes.append(gene) 
                                        expr.append(self.genes[gene].fc_pval[cond_fc][0])
                            except:
                                pass

                        if tot_genes != []:
                            if spacer not in spacer_sigma[sig].keys(): # init spacer in dict
                                spacer_sigma[sig][spacer] = {'genes_valid':[],'genes_tot':[],'expr':[]} # shape 15:{[expr]:[expr1,expr2],[genes_valid]:[g1,g2],[genes_tot]:[g1,g2,g3]} 
                            # genes_valid and genes_tot implemented in order to keep in memory all
                            # genes with that spacer length
                            spacer_sigma[sig]['genes_tot'].append(tot_genes) # gene in FC data + TSS data
                            spacer_sigma[sig][spacer]['genes_tot'].append(tot_genes)

                            if expr != [] and valid_genes != []:
                                spacer_sigma[sig]['expr'].append(expr)
                                spacer_sigma[sig]['genes_valid'].append(valid_genes)                            
                                spacer_sigma[sig][spacer]['expr'].append(expr)                               
                                spacer_sigma[sig][spacer]['genes_valid'].append(valid_genes)

                            print 'Success, FC',expr,'for genes',valid_genes,'added to sigma',sig,'with',spacer,'pb spacer'
                    except:
                        pass

        except Exception as e:
            print 'Error',e   

        if spacer_sigma == {}:
            print 'Unable to compute spacer_sigma'
        else:
            return spacer_sigma

###################### ANTOINE #############################

    def run_btssfinder(self,nom_liste_TSS,nameOut,freedom):
        '''
        running btsfinder with default parameters for E.Coli
        obj is Genome object
        list_TSS is the TSS list eg. biocyc
        nameOut is  the name of fasta file (whiout extension .fasta) or if don't exist will become file's name of all out
        freedom is number of additionnal base at the TSS region
        NEXT
        convert out-bTSSfinder.gff on basedir/data/[nom_liste_TSS]
        AND FINALY
        write the localization of new csv in file TSS.info to the next load.TSS()
        '''
        if nameOut == None:
            run_btssfinder(self,nom_liste_TSS,nom_liste_TSS,freedom)
        else:
            run_btssfinder(self,nom_liste_TSS,nameOut,freedom)
        gff2csv(self,nom_liste_TSS,nameOut,freedom)
        TSSinfo = basedir+"data/"+self.name+"/TSS/TSS.info"
        if os.path.exists(TSSinfo):
            exist = False
            f = open(TSSinfo,"r")
            for i in f.readlines():
                line = i.split('\t')
                if line[0] == nameOut:
                    exist = True
            f.close()
            if not exist:
                f = open(TSSinfo,"a")
                f.write(nameOut+'\t'+nameOut+'.csv'+'\t'+"2"+'\t'+"0"+'\t'+"2"+'\t'+"\\t"+'\t'+"1"+'\t'+"3"+'\t'+"4"+'\n')
                f.close()
        else:
            print "TSS info not found"
        print "Finished…"+'\n'+"Now, you can visualise file "+TSSinfo+" or you can just reload TSS list."
