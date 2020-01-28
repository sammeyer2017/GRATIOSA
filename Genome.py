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
from useful_functions import *
from itertools import groupby
from globvar import *
from Gene import Gene
from TSS import TSS
from TTS import TTS
from TU import TU
from datetime import datetime
from math import sqrt
from btssfinder import *FOpe
from scipy import stats

from BCBio import GFF
from Bio import SeqIO

#==============================================================================#

# -------------------
##### functions called by genome methods #####

def annotations_parser_general(annotations_filename,separator,tag_column,strand_column,left_column,right_column,start_line):
    ''' Called by load annotation, allows genes to be loaded from info file although most of the time, annotation is loaded
    from gff/gbk files
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
            try:
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

            except:
                pass

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

def annotations_gbk(file):
    ''' Called by load annotation, allows genes to be loaded from gbk
    '''
    genes_dict = {}
    gb_record = SeqIO.read(open(file,"r"), "genbank")
    for f in gb_record.features:
        try:
            if f.type == 'gene':
                left = f.location.start + 1
                right = int(f.location.end)
                strand = f.location.strand
                locus = f.qualifiers['locus_tag'][0]
                try:
                    name = f.qualifiers['gene'][0]
                except:
                    name = f.qualifiers['locus_tag'][0]

                genes_dict[locus]= Gene(annot_gbk=[locus,name,left,right,strand])

        except (IndexError, ValueError):
            print("Annotations : could not read line ")

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

def load_expr_cond(genes_dict, filename, condition, tag_col, nb_replicates, expr_col, separator, start_line):
    ''' Called by load_expr_cond, allows expression data to be loaded by specifying files, and where each
    information is (tag, expr...).
    '''
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

            st = int(expr_col) ; nb = int(nb_replicates)
            for i in range(st,st+nb):
                try:
                    genes_dict[line[tag_col]].add_expr_cond(float(line[i]),condition)
                except:
                    if line[tag_col] not in genes_dict.keys():
                        if line[tag_col] != '':
                            print(line[tag_col] + " not in annotation ")
                        else:
                            print("fc without locus")
    f.close()

def load_seq(filename):
    ''' Called by load_seq, allows genomic sequence to be loaded from .fasta file
    '''
    seq=str()
    my_file = open(filename, "r")
    for i in my_file.readlines():
        line=i.strip() #Removes \n
        if line != '':#Inspect if empty line
            if line[0]!=">":
                seq+=line
    my_file.close
    return seq


def load_TSS_cond(genes_dict, filename, TSS_column, start_line , separator, strandcol, genescol, sigcol, sitescol, scorecol,*args, **kwargs):
    ''' Called by load_TSS, allows TSS data to be loaded from .info file
    '''
    TSS_dict = {} # dict of TSS objects
    with open(filename, 'r') as f:
        i = 1
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
                strand = True if line[strandcol] in ["True","plus","+"] else False
                genes = line[genescol] if genescol != None else []
                sig = line[sigcol] if sigcol != None else None
                sites = line[sitescol] if sitescol != None else None
                
                # if TSS needs to be init
                if pos not in TSS_dict.keys():
                    # init object tss
                    TSS_dict[pos] = TSS(pos = pos)
                    TSS_dict[pos].add_strand(strand)
                    if genes != []:
                        TSS_dict[pos].add_genes(genes,genes_dict)
                        for gene in TSS_dict[pos].genes: # add TSS to gene attributes
                            genes_dict[gene].add_id_TSS(pos)

                # Add sigma factor and binding sites to the promoter dict
                if sig != None: # if sigma column
                    if sites != None:
                        TSS_dict[pos].add_promoter(sig, sites = sites)
                    else:
                        TSS_dict[pos].add_promoter(sig)

                if scorecol != None:
                    TSS_dict[pos].add_score(int(line[scorecol]))

            except Exception as e:
                print 'Error in line, wrong information type :',e

    return TSS_dict

def load_TTS_cond(filename, separator, start_line, leftcol, rightcol, strandcol, rhocol, seqcol, scorecol, genescol, *args, **kwargs):
    ''' Called by load_TTS, allows TTS data to be loaded from .info file
    '''
    TTS_dict = {} # dict of TSS objects
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
                left = int(line[leftcol])
                right = int(line[rightcol])
                strand = True if line[strandcol] in ["True","plus","+"] else False
                rho_dpdt = True if line[rhocol] in ["True","TRUE","1"] else False
                seq = line[seqcol] if seqcol != None else ""
                score = line[scorecol] if scorecol != None else None
                genes = line[genescol] if genescol != None else []

                newTTS = TTS(left = left, right = right, strand = strand, rho_dpdt = rho_dpdt, seq = seq, score =  score, genes = genes)
                TTS_dict[newTTS.start] = newTTS

            except Exception as e:
                print 'Error in line, wrong information type :',e

    return TTS_dict


def add_neighbour(dict_genes,list):
    for i in range(len(list)):
        if i != 0:
            dict_genes[list[i][1]].add_left_neighbour(list[i-1][1])
        if i != len(list)-1:
            dict_genes[list[i][1]].add_right_neighbour(list[i+1][1])
    return dict_genes

# ----------------------
def add_expression_to_genes(genes_dict, filename, tag_col, first_expression_col, is_log, separator):
    """ Adds expression data to Gene objects by parsing a file with as many
    columns as there are different conditions in the experiment, plus one for
    the gene names (first column).
    """        
    genes_valid = {} 
    with open(filename, 'r') as f:
        header=next(f)
        header=header.strip()
        if separator == '\\t':
            header = header.split('\t')
        else:
            header=header.split(separator)        
        header=header[first_expression_col:]
        genes_valid["conditions"] = header
        for line in f:
            line=line.strip()
            if separator == '\\t':
                line = line.split('\t')
            else:
                line = line.split(separator)
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
                    print(line[tag_col] + " not in annotation")
    return genes_valid

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

def load_TU_cond(filename, startcol, stopcol, strandcol, genescol, startline, separator, *args, **kwargs):
    ''' Called by load_TU, allows TU data to be loaded by specifying files, and where each information is
    '''
    TUs= {}
    with open(filename, 'r') as f:
        i=1
        while i < startline:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n')
            if separator == '\\t':
                line = line.split('\t')
            else:
                line=line.split(separator)
            try:
                TUs[int(float(line[startcol]))] = TU(start=int(float(line[startcol])), stop=int(float(line[stopcol])), orientation=line[strandcol], genes = line[genescol].split(","))   
            except Exception as e:
                pass
    f.close()
    return TUs



#####  #####

class Genome:

    def __init__(self, *args, **kwargs):
        """ Possible kwargs arguments: name, seq, length,
        """
        self.name = kwargs.get('name')
        self.length = kwargs.get('length')
        #self.genes=kwargs.get('genes')
        self.TSS_complete = {}
        self.TSS_plus={}
        self.TSS_minus={}

    def load_seq(self):
        self.seq=load_seq(basedir+"data/"+self.name+"/sequence.fasta").upper()
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


    def load_annotation(self):
        """ Load annotation. Two options : if gff file in directory -> load annotation from gff
        if no gff file in directory -> tries to load annotation.info (0 = file, 1 = separator ,2 =
        Locus column,3 = Strand column, 4,5 Left Rigth column, 6 start line)
        """
        custom = 0
        if custom:
            with open(basedir+"data/"+self.name+"/annotation/annotation.info","r") as f:
                for line in f:
                    line=line.strip()
                    line=line.split('\t')
                    self.genes=annotations_parser_general(basedir+"data/"+self.name+'/annotation/'+line[0],line[1],int(line[2]),int(line[3]),int(line[4]),int(line[5]),int(line[6]))
                f.close()

        else:     
            if os.path.exists(basedir+"data/"+self.name+"/annotation/sequence.gff3"):
                self.genes=annotations_parser_gff(basedir+"data/"+self.name+"/annotation/sequence.gff3")
            elif os.path.exists(basedir+"data/"+self.name+"/annotation/sequence.gff"):
                self.genes=annotations_parser_gff(basedir+"data/"+self.name+"/annotation/sequence.gff")

            else:
                try:
                    self.genes = annotations_gbk(basedir+"data/"+self.name+'/annotation/sequence.gbk')

                except Exception as e:
                    print e
                    print('No GFF file nor annotation.info, unable to load annotation')



    def load_TSS(self, *args, **kwargs):
        """ Load a TSS file info where indice 0 = condition, 1 = filename,
        2 = locus_tag, 3 = TSS_column, 4 = start_line, 5 = separator, 6 = strand column, 7 = Sig column
        8 = Sites column if much other condition give it in the seconde line of file and change TSS column """
        self.TSSs = {} # shape (dict of dict) : TSSs = {TSScond : {TSS:attr}}
        self.TSSs['all_TSS'] = {} # dict containing all TSS and where they appear (shape all_TSS = {pos:[conditions]})
        if not hasattr(self,"genes"):
            self.load_annotation()

        if os.path.exists(basedir+"data/"+self.name+"/TSS/TSS.info"):
            with open(basedir+"data/"+self.name+"/TSS/TSS.info","r") as f:
                skiphead = next(f) # skip head
                for line in f:
                    line = line.strip('\n')
                    line = line.split('\t')

                    filename = line[1] ; startline = int(line[4]) ; sep = line[5] ; strand = int(line[6]) ; TSScol = int(line[3])
                    genescol = int(line[2]) if line[2] != "" else None
                    sigcol = int(line[7]) if line[7] != "" else None
                    sitescol = int(line[8]) if line[8] != "" else None
                    scorecol = int(line[9]) if line[9] != "" else None
                    try: # successively try to load :
                        self.TSSs[line[0]] = load_TSS_cond(self.genes, basedir+"data/"+self.name+"/TSS/"+filename, TSScol, startline, sep,strand, genescol, sigcol, sitescol,scorecol)
                        # append all entries to all TSS dict
                        for entry in self.TSSs[line[0]].keys():
                            try: # works if entry already in dict
                                self.TSSs['all_TSS'][entry].append(line[0])
                            except: # init list of conditions for entry
                                self.TSSs['all_TSS'][entry] = []
                                self.TSSs['all_TSS'][entry].append(line[0])

                    except Exception as e:
                        print "Error loading",cond,e


        else:
            print("No TSS.info, unable to load TSS")


    def load_rpkm(self):
        """ Load a RPKM file information where indice 0 = Condition
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

###################### RAPH GENOME METHODS #############################

    def load_reads(self):
        '''
        Load paired end reads from .npz files that have been generated using process_bam_paired_end in useful_functions
        and which are described in reads.info
        New attribute reads : reads_pos & reads_neg, of shape {[condition] : .npy}, e.g. self.reads_pos[cond1]
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
        '''
        Load coverage either from .npz files which are described in cov.info
        or compute coverage itself from reads attribute
        New attribute cov : cov_pos & cov_neg, of shape {[condition] : .npy}, e.g. self.cov_pos[cond1]
        '''
        self.cov_pos = {} # cov on + strand
        self.cov_neg = {} # cov on - strand
        # function tries first to deal with cov_info and .npy files directly, if cov_info not available then
        # tries to open cov_txt.info, convert .txt files into .npy, create cov.info and load them
        if os.path.exists(basedir+"data/"+self.name+'/rnaseq_cov/cov.info'): # cov.info available, cov.info opening instead of cov_txt.info
            with open(basedir+"data/"+self.name+"/rnaseq_cov/cov.info","r") as f:
                header = next(f)
                for line in f: # for each condition
                    line=line.strip().split('\t')
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
                header = next(f)
                for line in f:
                    line = line.strip('\n').split('\t')
                    print 'Loading condition:',line[0]
                    # create .npy from .txt
                    cov_neg=np.loadtxt(basedir+"data/"+self.name+"/rnaseq_cov/"+line[1], usecols=[int(line[4])], skiprows= int(line[3])-1)
                    cov_pos=np.loadtxt(basedir+"data/"+self.name+"/rnaseq_cov/"+line[2], usecols=[int(line[4])], skiprows= int(line[3])-1)
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

    def load_cov_start_end(self):
        '''
        Load density of RNA fragment start / end from .npz files which are described in cov_start_stop.info
        New attribute cov : cov_start & cov_end for both + and - strands
        Corresponds to the density of RNA fragments start and ends (to locate TSS / TTS)
         & shape {[condition]:{pos:.npy, neg:.npy}}
        '''
        self.cov_start = {} # 
        self.cov_end = {} #
        # function tries first to deal with cov_info and .npy files directly, if cov_info not available then
        # tries to open cov_txt.info, convert .txt files into .npy, create cov.info and load them
        if os.path.exists(basedir+"data/"+self.name+'/rnaseq_cov/cov_start_stop.info'): # cov.info available, cov.info opening instead of cov_txt.info
            with open(basedir+"data/"+self.name+"/rnaseq_cov/cov_start_stop.info","r") as f:
                header = next(f)
                for line in f: # for each condition
                    line=line.strip().split('\t')
                    print 'Loading condition',line[0]
                    # load attributes
                    self.cov_start[line[0]] = {}
                    self.cov_end[line[0]] = {}

                    self.cov_start[line[0]][0] = np.load(basedir+"data/"+self.name+'/rnaseq_cov/'+line[1])["cov_start_neg"]
                    self.cov_start[line[0]][1] = np.load(basedir+"data/"+self.name+'/rnaseq_cov/'+line[1])["cov_start_pos"]
                    self.cov_end[line[0]][0] = np.load(basedir+"data/"+self.name+'/rnaseq_cov/'+line[1])["cov_end_neg"]
                    self.cov_end[line[0]][1] = np.load(basedir+"data/"+self.name+'/rnaseq_cov/'+line[1])["cov_end_pos"]

            f.close()

        else:
            print 'cov_start_stop.info not available please check /rnaseq_cov/ folder'

        self.cov_start_all = {0:np.sum([self.cov_start[x][0] for x in self.cov_start.keys()],axis=0), 1:np.sum([self.cov_start[x][1] for x in self.cov_start.keys()],axis=0)} # 
        self.cov_end_all = {0:np.sum([self.cov_end[x][0] for x in self.cov_end.keys()],axis=0), 1:np.sum([self.cov_end[x][1] for x in self.cov_end.keys()],axis=0)} #

        print 'Done'


    def compute_rpkm_from_cov(self, before=100):
        '''
        Adds rpkm values from coverage: along whole genes Before= number of bps to add before = to take into account
        DNA region upstream of the coding sequence of the gene
        '''
        if not self.genes: # if no genes loaded
            # try to load them
            self.load_annotation()
        if not hasattr(self,"cov_pos"):
            self.load_cov()

        try:
            for g in self.genes.keys(): # for each gene
                if self.genes[g].strand:
        # gene in + strand
                    for cond in self.cov_pos.keys(): # for each condition of cov
                        self.genes[g].add_single_rpkm(cond, np.mean(self.cov_pos[cond][(self.genes[g].left-before):self.genes[g].right]), np.sum(self.cov_pos[cond])+np.sum(self.cov_neg[cond]))
                elif not self.genes[g].strand:
        # gene in - strand
                    for cond in self.cov_neg.keys():
                        self.genes[g].add_single_rpkm(cond, np.mean(self.cov_neg[cond][self.genes[g].left:(self.genes[g].right+before)]), np.sum(self.cov_pos[cond])+np.sum(self.cov_neg[cond]))
        except:
            print("You need to load coverage pls")


    def load_fc_pval(self,*args, **kwargs):
        ''' Load fc and pval specified in fc.info. Fc info file columns (tab-separated) : condition, filename, tag column
        , fc column, file separator, startline, pvalcolumn (if not available leave empty)
        '''
        self.genes_valid = {} # for each condition, list of genes having a FC /pval value
        if not hasattr(self,'genes'):
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

    def load_expression(self,*args, **kwargs):
        ''' Load expression data specified in expression.info. Expression info file columns (tab-separated) : condition, filename, tag column
        , fc column, file separator, startline, pvalcolumn (if not available leave empty)
        '''
        if not hasattr(self,'genes_valid_expr'):
            self.genes_valid_expr = {}
        if not hasattr(self, "genes"): # if no genes loaded
            # try to load them
            self.load_annotation()

        if os.path.exists(basedir+"data/"+self.name+"/expression/expression.info"):
            with open(basedir+"data/"+self.name+"/expression/expression.info","r") as f:
                skiphead = next(f) # skip head
                for header in f:
                    header=header.strip()
                    header=header.split('\t')
                    self.genes_valid_expr[header[0]]=add_expression_to_genes(self.genes,basedir+"data/"+self.name+"/expression/"+header[0], int(header[1]), int(header[2]), header[3], header[4])
            f.close()
        else:
            print("No expression.info file, please create one")

    def compute_magic_prom(self,*arg,**kwargs):
        '''
        Extract sequences of the different promoter elements based on -35 and -10 coordinates, for all TSS conditions.
        '''
        if not hasattr(self, 'TSSs'): # if no TSS loaded
            self.load_TSS()
        if not hasattr(self, 'seq'): # if no seq loaded
            self.load_seq()
        shift = kwargs.get('shift',0) # number of nt to include beyond each region on either side, e.g. to compute angle
        prom_region = kwargs.get('prom',False)
        for cond_TSS in self.TSSs.keys():
            try:
                if cond_TSS != 'all_TSS':
                    for TSS in self.TSSs[cond_TSS].keys():
                        self.TSSs[cond_TSS][TSS].compute_magic_prom(self.seq,self.seqcompl,shift=shift,prom=prom_region)
            except:
                print 'Unable to compute magic prom :',cond_TSS


    def add_fake_expression(self,cond_fc):
        '''
        For a given FC condition cond_fc, add FC = 0 and p-value = 1 to all genes which are not in the condition 
        file but only in the annotation file e.g. Blot et al. transcriptomic dataset analysis
        '''
        self.load_fc_pval()
        for gene in self.genes.keys():               
            try: # if the gene has already an expression value for the condition, pass
                test = self.genes[gene].fc_pval[cond_fc]
            except:
                try: # if the gene has an expression value for another condition, only add key to dict and fake values
                    self.genes[gene].fc_pval[cond_fc] = (0,1)
                except: # otherwise init dict before adding fake values
                    self.genes[gene].fc_pval = {}
                    self.genes[gene].fc_pval[cond_fc] = (0,1)



    def compute_fc_from_expr(self, ctrls, conds, condname):
        ''' 
        Compute FC and p-values of a condition condname starting from a list of expression conditions : 
        control conditions ctrls, test conditions conds
        '''
        if not hasattr(self, 'genes_valid'):
            self.genes_valid = {}
        self.load_expression() # load expression values

        genes_val = [] # list containing all genes having valid expr
        for genename in self.genes.keys():
            gene = self.genes[genename]
            ctrlvals = [] # list of control values for that gene
            testvals = [] # list of test values for that gene
            try:
                for ctrl in ctrls:
                    ctrlvals.append(gene.expression[ctrl])
                for cond in conds:
                    testvals.append(gene.expression[cond])
                # add FC (meantest - meanctrl) and p-values (Student test)
                gene.add_fc_pval_cond(np.mean(testvals)-np.mean(ctrlvals),condname,stats.ttest_ind(ctrlvals,testvals,equal_var=False)[1])
                genes_val.append(genename)
            except:
                pass

        self.genes_valid[condname] = genes_val
        
    def load_neighbour_all(self):
        '''
        For each gene, find nearest neighbours (left and right) on genome, whatever their strand.
        '''
        if not self.genes:
            self.load_annotation()
        res={}
        # create dic with genes and position (res[position] = gene)
        for i in self.genes:
            res[int(self.genes[i].start)]=i
        # sort dic per position of genes
        l_res=sorted(list(res.items()), key=operator.itemgetter(0))
        # add neighbours to genes
        self.genes=add_neighbour(self.genes,l_res)


    def load_gene_orientation(self,*args,**kwargs):
        '''
        Compute gene orientation. If couple = 3, gene is considered divergent if left neighbour on - strand and right neighbour on + strand,
        convergent if left neighbour on + strand and right neighbour on - strand, tandem if left and right neighbours on same strand
        (whatever the strand of the given gene is). If couple = 2, gene is considered tandem if predecessor (left neighbour for gene on + strand,
        right neighbour for gene on - strand) is on same strand, divergent if the predecessor is on opposite strand.
        '''
        self.load_neighbour_all()
        bound = kwargs.get('bound',5000) # maximal distance for seeking neighbour, either left or right
        couple = kwargs.get('couple',3)
        res = {"tandem":0,"divergent":0,"convergent":0,"isolated":0}
        for gene in self.genes:
            try:               
                g = self.genes[gene]
                lg = self.genes[g.left_neighbour]
                rg = self.genes[g.right_neighbour]
                if couple == 3:
                    if (g.start - lg.start) < bound and (rg.start - g.start) < bound:
                        if not lg.strand and not rg.strand:
                            self.genes[gene].add_orientation('tandem')
                            res["tandem"] += 1
                        elif lg.strand and rg.strand:
                            self.genes[gene].add_orientation('tandem')
                            res["tandem"] += 1
                        elif lg.strand and not rg.strand:
                            self.genes[gene].add_orientation('convergent')
                            res["convergent"] += 1
                        elif not lg.strand and rg.strand:
                            self.genes[gene].add_orientation('divergent')
                            res["divergent"] += 1
                    else:
                        self.genes[gene].add_orientation('isolated')
                        res["isolated"] += 1
                
                elif couple == 2:
                    if g.strand:
                        if (g.start - lg.start) < bound:
                            if lg.strand:
                                self.genes[gene].add_orientation('tandem')
                                res["tandem"] += 1
                            elif not lg.strand:
                                self.genes[gene].add_orientation('divergent')
                                res["divergent"] += 1
                    if not g.strand:
                        if (rg.start - g.start) < bound:
                            if rg.strand:
                                self.genes[gene].add_orientation('divergent')
                                res["divergent"] += 1
                            elif not rg.strand:
                                self.genes[gene].add_orientation('tandem')
                                res["tandem"] += 1
            except:
                pass
                
        self.orientation = res

    def compute_state_from_fc(self,*args,**kwargs):
        '''
        Compute gene state from FC data. For each condition, below a given pvalue threshold, the gene is considered 
        significantly activated if its FC is above a given FC treshold, repressed below, and non affected either if its 
        pvalue is above the threshold, or if its FC is between + and - thesholds
        '''

        if not hasattr(self, 'genes_valid'):
            self.load_fc_pval() 
        thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
        thresh_fc = kwargs.get('thresh_fc', 0) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated
        for gene in self.genes.keys():
            g = self.genes[gene]
            for cond_fc in self.genes_valid.keys():
                try:               
                    if g.fc_pval[cond_fc][1] <= thresh_pval:
                        if g.fc_pval[cond_fc][0] <  0 - thresh_fc:
                            g.add_state_cond(cond_fc,'rep')
                        elif g.fc_pval[cond_fc][0] >  0 + thresh_fc:
                            g.add_state_cond(cond_fc,'act')
                        else:
                            g.add_state_cond(cond_fc,'non')
                    else:
                        g.add_state_cond(cond_fc,'non')
                except:
                    g.add_state_cond(cond_fc,'null')

    def load_TU(self,*args, **kwargs):
        ''' Load TUs specified in TU.info
        '''
        self.TUs = {}
        if not hasattr(self,"genes"): # if no genes loaded
            # try to load them
            self.load_annotation()

        if os.path.exists(basedir+"data/"+self.name+"/TU/TU.info"):
            with open(basedir+"data/"+self.name+"/TU/TU.info","r") as f:
                skiphead = next(f) # skip head
                for header in f:
                    header=header.strip()
                    header=header.split('\t')
                    try:                
                        self.TUs[header[0]]=load_TU_cond(basedir+"data/"+self.name+"/TU/"+header[1],int(header[2]),int(header[3]),int(header[4]),int(header[5]),int(header[6]), header[7])
                    except:
                        print("Error loading cond",header[0])
            f.close()
        else:
            print("No TU.info file, please create one")

    def load_TTS(self,*args, **kwargs):
        ''' Load TTS specified in TTS.info
        '''
        self.TTSs = {}
        if not hasattr(self,"genes"): # if no genes loaded
            # try to load them
            self.load_annotation()

        if os.path.exists(basedir+"data/"+self.name+"/TTS/TTS.info"):
            with open(basedir+"data/"+self.name+"/TTS/TTS.info","r") as f:
                skiphead = next(f) # skip head
                for header in f:
                    header=header.strip()
                    header=header.split('\t')
                    leftcol = int(header[2]) ; rightcol = int(header[3]) ; strandcol = int(header[4]) ; rhocol = int(header[10])
                    seqcol = int(header[7]) if header[7] != "" else None
                    scorecol = int(header[8]) if header[8] != "" else None
                    genescol = int(header[9]) if header[9] != "" else None
                    try:                
                        self.TTSs[header[0]]=load_TTS_cond(basedir+"data/"+self.name+"/TTS/"+header[1],header[6], int(header[5]), leftcol, rightcol, strandcol, rhocol, seqcol, scorecol, genescol)
                    except Exception as e:
                        print e
                        print("Error loading cond",header[0])
            f.close()
            self.TTSs["all"] = {}
            for x in self.TTSs.keys():
                self.TTSs["all"].update(self.TTSs[x])
        else:
            print("No TTS.info file, please create one")

###################### ANTOINE #############################

    def run_btssfinder(self,list_TSS,*args,**kwargs): #running bTSSfinder
        freedom = kwargs.get('free',0)
        nameOut = kwargs.get('out',list_TSS+'-btss-'+str(freedom))
        '''
        kwags possible : out & free
        bTSSfinder need to be install on this computer
        for more informations you can read btssfinder.py
        obj is Genome object
        list_TSS is the TSS list eg. biocyc
        out is  the name of fasta file (whiout extension .fasta) or if don't exist will become file's name of all out
        free is number of additionnal base at the TSS region
        NEXT
        convert out-bTSSfinder.gff on basedir/data/[nom_liste_TSS]
        AND FINALY
        write the localization of new csv in file TSS.info to the next load.TSS()
        '''
        try:
            test = self.TSSs[list_TSS]
            test = self.seq[0:4]
        except:
            self.load_TSS()
            self.load_seq()

        run_btssfinder(self,list_TSS,nameOut,freedom)

        gff2csv(self,list_TSS,nameOut,freedom)
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


