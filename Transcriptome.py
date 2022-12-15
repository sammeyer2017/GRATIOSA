#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from useful_functions_transcriptome_load import *
from useful_functions_transcriptome_compute import *
from globvar import *
from datetime import datetime
#from btssfinder import *
from scipy import stats
import Genome

#==============================================================================#


class Transcriptome:

    def __init__(self,name, *args, **kwargs):
        """ Possible kwargs arguments: name, seq, length,
        """
        self.name = name


    def load_cov_rnaseq(self,cond="all"):
        '''
        Load coverage either from .npz files which are described in cov.info
        or compute coverage itself from reads attribute
        New attribute cov : cov_rnaseq_pos & cov_rnaseq_neg, of shape {[condition] : .npy}, e.g. self.cov_rnaseq_pos[cond1]
        '''
        if not hasattr(self,"cov_rnaseq_pos"):
            self.cov_rnaseq_pos = {} # cov on + strand
            self.cov_rnaseq_neg = {} # cov on - strand

        # conversion from str to list (useful if only 1 condition is selected)
        if type(cond) == str:         
            cond = [cond]

        # function tries first to deal with cov_info and .npy files directly, if cov_info not available then
        # tries to open cov_txt.info, convert .txt files into .npy, create cov.info and load them
        list_cond_info = [] #list of the conditions described in the .info file
        path2dir = f"{basedir}data/{self.name}/cov_rnaseq/"
        if os.path.exists(f"{path2dir}cov.info"): # cov.info available, cov.info opening instead of cov_txt.info
            with open(f"{path2dir}cov.info","r") as f:
                header = next(f)
                for line in f: # for each condition
                    line=line.strip().split('\t')

                    list_cond_info.append(line[0])

                    # For each selected condition : 
                    if cond == ["all"] or line[0] in cond :
                        print('Loading condition',line[0])

                        # load attributes
                        self.cov_rnaseq_neg[line[0]]= np.load(path2dir+line[1])["cov_neg"]
                        self.cov_rnaseq_pos[line[0]]= np.load(path2dir+line[1])["cov_pos"]
            f.close()

        elif os.path.exists(f"{path2dir}cov_txt.info"):
            print('Unable to locate cov.info in /cov_rnaseq/')
            print('Working with .txt file (cov_txt.info)')
            file = open(f"{path2dir}cov.info",'w')
            file.write('Condition\tCov file\tDate\tReads file')
            file.close()
            with open(f"{path2dir}cov_txt.info","r") as f: # load cov_txt.info
                header = next(f)
                for line in f:
                    line = line.strip('\n').split('\t')
                    print('Loading condition:',line[0])
                    # create .npy from .txt
                    cov_rnaseq_neg=np.loadtxt(path2dir+line[1], usecols=[int(line[4])], skiprows= int(line[3])-1)
                    cov_rnaseq_pos=np.loadtxt(path2dir+line[2], usecols=[int(line[4])], skiprows= int(line[3])-1)
                    # load attributes
                    self.cov_rnaseq_neg[line[0]]= cov_rnaseq_neg
                    self.cov_rnaseq_pos[line[0]]= cov_rnaseq_pos
                    # save .npy into .npz
                    np.savez(f"{path2dir}{line[0]}_cov.npz", cov_rnaseq_pos=cov_rnaseq_pos, cov_rnaseq_neg=cov_rnaseq_neg)
                    # update cov.info
                    file=open(f"{path2dir}cov.info","a")
                    file.write('\n'+line[0]+'\t'+line[0]+'_cov.npz\t'+str(datetime.now())+'\tUnknown')
                    file.close()
            f.close()
        
        if not os.path.exists(f"{path2dir}cov.info") and not os.path.exists(f"{path2dir}cov_txt.info",):
            print('cov.info not available nor cov_txt.info, please check /cov_rnaseq/ folder')



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
        if os.path.exists(basedir+"data/"+self.name+'/cov_rnaseq/cov_start_stop.info'): # cov.info available, cov.info opening instead of cov_txt.info
            with open(basedir+"data/"+self.name+"/cov_rnaseq/cov_start_stop.info","r") as f:
                header = next(f)
                for line in f: # for each condition
                    line=line.strip().split('\t')
                    print('Loading condition',line[0])
                    # load attributes
                    self.cov_start[line[0]] = {}
                    self.cov_end[line[0]] = {}

                    self.cov_start[line[0]][0] = np.load(basedir+"data/"+self.name+'/cov_rnaseq/'+line[1])["cov_start_neg"]
                    self.cov_start[line[0]][1] = np.load(basedir+"data/"+self.name+'/cov_rnaseq/'+line[1])["cov_start_pos"]
                    self.cov_end[line[0]][0] = np.load(basedir+"data/"+self.name+'/cov_rnaseq/'+line[1])["cov_end_neg"]
                    self.cov_end[line[0]][1] = np.load(basedir+"data/"+self.name+'/cov_rnaseq/'+line[1])["cov_end_pos"]

            f.close()

        else:
            print('cov_start_stop.info not available please check /cov_rnaseq/ folder')

        self.cov_start_all = {0:np.sum([self.cov_start[x][0] for x in self.cov_start.keys()],axis=0), 1:np.sum([self.cov_start[x][1] for x in self.cov_start.keys()],axis=0)} # 
        self.cov_end_all = {0:np.sum([self.cov_end[x][0] for x in self.cov_end.keys()],axis=0), 1:np.sum([self.cov_end[x][1] for x in self.cov_end.keys()],axis=0)} #

        print('Done')


    def compute_rpkm_from_cov(self, before=100):
        '''
        Adds rpkm values from coverage: along whole genes Before= number of bps to add before = to take into account
        DNA region upstream of the coding sequence of the gene
        '''
        self.load_annotation()
        if not hasattr(self,"cov_rnaseq_pos"):
            self.load_cov_rnaseq()

        for g in self.genes.keys(): # for each gene
            try:
                if self.genes[g].strand:
        # gene in + strand
                    for cond in self.cov_rnaseq_pos.keys(): # for each condition of cov
                        self.genes[g].add_single_rpkm(cond, np.mean(self.cov_rnaseq_pos[cond][(self.genes[g].left-before):self.genes[g].right]), np.sum(self.cov_rnaseq_pos[cond])+np.sum(self.cov_rnaseq_neg[cond]))
                elif not self.genes[g].strand:
        # gene in - strand
                    for cond in self.cov_rnaseq_neg.keys():
                        self.genes[g].add_single_rpkm(cond, np.mean(self.cov_rnaseq_neg[cond][self.genes[g].left:(self.genes[g].right+before)]), np.sum(self.cov_rnaseq_pos[cond])+np.sum(self.cov_rnaseq_neg[cond]))
            except Exception as e:
                print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno), type(e).__name__, e)
                pass

   
    def load_fc_pval(self,*args, **kwargs):
        ''' Load fc and pval specified in fc.info. Fc info file columns (tab-separated) : condition, filename, tag column
        , fc column, file separator, startline, pvalcolumn (if not available leave empty)
        '''
        self.genes_valid = {} # for each condition, list of genes having a FC /pval value
        
        if not hasattr(self,'genes'):
            gen = Genome.Genome(self.name)            
            gen.load_annotation()
            self.genes = gen.genes
            
        if os.path.exists(basedir+"data/"+self.name+"/fold_changes/fc.info"):
            with open(basedir+"data/"+self.name+"/fold_changes/fc.info","r") as f:
                skiphead = next(f) # skip head
                for header in f:
                    header=header.strip()
                    header=header.split('\t')
                    print(f"Loading condition: {header[0]}")
                    try: # if p-value column specified works
                        self.genes_valid[header[0]]=load_fc_pval_cond(self.genes,basedir+"data/"+self.name+"/fold_changes/"+header[1],header[0],int(header[2]),int(header[3]),header[4],int(header[5]), p_value=int(header[6]))
                    except: # otherwise, set pvalue = 0
                        self.genes_valid[header[0]]=load_fc_pval_cond(self.genes,basedir+"data/"+self.name+"/fold_changes/"+header[1],header[0],int(header[2]),int(header[3]),header[4],int(header[5]))
            f.close()
        else:
            print("No fc.info file, please create one")

    def load_expression(self):
        ''' Load expression data specified in expression.info. Expression info file columns (tab-separated) : condition, filename, tag column
        , fc column, file separator, startline, pvalcolumn (if not available leave empty)
        '''
        if not hasattr(self,'genes_valid_expr'):
            self.genes_valid_expr = {}
        if not hasattr(self, "genes"): # if no genes loaded
            # try to load them
            gen = Genome(self.name)
            gen.load_annotation()
            self.genes = gen.genes

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

    def compute_state_from_fc(self,*args,**kwargs):
        '''
        Compute gene state from FC data. For each condition, below a given pvalue threshold, the gene is considered 
        significantly activated if its FC is above a given FC treshold, repressed below, and non affected either if its 
        pvalue is above the threshold, or if its FC is between + and - thesholds
        '''

        if not hasattr(self, 'genes_valid'):
            self.load_fc_pval() 
        if not hasattr(self, "genes"): # if no genes loaded
            # try to load them
            gen = Genome(self.name)
            gen.load_annotation()
            self.genes = gen.genes
        if not hasattr(self,"statesFC") :
            self.statesFC = {}

        thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
        thresh_fc = kwargs.get('thresh_fc', 0) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated
        

        for cond_fc in self.genes_valid.keys():
            self.statesFC[cond_fc] = {"rep" : [],"act":[],"non":[]} 
        for gene in self.genes.keys():
            g = self.genes[gene]
            for cond_fc in self.genes_valid.keys():
                try:               
                    if g.fc_pval[cond_fc][1] <= thresh_pval:
                        if g.fc_pval[cond_fc][0] <  0 - thresh_fc:
                            g.add_state_cond(cond_fc,'rep')
                            self.statesFC[cond_fc]['rep'].append(gene)
                        elif g.fc_pval[cond_fc][0] >  0 + thresh_fc:
                            g.add_state_cond(cond_fc,'act')
                            self.statesFC[cond_fc]['act'].append(gene)
                        else:
                            g.add_state_cond(cond_fc,'non')
                            self.statesFC[cond_fc]['non'].append(gene)
                    else:
                        g.add_state_cond(cond_fc,'non')
                        self.statesFC[cond_fc]['non'].append(gene)
                except:
                    if not hasattr(self.statesFC[cond_fc],'null') :
                        self.statesFC[cond_fc]['null'] = [gene] 
                    else :
                        self.statesFC[cond_fc]['null'].append(gene)
                    g.add_state_cond(cond_fc,'null')
