#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from globvar import *
from Gene import Gene
from TSS_TTS_TU import TSS, TTS, TU
from Bio import SeqIO

#==============================================================================#

##### functions called by genome methods #####

def load_seq(filename):
    ''' Called by load_seq, allows genomic sequence to be loaded from .fasta file
    '''
    seq=str()
    seq_file = open(filename, "r")
    for i in seq_file.readlines():
        line=i.strip() #Removes \n
        if line != '':#Inspect if empty line
            if line[0]!=">":
                seq+=line
    seq_file.close
    return seq

def load_annot_general(annotations_filename,separator,tag_column,name_column,ID_column,strand_column,left_column,right_column,ASAP_column,start_line):
    ''' Called by load annotation, allows genes to be loaded from info file although most of the time, annotation is loaded
    from gff files
    '''
    genes_dict = {}
    with open(annotations_filename, 'r') as f:
        head=next(f)
        headl=head.strip()
        if separator == '\\t': separator = '\t'
        headl=headl.split(separator)
        i=1
        while i < start_line:
            header=next(f)
            i+=1
        for line in f:
            try:
                line=line.strip()
                line=line.split(separator)
                locus = line[tag_column]
                name = line[name_column]
                ID = line[ID_column]
                ASAP = line[ASAP_column]
                if line[strand_column] in ["complement","-1","-",False]:
                    strand = False
                elif line[strand_column] in ["forward","1","+",True]:
                    strand = True
                else:
                    print("Unrecognized strand type")
                    sys.exit("Unrecognized strand type")
                left = line[left_column]
                right = line[right_column]            
                genes_dict[locus]=Gene(locus,name,ID,left,right,strand)
            except Exception as e :
                if line != [''] :
                    print(f"WARNING: Error with {line}")
    return genes_dict


def load_gff(annotations_filename):
    ''' Called by load annotation, allows genes to be loaded from gff
        ASAP name will work only if ASAP name is in the "name" line 
        or the line next to it
    '''
    genes_dict = {}
    annot_file = open(annotations_filename, "r")
    count_ASAP = 0 
    for line in annot_file.readlines():
        ASAP_name=''
        if line[0]!= '#':
            if line != '\n':
                line=line.split('\t')
                gene_info={}
                for x in line[8].split(';'):
                    x=x.strip().split('=')
                    gene_info[x[0]]=x[1]            
                if('locus_tag' in gene_info):
                    locus = gene_info['locus_tag']
                    if "ID" in gene_info.keys() :
                        ID = gene_info["ID"] 
                    else :
                        ID = locus
                    if "gene" in gene_info.keys() : 
                        name = gene_info["gene"] 
                    else :
                        name = locus
                
                    if "Dbxref" in gene_info.keys() : 
                        db_ref = "Dbxref"
                    elif "db_xref" in gene_info.keys() :
                        db_ref = "db_xref"
                    else : 
                        db_ref = None
                    if db_ref is not None : 
                        dbxref = gene_info[db_ref].split(",")
                        for db in dbxref :
                            if "ASAP" in db :
                                ASAP_name = db.split(":")[1]
                                count_ASAP += 1
                    left, right = map(int,line[3:5])
                    strand = True if line[6] == "+" else False
                    genes_dict[locus]= Gene(locus,name,ID,left,right,strand,ASAP_name)
                
                elif line[2] == "CDS" and ASAP_name =='':
                    tleft, tright = map(int,line[3:5])
                    if tleft == left and tright == right : 
                        if "Dbxref" in gene_info.keys() : 
                            db_ref = "Dbxref"
                        elif "db_xref" in gene_info.keys() :
                            db_ref = "db_xref"
                        else : 
                            db_ref = None
                        if db_ref is not None : 
                            dbxref = gene_info[db_ref].split(",")
                            for db in dbxref :
                                if "ASAP" in db :
                                    ASAP_name = db.split(":")[1]
                                    count_ASAP += 1
                        genes_dict[locus]= Gene(locus,name,ID,left,right,strand,ASAP_name)
    if count_ASAP < 100 :
        print(f"Only {count_ASAP} genes have a ASAP annotation. If your annotation file should"
              " have more ASAP annotations, please verify that ASAP annotations are either"
              " on the same line than the locus_tags or on a \"CDS\" line right after it")
    return genes_dict              


def load_TSS_cond(genes_dict, filename, TSS_column, start_line , separator, strandcol, genescol, sigcol, sitescol, scorecol,*args, **kwargs):
    ''' Called by load_TSS, allows TSS data to be loaded from .info file
    '''
    TSS_dict = {} # dict of TSS objects
    with open(filename, 'r') as f:
        i = 0
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
                strand = True if line[strandcol] in ["True","plus","+","1"] else False
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
                    TSS_dict[pos].add_score(int(float(line[scorecol])))

            except Exception as e:
                print('Error in line, wrong information type :',e)

    return TSS_dict

def load_TTS_cond(filename, separator, start_line, leftcol, rightcol, strandcol, rhocol, seqcol, scorecol, genescol, *args, **kwargs):
    ''' Called by load_TTS, allows TTS data to be loaded from .info file
    '''
    TTS_dict = {} # dict of TTS objects
    with open(filename, 'r') as f:
        i=0
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
                print('Error in line, wrong information type :',e)

    return TTS_dict

def load_TU_cond(filename, startcol, stopcol, strandcol, genescol, startline, separator):
    ''' Called by load_TU, allows TU data to be loaded by specifying files, and where each information is
    '''
    TUs= {}
    with open(filename, 'r') as f:
        i=0
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
                TUs[int(float(line[startcol]))] = TU(start=int(float(line[startcol])), stop=int(float(line[stopcol])), strand=line[strandcol], genes = line[genescol].split(","))   
            except Exception as e:
                print(e)
    f.close()
    return TUs