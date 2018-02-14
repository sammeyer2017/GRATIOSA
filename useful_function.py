# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd
import inspect as inspect
import math
import matplotlib.pyplot as plt
import operator
from globvar import *
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import motifs



#==============================================================================#

# -------------------
# useful function for weight matrix

def create_matrix_from_file(filename, factor):
    i =1
    with open(filename,"r") as f:
        for line in f:
            if 'Transcription Factor Name: '+factor in line:
                while i <5:
                    i+=1
                    header=next(f)
                m = motifs.Motif()
                a = IUPAC.unambiguous_dna
                m.add_instance(Seq(header.strip(),a))
                l=('*')
                j=1
                while j != len(header.strip()):
                    j+=1
                    l+='*'
                while header[0] != '\n':
                    header=next(f)
                    if header[0] != '\n':
                        m.add_instance(Seq(header.strip(),a))
                m.set_mask(l)
                if header[0]=='\n':
                    return m

def create_matrix_from_file_2(filename, factor):
    i =1
    with open(filename,"r") as f:
        for line in f:
            if 'Transcription Factor Name: '+factor in line:
                while i <5:
                    i+=1
                    header=next(f)
                instances=[Seq(header.strip())]
                while header[0] != '\n':
                    header=next(f)
                    if header[0] != '\n':
                        instances.append(Seq(header.strip()))
                m = motifs.create(instances)
                if header[0]=='\n':
                    return m


#==============================================================================#

# -------------------
# useful function for SIST


def create_fasta(filename, start, end, sequence):
    with open(filename,"r") as f:
        header = next(f)
        fw=open("SIST/sequence.fasta","w")
        fw.write(header)
        f.close()
        fw.write(sequence[start-1:end])
        fw.close()

def load_profile(filename, start, end, sequence,*args, **kwargs):
    create_fasta(filename, start, end, sequence)
    profile={}
    option = kwargs.get('option')
    if option == 'A':
        profile['cruciform']={}
        profile['Z-DNA']={}
        profile['melt']={}
        os.system("perl "+basedir+"code/SIST/master.pl -a A -f "+basedir+"code/SIST/sequence.fasta -o "+basedir+"code/SIST/new.txt")
        os.system("rm "+basedir+"code/SIST/*.html")
    elif option == 'Z':
        profile['Z-DNA']={}
        os.system("perl "+basedir+"code/SIST/master.pl -a Z -f "+basedir+"code/SIST/sequence.fasta -o "+basedir+"code/SIST/new.txt")
    elif option == 'C':
        profile['cruciform']={}
        os.system("perl "+basedir+"code/SIST/master.pl -a C -f "+basedir+"code/SIST/sequence.fasta -o "+basedir+"code/SIST/new.txt")
        os.system("rm "+basedir+"code/SIST/*.html")
    else:
        profile['melt']={}
        os.system("perl "+basedir+"code/SIST/master.pl -a M -f "+basedir+"code/SIST/sequence.fasta -o "+basedir+"code/SIST/new.txt")
    with open(basedir+"code/SIST/new.txt","r") as f:
        header=next(f)
        while header[0] !='P':
            header=next(f)
        for line in f:
            line=line.strip()
            line=line.split('\t')
            if option == 'Z':
                profile['Z-DNA'][start]=float(line[1])
            elif option == 'C':
                profile['cruciform'][start]=float(line[1])
            elif option =='A':
                profile['melt'][start]=float(line[1])
                profile['Z-DNA'][start]=float(line[2])
                profile['cruciform'][start]=float(line[3])
            else:
                profile['melt'][start]=float(line[1])
            start+=1
    os.system("rm "+basedir+"code/SIST/*.txt")
    os.system("rm "+basedir+"code/SIST/*.fasta")
    return profile



#==============================================================================#

# -------------------
# useful function for get_coverage


def download_pair(filename,name):
    with open(filename, 'r') as f:
        header=next(f)
        pair=0
        os.system("bowtie-build "+basedir+"data/"+name+"/sequence.fasta "+basedir+"data/"+name+"/rnaseq_cov/index")
        list_cov=[]
        i=0
        for line in f:
            line=line.strip()
            line=line.split('\t')
            if i ==0:
                i=search_link(line)
            new=line[i].split('/')[-1]
            new=new.split('.')[0]
            condition=line[0]
            if pair == 0:
                first=new+'.fastq'
            if pair == 1:
                second=new+'.fastq'
            pair +=1
            if not(os.path.exists(basedir+"data/"+name+"/rnaseq_cov/"+condition+"_plus.txt")):
                os.system("wget " + line[i] + " -P "+basedir+"data/"+name+"/rnaseq_cov")
                os.system("gunzip "+basedir+"data/"+name+"/rnaseq_cov/"+new+'.fastq')
                if pair ==2:
                    pair =0
                    os.system("bowtie -S "+basedir+"data/"+name+"/rnaseq_cov/index -1 "+basedir+"data/"+name+"/rnaseq_cov/"+first+" -2 "+basedir+"data/"+name+"/rnaseq_cov/"+second+" "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".sam")
                    os.system("rm "+basedir+"data/"+name+"/rnaseq_cov/"+first+"*")
                    os.system("rm "+basedir+"data/"+name+"/rnaseq_cov/"+second+"*")
                    os.system("samtools view -S -b "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".sam > "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".bam")
                    os.system("rm "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".sam")
                    os.system("bedtools genomecov -ibam "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".bam -strand + -d > "+basedir+"data/"+name+"/rnaseq_cov/"+condition+"_plus.txt")
                    os.system("bedtools genomecov -ibam "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".bam -strand - -d > "+basedir+"data/"+name+"/rnaseq_cov/"+condition+"_minus.txt")
                    os.system("rm "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".bam")
                    list_cov.append(condition)
                    list_cov.append(condition+"_plus.txt")
                    list_cov.append(condition+"_minus.txt")
            if pair == 2:
                pair =0
        os.system("rm "+basedir+"data/"+name+"/rnaseq_cov/index*")
    f.close()
    return list_cov

def download_single(filename,name):
    with open(filename, 'r') as f:
        header=next(f)
        pair=0
        os.system("bowtie-build "+basedir+"data/"+name+"/sequence.fasta "+basedir+"data/"+name+"/rnaseq_cov/index")
        list_cov=[]
        i=0
        for line in f:
            line=line.strip()
            line=line.split('\t')
            if i ==0:
                i=search_link(line)
            new=line[i].split('/')[-1]
            new=new.split('.')[0]
            condition=line[0]
            single=new+'.fastq'
            if not(os.path.exists(basedir+"data/"+name+"/rnaseq_cov/"+condition+"_plus.txt")):
                os.system("wget " + line[i] + " -P "+basedir+"data/"+name+"/rnaseq_cov")
                os.system("gunzip "+basedir+"data/"+name+"/rnaseq_cov/"+new+'.fastq')
                os.system("bowtie -S "+basedir+"data/"+name+"/rnaseq_cov/index "+basedir+"data/"+name+"/rnaseq_cov/"+single+" "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".sam")
                os.system("rm "+basedir+"data/"+name+"/rnaseq_cov/"+single+"*")
                os.system("samtools view -S -b "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".sam > "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".bam")
                os.system("rm "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".sam")
                os.system("bedtools genomecov -ibam "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".bam -strand + -d > "+basedir+"data/"+name+"/rnaseq_cov/"+condition+"_plus.txt")
                os.system("bedtools genomecov -ibam "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".bam -strand - -d > "+basedir+"data/"+name+"/rnaseq_cov/"+condition+"_minus.txt")
                os.system("rm "+basedir+"data/"+name+"/rnaseq_cov/"+condition+".bam")
                list_cov.append(condition)
                list_cov.append(condition+"_plus.txt")
                list_cov.append(condition+"_minus.txt")
        os.system("rm "+basedir+"data/"+name+"/rnaseq_cov/index*")
    f.close()
    return list_cov



def create_cov_info(list_cov,name):
    l=list_cov
    pair=0
    fw=open(basedir+"data/"+name+"/rnaseq_cov/cov.info","a")
    for i in l:
        fw.write(i)
        pair+=1
        if pair == 1:
            fw.write('\t')
        elif pair == 2:
            fw.write('\t2\t')
        elif pair == 3:
            pair =0
            fw.write('\t2\n')
    fw.close()

def search_link(list):
    i=0
    for j in list:
        if 'ftp' in j:
            return i
        i+=1
