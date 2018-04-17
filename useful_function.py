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
from TSS import TSS



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


def compare_TSS_lists(TSSd1,TSSd2,*arg,**kwargs):
    match = {'len':[len(TSSd1.keys()),len(TSSd2.keys())]}
    bound = kwargs.get('bound',5)
    sig = kwargs.get('sig',sigfactor)

    match[0] = {'tot':0,'tot_sig':0,'m1035':0,'m10!35':0,'m35!10':0,'discr':0,'discrm1':0}
    for b in range(0,bound+1):
    # tot : total TSS matching
    # tot_sig : among TSS matching, those which have a spacer for sigfactor
        #match[b] = {'tot':0,'tot_sig':0,'m10':0,'m35':0}
        match[0][b] = 0
        # for bb in range(0,bound+1):
        #     match[b][bb] = 0
    dist = [[],[]]

    for TSS1 in TSSd1.keys():
        for TSS2 in TSSd2.keys():
            TSSdiff = abs(TSS1-TSS2) # shift between TSS
            try: # if TSSdiff between 0 (perfect match) and bound
                match[TSSdiff]['tot'] += 1
                SPdiff = abs(len(TSSd1[TSS1].promoter[sig]['spacer'])-len(TSSd2[TSS2].promoter[sig]['spacer'])) #shift between spacer lengths
                match[TSSdiff]['tot_sig'] += 1 # both TSS have spacer for sigfactor
                
                if TSSd1[TSS1].strand == True:
                    d1 = TSSd1[TSS1].pos - TSSd1[TSS1].promoter[sig]['sites'][1]-1
                    d2 = TSSd2[TSS2].pos - TSSd2[TSS2].promoter[sig]['sites'][1]-1
                
                elif TSSd2[TSS2].strand == False:
                    d1 = TSSd1[TSS1].promoter[sig]['sites'][0] - TSSd1[TSS1].pos -1
                    d2 = TSSd2[TSS2].promoter[sig]['sites'][0] - TSSd2[TSS2].pos -1

                dist[0].append(d1)
                dist[1].append(d2)

                if d == 0:
                    match[TSSdiff]['discr'] += 1
                elif d == 1:
                    match[TSSdiff]['discrm1'] += 1

                if TSSd1[TSS1].promoter[sig]['minus10'] == TSSd2[TSS2].promoter[sig]['minus10'] and TSSd1[TSS1].promoter[sig]['minus35'] == TSSd2[TSS2].promoter[sig]['minus35']:
                    match[TSSdiff]['m1035'] += 1
                if TSSd1[TSS1].promoter[sig]['minus10'] == TSSd2[TSS2].promoter[sig]['minus10'] and TSSd1[TSS1].promoter[sig]['minus35'] != TSSd2[TSS2].promoter[sig]['minus35']:
                    match[TSSdiff]['m10!35'] += 1
                if TSSd1[TSS1].promoter[sig]['minus10'] != TSSd2[TSS2].promoter[sig]['minus10'] and TSSd1[TSS1].promoter[sig]['minus35'] == TSSd2[TSS2].promoter[sig]['minus35']:
                    match[TSSdiff]['m35!10'] += 1

                match[TSSdiff][SPdiff] += 1
            except: # key error
                pass

    return match,dist

def draw_d(d):
    plt.hist([d[0],d[1]], label=['Biocyc','bTSS'], range = (4,9), bins = 6,normed=True)
    plt.ylabel('Frequency',fontweight='bold')
    plt.xlabel('Discriminator length',fontweight='bold')
    plt.legend()
    plt.show()
    # plt.title('Discriminator')



def valid_TSS(gen,cond_fc,cond_tss,thresh_fc,*arg,**kwargs):
    '''
    For given TSS cond / FC cond, extract valid TSS (TSS with genes having a valid FC value).
    '''
    sig = kwargs.get('sig',sigfactor)
    if not hasattr(gen, 'genes_valid'): # if no FC loaded 
        gen.load_fc_pval()
    if not hasattr(gen, 'TSSs'): # if no TSS loaded
        gen.load_TSS() 
    if not hasattr(gen, 'seq'): # if no seq loaded
        gen.load_seq() 

    validTSS = {'act':{'TSS':{},'values':[]},'rep':{'TSS':{},'values':[]}, 'all':{'TSS':{},'values':[]}}
    for TSSpos in gen.TSSs[cond_tss].keys(): # for all TSS in cond_tss
        TSSu = gen.TSSs[cond_tss][TSSpos] # single TS
        expr = []
        try:
            if sig in TSSu.promoter.keys():
                try:
                    TSSu.compute_magic_prom(gen.seq,gen.seqcompl)
                    spacer = len(TSSu.promoter[sig]['spacer'])
                except: # if spacer length instead of sites coordinates
                    spacer = TSSu.promoter[sig]['sites'][0]

                if spacer in spacers:
                    for gene in TSSu.genes:
                        if gene in gen.genes_valid[cond_fc]:
                            expr.append(gen.genes[gene].fc_pval[cond_fc][0])
                    if expr != []:
                        validTSS['all']['TSS'][TSSpos] = TSSu
                        val = np.mean(expr)
                        validTSS['all']['values'].append(val)
                        if val > float(thresh_fc):
                            validTSS['act']['TSS'][TSSpos] = TSSu
                            validTSS['act']['values'].append(val)
                        else:
                            validTSS['rep']['TSS'][TSSpos] = TSSu
                            validTSS['rep']['values'].append(val)
        except Exception as e:
            print e
    
    return validTSS
