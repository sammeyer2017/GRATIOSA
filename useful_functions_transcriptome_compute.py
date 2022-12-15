#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from globvar import *
import os
import pysam
from datetime import datetime

#==============================================================================#

##### functions called by transcriptome methods #####

########## BAM2NPZ RNA_SEQ ANALYSIS, RAPHAEL FORQUET ##########

def process_bam_paired_end(bam_file, gen): # /!\ paired-end only /!\ -> return fragments for each strand
    if not os.path.exists(basedir+"data/"+gen.name+'/rnaseq_reads/'+bam_file+".bai"): # if index needed, created using samtools (.bai file)
        os.system("samtools index -b %s"%bam_file)
    
    print(gen.name)
    bamfile = pysam.AlignmentFile(basedir+"data/"+gen.name+'/rnaseq_reads/'+bam_file, "rb") # BAM opening, alignment file object
    print('Header :',bamfile.header)
    print('Reads mapped :',bamfile.mapped)
    print('Reads without coordinate :',bamfile.nocoordinate)
    print('Reads not mapped :',bamfile.unmapped)
    # /!\ 0-based coordinate system /!\
    # fastest way to build a numpy matrix -> store every coordinate of interest into lists, then merge into numpy array
    # lists storing coordinates of paired-end fragments for +/- strands
    Rpos_start = []
    Rpos_end = []
    Rneg_start = []
    Rneg_end = []


    for read in bamfile.fetch(): 
        # if correctly paired and R1/R2 in different strands and fragment size < 1000 kb
        if read.is_read1 and read.is_paired and not read.mate_is_unmapped and read.is_reverse != read.mate_is_reverse and abs(read.template_length) < 1000 :
            if read.is_reverse: 
                Rneg_start.append(read.reference_end)
                Rneg_end.append(read.reference_end + abs(read.template_length))
            elif not read.is_reverse:
                Rpos_start.append(read.reference_start)
                Rpos_end.append(read.reference_start + read.template_length)
    bamfile.close();
 
    # conversion step into numpy array
    Rpos_start = np.array(Rpos_start,dtype=int)
    Rpos_end = np.array(Rpos_end,dtype=int)
    Rneg_start = np.array(Rneg_start,dtype=int)
    Rneg_end = np.array(Rneg_end,dtype=int)

    Rpos = np.column_stack((Rpos_start,Rpos_end))
    Rneg = np.column_stack((Rneg_start,Rneg_end))

    # delete rows where one coordinate is missing
    Rpos = Rpos[np.isfinite(Rpos).all(axis=1)] 
    Rneg = Rneg[np.isfinite(Rneg).all(axis=1)]

    # if reads.info not exist
    if not os.path.exists(basedir+"data/"+gen+'/rnaseq_reads/reads.info'):
        file = open(basedir+"data/"+gen+'/rnaseq_reads/reads.info','w') 
        file.write('Condition\tReads file\tDate\tBAM file')
        file.close() 
    # save results ; .npz contains two .npy Rpos and Rneg
    file = open(basedir+"data/"+gen+'/rnaseq_reads/reads.info','a')
    file.write('\n'+bam_file[:-4]+'\t'+bam_file[:-4]+'_reads.npz\t'+str(datetime.now())+'\t'+bam_file)  
    file.close()    
    np.savez(basedir+"data/"+gen+'/rnaseq_reads/'+bam_file[:-4]+'_reads.npz', Rpos=Rpos, Rneg=Rneg)


# compute coverage from reads (.npz, one .npy per strand)
def cov_from_reads(npz_file, gen, genome_length):
    # load npz file
    npzfile = np.load(basedir+"data/"+gen+'/rnaseq_reads/'+npz_file)
    Rpos = npzfile["Rpos"]
    Rneg = npzfile["Rneg"]
    # init cov
    cov_pos = np.zeros(genome_length, dtype=int)
    cov_neg = np.zeros(genome_length, dtype=int)
    # compute cov
    # on positive strand
    for start,end in Rpos:
        cov_pos[start:end+1] += 1
    # on negative strand
    for start,end in Rneg:
        cov_neg[start:end+1] += 1
    # if rnaseq_cov folder not exist        
    # if cov.info not exist
    if not os.path.exists(basedir+"data/"+gen+'/rnaseq_cov/cov.info'):
        file = open(basedir+"data/"+gen+'/rnaseq_cov/cov.info','w') 
        file.write('Condition\tCov file\tDate\tReads file')
        file.close() 
    # save results
    file = open(basedir+"data/"+gen+'/rnaseq_cov/cov.info','a')
    file.write('\n'+npz_file[0:-10]+'\t'+npz_file[:-4]+'_cov.npz\t'+str(datetime.now())+'\t'+npz_file)  
    file.close()
    # save results ; .npz contains two .npy cov_pos and cov_neg
    np.savez(basedir+"data/"+gen+'/rnaseq_cov/'+npz_file[:-4]+'_cov.npz', cov_pos=cov_pos, cov_neg=cov_neg)

def paf_to_npz(gen):
# PAF format: col 1 = read name, 2 = read length, 3 = read start, 4 = read end, 5 = strand, 6 = query name
# 7 = mapping length, 8 = mapping start, 9 = mapping end, 10 = nb of matching bases in the mapping
# 11 = nb of bases including gaps in the mapping, 12 = mapping quality (0-255 with 255=missing)
    d = {"+":[[],[]], "-":[[],[]]}
    fname = "RNAmRNA.paf"
    path = "/home/raphael/Documents/topo/RNAseq_data/2021_nanopore/raw_reads/" + fname
    with open(path,"r") as f:
        for line in f:
            line=line.strip().split('\t')
            d[line[4]][0].append(int(line[7]))
            d[line[4]][1].append(int(line[8]))
    f.close()   

    Rpos_start = np.array(d["+"][0],dtype=int)
    Rpos_end = np.array(d["+"][1],dtype=int)
    Rneg_start = np.array(d["-"][1],dtype=int)
    Rneg_end = np.array(d["-"][0],dtype=int)

    Rpos = np.column_stack((Rpos_start,Rpos_end))
    Rneg = np.column_stack((Rneg_start,Rneg_end))

    # delete rows where one coordinate is missing
    # Rpos = Rpos[np.isfinite(Rpos).all(axis=1)] 
    # Rneg = Rneg[np.isfinite(Rneg).all(axis=1)]
    print(np.shape(Rpos),np.shape(Rneg))

    # if reads.info not exist
    if not os.path.exists(basedir+"data/"+gen.name+'/rnaseq_reads/reads.info'):
        file = open(basedir+"data/"+gen.name+'/rnaseq_reads/reads.info','w') 
        file.write('Condition\tReads file\tDate\tBAM file')
        file.close() 
    # save results ; .npz contains two .npy Rpos and Rneg
    file = open(basedir+"data/"+gen.name+'/rnaseq_reads/reads.info','a')
    file.write('\n'+fname[:-4]+'\t'+fname[:-4]+'_reads.npz\t'+str(datetime.now())+'\t'+fname)  
    file.close()    
    np.savez(basedir+"data/"+gen.name+'/rnaseq_reads/'+fname[:-4]+'_reads.npz', Rpos=Rpos, Rneg=Rneg)

##### MAIN #####
# for a list of .bam
# samples=["E%d"%x for x in range(1,15)]+["F%d"%x for x in range(1,9)+[13,14]]

# for s in samples:
#     process_bam_paired_end('%s.bam'%s)
#     cov_from_reads('%s_reads.npz'%s, 4922802) # dickeya genome length    
########## BAM2NPZ RNA_SEQ ANALYSIS, RAPHAEL ##########

# compute coverage from reads (.npz, one .npy per strand)
def cov_start_stop_from_reads(gen):
    if not hasattr(gen, 'seq'):
        gen.load_seq()
    genome_length = len(gen.seq)

    with open(basedir+"data/"+gen.name+"/rnaseq_reads/reads.info","r") as f:
        header = next(f)
        for line in f: # for each condition
            line=line.strip().split('\t')
            print('Loading condition',line[0])   
            # load npz file corresponding to condition
            npzfile = np.load(basedir+"data/"+gen.name+'/rnaseq_reads/'+line[1])
            Rpos = npzfile["Rpos"]
            Rneg = npzfile["Rneg"]
            # init cov
            cov_start = {0:np.zeros(genome_length, dtype=int), 1:np.zeros(genome_length, dtype=int)}
            cov_end = {0:np.zeros(genome_length, dtype=int), 1:np.zeros(genome_length, dtype=int)}
            # compute cov
            # on positive strand
            for start,end in Rpos:
                try:
                    cov_start[1][start-1] += 1
                    cov_end[1][end-1] += 1
                except:
                    pass
            # on negative strand
            for start,end in Rneg:
                try:
                    cov_start[0][start-1] += 1
                    cov_end[0][end-1] += 1
                except:
                    pass
            # if rnaseq_cov folder not exist        
            # if cov.info not exist
            if not os.path.exists(basedir+"data/"+gen.name+'/rnaseq_cov/cov_start_stop.info'):
                file = open(basedir+"data/"+gen.name+'/rnaseq_cov/cov_start_stop.info','w') 
                file.write('Condition\tCov file\tDate\tReads file')
                file.close() 
            # save results
            file = open(basedir+"data/"+gen.name+'/rnaseq_cov/cov_start_stop.info','a')
            fname = line[0] + "_cov_start_end"
            file.write('\n{}\t{}.npz\t{}\t{}'.format(line[0],fname,str(datetime.now()),line[1]))  
            file.close()
            # save results ; .npz contains four .npy: start on pos, stop on pos, start on neg, stop on neg
            np.savez(basedir+"data/"+gen.name+'/rnaseq_cov/'+fname+'.npz', cov_start_pos= cov_start[1],
                cov_end_pos= cov_end[1], 
                cov_start_neg= cov_start[0], 
                cov_end_neg= cov_end[0])