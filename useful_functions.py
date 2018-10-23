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


########## BAM2NPZ RNA_SEQ ANALYSIS, RAPHAEL ##########
import os,pysam, numpy as np, pandas as pd
from datetime import datetime

def process_bam_paired_end(bam_file): # /!\ paired-end only /!\ -> return fragments for each strand
    if not os.path.exists(bam_file+".bai"): # if index needed, created using samtools (.bai file)
        os.system("samtools index -b %s"%bam_file)
    
    bamfile = pysam.AlignmentFile(bam_file, "rb") # BAM opening, alignment file object
    print 'Header :',bamfile.header
    print 'Reads mapped :',bamfile.mapped
    print 'Reads without coordinate :',bamfile.nocoordinate
    print 'Reads not mapped :',bamfile.unmapped
    # /!\ 0-based coordinate system /!\
    # fastest way to build a numpy matrix -> store every coordinate of interest into lists, then merge into numpy array
    # lists storing coordinates of paired-end fragments for +/- strands
    Rpos_start = []
    Rpos_end = []
    Rneg_start = []
    Rneg_end = []
    # per default, only mapped reads to reference fetched
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

    # if rnaseq_reads folder not exist
    if not os.path.exists(os.getcwd()+'/rnaseq_reads'):
        os.makedirs('rnaseq_reads')
    # if reads.info not exist
    if not os.path.exists(os.getcwd()+'/rnaseq_reads/reads.info'):
        file = open(os.getcwd()+'/rnaseq_reads/reads.info','w') 
        file.write('Condition\tReads file\tDate\tBAM file')
        file.close() 
    # save results ; .npz contains two .npy Rpos and Rneg
    file = open(os.getcwd()+'/rnaseq_reads/reads.info','a')
    file.write('\n'+bam_file[:-4]+'\t'+bam_file[:-4]+'_reads.npz\t'+str(datetime.now())+'\t'+bam_file)  
    file.close()    
    np.savez(os.getcwd()+'/rnaseq_reads/'+bam_file[:-4]+'_reads.npz', Rpos=Rpos, Rneg=Rneg)


# compute coverage from reads (.npz, one .npy per strand)
def cov_from_reads(npz_file, genome_length):
    # load npz file
    npzfile = np.load(os.getcwd()+'/rnaseq_reads/'+npz_file)
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
    if not os.path.exists(os.getcwd()+'/rnaseq_cov'):
        os.makedirs('rnaseq_cov')
    # if cov.info not exist
    if not os.path.exists(os.getcwd()+'/rnaseq_cov/cov.info'):
        file = open(os.getcwd()+'/rnaseq_cov/cov.info','w') 
        file.write('Condition\tCov file\tDate\tReads file')
        file.close() 
    # save results
    file = open(os.getcwd()+'/rnaseq_cov/cov.info','a')
    file.write('\n'+npz_file[0:-10]+'\t'+npz_file[:-4]+'_cov.npz\t'+str(datetime.now())+'\t'+npz_file)  
    file.close()
    # save results ; .npz contains two .npy cov_pos and cov_neg
    np.savez(os.getcwd()+'/rnaseq_cov/'+npz_file[:-4]+'_cov.npz', cov_pos=cov_pos, cov_neg=cov_neg)

##### MAIN #####
# for a list of .bam
# samples=["E%d"%x for x in range(1,15)]+["F%d"%x for x in range(1,9)+[13,14]]

# for s in samples:
#     process_bam_paired_end('%s.bam'%s)
#     cov_from_reads('%s_reads.npz'%s, 4922802) # dickeya genome length    
########## BAM2NPZ RNA_SEQ ANALYSIS, RAPHAEL ##########


"""
Vincent CABELI
"""

import sys
import os
import numpy as np
import scipy.spatial.distance
from scipy.optimize import curve_fit
from numpy.lib.stride_tricks import as_strided
from random import shuffle
from Domain import Domain


def periodic_dist(genes, max_size):
    """ Computes the periodic distance matrix of points on a circle of size
    max_size. 
    """
    pos_left = np.transpose([[x.left for x in genes]])
    pos_right = np.transpose([[x.right for x in genes]])
    pos_middle = np.transpose([[x.middle for x in genes]])
    distance_matrix = scipy.spatial.distance_matrix(pos_middle, pos_middle)
    #lambda u,v
    distance_matrix = np.where(distance_matrix > 0.5*max_size, 
        max_size-distance_matrix, distance_matrix)
    return distance_matrix


def change_dict_keys(genes_dict):
    """ Change a genes dictionnary's keys to the gene names.
    """
    for tag in list(genes_dict.keys()):
        gene = genes_dict[tag]
        genes_dict[gene.name] = gene
        del genes_dict[tag]
    return genes_dict


def gen_successive_domains(genes_dict, k):
    """ Generates overlapping domains of adjacent genes on the whole genome. It
    will return a list of n_genes-k domains.
    """
    sorted_genes = list(genes_dict.values())
    sorted_genes.sort(key = lambda x: x.start)
    domains = []
    for i in range(len(sorted_genes)-k):
        if all([hasattr(gene, 'mean_expression') for gene in sorted_genes[i:i+k]]):
            domains.append(Domain(sorted_genes[i:i+k]))
    return domains


def gen_random_domains(sorted_genes, n, domains=None, b=0.575):
    """ Generates a list of n random domains of adjacent genes. Uses an
    exponential law to determine the length of each domain, with the formula 
    for exponential decay exp(-b * x). Default = 0.575 as observed in the 
    synteny segments size distribution. Will determine a new b parameter by
    fitting an exponential curve if a domains list is specified.
    """
    n_genes = len(sorted_genes)
    # consider max 30 genes domains
    max_size = 30
    if domains:
        t = [len([domain for domain in domains if 
                  domain.number_of_genes==i]) for i in range(2, max_size)]
        popt, pcov = curve_fit(exp_to_fit, np.asarray(list(range(2, max_size))), t)
        b = popt[1]
    probs = np.exp(-b * np.asarray(list(range(max_size))))
    rand_domains = []
    rand_domains_ext = []
    for i in range(n):
        gene_pos = int(np.random.uniform(1,n_genes-max_size))
        size = np.max(np.where(np.random.rand() < probs)) + 2
        rand_domains.append(Domain(sorted_genes[gene_pos:gene_pos+size]))
        rand_domains_ext.append(Domain(sorted_genes[gene_pos-1:gene_pos+size+1]))
    return (rand_domains, rand_domains_ext)


def gen_random_segment_like_domains(sorted_genes, n, maps, domains,
                                    successions_types=None, 
                                    successions_types_counts=None):
    """ Generates a list of len(domains) random domains with the same number of
    convergent, divergent and same sense gene successions as those observed in
    the given list of domains. 
    """
    rand_domains = []
    rand_domains_ext = []
    if successions_types and successions_types_counts:
        for i, successions_type in enumerate(successions_types):
            n_genes = sum(successions_type)+1
            candidates = np.where([(list(row)==successions_type) for row in maps[n_genes-2]])[0]
            positions = np.random.choice(candidates, size=successions_types_counts[i])
            for pos in positions:
                rand_domains.append(Domain(sorted_genes[pos:pos+n_genes]))
                rand_domains_ext.append(Domain(sorted_genes[pos-1:pos+n_genes+1]))
        return (rand_domains, rand_domains_ext)
    else:
        for domain in domains:
            domain_successions = [count_successions([domain], np.asarray([True,True])) +
                                  count_successions([domain], np.asarray([False,False])),
                                  count_successions([domain], np.asarray([False,True])),
                                  count_successions([domain], np.asarray([True,False]))]
            candidates = np.where([(list(row)==domain_successions) for row in maps[domain.number_of_genes-2]])[0]
            pos = np.random.choice(candidates)
            rand_domains.append(Domain(sorted_genes[pos:pos+domain.number_of_genes]))
            rand_domains_ext.append(Domain(sorted_genes[pos-1:pos+domain.number_of_genes+1]))
        return (rand_domains, rand_domains_ext)


def count_successions_types(domains):
    seen_successions = []
    seen_successions_count = []
    for domain in domains:
        domain_successions = [count_successions([domain], np.asarray([True,True])) +
                              count_successions([domain], np.asarray([False,False])),
                              count_successions([domain], np.asarray([False,True])),
                              count_successions([domain], np.asarray([True,False]))]
        if domain_successions not in seen_successions:
            seen_successions.append(domain_successions)
            seen_successions_count.append(1)
        else:
            seen_successions_count[seen_successions.index(domain_successions)] += 1
    return (seen_successions, seen_successions_count)


def convergent_significance(genes_dict, domains, n):
    """ Generates n times random domains list of the same size as the list
    passed as argument. Returns a tuple with the number of observed convergences
    as well as the number of convergences observed only in direct domain 
    neighbours.
    """
    sorted_genes = list(genes_dict.values())
    sorted_genes.sort(key = lambda x: x.left)
    n_convergents = []
    diff = []
    # fit exp and get parameter b
    max_size = max([x.number_of_genes for x in domains])
    print('Fitting exponential...')
    domains_sizes = [len([domain for domain in domains if 
                     domain.number_of_genes==i]) for i in range(2, max_size)]
    popt, pcov = curve_fit(exp_to_fit, np.asarray(list(range(2, max_size))), domains_sizes)
    b = popt[1]
    print('Generating maps...')
    maps = gen_maps(sorted_genes, max_size)
    n_domains = len(domains)
    print('Generating domains...')
    (successions_types, successions_types_counts) = count_successions_types(domains)
    for i in range(n):
        print(i)
        #(d, d_ext) = gen_random_domains(sorted_genes, n_domains, b=b)
        (d, d_ext) = gen_random_segment_like_domains(sorted_genes, 1, maps, 
                                                     domains, 
                                                     successions_types,
                                                     successions_types_counts)
        n_convergents.append(count_successions(d, np.array([False, True])))
        diff.append(count_successions(d_ext, np.array([False,True])) - n_convergents[i])
    return (n_convergents, diff)


def count_successions(domains, succession):
    """ Counts the successions of 
    """
    n = 0
    b = np.size(succession)
    for domain in domains:
        for i in range(domain.number_of_genes-1):
            if np.all(succession == domain.genes_orientations[i:i+b]):
                n += 1
    return n
        

def exp_to_fit(x, a, b):     
    return a * np.exp(-b * x)


def gen_maps(sorted_genes, max_size):
    """ Generates a [max_size-1, len(sorted_genes)-max_size, 3]-shaped matrix 
    that contains a map of
    """
    maps = np.zeros((max_size, len(sorted_genes)-max_size, 3), dtype='int')
    for size in range(max_size):
        # allow for a 1-gene displacement for extended domains
        for pos in range(1, len(sorted_genes)-max_size-1):
            domain = Domain(sorted_genes[pos:pos+size+2])
            maps[size,pos,:] = [count_successions([domain], np.asarray([True,True])) + 
                                count_successions([domain], np.asarray([False,False])), 
                                count_successions([domain], np.asarray([False,True])),
                                count_successions([domain], np.asarray([True,False]))]
    return maps


def compare_domains_lists(domains1, domains2):
    """ Compares two domains lists and returns
    """
    max_size = np.max([x.number_of_genes for x in domains1] + 
                      [x.number_of_genes for x in domains2])
    equalities = np.zeros(max_size, dtype='int')
    for domain1 in domains1:
        for domain2 in domains2:
            if domain1 == domain2 or domain1.includes(domain2) or domain2.includes(domain1):
                equalities[domain1.number_of_genes-1] += 1
                break
    return equalities


def build_domains_list_from_matrix(genes_dict, domains_matrix):
    """ Builds a list of Domain objects from a square matrix that is returned
    by the Dickeya_parser.domains_matrix() function. 
    """
    if (len(np.shape(domains_matrix)) != 2) or (
        np.shape(domains_matrix)[0] != np.shape(domains_matrix)[1]):
        print('The matrix is not a square matrix, cannot build list of domains')
        return

    genes = list(genes_dict.values())
    domains = []
    i = 0
    j = 0
    current_domain = []
    seen_indices = []
    for i in range(np.shape(domains_matrix)[1]):
        if i in seen_indices:
            continue
        for j in range(i, np.shape(domains_matrix)[1]):
            if domains_matrix[i,j]:
                current_domain.append(genes[j])
                seen_indices.append(j)
        if current_domain:
            domains.append(Domain(current_domain))
            current_domain = []
    return domains


# ------------


from itertools import groupby

def cov_is_zero(cov,winlen=20,maxval=5):
    """ cov analysis: extract positions where cov is zero
    criterium: maximum value of winlen consecutive basepairs
    """
    wherezero = cov<=maxval
    a=[(k,len(list(g))) for k,g in groupby(wherezero)]
    q=np.zeros(len(wherezero),dtype=bool)
    index=0
    for (k,le) in a:
        if k==True and le>=winlen:
            q[index:index+le]=True
        index+=le
    return q


"""
Vincent CABELI
"""
    