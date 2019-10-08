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
from matplotlib_venn import venn2

from scipy import stats
from statsmodels.stats import weightstats

plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})
plt.rcParams.update({'font.size': 11})
plt.rcParams.update({'legend.fontsize': 9})
plt.rcParams.update({'font.family': "Arial"})

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





########## BAM2NPZ RNA_SEQ ANALYSIS, RAPHAEL ##########
import os,pysam, numpy as np, pandas as pd
from datetime import datetime

def process_bam_paired_end(bam_file, gen): # /!\ paired-end only /!\ -> return fragments for each strand
    if not os.path.exists(basedir+"data/"+gen+'/rnaseq_reads/'+bam_file+".bai"): # if index needed, created using samtools (.bai file)
        os.system("samtools index -b %s"%bam_file)
    
    bamfile = pysam.AlignmentFile(basedir+"data/"+gen+'/rnaseq_reads/'+bam_file, "rb") # BAM opening, alignment file object
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







def compare_fc_conds(self,c1,c2,*args,**kwargs):
    self.load_fc_pval()
    thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
    thresh_fc = kwargs.get('thresh_fc', 0) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated
    res = {}
    states = ['+*','+','-*','-']
    for st1 in states:
        res[st1] = {}
        for st2 in states:
            res[st1][st2] = 0

    for gene in self.genes.keys():
        try:
            g = self.genes[gene]
            if g.fc_pval[c1][1] <= thresh_pval:
                stat1 = '*'
            else:
                stat1 = ''
            if g.fc_pval[c1][0] > thresh_fc:
                sign1 = '+'
            else:
                sign1 = '-'


            if g.fc_pval[c2][1] <= thresh_pval:
                stat2 = '*'
            else:
                stat2 = ''
            if g.fc_pval[c2][0] > thresh_fc:
                sign2 = '+'
            else:
                sign2 = '-'

            st1 = sign1+stat1 ; st2 = sign2+stat2
            res[st1][st2] += 1
        except:
            pass

    df = pd.DataFrame.from_dict(res,orient='index')
    df.sort_index(axis=1,inplace=True) ; df.sort_index(axis=0,inplace=True)
    df.to_csv('/home/raphael/Documents/test.csv')
    

def compare_genomes(gen1,gen2,*args,**kwargs):
    '''
    Compare genes in common among two genomes based on name, and compare roughly promoter sequences
    '''
    gen1.load_annotation()
    gen1.load_seq()

    gen2.load_annotation()
    gen2.load_seq()

    prom = kwargs.get('prom',250)
    res1 = {}
    res2 = {}
    for gene in gen1.genes.keys():
        g = gen1.genes[gene]
        if g.strand:
            seq = gen1.seq[g.start-prom:g.start]
        else:
            seq = gen1.seqcompl[g.start:g.start+prom]
        res1[g.name] = seq
          
    for gene in gen2.genes.keys():
        g = gen2.genes[gene]
        if g.strand:
            seq = gen2.seq[g.start-prom:g.start]
        else:
            seq = gen2.seqcompl[g.start:g.start+prom]
        res2[g.name] = seq
    print res1
    print res2
    print len(res1.keys()),len(res2.keys())

    res = {}
    for g1 in res1.keys():
        if g1 in res2.keys():
            mismatches = sum(c1!=c2 for c1,c2 in zip(res1[g1],res2[g1]))
            if mismatches not in res.keys():
                res[mismatches] = 1
            else:
                res[mismatches] += 1
    print res
    print sum(res[c] for c in res.keys())


def load_TSS_from_common_genome(gen,genref,*arg,**kwargs):
    prom_region = kwargs.get('prom',False)
    genref.compute_magic_prom(prom=prom_region)
    gen.load_annotation()
    gen.load_seq()
    check = kwargs.get('check',200)
    thresh = kwargs.get('thresh',200)
    gnames = {} # dict {gene_name:gene_id} to map names to ID
    gen.TSSs = {} # init TSS
    for g in gen.genes.keys():
        gnames[gen.genes[g].name] = g

    for TSScond in genref.TSSs.keys():
        try:

            gen.TSSs[TSScond] = genref.TSSs[TSScond] # init TSS cond
            for pos in gen.TSSs[TSScond].keys(): # for each TSS
                TSS = gen.TSSs[TSScond][pos]
                newgenes = []
                for gene in TSS.genes:
                    try:
                        gref = genref.genes[gene]
                        g = gen.genes[gnames[gref.name]]
                        if g.strand and gref.strand:
                            s = gen.seq[g.start-check-1:g.start]
                            sref = genref.seq[gref.start-check-1:gref.start]
                        elif not g.strand and not gref.strand:
                            s = gen.seqcompl[g.start-1:g.start+check][::-1]
                            sref = genref.seqcompl[gref.start-1:gref.start+check][::-1]
                        if sum(c1!=c2 for c1,c2 in zip(s,sref)) <= thresh:
                            newgenes.append(gnames[gref.name])

                    except:
                        pass
                TSS.genes = newgenes
        except:
            pass

def load_fc_from_common_genome(genTSS,genFC,*arg,**kwargs):
    genFC.load_fc_pval() # genome from which we want to extrapolate FC (ecoli_b)
    genTSS.load_fc_pval() # genome where we want to extrapolate FC (ecoli)
    check = kwargs.get('check',200) # nb of bp to test upstream gene start to evaluate gene similarity in genomes
    thresh = kwargs.get('thresh',10) # if nb of different nucleotides in the genome > thresh for a given gene, no extrapolation (5% differences max)
    
    if not hasattr(genTSS, 'seq'):
        genTSS.load_seq()
    if not hasattr(genFC, 'seq'):
        genFC.load_seq()

    gnames = {} # dict {gene_name:gene_id} to map names to ID for genome
    for g in genFC.genes.keys():
        gnames[genFC.genes[g].name] = g

    for cond in genFC.genes_valid.keys(): # prepare conditions to be loaded in reference genome
        genTSS.genes_valid[cond] = []

    d = {}
    for idref in genTSS.genes.keys(): # for each gene of reference genome
    # test whether or not upstream regions of gene start are similar between the two genomes
        try:
            gref = genTSS.genes[idref]
            g = genFC.genes[gnames[gref.name]] # works if the gene name in genTSS exists in gen 
            if g.strand and gref.strand:
                s = genFC.seq[g.start-check-1:g.start]
                sref = genTSS.seq[gref.start-check-1:gref.start]
            elif not g.strand and not gref.strand:
                s = genFC.seqcompl[g.start-1:g.start+check][::-1]
                sref = genTSS.seqcompl[gref.start-1:gref.start+check][::-1]

            s = sum(c1!=c2 for c1,c2 in zip(s,sref))
            try:
                d[s] += 1
            except:
                d[s] = 1
            if s <= thresh:
                gref.fc_pval = g.fc_pval
                for cond in gref.fc_pval.keys():
                    genTSS.genes_valid[cond].append(idref)
        except:
            pass
    print d

def barplot_annotate_brackets(num1, num2, text, center, height, yerr=None, dh=.05, barh=.05, fs=12, maxasterix=None, dt=0, bold = False):
    """ 
    Annotate barplot with p-values.

    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """
    lx, ly = center[num1], height[num1]
    rx, ry = center[num2], height[num2]
    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh *= (ax_y1 - ax_y0)
    barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh)

    plt.plot(barx, bary, c='black')
    if bold:
        plt.text(mid[0],mid[1]+dt, text, fontsize=fs, ha='center', fontweight = "bold")
    else:
        plt.text(mid[0],mid[1]+dt, text, fontsize=fs, ha='center')



def significance(pval):
    if pval <= 0.001:
        s = '***'
    elif pval <= 0.01:
        s = '**' 
    elif pval <= 0.05:
        s = '*' 
    else:
        s = 'ns'
    return s


def fill_cov_file(gen):
    gen.load_seq()
    with open(basedir+"data/"+gen.name+"/rnaseq_cov/cov_txt.info","r") as file:
        header = next(file)
        for line in file:
            line = line.strip('\n').split('\t')
            print 'Loading condition:',line[0]
            cov_neg = {}
            cov_pos = {}
            with open(basedir+"data/"+gen.name+"/rnaseq_cov/"+line[1], 'r') as f:
                i = 1
                while i < int(line[3]):
                    header=next(f)
                    i+=1
                for l in f:
                    l = l.strip('\n').split('\t')
                    cov_neg[int(l[0])] = float(l[1])
            f.close()
            fnew = open(basedir+"data/"+gen.name+"/rnaseq_cov/"+line[1][:-4]+'_treated.wig', 'w')
            for i in range(1,len(gen.seq)+1):
                try:
                    cov = cov_neg[i]
                except:
                    cov = 0
                fnew.write('{}\t{}\n'.format(str(i),str(cov)))
            fnew.close()

            with open(basedir+"data/"+gen.name+"/rnaseq_cov/"+line[2], 'r') as f:
                i = 1
                while i < int(line[3]):
                    header=next(f)
                    i+=1
                for l in f:
                    l = l.strip('\n').split('\t')
                    cov_pos[int(l[0])] = float(l[1])
            f.close()
            fnew = open(basedir+"data/"+gen.name+"/rnaseq_cov/"+line[2][:-4]+'_treated.wig', 'w')
            for i in range(1,len(gen.seq)+1):
                try:
                    cov = cov_pos[i]
                except:
                    cov = 0
                fnew.write('{}\t{}\n'.format(str(i),str(cov)))
            fnew.close()

def export_rpkm(gen):
    gen.com
    res = {}
    for gene in gen.genes.keys():
        try:
            res[gene] = gen.genes[gene].rpkm
        except:
            pass

        df = pd.DataFrame.from_dict(res,orient='index', dtype=float) ; df.sort_index(axis=1,inplace=True)
        df = df.applymap(np.log2)
        df.to_csv(basedir+"data/"+gen.name+"/expression/log2rpkm.csv") 


def write_FC_from_rpkm(gen):
    gen.load_cov()
    gen.compute_rpkm_from_cov()
    f = open(basedir+"data/"+gen.name+"/expression/seco_rnaseq_clean_FC.csv", 'w')
    f.write('gene\tcold_shock_10mn\tcold_shock_30mn\n')
    for gene in gen.genes.keys():
        g = gen.genes[gene]
        try:
            c1 = g.rpkm['WT'] 
            c2 = g.rpkm['WT2'] 
            c = np.mean([c1,c2])

            t1 = g.rpkm['cold_shock_10mn']
            fc1 = np.log2(t1/c)
            
            t2 = g.rpkm['cold_shock_30mn']
            fc2 = np.log2(t2/c)

            f.write('{}\t{}\t{}\n'.format(gene,str(fc1),str(fc2)))
        except Exception as e:

            pass
    f.close()


def load_expression_bartholomaus(gen,*arg,**kwargs):
    gen.load_annotation()
    gnames = {} # dict {gene_name:gene_id} to map names to ID
    for g in gen.genes.keys():
        gnames[gen.genes[g].name] = g

    path = basedir+"data/"+gen.name+"/expression/"+"jozefczuk_cold-heat-ox/"
    for file in os.listdir(path):
        df = pd.read_table(path+file,sep="\t")
        def replace_name(df,d):
            '''
            Clean DF : for each row, convert gene names to bnames
            '''
            def cleandf(row,row_accumulator,d):
                old_row = row.to_dict()
                new_row = row.to_dict()
                try:           
                    g = d[new_row['ORF']]
                    new_row['ORF'] = g
                    row_accumulator.append(new_row)
                except:
                    pass
            
            new_rows = []
            df.apply(cleandf,axis=1,args=(new_rows,d))
            new_df = pd.DataFrame(new_rows)
            return new_df

        newdf = replace_name(df,gnames)
        newdf.to_csv(path+file[:-4]+"_clean.csv",sep='\t',header=True,index=False)
    #     with open(path+file, 'r') as f:
    #         header = next(f).strip().split('\t')
    #         gen.genes_valid[header[2]] = []
    #         for line in f:
    #             line=line.strip().split('\t')
    #             try:
    #                 if not hasattr(gen.genes[gnames[line[0]]],'fc_pval'):
    #                     gen.genes[gnames[line[0]]].fc_pval={}
    #                 gen.genes[gnames[line[0]]].fc_pval[header[2]] = (math.log(float(line[2])/float(line[1]),2),0)
    #                 gen
    #                 gen.genes[gnames[line[0]]].fc_pval[header[3]] = (math.log(float(line[3])/float(line[1]),2),0)
    #             except Exception as e:
    #                 print e
    #                 pass
    # gen.load_fc_pval()
    # gnames = {} # dict {gene_name:gene_id} to map names to ID
    # for g in gen.genes.keys():
    #     gnames[gen.genes[g].name] = g
    # file = kwargs.get("file","stress_osmotic_heat.csv")


def rpkm_from_raw_counts(gen,*arg,**kwargs):
    file = kwargs.get('file','seco_rnaseq.csv')
    gen.load_annotation()
    df = pd.read_table(basedir+"data/"+gen.name+"/expression/"+file,sep="\t")
    def replace_name(df,d):
        '''
        Clean DF : for each row, convert gene names to bnames
        '''
        def cleandf(row,row_accumulator,d):
            new_row = row.to_dict()
            try:
                new_row['length'] = gen.genes[new_row['Gene']].length
                row_accumulator.append(new_row)
            except:
                pass
        new_rows = []
        df.apply(cleandf,axis=1,args=(new_rows,d))
        new_df = pd.DataFrame(new_rows)
        return new_df

    newdf = replace_name(df,{})
    newdf.to_csv(basedir+"data/"+gen.name+"/expression/"+file[:-4]+"_clean.csv",sep='\t',header=True,index=False)


def check_TSS_strand_add_genes(gen,cond):
    f = open(basedir+"data/"+gen.name+"/TSS/checked.txt",'w')
    f.write("TSS\tStrand\tGenes\n")
    for TSSpos in gen.TSSs[cond].keys(): # for all TSS in cond_tss
        TSSu = gen.TSSs[cond][TSSpos] # single TS
        genes = []
        for gene in gen.genes.keys():
            g = gen.genes[gene]
            if TSSu.strand:
                if g.start - TSSpos <= 200 and  g.start - TSSpos > 1:
                    genes.append(gene)
            else:
                if TSSpos - g.start <= 200 and  TSSpos - g.start > 1:
                    genes.append(gene)             

        f.write("{}\t{}\t{}\n".format(TSSpos,TSSu.strand,','.join(genes)))
    f.close()

def extract_prom(gen,cond):
    gen.compute_magic_prom()
    if not hasattr(gen,'genes_valid'):
        gen.load_fc_pval()
    f = open(basedir+"data/"+gen.name+"/TSS/promoters.txt",'w')
    f.write("TSS\tPromoter\n")
    for TSS in gen.TSSs[cond].keys():
        TSSu = gen.TSSs[cond][TSS]
        try:
            f.write("{}\t{}\n".format(str(TSS),TSSu.promoter[sigfactor]['discr_model']))
        except:
            pass
    f.close()








def ci_mean_binomial(obs,n):
    p = float(obs) / n
    sd = float(stats.binom.std(n, p, loc=0))
    sdmean = sd/np.sqrt(n)
    cimean = obs - weightstats._tconfint_generic(obs,sdmean,n-1,0.05,"two-sided")[0]
    return p,sd,sdmean, cimean

def ci_mean_std(mean,std,n):
    sdmean = std/np.sqrt(n)  
    cimean = n - weightstats._tconfint_generic(mean,sdmean,n-1,0.05,"two-sided")[0]
    return cimean

# function for calculating the t-test for two independent samples
def independent_ttest(m1,m2,std1,std2,n1,n2,alpha):
    # standard error on the difference between the samples
    sed = np.sqrt(std1**2.0 + std2**2.0)
    # calculate the t statistic
    t_stat = (m1 - m2) / sed
    # degrees of freedom
    df = n1 + n2 - 2
    # calculate the critical value
    cv = stats.t.ppf(1.0 - alpha, df)
    # calculate the p-value
    p = (1.0 - stats.t.cdf(abs(t_stat), df)) * 2.0
    # return everything
    return t_stat, p






def student_mean_ci(data):
    mean = np.mean(data)
    # evaluate sample variance by setting delta degrees of freedom (ddof) to
    # 1. The degree used in calculations is N - ddof
    stddev = std(data, ddof=1)
    # Get the endpoints of the range that contains 95% of the distribution
    t_bounds = stats.t.interval(0.95, len(data) - 1)
    # sum mean to the confidence interval
    ci = [mean + critval * stddev / sqrt(len(data)) for critval in t_bounds]
    return mean, (ci[0],ci[1])



def clean_operon_file(gen):
    gen.load_annotation()
    path = basedir+"data/"+gen.name+"/operon_raph_raw"
    d = {} ; operons = [["Start","Stop","Strand","Genes"]]
    for g in gen.genes.keys():
        d[gen.genes[g].start] = g

    if os.path.exists(path):
        with open(path,"r") as f:
            header = next(f)
            for line in f:
                try:
                    line = line.strip('\n')
                    line=line.split('\t')
                    genes = line[1].split(",")
                    starts = [] ; stops = [] ; strands = []
                    for gene in genes:
                        g = gen.genes[gene]
                        starts.append(g.start) ; stops.append(g.end) ; strands.append(g.strand)
                    starts = sorted(starts) ; stops = sorted(stops)
                    if all(item == True for item in strands):
                        operons.append([starts[0],stops[-1],True,[d[st] for st in starts]])
                    elif all(item == False for item in strands):
                        operons.append([starts[-1],stops[0],False,[d[st] for st in starts][::-1]])
                except Exception as e:
                    print e
                    pass
    df = pd.DataFrame(operons)
    df.to_csv(basedir+"data/"+gen.name+"/operon_raph_clean.csv",header=True,sep="\t")    


def correct_fasta(gen,res):
    print gen.name
    for gene in gen.genes.keys():
        g = gen.genes[gene]
        if g.name == "mioC" or g.name == "mnmG":
            print gene, g.name, g.start
            res[gene] = (g.name,g.start)

    if res["mioC"][1] > 1000 and res["mioC"][1] > res["mnmG"][1]:
        newstart = res["mioC"][1] ; newstop = res["mnmG"][1]


def load_annot_detailed(gen, *arg,**kwargs):
    if not hasattr(gen,'genes'):
        gen.load_annotation()

    with open(basedir+"data/"+gen.name+"/annotation/annot_detailed.info","r") as f:
        skiphead = next(f) # skip head        
        for line in f:
            line=line.strip()
            line=line.split('\t')
            locuscol = int(line[0]) ; ordercol = int(line[1]) ; operoncol = int(line[2])
            snamecol = int(line[3]) ; domaincol = int(line[4]) ; LScol = int(line[5])
            replichorecol = int(line[6]) ; myclassifcol = int(line[7]) ; mysubclasscol = int(line[8])
            mylvlcol = int(line[9])
        f.close()    
    
    with open(basedir+"data/"+gen.name+"/annotation/annot_detailed.csv","r") as f:
        skiphead = next(f) # skip head        
        for line in f:
            line=line.strip()
            line=line.split('\t')
            try:
                gen.genes[line[locuscol]].order = int(line[ordercol])
                gen.genes[line[locuscol]].operon = int(line[operoncol])
                gen.genes[line[locuscol]].sname = line[snamecol]
                gen.genes[line[locuscol]].domain = int(line[domaincol])
                gen.genes[line[locuscol]].leading_strand = int(line[LScol])
                gen.genes[line[locuscol]].replichore = line[replichorecol]
                gen.genes[line[locuscol]].myclassif = line[myclassifcol]
                gen.genes[line[locuscol]].mysubclass = line[mysubclasscol]
                gen.genes[line[locuscol]].mylvl = line[mylvlcol]
                gen.genes[line[locuscol]].locus_tag = line[locuscol]
            except Exception as e:
                pass

        f.close()  
    gen.annot_detailed = True

def export_annot_detailed(gen, l, name, *arg,**kwargs):
    res = {}
    for g in l:
        res[g] = gen.genes[g].__dict__

    df = pd.DataFrame.from_dict(res, orient='index')
    df.sort_index(inplace=True)
    df = df[["locus_tag","sname","left","right","strand","replichore","leading_strand","order","domain","myclassif","mysubclass","mylvl"]]
    df.to_csv(basedir+"data/"+gen.name+"/annotation/"+name+".csv",sep='\t',encoding='utf-8', index=False)






def write_genes_fc(gen,*args,**kwargs):
    gen.load_fc_pval()
    thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
    thresh_fc = kwargs.get('thresh_fc', 1) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated 
    # pathdb = "/home/raphael/Documents/GO_ihfA/"
    # if not os.path.exists(pathdb):
    #   os.makedirs(pathdb)

    res = {}
    for cond_fc in gen.genes_valid.keys():
        res[cond_fc] = {'act':[],'rep':[],'non':[]}
        for gene in gen.genes_valid[cond_fc]:
            g = gen.genes[gene]
            if g.fc_pval[cond_fc][1] <= thresh_pval:
                if g.fc_pval[cond_fc][0] <  0 - thresh_fc:
                    res[cond_fc]['rep'].append(gene)
                elif g.fc_pval[cond_fc][0] >  0 + thresh_fc:
                    res[cond_fc]['act'].append(gene)
                else:
                    res[cond_fc]['non'].append(gene)
            else:
                    res[cond_fc]['non'].append(gene)

    groups = [["WT_expo_vs_ihfA_expo", "WT_stat_vs_ihfA_stat", "WT_PGA_stat_vs_ihfA_PGA_stat"],
    ["WT_nov_expo_vs_ihfA_nov_expo","WT_nov_stat_vs_ihfA_nov_stat","WT_PGA_nov_stat_vs_ihfA_PGA_nov_stat"],
    ["WT_nov_expo_vs_WT_expo","WT_nov_stat_vs_WT_stat","WT_PGA_nov_stat_vs_WT_PGA_stat"],
    ["ihfA_nov_expo_vs_ihfA_expo","ihfA_nov_stat_vs_ihfA_stat","ihfA_PGA_nov_stat_vs_ihfA_PGA_stat"]]

    lgd = ["ihfA responsive","ihfA responsive under relaxation","relaxation responsive","relaxation responsive without ihfA"]

    reslist = {}
    i = 0
    for group in groups:
        A = res[group[0]] ; B = res[group[1]] ; C = res[group[2]]
        Aact = A["act"] ; Arep = A["rep"]
        Bact = B["act"] ; Brep = B["rep"]
        Cact = C["act"] ; Crep = C["rep"]

        allact = set(Aact) | set(Bact) | set(Cact) 
        allrep = set(Arep) | set(Brep) | set(Crep)

        reslist[lgd[i]] = {"act":list(allact),"rep":list(allrep)}
        i += 1
        print group
        print len(Aact),len(Bact),len(Cact),len(allact),len(Aact)+len(Bact)+len(Cact)
        print len(Arep),len(Brep),len(Crep),len(allrep),len(Arep)+len(Brep)+len(Crep)

    for group in reslist.keys():
        for reg in reslist[group].keys():
            pathDB = basedir+"data/"+gen.name+"/GO_analysis/res_lists_DB/{}_{}.txt".format(group,reg)           
            file = open(pathDB,"w")
            for g in reslist[group][reg]:
                file.write(g+"\n")
            file.close()


    # ihf = kwargs.get('ihfa', False)
    # if not ihf:
    #     for cond_fc in res.keys():
    #         res[cond_fc]['act'] = ','.join(res[cond_fc]['act'])
    #         res[cond_fc]['rep'] = ','.join(res[cond_fc]['rep'])
    #         res[cond_fc]['non'] = ','.join(res[cond_fc]['non'])
    #     df = pd.DataFrame.from_dict(res,orient='index')
    #     df.drop(['non'],axis=1,inplace=True)
    #     df.to_csv(pathdb+'list_genes.csv')

    # if ihf:
    #     for c in res.keys():
    #         titl = "{}_act_{}.txt".format(len(res[c]['act']),c)           
    #         file = open(pathdb+titl,"w")
    #         for g in res[c]['act']:
    #             file.write(g+"\n")
    #         file.close()

    #         titl = "{}_rep_{}.txt".format(len(res[c]['rep']),c)           
    #         file = open(pathdb+titl,"w")
    #         for g in res[c]['rep']:
    #             file.write(g+"\n")
    #         file.close()



    #    pairs = [["WT_expo_vs_ihfA_expo","WT_stat_vs_ihfA_stat"],["WT_nov_expo_vs_ihfA_nov_expo","WT_nov_stat_vs_ihfA_nov_stat"],["ihfA_nov_expo_vs_ihfA_expo","ihfA_nov_stat_vs_ihfA_stat"],["WT_nov_expo_vs_WT_expo","WT_nov_stat_vs_WT_stat"], ["ihfA_PGA_stat_vs_ihfA_stat","ihfA_PGA_nov_stat_vs_ihfA_nov_stat"], ["WT_PGA_stat_vs_WT_stat","WT_PGA_nov_stat_vs_WT_nov_stat"],["WT_stat_vs_WT_expo","ihfA_stat_vs_ihfA_expo"],["WT_nov_stat_vs_WT_nov_expo","ihfA_nov_stat_vs_ihfA_nov_expo"]]


    
    # pairs = [["WT_expo_vs_ihfA_expo","WT_nov_expo_vs_ihfA_nov_expo"],["WT_stat_vs_ihfA_stat","WT_nov_stat_vs_ihfA_nov_stat"]]
    # pairs = [["WT_nov_expo_vs_WT_expo","ihfA_nov_expo_vs_ihfA_expo"], ["WT_nov_stat_vs_WT_stat","ihfA_nov_stat_vs_ihfA_stat"]]
    # for pair in pairs:
    #     A = res[pair[0]] ; B = res[pair[1]]

        # Aact = set(A["act"]) ; Arep = set(A["rep"]) ; Bact = set(B["act"]) ; Brep = set(B["rep"])

        # titl = "{}_{}.txt".format(pair[0],pair[1])           
        # file = open(pathdb+titl,"w")
        # for cond in [(Aact,pair[0],"act"),(Arep,pair[0],"rep"),(Bact,pair[1],"act"),(Brep,pair[1],"rep")]:
        #     file.write(cond[2]+" "+cond[1]+"\n")
        #     for g in cond[0]:
        #         file.write(g+"\n")
        # file.close()


        # for reg in ["act","rep"]:
        #     Areg = set(A[reg]) ; Breg = set(B[reg])
        #     inter = Areg & Breg ; Aonly = Areg - Breg ; Bonly = Breg - Areg
        #     export_annot_detailed(gen, inter, "{}_{}_I_{}".format(reg,pair[0],pair[1]))
        #     export_annot_detailed(gen, Aonly, "{}_{}_only".format(reg,pair[0]))
        #     export_annot_detailed(gen, Bonly, "{}_{}_only".format(reg,pair[1]))

        # fig_width = 3.5 ; fig_height = 3.5
        # fig = plt.figure(figsize=(fig_width,fig_height))
        # d = venn2([Aact, Bact], set_labels = (pair[0], pair[1]))
        # d.get_label_by_id('A').set_fontsize(5)
        # d.get_label_by_id('B').set_fontsize(5)
        # plt.title("Activated")
        # plt.tight_layout()
        # plt.savefig(pathdb+titl[0:-4]+'.png',transparent=False, dpi=300)
        # plt.close('all')

        
        # inter = Arep & Brep ; Aonly = Arep - Brep ; Bonly = Brep - Arep

        # titl = "{}_rep_{}_{}.txt".format(str(len(inter)),pair[0],pair[1])           
        # file = open(pathdb+titl,"w")
        # for g in inter:
        #     file.write(g+"\n")
        # file.close()

        # fig_width = 3.5 ; fig_height = 3.5
        # fig = plt.figure(figsize=(fig_width,fig_height))
        # d = venn2([Arep, Brep], set_labels = (pair[0], pair[1]))
        # d.get_label_by_id('A').set_fontsize(5)
        # d.get_label_by_id('B').set_fontsize(5)
        # plt.title("Repressed")
        # plt.tight_layout()
        # plt.savefig(pathdb+titl[0:-4]+'.png',dpi=300,transparent=False)
        # plt.close('all')

def compute_zscore_binomial(nbtot, nbobs, pexp):
    '''
    Approximation of binomial law by normal law and zscore computation
    '''
    return (nbobs - nbtot*pexp) / (np.sqrt((nbtot)*pexp*(1-pexp)))

def compute_zscore_2means(X1, X2, mudiff, sd1, sd2, n1, n2):
    pooledSE = np.sqrt(sd1**2/n1 + sd2**2/n2)
    z = ((X1 - X2) - mudiff)/pooledSE
    return z


def correct_sequence(gen,orient,mioC):
    print gen.name
    gen.load_seq()
    print len(gen.seq)

    if orient:
        newseq = gen.seq[mioC-1:] + gen.seq[0:mioC-1]
    else:
        newseq = gen.seqcompl[mioC-1:] + gen.seqcompl[0:mioC-1]     
        newseq = newseq[::-1
        ]
    print len(newseq)
    f = open(basedir+"data/"+gen.name+"/sequence.fasta","w")
    f.write("> "+gen.name+"\n")
    f.write(newseq)
    f.close()


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