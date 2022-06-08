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
#from Bio.Alphabet import IUPAC
from Bio import motifs
from TSS import TSS
from matplotlib_venn import venn2

from scipy import stats
from statsmodels.stats import weightstats
from operator import itemgetter
import itertools

plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})
plt.rcParams.update({'font.size': 11})
plt.rcParams.update({'legend.fontsize': 9})
plt.rcParams.update({'font.family': "Arial"})

#==============================================================================#

# -------------------
# useful function for weight matrix

#def create_matrix_from_file(filename, factor):
#    i =1
#    with open(filename,"r") as f:
#        for line in f:
#            if 'Transcription Factor Name: '+factor in line:
#                while i <5:
#                    i+=1
#                    header=next(f)
#                m = motifs.Motif()
#                a = IUPAC.unambiguous_dna
#                m.add_instance(Seq(header.strip(),a))
#                l=('*')
#                j=1
#                while j != len(header.strip()):
#                    j+=1
#                    l+='*'
#                while header[0] != '\n':
#                    header=next(f)
#                    if header[0] != '\n':
#                        m.add_instance(Seq(header.strip(),a))
#               m.set_mask(l)
#                if header[0]=='\n':
#                    return m

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
    print(res1)
    print(res2)
    print(len(res1.keys()),len(res2.keys()))

    res = {}
    for g1 in res1.keys():
        if g1 in res2.keys():
            mismatches = sum(c1!=c2 for c1,c2 in zip(res1[g1],res2[g1]))
            if mismatches not in res.keys():
                res[mismatches] = 1
            else:
                res[mismatches] += 1
    print(res)
    print(sum(res[c] for c in res.keys()))

# It happens often that expression and TSS data come from different subspecies with different gene IDs
# and slight structural differences. In that case, you can use the following function, that test whether
# or not two genes with the same name have similar sequences. If it is the case, expression / TSS data are
# extrapolated from the other subspecie (see frequent_cmd)

def load_fc_from_common_genome(genTSS,genFC,*arg,**kwargs):
    genFC.load_fc_pval() # genome from which we want to extrapolate FC (e.g. ecoli_b)
    genTSS.load_fc_pval() # genome where we want to extrapolate FC (e.g. ecoli K12)
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
    print(d)


def barplot_annotate_brackets(num1, num2, text, center, height, yerr=None, dh=.05, barh=.05, fs=12, maxasterix=None, dt=0, bold = False, lw=1):
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
    :param lw: bar linewidth
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

    plt.plot(barx, bary, c='black', linewidth=lw)
    if bold:
        plt.text(mid[0],mid[1]-dt, text, fontsize=fs, ha='center', fontweight = "bold")
    else:
        plt.text(mid[0],mid[1]-dt, text, fontsize=fs, ha='center')



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


########## BAM2NPZ RNA_SEQ ANALYSIS, RAPHAEL FORQUET ##########
import os,pysam, numpy as np, pandas as pd
from datetime import datetime

def process_bam_paired_end(bam_file, gen): # /!\ paired-end only /!\ -> return fragments for each strand
    if not os.path.exists(basedir+"data/"+gen+'/rnaseq_reads/'+bam_file+".bai"): # if index needed, created using samtools (.bai file)
        os.system("samtools index -b %s"%bam_file)
    
    bamfile = pysam.AlignmentFile(basedir+"data/"+gen+'/rnaseq_reads/'+bam_file, "rb") # BAM opening, alignment file object
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

def hierarchical_clustering(genes, corr, thresh):
    '''
    genes = list of genes, corr = matrix of correlations among genes, thresh = correlation cutoff to end clustering
    hierarchical clustering with constraint on neighbours

    '''
    stop = False # stop = True when clusters can't fuse anymore because of thresh_corr reached
    # init, each is cluster is a gene
    clusters = [[x] for x in genes]
    while not stop: # when clusters can still fuse 
        if len(clusters) > 1:
            corrs = [] # list of correlations among clusters
            idx = range(len(clusters)) # idx of clusters
            # combinations of neighbours
            combinations = [x for x in list(itertools.combinations(idx, 2)) if x[1] - x[0] == 1]
            for c in combinations:
            # correlation among neighbours clusters
                l = []
                for g1 in clusters[c[0]]:
                    for g2 in clusters[c[1]]:
                        try:
                            l.append(corr[g1][g2])
                        except Exception as e: # no correlation data for the pair of genes
                            pass
                
                # median = np.nan if no correlation data ; in that case, we still fuse the clusters
                corrs.append((c[0], c[1], np.median(l)))

            higher = max(corrs, key=itemgetter(2))
            if higher[2] >= thresh or np.isnan(higher[2]):
                idx_before = [x for x in idx if x < higher[0]]
                idx_after = [x for x in idx if x > higher[1]]
                clusters = [clusters[x] for x in idx_before] + [clusters[higher[0]] + clusters[higher[1]]] + [clusters[x] for x in idx_after] 
            else:
                stop = True
        else:
            stop = True

    return clusters
