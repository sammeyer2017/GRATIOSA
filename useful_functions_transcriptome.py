#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import math
from globvar import *
import os
import pysam
import Genome
#==============================================================================#


##### functions called by transcriptome methods #####
def add_single_rpkm_to_genes(genes_dict, expression_filename, condition, TSS_column, start_line, separator,tag_column):
    """ Adds rpkm data to Gene objects by parsing a file with two columns:
    gene name and value
    """
    if separator == '\\t' : separator = '\t' 
    with open(expression_filename, 'r') as f:
        i=1
        while i != start_line:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n').split(separator)
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


def add_expression_to_genes(genes_dict, cond, filename, tag_col, expr_col, is_log, separator):
    """ Adds expression data to Gene objects by parsing a file with as many
    columns as there are different conditions in the experiment, plus one for
    the gene names (first column).
    """        
    genes_valid = {} 
    except_locus = []
    expr0 = 0
    if separator == '\\t' : separator = '\t' 
    with open(filename, 'r') as f:
        header=next(f)
        header=header.strip().split(separator)        
        header=header[expr_col:]
        genes_valid["conditions"] = header
        for line in f:
            line=line.strip().split(separator)
            try:
                if is_log == 'no':
                    genes_dict[line[tag_col]].add_expression(cond,math.log(float(line[expr_col]),2))
                elif is_log == 'yes':
                    genes_dict[line[tag_col]].add_expression(cond,float(line[expr_col]))
                else : 
                    print("Please set 'is_log' to 'yes' or 'no' in the expression.info file")
                genes_valid[line[tag_col]] = genes_dict[line[tag_col]]

            except KeyError:
                if line[tag_col] == 'none':
                    print("expressions without locus tag")
                else:
                    except_locus.append(line[tag_col])
            except Exception as e: 
                    if int(float(line[expr_col])) == 0:
                        expr0 += 1
                    else : 
                        print(e)

    print(f"{len(genes_valid.keys())} expressions were successfully loaded")
    if except_locus :
        if len(except_locus) > 20 :
            print(f"{len(except_locus)} locus are not in annotation")
        else : 
            print(f"{except_locus} are not in annotation" )

    if expr0: 
        print(f"{expr0} expressions are null and led to math domain error")
    return genes_valid


def load_fc_pval_cond(genes_dict, filename, condition, tag_col, fc_col, separator, start_line, *args, **kwargs):
    ''' Called by load_fc_pval, allows expression data to be loaded by specifying files, and where each
    information is (tag, fc, pval...). If no p-value column, assigns pval = 0 to each gene
    '''
    genes_valid = [] # list containing all genes having valid FC / pval
    p_val_col= kwargs.get('p_value')
    if separator == '\\t' : separator = '\t' 
    except_locus = []
    except_fc = 0 
    with open(filename, 'r') as f:
        i=0
        while i < start_line:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n').split(separator)
            try:
                if p_val_col:
                    genes_dict[line[tag_col]].add_fc_pval_cond(float(line[fc_col]),condition, float(line[p_val_col]))
                else:
                    genes_dict[line[tag_col]].add_fc_pval_cond(float(line[fc_col]),condition, float(0))
                genes_valid.append(line[tag_col])
            except:
                if line[tag_col] not in genes_dict.keys():
                    if line[tag_col] != '':
                        except_locus.append(line[tag_col])
                    else:
                        except_fc += 1
    f.close()
    if except_locus :
        if len(except_locus) > 20 :
            print(f"\t{len(except_locus)} locus are not in annotation")
        else : 
            print(f"\t{except_locus} are not in annotation" )
    if except_fc :
        print(f"\t{except_fc} fc without locus")
    return genes_valid

def load_expr_cond(genes_dict, filename, condition, tag_col, nb_replicates, expr_col, separator, start_line):
    ''' Called by load_expr_cond, allows expression data to be loaded by specifying files, and where each
    information is (tag, expr...).
    '''
    if separator == '\\t' : separator = '\t'
    except_locus = []
    with open(filename, 'r') as f:
        i=1
        while i < start_line:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n').split(separator)
            st = int(expr_col) ; nb = int(nb_replicates)
            for i in range(st,st+nb):
                try:
                    genes_dict[line[tag_col]].add_expr_cond(float(line[i]),condition)
                except:
                    if line[tag_col] not in genes_dict.keys():
                        if line[tag_col] != '':
                            except_locus.append(line[tag_col])
                        else:
                            print("fc without locus")
    f.close()
    if except_locus :
        if len(except_locus) > 20 :
            print(f"{len(except_locus)} locus are not in annotation")
        else : 
            print(f"{except_locus} are not in annotation" )

def process_bam_paired_end(gen): # /!\ paired-end only /!\ -> return fragments for each strand
    """
    tophat2 -o tophat_out/WT genome/MG1655 SRR12057814_1.fastq.gz SRR12057814_2.fastq.gz
    """
    exist_cond = []
    path2dir = f"{basedir}data/{gen.name}/rnaseq_reads/"
    if not os.path.exists(path2dir+'reads.info'):
        file = open(path2dir+'reads.info','w') 
        file.write('Condition\tReads file\tBAM file\n')
    else : 
        file = open(path2dir+'reads.info','r')
        header = next(file)
        for line in file : 
            exist_cond.append(line.split("\t")[0])
    file.close()
    with open(path2dir+"bam_files.info","r") as f:
        header = next(f)
        for line in f: # for each condition
            line=line.strip().split('\t')
            bam_cond = line[0]
            if bam_cond not in exist_cond :
                bam_file=line[1]
                print(f"Converting .bam file for: {bam_cond}")
                if not os.path.exists(path2dir+bam_cond+".bai"): # if index needed, created using samtools (.bai file)
                    os.system("samtools index -b %s"%path2dir+bam_file)
                bamfile = pysam.AlignmentFile(path2dir+bam_file, "rb") # BAM opening, alignment file object

                # /!\ 0-based coordinate system /!\
                # lists storing coordinates of paired-end fragments for +/- strands
                Rpos_start, Rpos_end, Rneg_start, Rneg_end = [], [], [], []

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

                # save results ; .npz contains two .npy Rpos and Rneg
                file = open(path2dir+'reads.info','a')
                file.write(bam_cond+'\t'+bam_cond+'_reads.npz\t'+bam_file+'\n')  
                file.close()    
                np.savez(path2dir+bam_cond+'_reads.npz', Rpos=Rpos, Rneg=Rneg)


# compute coverage from reads (.npz, one .npy per strand)
def cov_from_reads(tr):
    # load npz file
    gen = Genome.Genome(tr.name)
    if not hasattr(gen, "seq"):
        gen.load_seq()
    genome_length = gen.length
    exist_cond = []
    path2cov = f"{basedir}data/{gen.name}/rnaseq_cov/"
    path2reads = f"{basedir}data/{gen.name}/rnaseq_reads/"
    if not os.path.exists(path2cov+'cov.info'):
        file = open(path2cov+'cov.info','w') 
        file.write('Condition\tCov file\n')
    else : 
        file = open(path2cov+'cov.info','r')
        header = next(file)
        for line in file : 
            exist_cond.append(line.split("\t")[0])
    file.close()
    with open(path2reads+"reads.info","r") as f:
        header = next(f)
        for line in f: # for each condition
            line=line.strip().split('\t')

            if line[0] not in exist_cond :
                print('Converting reads in coverage for:',line[0])   
                # load npz file corresponding to condition
                npzfile = np.load(path2reads+line[1])
                Rpos = npzfile["Rpos"]
                Rneg = npzfile["Rneg"]
                # init cov
                cov_pos = np.zeros(genome_length, dtype=int)
                cov_neg = np.zeros(genome_length, dtype=int)
                # compute cov on positive strand
                for start,end in Rpos:
                    cov_pos[start:end+1] += 1
                # on negative strand
                for start,end in Rneg:
                    cov_neg[start:end+1] += 1

                # save results
                file = open(path2cov+'cov.info','a')
                file.write(f"{line[0]}\t{line[0]}_cov.npz\n")  
                file.close()
                # save results ; .npz contains two .npy cov_pos and cov_neg
                np.savez(path2cov+line[0]+'_cov.npz', cov_pos=cov_pos, cov_neg=cov_neg)

# compute coverage from reads (.npz, one .npy per strand)
def cov_start_stop_from_reads(tr):
    gen = Genome.Genome(tr.name)
    if not hasattr(gen, "seq"):
        gen.load_seq()
    genome_length = gen.length
    exist_cond = []
    path2cov = f"{basedir}data/{tr.name}/rnaseq_cov/"
    path2reads = f"{basedir}data/{tr.name}/rnaseq_reads/"
    if not os.path.exists(path2cov+'cov_start_stop.info'):
        file = open(path2cov+'cov_start_stop.info','w') 
        file.write('Condition\tCov file\n')
    else : 
        with open(path2cov+'cov_start_stop.info','r') as file:
            header = next(file)
            for line in file : 
                exist_cond.append(line[0])
    file.close()
    with open(path2reads+"reads.info","r") as f:
        header = next(f)
        for line in f: # for each condition
            line=line.strip().split('\t')
            if line[0] not in exist_cond :
                print('Searching for start and stop from reads in:',line[0])   
                # load npz file corresponding to condition
                npzfile = np.load(path2reads+line[1])
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
                # save results
                file = open(path2cov+'cov_start_stop.info','a')
                fname = line[0] + "_cov_start_end"
                file.write(f'{line[0]}\t{fname}.npz\n')  
                file.close()
                # save results ; .npz contains four .npy: start on pos, stop on pos, start on neg, stop on neg
                np.savez(path2cov+fname+'.npz', cov_start_pos= cov_start[1],
                    cov_end_pos= cov_end[1], 
                    cov_start_neg= cov_start[0], 
                    cov_end_neg= cov_end[0])