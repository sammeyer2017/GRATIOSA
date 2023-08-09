#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions called by Transcriptome methods
"""
import numpy as np
from GRATIOSA.globvar import *
from GRATIOSA import Genome
import os
import pysam
import math

def add_expression_to_genes(genes_dict, 
                            cond, 
                            filename, 
                            tag_col, 
                            expr_col, 
                            is_log, 
                            separator):
    """
    Called by load_expression, adds expression data to Gene objects by 
    parsing a file with one column containing the locus tags of the genes and 
    one containing the expression data. 

    Creates a new attribute "expression" of the Gene objects contained in 
    the genes_dict given as argument. This attribute is a dictionary of shape 
    {condition: expression level [float.]}

    Args:
        genes_dict (dict.): dictionary of shape {locus: Gene object}
        cond (str.): condition name
        filename (str.): path to the file containing the expression data
        tag_col (int.): index of the column containing the locus tags in the 
                file
        expr_col (int.): index of the column containing the expression data
        is_log (bool.): True if the data are already in log2, such as log2RPKM.
                        False if the data are RPKM.
        separartor (str.): file separator

    Returns:
        Dictionary: Dict. of shape {'locus':Gene object} with in key the locus 
                tags associated to an expression data and in value the gene 
                object associated to each of these locus tags.
    Note: 
        Column numbering starts at 0.
    """
    genes_valid = {}
    except_locus = []
    expr0 = 0
    if separator == '\\t':
        separator = '\t'
    with open(filename, 'r') as f:
        header = next(f)
        header = header.strip().split(separator)
        header = header[expr_col:]
        genes_valid["conditions"] = header
        for line in f:
            line = line.strip().split(separator)
            try:
                if is_log == 'no':
                    genes_dict[line[tag_col]].add_expression(
                        cond, math.log(float(line[expr_col]), 2))
                elif is_log == 'yes':
                    genes_dict[line[tag_col]].add_expression(
                        cond, float(line[expr_col]))
                else:
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
                else:
                    print(e)

    print(f"{len(genes_valid.keys())} expressions were successfully loaded")
    if except_locus:
        if len(except_locus) > 20:
            print(f"{len(except_locus)} locus are not in annotation")
        else:
            print(f"{except_locus} are not in annotation")

    if expr0:
        print(f"{expr0} expressions are null and led to math domain error")
    return genes_valid


def load_fc_pval_cond(genes_dict, 
                      filename, 
                      condition,
                      tag_col, 
                      fc_col, 
                      separator, 
                      start_line, 
                      p_val_col=None):
    '''
    Called by load_fc_pval, allows expression data to be loaded by specifying
    a file containing the following information (one column for each type
    of information): locus tags of the genes, log2(Fold-Changes) and p-value.

    Args:
        genes_dict (dict.): dictionary of shape {locus: Gene object}
        filename (str.): path to the file containing the expression data
        cond (str.): condition name
        tag_col (int.): index of the column containing the locus tags in the 
                file
        fc_col (int.): index of the column containing the log2FC
        separartor (str.): file separator
        start_line (int.): file start line
        p_val_col (int.): index of the column containing the p-values in the
                file. If no index is indicated (default: None), the p-values 
                will be assigned to 0 for all genes.
    
    Returns:
        List: locus tags of all genes having a valid log2FC

    Note: 
        Column numbering starts at 0.
    '''
    genes_valid = []
    if separator == '\\t':
        separator = '\t'
    except_locus = []
    except_fc = 0
    with open(filename, 'r') as f:
        i = 0
        while i < start_line:
            header = next(f)
            i += 1
        for line in f:
            line = line.strip('\n').split(separator)
            try:
                if p_val_col:
                    genes_dict[line[tag_col]].add_fc_pval_cond(
                        float(line[fc_col]), condition, float(line[p_val_col]))
                else:
                    genes_dict[line[tag_col]].add_fc_pval_cond(
                        float(line[fc_col]), condition, float(0))
                genes_valid.append(line[tag_col])
            except BaseException:
                if line[tag_col] not in genes_dict.keys():
                    if line[tag_col] != '':
                        except_locus.append(line[tag_col])
                    else:
                        except_fc += 1
    f.close()
    if except_locus:
        if len(except_locus) > 20:
            print(f"\t{len(except_locus)} locus are not in annotation")
        else:
            print(f"\t{except_locus} are not in annotation")
    if except_fc:
        print(f"\t{except_fc} fc without locus")
    return genes_valid


def process_bam_paired_end(tr):
    """
    Called by load_cov_start_end and load_rnaseq_cov, process_bam_paired_end
    converts paired-end .bam files in .npy files containing the fragments.
    These files will enable a faster data loading for the next data
    importation with load_cov_start_end or load_rnaseq_cov.
    
    Creates:
        * For each condition listed in the bam_files.info file, this function
          creates a npz file containing the coverage data on both DNA strands
          (Rpos.npy and Rneg.npy). These files are saved in the /rnaseq_reads/
          directory
        *  A reads.info file containing the information required to associate
           each .npz  file to a condition name is also created (or completed if
           this file already exists) in the /rnaseq_reads/ directory. This file
           contains the following columns:
           [0] Condition [1] Reads filename [2] BAM filename

    Args:
        tr: Transcriptome instance

    Note: 
        The paired-end .bam reads files are, with a bam_files.info file, in the
        /rnaseq_reads/ directory. The bam_files.info file contains the following
        information:  [0] Condition [1] Reads filename for this condition
    
    Warnings: 
        Only for paired-end reads !
    """

    exist_cond = []
    path2dir = f"{basedir}data/{tr.name}/rnaseq_reads/"
    if not os.path.exists(path2dir + 'reads.info'):
        file = open(path2dir + 'reads.info', 'w')
        file.write('Condition\tReads filename\tBAM filename\n')
    else:
        file = open(path2dir + 'reads.info', 'r')
        header = next(file)
        for line in file:
            exist_cond.append(line.split("\t")[0])
    file.close()
    with open(path2dir + "bam_files.info", "r") as f:
        header = next(f)
        for line in f:  # for each condition
            line = line.strip().split('\t')
            bam_cond = line[0]
            if bam_cond not in exist_cond:
                bam_file = line[1]
                print(f"Converting .bam file for: {bam_cond}")
                if not os.path.exists(f"{path2dir}{bam_cond}.bai"):  
                        # if index needed, created using samtools (.bai file)
                    os.system("samtools index -b %s" % path2dir + bam_file)
                # BAM opening, alignment file object
                bamfile = pysam.AlignmentFile(path2dir + bam_file, "rb")

                # /!\ 0-based coordinate system /!\
                # lists storing coordinates of paired-end fragments for +/-
                # strands
                Rpos_start, Rpos_end, Rneg_start, Rneg_end = [], [], [], []

                for read in bamfile.fetch():
                    # if correctly paired and R1/R2 in different strands and
                    # fragment size < 1000 kb
                    if read.is_read1 and read.is_paired and not read.mate_is_unmapped and read.is_reverse != read.mate_is_reverse and abs(
                            read.template_length) < 1000:
                        if read.is_reverse:
                            Rneg_start.append(read.reference_end)
                            Rneg_end.append(
                                read.reference_end + abs(read.template_length))
                        elif not read.is_reverse:
                            Rpos_start.append(read.reference_start)
                            Rpos_end.append(
                                read.reference_start + read.template_length)
                bamfile.close()

                # conversion step into numpy array
                Rpos_start = np.array(Rpos_start, dtype=int)
                Rpos_end = np.array(Rpos_end, dtype=int)
                Rneg_start = np.array(Rneg_start, dtype=int)
                Rneg_end = np.array(Rneg_end, dtype=int)

                Rpos = np.column_stack((Rpos_start, Rpos_end))
                Rneg = np.column_stack((Rneg_start, Rneg_end))

                # delete rows where one coordinate is missing
                Rpos = Rpos[np.isfinite(Rpos).all(axis=1)]
                Rneg = Rneg[np.isfinite(Rneg).all(axis=1)]

                # save results ; .npz contains two .npy Rpos and Rneg
                file = open(path2dir + 'reads.info', 'a')
                file.write(f"{bam_cond}\t{bam_cond}_reads.npz\t{bam_file}\n")
                file.close()
                np.savez(f"{path2dir}{bam_cond}_reads.npz", 
                         Rpos=Rpos, 
                         Rneg=Rneg)


def cov_from_reads(tr):
    """
    Called by load_rnaseq_cov, cov_from_reads computes the coverage from 
    paired-end reads previously converted to .npz format (containing one .npy 
    per strand: Rpos.npy and Rneg.npy) with the process_bam_paired_end 
    function. 

    Creates:
        * For each condition listed in the reads.info file, this function
          creates a npz file containing the coverage data on both DNA strands:
          cov_pos.npy and cov_neg.npy. These files are saved in the 
          /rnaseq_cov/ directory and will enable a faster data loading for the 
          next data importation with load_rnaseq_cov.

        * A cov.info file containing the information required to associate
          each .npz  file to a condition name is also created (or completed if
          this file already exists) in the /rnaseq_cov/ directory. This file
          contains the following columns:  [0] Condition [1] Coverage filename

    Args:
        tr: Transcriptome instance

    Note: 
        The .npz files have to be, with a reads.info file (also created 
        by the process_bam_paired_end function), in the /rnaseq_reads/ directory. 
        This reads.info file contains at least the following information: 
        [0] Condition [1] Reads filename
    """
    # load npz file
    gen = Genome.Genome(tr.name)
    if not hasattr(gen, "seq"):
        gen.load_seq()
    genome_length = gen.length
    exist_cond = []
    path2cov = f"{basedir}data/{gen.name}/rnaseq_cov/"
    path2reads = f"{basedir}data/{gen.name}/rnaseq_reads/"
    if not os.path.exists(path2cov + 'cov.info'):
        file = open(path2cov + 'cov.info', 'w')
        file.write('Condition\tCoverage filename\n')
    else:
        file = open(path2cov + 'cov.info', 'r')
        header = next(file)
        for line in file:
            exist_cond.append(line.split("\t")[0])
    file.close()
    with open(path2reads + "reads.info", "r") as f:
        header = next(f)
        for line in f:  # for each condition
            line = line.strip().split('\t')

            if line[0] not in exist_cond:
                print('Converting reads in coverage for:', line[0])
                # load npz file corresponding to condition
                npzfile = np.load(path2reads + line[1])
                Rpos = npzfile["Rpos"]
                Rneg = npzfile["Rneg"]
                # init cov
                cov_pos = np.zeros(genome_length, dtype=int)
                cov_neg = np.zeros(genome_length, dtype=int)
                # compute cov on positive strand
                for start, end in Rpos:
                    cov_pos[start:end + 1] += 1
                # on negative strand
                for start, end in Rneg:
                    cov_neg[start:end + 1] += 1

                # save results
                file = open(path2cov + 'cov.info', 'a')
                file.write(f"{line[0]}\t{line[0]}_cov.npz\n")
                file.close()
                # save results ; .npz contains two .npy cov_pos and cov_neg
                np.savez(f"{path2cov}{line[0]}_cov.npz", 
                         cov_pos=cov_pos,
                         cov_neg=cov_neg)


def cov_start_stop_from_reads(tr):
    """
    Called by load_cov_start_end, cov_start_stop_from_reads computes the 
    density of RNA fragment starts and ends from paired-end reads previously 
    converted to .npz format (containing one .npy per strand: Rpos.npy and 
    Rneg.npy) with the process_bam_paired_end function.
    
    Creates:
        * For each condition listed in the reads.info file, this function
          creates a npz file containing the density of RNA fragment starts and
          ends on both DNA strands: cov_start_pos.npy, cov_start_neg.npy,
          cov_end_pos.npy and cov_end_neg.npy.
          These files are saved in the /rnaseq_cov/ directory and will enable
          a faster data loading for the next data importation with the
          cov_start_stop_from_reads function.
        * A cov_start_stop.info file containing the information required to 
          associate each .npz  file to a condition name is also created (or 
          completed if this file already exists) in the /rnaseq_cov/ directory. 
          This file contains the following columns:  
          [0] Condition [1] Coverage filename

    Args:
        tr: Transcriptome instance

    Note: 
        The .npz files have to be, with a reads.info file (also created by 
        process_bam_paired_end function), in the /rnaseq_reads/ directory. 
        This reads.info file contains at least the following information:
        [0] Condition [1] Reads filename
    """
    gen = Genome.Genome(tr.name)
    if not hasattr(gen, "seq"):
        gen.load_seq()
    genome_length = gen.length
    exist_cond = []
    path2cov = f"{basedir}data/{tr.name}/rnaseq_cov/"
    path2reads = f"{basedir}data/{tr.name}/rnaseq_reads/"
    if not os.path.exists(path2cov + 'cov_start_stop.info'):
        file = open(path2cov + 'cov_start_stop.info', 'w')
        file.write('Condition\tCov file\n')
    else:
        with open(path2cov + 'cov_start_stop.info', 'r') as file:
            header = next(file)
            for line in file:
                exist_cond.append(line[0])
    file.close()
    with open(path2reads + "reads.info", "r") as f:
        header = next(f)
        for line in f:  # for each condition
            line = line.strip().split('\t')
            if line[0] not in exist_cond:
                print('Searching for start and stop from reads in:', line[0])
                # load npz file corresponding to condition
                npzfile = np.load(path2reads + line[1])
                Rpos = npzfile["Rpos"]
                Rneg = npzfile["Rneg"]
                # init cov
                cov_start = {0: np.zeros(genome_length, dtype=int), 
                             1: np.zeros(genome_length, dtype=int)}
                cov_end = {0: np.zeros(genome_length, dtype=int), 
                           1: np.zeros(genome_length, dtype=int)}
                # compute cov
                # on positive strand
                for start, end in Rpos:
                    try:
                        cov_start[1][start - 1] += 1
                        cov_end[1][end - 1] += 1
                    except BaseException:
                        pass
                # on negative strand
                for start, end in Rneg:
                    try:
                        cov_start[0][start - 1] += 1
                        cov_end[0][end - 1] += 1
                    except BaseException:
                        pass
                # save results
                file = open(path2cov + 'cov_start_stop.info', 'a')
                fname = line[0] + "_cov_start_end"
                file.write(f'{line[0]}\t{fname}.npz\n')
                file.close()
                np.savez(f"{path2cov}{fname}.npz", 
                         cov_start_pos=cov_start[1],
                         cov_end_pos=cov_end[1],
                         cov_start_neg=cov_start[0],
                         cov_end_neg=cov_end[0])
