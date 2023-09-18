#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from datetime import datetime
from GRATIOSA.useful_functions_transcriptome import *
from GRATIOSA.globvar import *
from GRATIOSA.Genome import Genome

class Transcriptome:
    '''
    The Transcriptome class is the key class for loading expression data and
    differential expression data for analysis and comparison. Typical data are
    data obtained with high throughput methods such as RNASeq and processed with
    common tools such as DeSeq2.

    Each Transcriptome instance has to be initialized with an organism name
    
        >>> tr = Transcriptome.Transcriptome("dickeya")
    '''

    def __init__(self, name):
        self.name = name

    def load_expression(self):
        '''
        loads expression data (such as log2rpkm) from files that are in the
        /expression/ directory.

        Creates:
            self.genes[locus].expression (dict.): new attribute of Gene 
                    instances related to the Transcriptome instance given as 
                    argument. Dictionary of shape 
                    {condition: expression level (float.)}
            self.genes_valid_expr(dict.): Dictionary of shape 
                    {condition: list of genes having an expression value}

        Note:
            The data importation requires an expression.info file, containing
            column indices of each information in the data file and some 
            additional information, in the following order:
            [0] Condition [1] Filename [2] Locus_tag column
            [3] Expression column [4] is Log ? (boolean) [5] Separator

        Warning:
            This method needs a genomic annotation. If no annotation is
            loaded, the load_annotation method with the default "sequence.gff3"
            file is computed. To use another annotation, please load an 
            annotation to your Transcriptome instance with the following commands 
            before using this method:
                >>> from GRATIOSA import Genome, Transcriptome
                >>> tr = Transcriptome.Transcriptome("ecoli")
                >>> g = Genome.Genome(tr.name)
                >>> g.load_annotation(annot_file=chosen_file)
                >>> tr.genes = g.genes

        Example:
            >>> tr = Transcriptome.Transcriptome("dickeya")
            >>> tr.load_expression()
            >>> tr.genes["Dda3937_00004"].expression
            {'WT_nov_stat(E10)_rpkm': 66.6413789332929,
            'WT_stat(E3)': 6.980245227157392,
            'WT_PGA_stat(E13)_rpkm': 13.9428053948966}
        '''
        if not hasattr(self, 'genes_valid_expr'):
            self.genes_valid_expr = {}
        if not hasattr(self, "genes"):
            gen = Genome(self.name)
            gen.load_annotation()
            self.genes = gen.genes

        path2dir = f"{basedir}data/{self.name}/expression/"
        if os.path.exists(path2dir + "expression.info"):
            with open(path2dir + "expression.info", "r") as f:
                skiphead = next(f)  # skip head
                for header in f:
                    header = header.strip().split('\t')
                    path2file = path2dir + header[1]
                    print(f"Loading condition: {header[0]}")
                    self.genes_valid_expr[header[0]] = add_expression_to_genes(
                        genes_dict=self.genes,
                        cond=header[0],
                        filename=path2file,
                        tag_col=int(header[2]),
                        expr_col=int(header[3]),
                        is_log=header[4],
                        separator=header[5])
            f.close()
        else:
            print("No expression.info file, please create one")

    def compute_fc_from_expr(self, ctrls, conds, condname):
        '''
        Compute_fc_from_expr computes, for each gene:
            * the mean of expression values (such as log2rpkm) in a list of 
              controls replicates (ctrls)
            * the mean of expression values in a list of tests replicates (conds)
            * the fold-change ie the difference between the mean control 
              and the mean test values
            * the pvalue of the Student-test for these two means

        Creates:
            * self.genes[locus].fc_pval[condname] (tuple)
                    new attribute of Gene instances related to the 
                    Transcriptome instance given as argument.Tuple of shape 
                    (fold-change,pvalue).
            * self.genes_valid_fc[condname] (list.)
                    list of genes having a fold-change and a pvalue corresponding 
                    to the new condition

        Args:
            conds (list of str.): list of expression conditions to use as test
            ctrls (list of str.): list of expression conditions to use as control
            condname (str.): condition name to use for the new fold-changes
                             and p-values
       
        Note:
            The t-test is computed with scipy.stats.ttest_ind

        Warning:
            This method needs a genomic annotation. If no annotation is
            loaded, the load_annotation method with the default "sequence.gff3"
            file is computed. To use another annotation, please load an 
            annotation to your Transcriptome instance with the following commands 
            before using this method:
                >>> from GRATIOSA import Genome, Transcriptome
                >>> tr = Transcriptome.Transcriptome("ecoli")
                >>> g = Genome.Genome(tr.name)
                >>> g.load_annotation(annot_file=chosen_file)
                >>> tr.genes = g.genes

        Example:
            >>> tr = Transcriptome.Transcriptome("ecoli")
            >>> tr.compute_fc_from_expr(["Ctrl1","Ctrl2"],
                                        ["Test1","Test2"],
                                        ["FC_Test12_Ctrl12"])
            >>> tr.genes["b0002"].fc_pval["FC_Test12_Ctrl12"]
            (-1.786, 0.003125)
        '''
        if not hasattr(self, 'genes_valid'):
            self.genes_valid = {}
        self.load_expression()  # load expression values

        genes_val = []  # list containing all genes having valid expr
        for genename in self.genes.keys():
            gene = self.genes[genename]
            ctrlvals = []  # list of control values for that gene
            testvals = []  # list of test values for that gene
            try:
                for ctrl in ctrls:
                    ctrlvals.append(gene.expression[ctrl])
                for cond in conds:
                    testvals.append(gene.expression[cond])
                # add FC (meantest - meanctrl) and p-values (Student test)
                gene.add_fc_pval_cond(
                    np.mean(testvals) -
                    np.mean(ctrlvals),
                    condname,
                    stats.ttest_ind(
                        ctrlvals,
                        testvals,
                        equal_var=False)[1])
                genes_val.append(genename)
            except BaseException:
                pass

        self.genes_valid[condname] = genes_val

    def load_fc_pval(self):
        '''
        loads Fold-changes and p-values (if available) from data files that 
        are in the /fold_changes/ directory.

        Creates:
            * self.genes[locus].fc_pval (dict.): new attribute of Gene instances
                    related to the Transcriptome instance given as argument.
                    Dictionary of shape {condition: (fold-change,pvalue)}.
                    If no p-value is available in the data file, 0 is set as
                    p-value.
            * self.genes_valid_fc (dict.): Dictionary of shape 
                    {condition: list of genes having at least a fold-change}

        Note:
            The data importation requires a fc.info file, containing
            column indices of each information in the data file and some 
            additional information, in the following order:
            [0] Condition [1] Filename [2] Locus_tag column
            [3] FC column [4] Separator [5] Startline [6] P-value column

        Warning:
            This method needs a genomic annotation. If no annotation is
            loaded, the load_annotation method with the default "sequence.gff3"
            file is computed. To use another annotation, please load an 
            annotation to your Transcriptome instance with the following commands 
            before using this method:
                >>> from GRATIOSA import Genome, Transcriptome
                >>> tr = Transcriptome.Transcriptome("ecoli")
                >>> g = Genome.Genome(tr.name)
                >>> g.load_annotation(annot_file=chosen_file)
                >>> tr.genes = g.genes

        Example:
            >>> tr = Transcriptome.Transcriptome("ecoli")
            >>> tr.load_fc_pval()
            >>> tr.genes["b0002"].fc_pval
            {'WT': (-1.786, 0.003125),
             'osmotic': (-1.81603942323042, 0.0),
             'heat': (1.62203513077067, 0.0)}
        '''
        self.genes_valid_fc = {}

        if not hasattr(self, "genes"):
            gen = Genome(self.name)
            gen.load_annotation()
            self.genes = gen.genes

        path2dir = f"{basedir}data/{self.name}/fold_changes/"
        if os.path.exists(path2dir + "fc.info"):
            with open(path2dir + "fc.info", "r") as f:
                skiphead = next(f)  # skip head
                for header in f:
                    header = header.strip()
                    header = header.split('\t')
                    print(f"Loading condition: {header[0]}")
                    try:  # if p-value column specified works
                        self.genes_valid_fc[header[0]] = load_fc_pval_cond(
                            genes_dict=self.genes,
                            filename=path2dir + header[1],
                            condition=header[0],
                            tag_col=int(header[2]),
                            fc_col=int(header[3]),
                            separator=header[4],
                            start_line=int(header[5]),
                            p_value=int(header[6]))
                    except BaseException:  # otherwise, set pvalue = 0
                        self.genes_valid_fc[header[0]] = load_fc_pval_cond(
                            genes_dict=self.genes,
                            filename=path2dir + header[1],
                            condition=header[0],
                            tag_col=int(header[2]),
                            fc_col=int(header[3]),
                            separator=header[4],
                            start_line=int(header[5]))
            f.close()
        else:
            print("No fc.info file, please create one")

    def compute_state_from_fc(self, thresh_pval=0.05, thresh_fc=0):
        '''
        Loads Fold-changes and p-values data and computes genes state from 
        these data. If no p-values is available in the data file, the default 
        value of 0 is assigned to each gene's pvalue.

        A gene is considered:
            * `activated` if its FC is above the FC threshold given as argument
              and its pvalue is below the pvalue threshold given as 
              argument ie. FC > thresh_FC and pval < thresh_pval
            * `repressed` if its FC is below the opposite of the FC threshold
              and its pvalue is below the pvalue threshold
              ie. FC < -thresh_FC and pval < thresh_pval
            * `not affected` either if its pvalue is above the threshold,
              or if its FC is between the - thresh_FC and + thresh_FC

        Creates:
           
            * self.statesFC (dict. of dict.)
                    New attribute of the  Transcriptome instance. 
                    Dictionary containing one subdictionary per condition 
                    listed in fc.info. Each subdictionary contains the list of 
                    genes corresponding to  each state ('act', 'rep', 'non' or 
                    'null').
            * self.genes[locus].fc_pval (dict.) 
                    new attribute of Gene instances
                    related to the Transcriptome instance given as argument.
                    Dictionary of shape {condition: (FC,pvalue)}.
            * self.genes[locus].state (dict.) 
                    new attribute of Gene instances
                    related to the Transcriptome instance given as argument.
                    Dictionary of shape {condition: state} with state either
                        * `act` if the gene is activated
                        * `rep` if the gene is repressed
                        * `non` if the gene is not affected
                        * `null` if the gene is not present in the data


        Args:
            thresh_pval (Float): pvalue threshold used for the genes 
                    classification
            thresh_fc (Float): fold-change threshold used for the genes 
                    classification
    
        Note:
            The FC and pvalues data importation requires an FC.info file, in the
            fold_changes directory, containing column indices of each information 
            in the data file, in the following order:
            [0] Condition [1] Filename [2] Locus_tag column [3] FC columns
            [4] Separator [5] File start line [6] P-value column
            See load_fc_pval for more details about the data importation.

        Warning:
            This method needs a genomic annotation. If no annotation is
            loaded, the load_annotation method with the default "sequence.gff3"
            file is computed. To use another annotation, please load an 
            annotation to your Transcriptome instance with the following commands 
            before using this method:
                >>> from GRATIOSA import Transcriptome, Genome
                >>> tr = Transcriptome.Transcriptome("ecoli")
                >>> g = Genome.Genome(tr.name)
                >>> g.load_annotation(annot_file=chosen_file)
                >>> tr.genes = g.genes

        Example:
            >>> tr = Transcriptome.Transcriptome("ecoli")
            >>> tr.compute_state_from_fc()
            >>> tr.genes['b0001'].state
            {'osmotic': 'rep', 'acidic_1mn': 'act'}
            >>> tr.genes['b0001'].fc_pval
            {'osmotic': (-1.73717046009437, 0.0), 'acidic_1mn': (1.73, 0.0)}
            >>> tr.statesFC['osmotic']['rep']
            ['b0001', 'b0002', 'b0003', 'b0004',...]
        '''

        if not hasattr(self, 'genes_valid_fc'):
            self.load_fc_pval()
        if not hasattr(self, "genes"):  # if no genes loaded
            # try to load them
            gen = Genome(self.name)
            gen.load_annotation()
            self.genes = gen.genes
        if not hasattr(self, "statesFC"):
            self.statesFC = {}

        for cond_fc in self.genes_valid_fc.keys():
            self.statesFC[cond_fc] = {"rep": [], "act": [], "non": []}
        for gene in self.genes.keys():
            g = self.genes[gene]
            for cond_fc in self.genes_valid_fc.keys():
                try:
                    if g.fc_pval[cond_fc][1] <= thresh_pval:
                        if g.fc_pval[cond_fc][0] < 0 - thresh_fc:
                            g.add_state_cond(cond_fc, 'rep')
                            self.statesFC[cond_fc]['rep'].append(gene)
                        elif g.fc_pval[cond_fc][0] > 0 + thresh_fc:
                            g.add_state_cond(cond_fc, 'act')
                            self.statesFC[cond_fc]['act'].append(gene)
                        else:
                            g.add_state_cond(cond_fc, 'non')
                            self.statesFC[cond_fc]['non'].append(gene)
                    else:
                        g.add_state_cond(cond_fc, 'non')
                        self.statesFC[cond_fc]['non'].append(gene)
                except BaseException:
                    if not hasattr(self.statesFC[cond_fc], 'null'):
                        self.statesFC[cond_fc]['null'] = [gene]
                    else:
                        self.statesFC[cond_fc]['null'].append(gene)
                    g.add_state_cond(cond_fc, 'null')

    def load_rnaseq_cov(self, cond="all", compute_from_bam=False):
        '''
        load_rnaseq_cov loads a RNASeq coverage to a Transcriptome instance
        either from:
            * .npz files (described in a cov.info file and computed when new 
              data are loaded for the first time)  which are in the /rnaseq_cov/ 
              directory,
            * coverage files which are in the /rnaseq_cov/ directory and are 
              described in cov_txt.info file,
            * paired-end .bam reads files which are in the /rnaseq_reads/ 
              directory  and are treated with 
              useful_functions_transcriptome.cov_from_reads function.

        Creates 2 new attributes of the Transcriptome instance:
            * self.rnaseq_cov_pos: 
                    dictionary of shape {condition: cov+}
                    with cov+ an array containing one signal data per genomic
                    position for the + strand (forward coverage)
            * self.rnaseq_cov_neg: 
                    idem for the - strand (reverse coverage)

        

        Args:
            cond (Optional [list of str.]): selection of one or several 
                    conditions (1 condition corresponds to 1 data file).
                    By default: cond ='all' ie all available coverages are 
                    loaded.
            compute_from_bam (Boolean.): if True, computes, using 
                    useful_functions_transcriptome.cov_from_reads, coverage 
                    from paired-end .bam reads files that are, with a 
                    bam_files.info file, in the /rnaseq_reads/ directory.
                    WARNING: Works only with paired-end files

        Note:
            To use directly new coverages data, both forward and reverse coverage 
            data files are required. These files contain the coverage for each 
            genomic position (one line = one genomic position, no position can be 
            ommited), therefore if you use bedtools genomecov, please use the -d 
            option.
            The importation of these new data requires a cov_txt.info 
            file, in the /rnaseq_cov/ directory, containing the following 
            information:
            [0] Condition
            [1] Reverse coverage filename
            [2] Forward coverage filename
            [3] Startline (has to be the same for both coverage files)
            [4] Column containing the coverage data

            Coverages from the reverse file will be added to the Transcriptome
            instance as "rnaseq_cov_neg" attribute and coverages from the forward
            file will be the "rnaseq_cov_pos"attribute.

        Note:
            To compute coverages data from reads files, the data files and a
            bam_files.info file have to be in the /rnaseq_reads/ directory.
            The bam_files.info files contains:
            [0] Condition [1] Reads filename

        Note:
            During the first importation of new coverages data, a .npy file
            containing the data is created and will enable a faster data loading
            for the next data importation. A cov.info file containing information
            for this fast importation is also in the /rnaseq_cov/ directory and
            contains in the following columns:
            [0] Condition [1] Coverage filename
            This cov.info file automatically completes itself when the data are
            loaded for the first time.

        Example:
            >>> tr = Transcriptome.Transcriptome("ecoli")
            # To load only the coverage for one condition named "WT"
            # in the cov.info or cov_txt.info file
            >>> tr.load_rnaseq_cov(["WT"]) 
            >>> tr.rnaseq_cov_neg["WT"]
            array([0., 0., 0., ..., 0., 0., 0.])
            >>> tr.rnaseq_cov_pos["WT"]
            array([0., 0., 0., ..., 0., 0., 0.])
            # To compute the coverage from a new .bam files obtained from 
            # paired-end data for example with: tophat2 genome/MG1655 
            # Test_1.fastq.gz Test_2.fastq.gz
            >>> tr = Transcriptome.Transcriptome("ecoli")
            >>> tr.load_rnaseq_cov(compute_from_bam=True)
            >>> tr.rnaseq_cov_pos["Test"]
            array([10., 10., 10., ..., 0., 0., 0.])
            >>> tr.rnaseq_cov_pos["Test"]
            array([0., 0., 0., ..., 20., 20., 20.])
        '''
        if not hasattr(self, "rnaseq_cov_pos"):
            self.rnaseq_cov_pos = {}  # cov on + strand
            self.rnaseq_cov_neg = {}  # cov on - strand

        # conversion from str to list (useful if only 1 condition is selected)
        if isinstance(cond, str):
            cond = [cond]

        path2dir = f"{basedir}data/{self.name}/rnaseq_cov/"

        if compute_from_bam:
            process_bam_paired_end(self)
            cov_from_reads(self)

        if os.path.exists(f"{path2dir}cov_txt.info"):
            if not os.path.exists(f"{path2dir}cov.info"):
                file = open(f"{path2dir}cov.info", 'w')
                file.write('Condition\tCov file\n')
            else:
                file = open(f"{path2dir}cov.info", 'a')
            with open(f"{path2dir}cov_txt.info", "r") as f:
                header = next(f)
                for line in f:
                    line = line.strip('\n').split('\t')
                    print(f'Updating cov.info with the new data: {line[0]}')
                    # create .npy from .txt
                    rnaseq_cov_neg = np.loadtxt(
                        path2dir + line[1], usecols=[int(line[4])], skiprows=int(line[3]))
                    rnaseq_cov_pos = np.loadtxt(
                        path2dir + line[2], usecols=[int(line[4])], skiprows=int(line[3]))
                    # save .npy into .npz
                    np.savez(
                        f"{path2dir}{line[0]}_cov.npz",
                        cov_pos=rnaseq_cov_pos,
                        cov_neg=rnaseq_cov_neg)
                    # update cov.info
                    file.write(f"{line[0]}\t{line[0]}_cov.npz\n")
            with open(f"{path2dir}cov_txt.info", "w") as f:
                f.write(header)
            f.close()
            file.close()

        if os.path.exists(f"{path2dir}cov.info"):
            with open(f"{path2dir}cov.info", "r") as f:
                header = next(f)
                loaded_cond = []
                for line in f:
                    line = line.strip().split('\t')
                    # For each selected condition:
                    if cond == ["all"] or line[0] in cond:
                        print('Loading condition', line[0])
                        loaded_cond.append(line[0])
                        # load attributes
                        self.rnaseq_cov_neg[line[0]] = np.load(
                            path2dir + line[1])["cov_neg"]
                        self.rnaseq_cov_pos[line[0]] = np.load(
                            path2dir + line[1])["cov_pos"]
                if cond != ["all"]:
                    unloaded = set(loaded_cond) ^ set(cond)
                    if len(unloaded) != 0:
                        print(f"{unloaded} data were not found, please check .info files")
            f.close()

        else:
            print(
                'cov.info not available nor cov_txt.info, please check /rnaseq_cov/ folder')

    def load_cov_start_end(self, compute_from_bam=False):
        '''
        Useful to locate TSS and TTS, load_cov_start_end loads the density of
        RNA fragment starts and ends from .npz files described in
        cov_start_stop.info, in the /rnaseq_cov/ directory.
        .npz files and .info files are obtained using the
        useful_functions_transcriptome.cov_start_stop_from_reads

        Creates 2 new attributes of the Transcriptome instance:
            * self.cov_start (Dict. of dict.)
                    dictionary containing one subdictionary per condition.
                    Each subdictionary has 2 keys:
                    "0" for the start positions on the - strand and
                    "1" for the start positions on the + strand
                    and the subdictionary values are arrays containing
                    all start positions on the corresponding strand
            * self.cov_end (Dict. of dict.)
                    idem for the end positions

        Args:
            compute_from_bam (Boolean): if True, computes, using 
                    useful_functions_transcriptome.cov_start_stop_from_reads,
                    coverage from paired-end .bam reads files that are, with 
                    a bam_files.info file, in the /rnaseq_reads/ directory.

        Note:
            To compute new coverages data from paired-end .bam reads files, the 
            data files and a bam_files.info file have to be in the /rnaseq_reads/ 
            directory and the input argument "compute_from_bam" has to be set to 
            "True".The bam_files.info files contains:
            [0] Condition [1] Reads filename
            During the importation of these new data, a corresponding .npz files
            is created and the cov_start_stop.info file automatically completes 
            itself. For the next importations of these data, only the npz. and 
            the cov_start_stop.info will be loaded, thus allowing faster loading.

        Warning:
            compute_from_bam works only with paired-end files

        Example:
            >>> tr = Transcriptome.Transcriptome("ecoli")
            >>> tr.load_cov_start_end(compute_from_bam=True)
            >>> tr.cov_end["WT_paired"]
            {0: array([0, 0, 0, ..., 0, 0, 0]),
             1: array([0, 0, 0, ..., 0, 0, 0])}
        '''

        if compute_from_bam:
            process_bam_paired_end(self)
            cov_start_stop_from_reads(self)

        self.cov_start = {}
        self.cov_end = {}
        path2dir = f"{basedir}data/{self.name}/rnaseq_cov/"
        if os.path.exists(f"{path2dir}cov_start_stop.info"):
            with open(f"{path2dir}cov_start_stop.info", "r") as f:
                header = next(f)
                for line in f:  # for each condition
                    line = line.strip().split('\t')
                    print('Loading condition', line[0])
                    # load attributes
                    self.cov_start[line[0]] = {}
                    self.cov_end[line[0]] = {}
                    self.cov_start[line[0]][0] = np.load(
                        path2dir + line[1])["cov_start_neg"]
                    self.cov_start[line[0]][1] = np.load(
                        path2dir + line[1])["cov_start_pos"]
                    self.cov_end[line[0]][0] = np.load(
                        path2dir + line[1])["cov_end_neg"]
                    self.cov_end[line[0]][1] = np.load(
                        path2dir + line[1])["cov_end_pos"]
            f.close()

        else:
            print('cov_start_stop.info not available please check /rnaseq_cov/ folder')

    def compute_log2rpkm_from_cov(self, cond='all', before=100):
        '''
        Compute_log2rpkm_from_cov method:
            * loads the RNASeq coverage (using load_rnaseq_cov method)
            * computes and loads RPKM value on the Gene instances
              (using Gene.add_single_rpkm )
            * saves the log2(RPKM) data in a new .csv file
              (1 file for all conditions, with 1 column per condition)
            * adds the new .csv file informations in the expression.info file 
              in the /expression/ directory

        Creates:
            * self.genes[locus].rpkm (dict.):
                    new attribute of Gene instances related to the 
                    Transcriptome instance given as argument.
                    Dictionary of shape {condition: RPKM value}, containing 
                    the RPKM for each computed condition.

            * log2rpkm_from_cov_{datetime.now()()}.csv: csv file containing the
                    log2(RPKM) values (1 row per gene, 1 column per condition).
                    If a RPKM is equal to 0, log2RPKM is set to 1.
                    This file is saved in the /expression/ directory and is
                    automatically added in the expression.info file.

        Args: 
            cond (list of str.):
                    selection of one or several RNASeq coverage condition
                    By default: cond ='all' ie all available coverages are
                    converted. All selected conditions have to be listed in
                    cov.info file. If one of the chosen condition is already
                    registered in the expression.info file in the expression/
                    directory, this condition is ignored.
            before (int.):
                    number of bps before the gene start to take into account

        Warning: 
            This method needs a genomic annotation. If no annotation is
            loaded, the load_annotation method with the default "sequence.gff3"
            file is computed. To use another annotation, please load an 
            annotation to your Transcriptome instance with the following commands 
            before using this method:
                >>> from GRATIOSA import Genome, Transcriptome
                >>> tr = Transcriptome.Transcriptome("ecoli")
                >>> g = Genome.Genome(tr.name)
                >>> g.load_annotation(annot_file=chosen_file)
                >>> .tr.genes = g.genes

        Example:
            >>> tr = Transcriptome.Transcriptome("ecoli")
            >>> tr.compute_log2rpkm_from_cov(cond="WT")
            >>> tr.genes["b0002"].rpkm["WT"]
            142
        '''

        # Loading the annotation and the RNASeq coverage
        if not hasattr(self, "genes"):
            gen = Genome(self.name)
            gen.load_annotation()
            self.genes = gen.genes
        if not hasattr(self, "rnaseq_cov_pos"):
            self.load_rnaseq_cov()

        # Checking if each condition is in cov.info and not in expression.info
        path2dir = f"{basedir}data/{self.name}/expression/"
        if isinstance(cond, str):
            cond = [cond]
        if cond == ["all"]:
            cond = self.rnaseq_cov_pos.keys()
        else:
            condpb = set(cond) - set(self.rnaseq_cov_pos.keys())
        if condpb:
            print(f"Following conditions are not in cov.info: {condpb}")

        if os.path.exists(f"{path2dir}expression.info"):
            existing_cond = []
            with open(f"{path2dir}expression.info", 'r')as f:
                header = next(f)
                for line in f:
                    line = line.strip('\n').split('\t')
                    existing_cond.append(line[0])
            f.close()

        cond2convert = list(set(cond) - set(existing_cond) - condpb)
        ignoredcond = set(cond) & set(existing_cond)
        if ignoredcond:
            print(f"Following conditions are already in expression.info: {ignoredcond})")

        # Computing the RPKM
        if cond2convert:
            for g in self.genes.keys():  # for each gene
                try:
                    if self.genes[g].strand:
                        # gene in + strand
                        for c in cond2convert:  # for each condition of cov
                            self.genes[g].add_single_rpkm(c, np.mean(self.rnaseq_cov_pos[c][(
                                self.genes[g].left - before):self.genes[g].right]), np.sum(self.rnaseq_cov_pos[c]) + np.sum(self.rnaseq_cov_neg[c]))
                    elif not self.genes[g].strand:
                        # gene in - strand
                        for c in cond2convert:
                            self.genes[g].add_single_rpkm(c, np.mean(self.rnaseq_cov_neg[c][self.genes[g].left:(
                                self.genes[g].right + before)]), np.sum(self.rnaseq_cov_pos[c]) + np.sum(self.rnaseq_cov_neg[c]))
                except Exception as e:
                    print(f'WARNING: {g} , {cond}, {e}')
                    pass

        # Saving log2RPKM in a .csv file and completing the expression.info
        # file
            file = []
            for g in self.genes.keys():
                g_log2rpkm = []
                for c in cond2convert:
                    if self.genes[g].rpkm[c] != 0:
                        g_log2rpkm.append(np.log2(self.genes[g].rpkm[c]))
                    else:
                        g_log2rpkm.append(1)
                file.append(g_log2rpkm)
            df = pd.DataFrame(data=file, columns=cond2convert,
                              index=self.genes.keys())
            date = str(datetime.now())[:-10].replace(" ", "_")
            filename = f"log2rpkm_from_cov_{date}.csv"
            df.to_csv(path2dir + filename, sep=',', index=True)

            if not os.path.exists(f"{path2dir}expression.info"):
                f = open(f"{path2dir}expression.info", 'w')
                f.write('Condition\tFilename\tTagColumn\tExprColumn\tIslog(yes/no)\tSeparator\n')
            else:
                f = open(f"{path2dir}expression.info", 'a')

            for n, c in enumerate(cond2convert):
                f.write(f'{c}\t{filename}\t0\t{n+1}\tyes\t,\n')
        else:
            print("All selected conditions are already in expression.info")