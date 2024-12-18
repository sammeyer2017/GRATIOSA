#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions called by Genome methods
"""
import os
import numpy as np
from GRATIOSA.globvar import *
from GRATIOSA.Gene import Gene
from GRATIOSA.TSS_TTS_TU import TSS, TTS, TU


def read_seq(filename):
    '''
    Called by load_seq, read_seq allows name(s) and genomic sequence(s) to be loaded from a
    .fasta file

    Args:
        filename (str.): name of the file containing the DNA sequence in
                         FASTA format.

    Returns:
        tuple with list of names (str.) and list of genomic sequences (str.)

    '''
    names=[]
    seqs=[]
    seq_file = open(filename, "r")
    seq=''
    for i in seq_file.readlines():
        line = i.strip()  # Removes \n
        if line != '':  # Inspect if empty line
            if line[0] == ">":
                names.append(line[1:])
                seqs.append(seq)
                seq=""
            else:
                seq += line
    seq_file.close
    seqs.append(seq)
    return names,seqs[1:]

def seq_compl(seq):
    l=seq
    l = l.replace('A', 't').replace('T', 'a')
    l = l.replace('C', 'g').replace('G', 'c')
    l = l.replace('a', 'A').replace('t', 'T')
    l = l.replace('c', 'C').replace('g', 'G')
    return l

    
def update_NCBI_genomes():
    """
    Tries to download the current list of NCBI complete bacterial genomes, and store it in the file "accession_number_complete.txt" in the base directory. This is useful if you want to create a new organism and download the annotation automatically from NCBI, and you do not have the NCBI command-line tool installed... Caution, this requires to download a large file (~ 1 Gb). The stored file contains only the complete genomes and is much lighter (~ 20 Mb). 
    """
    print("Downloading the current version of the NCBI complete genome list to accession_number_complete.txt file in the base directory")
    print("It may take a while, the complete file is ~1Gb to download")
    print("Check the large temporary file at location %s/assembly_summary.txt"%os.getcwd())
    #os.system("wget -nv ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt >/dev/null 2>&1")
    os.system("curl -o assembly_summary.txt ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/assembly_summary.txt")
    if not os.path.exists("assembly_summary.txt"):
        print("Error downloading the file!!")
        return 1
    else:
        os.system(f"grep Complete assembly_summary.txt >  {basedir}data/assembly_summary_complete.txt")
        os.system("rm assembly_summary.txt")
        return 0

def load_gff(annotations_filename, genome_dict, features=["gene"]):
    '''
    Called by load_annotation, load_gff allows genes to be loaded from .gff 
    or .gff3 annotation file, using the defined list of features and using 
    the "locus_tag" as primary field to name the genes, with "gene" being the 
    secondary field if locus_tag does not exist. For each gene of the
    annotation, one Gene object is created and intialized with the following attributes: 
    locus_tag, ID, name, left (left coordinate), right (right coordinate), 
    middle, length, strand, start, end, ASAP ID (legacy), and genome.

    Args:
        filename (str.): name of the annotation file.
        genome_dict (dict.): dictionary giving the genome object associated to 
        each chromosome name found in the gff file
        features (list of str.): list of attributes (column 3) to be imported: 
        usually ["gene"] (default), but can be adapted, e.g., for "exon", "rRNA", 
    etc. A feature will be retained if it contains a 

    Returns:
        Dictionary: Dictionary of shape {locus: Gene object} with each Gene 
        object initialized with the following attributes:

            * locus_tag: "locus_tag" recorded in the annotation,
            * ID: "ID" or "locus_tag" if no "ID" is associated this gene in 
              the annotation,
            * name: "gene" or "locus_tag" if no "gene" is associated to this
              gene in the annotation,
            * ASAP: "ASAP" name (By default: '').
            * strand: gene strand,
            * left, right: gene coordinates (does not take into account the
              strand, ie right > left)
            * start, end, middle: positions of the beginning, the middle and
              the end of the gene
            * length: gene length (=right-left)
            * genome: genome object where this gene belongs (useful only for 
              genomes with multiple chromosomes)

    Warnings: 
        ASAP name  is found only if it is noted on the same line as the
        "name" line or the line next to it (legacy of an annotation file 
        of Dickeya dadantii)
    '''
    genes_dict = {}
    annot_file = open(annotations_filename, "r")
    for line in annot_file.readlines():
        ASAP_name = ''
        if line[0] != '#':
            if line != '\n':
                line = line.split('\t')

                if line[2] in features:
                    gene_info = {}
                    for x in line[8].split(';'):
                        try : 
                            x = x.strip().split('=')
                            gene_info[x[0]] = x[1]
                        except Exception as e: 
                            if  x != ['']:
                                print(f"Error {e} in '{line[8]}'")

                    if ('locus_tag' in gene_info.keys()):
                        locus = gene_info['locus_tag']
                    elif ("gene" in gene_info.keys()):
                        locus = gene_info["gene"]
                    else:
                        try:
                            del locus
                        except:
                            pass

                    if "locus" in locals():# isinstance(locus,str):
                        # the gene will be recorded
                        try:
                            ge=genome_dict[line[0]]
                        except:
                            print("problem with gene %s: chromosome name (column 1) does not match any known chromosome name for the organism. Check consistency between the reference sequence and annotation file"%locus)
                            return 1
                        if "ID" in gene_info.keys():
                            ID = gene_info["ID"]
                        else:
                            ID = locus
                        if "gene" in gene_info.keys():
                            name = gene_info["gene"]
                        elif "gene_synonym" in gene_info.keys():
                            name = gene_info["gene_synonym"]
                        else:
                            name = locus
                        if "product" in gene_info.keys():
                            product = gene_info["product"]
                        else:
                            product = ""
                        if "Dbxref" in gene_info.keys():
                            db_ref = "Dbxref"
                        elif "db_xref" in gene_info.keys():
                            db_ref = "db_xref"
                        else:
                            db_ref = None
                        if db_ref is not None:
                            dbxref = gene_info[db_ref].split(",")
                            for db in dbxref:
                                if "ASAP" in db:
                                    ASAP_name = db.split(":")[1]
                        left, right = map(int, line[3:5])
                        strand = True if line[6] == "+" else False
                        notes = line[8]
                        feature=line[2]
                        genes_dict[locus] = Gene(locus, feature, name, ID, ge, left, right, strand, ASAP_name, product, notes)
                    """
                    elif line[2] == "CDS" and ASAP_name == '':
                        tleft, tright = map(int, line[3:5])
                        if tleft == left and tright == right:
                            if "Dbxref" in gene_info.keys():
                                db_ref = "Dbxref"
                            elif "db_xref" in gene_info.keys():
                                db_ref = "db_xref"
                            else:
                                db_ref = None
                            if db_ref is not None:
                                dbxref = gene_info[db_ref].split(",")
                                for db in dbxref:
                                    if "ASAP" in db:
                                        ASAP_name = db.split(":")[1]
                            genes_dict[locus] = Gene(
                                locus, name, ID, ge, left, right, strand, ASAP_name)
                    """
    return genes_dict


def load_annot_general(annotations_filename, 
                       separator, 
                       tag_column, 
                       name_column,
                       ID_column,
                       genome_column,
                       strand_column, 
                       left_column, 
                       right_column, 
                       ASAP_column, 
                       start_line,
                       features=["gene"]):
    '''
    Legacy function, not very useful. Called by load annotation, load_annot_general 
    allows genes to be loaded from any csv annotation file, not following gff format. 
    This is almost never useful, as annotation is usually loaded
    from gff files (using the load_gff function). The annotation file has to
    contains one row per genes and the following columns: locus_tag, ID, 
    name, ASAP, left coordinate, right coordinate and strand.
    For each gene of the annotation, one Gene object is created and 
    intialized with the following attributes: locus_tag, ID, name, left 
    (left coordinate), right (right coordinate), middle, length, strand, 
    start, end and ASAP.

    Args:
        filename (str.): name of the annotation file
        separator (str.): file separator
        tag_col (int.): index of the column containing the locus tags 
        name_column (int.):index of the column containing the names 
        ID_column (int.): index of the column containing the ID
        strand_column (int.): index of the column containing the strand 
        left_column (int.): index of the column containing the left 
                coordinates
        right_column (int.): index of the column containing the right 
                coordinates
        ASAP_column (int.):index of the column containing the ASAP names
        start_line (int.): file start line
    
    Returns:
        Dictionary: Dict. of shape {locus: Gene object} with each Gene object
        initialized with the following attributes: locus_tag, ID, name,
        left (left coordinate), right (right coordinate), middle, length, 
        strand, start, end and ASAP.

    Note: 
        Column numbering starts at 0.

        Strands can be noted with one of the following writings:

            * Forward: "forward","1","+","true","plus" (case insensitive)
            * Reverse: "complement","-1","-","false","minus" (case insensitive)
    '''
    genes_dict = {}
    with open(annotations_filename, 'r') as f:
        head = next(f)
        headl = head.strip()
        if separator == '\\t':
            separator = '\t'
        headl = headl.split(separator)
        i = 0
        while i < start_line:
            header = next(f)
            i += 1
        for line in f:
            try:
                line = line.strip()
                line = line.split(separator)
                locus = line[tag_column]
                name = line[name_column]
                ID = line[ID_column]
                ASAP = line[ASAP_column]
                if line[strand_column].lower() in ["complement", "-1", "-", "false", "minus"]:
                    strand = False
                elif line[strand_column].lower() in ["forward", "1", "+", "true", "plus"]:
                    strand = True
                else:
                    print("Unrecognized strand type")
                    sys.exit("Unrecognized strand type")
                left = line[left_column]
                right = line[right_column]
                genes_dict[locus] = Gene(
                    locus, name, ID, left, right, strand, ASAP)
            except Exception as e:
                if line != ['']:
                    print(f"WARNING: Error with {line}")
    return genes_dict


def load_TSS_cond(genes_dict, 
                  filename, 
                  TSS_column, 
                  start_line, 
                  separator,
                  strandcol, 
                  genescol=None, 
                  sigcol=None, 
                  sitescol=None, 
                  scorecol=None,):
    '''
    Called by load_TSS, load_TSS_cond allows TSS data to be loaded from any 
    file with one row per TSS and at least the following columns (separated 
    by a separator that is not a comma): TSS position and DNA strand. 
    Additional columns can be locus tags, sigma factors, sites positions and 
    score. The sites positions columns is divided in 4 part separated by 
    commas: -10 element left coordinate, -10 element right coordinate, 
    -35 element left coordinate and -35 element right coordinate.
    For each TSS of the annotation, one TSS object is created and intialized 
    with the attributes available in the data file (at least pos and strand).

    Args:
        genes_dict (dict.): dictionary of shape {locus tag: Gene object}
        filename (str.): name of the annotation file
        TSS_column (int.): index of the column containing the TSS positions
        start_line (int.): file start line
        separator (str.): file separator (other than commas !)
        strandcol (int.): index of the column containing the strands in the 
                file
        genescol (int.): index of the column containing the locus tags of the
                genes associated to each TSS in the file. By default: None 
                (ie not in the file)
        sigcol (int.): index of the column containing the name of the sigma 
                factor associated to each TSS in the file
                By default: None (ie not on file)
        sitescol (int.): index of the column containing the sites positions
                         (ie the coordinates of the -10 and the -35 elements)
                         associated to each TSS in the file
                         By default: None (ie not on file)
        scorecol (int.): index of the TSS scores column in the file
                         By default: None (ie not on file)


    Returns:
        Dictionary: dict. of shape {TSS position: TSS object} with each TSS object
        initialized with, at least, the following attributes: pos and strand.
        Depending on the data present in the file and passed as arguments,
        the following attributes may also have been added: genes, score and
        promoter. The promoter attribute is a dictionary containing, for each
        sigma factor (keys) a subdictionary (value). The first created key of
        this subdictionary is "sites" and the associated value is a tuple
        containing the positions of promoter elements (-10l,-10r,-35l,-35r)
        with l = left coordinate and r = right coordinate.

    Note: 
        Column numbering starts at 0.

        Strands can be noted with one of the following writings:
        
            * Forward: "forward","1","+","true","plus" (case insensitive)
            * Reverse: "complement","-1","-","false","minus" (case insensitive)
    '''
    TSS_dict = {}  # dict of TSS objects
    with open(filename, 'r') as f:
        i = 0
        while i < start_line:
            header = next(f)
            i += 1
        for line in f:
            line = line.strip('\n')
            if separator == '\\t':
                line = line.split('\t')
            else:
                line = line.split(separator)
            try:
                pos = int(line[TSS_column])
                if line[strandcol].lower() in ["complement", "-1",
                                               "-", "false", "minus"]:
                    strand = False
                elif line[strandcol].lower() in ["forward", "1", "+", "true", "plus"]:
                    strand = True
                else:
                    strand = None
                genes = line[genescol] if genescol is not None else []
                sig = line[sigcol] if sigcol is not None else None
                sites = line[sitescol] if sitescol is not None else None
                score = line[scorecol] if scorecol is not None else None

                # if TSS needs to be init
                if pos not in TSS_dict.keys():
                    # init object tss
                    TSS_dict[pos] = TSS(pos=pos, strand=strand, score=score)
                    if genes != []:
                        TSS_dict[pos].add_genes(genes, genes_dict)
                        # add TSS to gene attributes
                        for gene in TSS_dict[pos].genes:
                            genes_dict[gene].add_id_TSS(pos)

                # Add sigma factor and binding sites to the promoter dict
                if sig is not None:  # if sigma column
                    if sites is not None:
                        TSS_dict[pos].add_promoter(sig, sites=sites)
                    else:
                        TSS_dict[pos].add_promoter(sig)

            except Exception as e:
                print('Error in line, wrong information type:', e)

    return TSS_dict


def load_TTS_cond(filename, 
                  separator, 
                  start_line, 
                  leftcol, 
                  rightcol, 
                  strandcol,
                  rhocol, 
                  seqcol=None, 
                  scorecol=None, 
                  genescol=None, 
                  *args, **kwargs):
    '''
    Called by load_TTS, load_TTS_cond allows TTS data to be loaded from any 
    file with one row per TTS and at least the following columns:
    left coordinate, right coordinate, strand and rho_dpdt. Additional 
    columns can be genes locus tags, TTS sequence and TTS score.
    For each TTS of the annotation, one TTS object is created and intialized 
    with the attributes available in the data file (at least left, right, 
    rho_dpdt and strand).

    Args:
        filename (str.): name of the annotation file
        separator (str.): file separator
        start_line (int.): file start line
        leftcol (int.): index of the column containing the TTS left 
                coordinates
        rightcol (int.): index of the column containing the TTS right 
                coordinates
        strandcol (int.): index of the column containing the TTS strands
        rhocol (int.): index of the column containing the information about
                the rho dependency of the TTS
        seqcol (int.): index of the column containing the TTS sequences
                By default: None (ie not on file)
        scorecol (int.): index of the TTS scores column in the file
                By default: None (ie not on file)
        genescol (int.): index of the column containing the locus tags of the
                genes associated to each TTS in the file.
                By default: None (ie not on file)

    Returns:
        Dictionary: Dict. of shape {TTS position: TTS object} with each TTS object
        initialized with, at least, the following attributes: left,
        right, start, end, strand and rho_dpdt. Depending on the data 
        present in the file and passed as arguments, the following attributes 
        may also have been added: genes, score and seq.

    Note: 
        Column numbering starts at 0.

        Strands can be noted with one of the following writings:
        
            * Forward: "forward","1","+","true","plus" (case insensitive)
            * Reverse: "complement","-1","-","false","minus" (case insensitive)

        Rho dependent TTS can be noted with one of the following writings: 
        "True","1" (case insensitive)
    '''

    TTS_dict = {}  # dict of TTS objects
    with open(filename, 'r') as f:
        i = 0
        while i < start_line:
            header = next(f)
            i += 1
        for line in f:
            line = line.strip('\n')
            if separator == '\\t':
                line = line.split('\t')
            else:
                line = line.split(separator)
            try:
                left = int(line[leftcol])
                right = int(line[rightcol])
                if line[strandcol].lower() in ["complement", "-1", "-", "false", "minus"]:
                    strand = False
                elif line[strandcol].lower() in ["forward", "1", "+", "true", "plus"]:
                    strand = True
                else:
                    strand = None
                rho_dpdt = True if line[rhocol].lower() in [
                    "true", "1"] else False
                seq = line[seqcol] if seqcol is not None else ""
                score = line[scorecol] if scorecol is not None else None
                genes = line[genescol] if genescol is not None else []

                newTTS = TTS(
                    left=left,
                    right=right,
                    strand=strand,
                    rho_dpdt=rho_dpdt,
                    seq=seq,
                    score=score,
                    genes=genes)
                TTS_dict[newTTS.start] = newTTS

            except Exception as e:
                print('Error in line, wrong information type:', e)

    return TTS_dict


def load_TU_cond(filename, 
                 IDcol,
                 startcol, 
                 endcol, 
                 strandcol,
                 startline, 
                 separator,
                 genescol=None, 
                 exprcol=None,
                 TSScol=None,
                 TTScol=None):
    '''
    Called by load_TU, load_TU_cond allows TU data to be loaded from any file
    with one row per TU and the following columns: start position, end 
    position, strand and genes locus tags. For each TU of the annotation, one 
    TU object is created and intialized with the attributes start, end, 
    strand, genes and expression.

    Args:
        filename (str.): name of the annotation file
        IDcol (int.): index of the column containing the TU ID
        startcol (int.): index of the column containing the TU start positions
        endcol (int.): index of the column containing the TU end positions
        strandcol (int.): index of the column containing the TU strands
        start_line (int.): file start line
        separator (str.): file separator
        genescol (int.): index of the column containing the locus tags of the 
                genes associated to each TU in the file, separated by commas. 
                By default: None (ie not on file)
        exprcol (int.): index of the column containing the TU expression (def. None)
        TSScol (int.): index of the column with TSS position
        TTScol (int.): index of the column with TTS position

    Returns:
        Dictionary: dict. of shape {TU start: TU object} with each TU object
        initialized with the attributes start, end, strand and genes.

    Note: 
        Column numbering starts at 0.

        Strands can be noted with one of the following writings:
        
            * Forward: "forward","1","+","true","plus" (case insensitive)
            * Reverse: "complement","-1","-","false","minus" (case insensitive)
    '''

    TUs = {}
    with open(filename, 'r') as f:
        i = 0
        while i < startline:
            header = next(f)
            i += 1
        for line in f:
            line = line.strip('\n')
            if separator == '\\t':
                line = line.split('\t')
            else:
                line = line.split(separator)
            try:
                if line[strandcol].lower() in ["complement", "-1", "-", "false", "minus"]:
                    strand = False
                elif line[strandcol].lower() in ["forward", "1", "+", "true", "plus"]:
                    strand = True
                else:
                    strand = None
                if genescol != None:
                    genesnames=line[genescol].split(",")
                else:
                    genesnames = None
                if exprcol != None:
                    expr = float(line[exprcol])
                else :
                    expr = None
                if TSScol != None:
                    # we take out everything that is not a number
                    TSSpos = int("".join([ele for ele in line[TSScol] if ele.isdigit()]))
                    #TSSpos = int(line[TSS])
                else:
                    TSSpos = None
                if TTScol != None:
                    TTSpos = int("".join([ele for ele in line[TTScol] if ele.isdigit()]))
                    #TTSpos = int(line[TTS])
                else:
                    TTSpos = None
                TUs[int(line[IDcol])] = TU(start=int(line[startcol]), 
                                            end=int(line[endcol]), 
                                            strand=strand, 
                                            genes=genesnames,
                                            expression=expr,
                                            TSS=TSSpos,
                                            TTS=TTSpos)
            except Exception as e:
                print(e)
    f.close()
    return TUs
