#! /usr/bin/env python
# -*- coding: utf-8 -*-
from GRATIOSA.globvar import *
from datetime import datetime
from GRATIOSA import stat_analysis
import os
import re
import pandas as pd

class GO:
    '''
    The GO (Gene Ontology) class contains the information (names, definitions) 
    of the GO terms of the different categories (biological process, molecular 
    function and cellular component). The associated function 
    (GO_enrichment_test) allows the study of functional enrichment.
    '''

    def __init__(self,
                 filename="go-basic.obo",
                 path2dir=f"{basedir}data/",
                 obo_reload=False):
        """
        Args:
            filename (str.): name of the Gene Ontology .obo file
            path2dir (str.): path to the directory containing the .obo file
            obo_reload (Bool.) : if true, the Gene Ontology basic.obo file is
                    reloaded. This file allows each GO term to be associated 
                    with a category (GOp - biological_process, GOf - 
                    molecular_function or GOc - cellular_component), a name 
                    and a description

        Initializes GO object with the following attributes:
            dict_GO (dict. of dict.) 
                    dictionary containing one subdictionary
                    per GO term. Each subdictionary contains the "Name", the 
                    "Namespace" (biological_process, molecular_function or 
                    cellular_component) and the definition corresponding to 
                    the GO term.
            dict_GOc (dict. of dict.) 
                    dictionary containing only the GOterms
                    from "cellular_component" namespace. This dictionary 
                    contains one subdictionary per GO term, with the 
                    following shape self.dict_GOc["GOterm1"] = {"Name": GO 
                    name, "Definition": associated definition}
            dict_GOf (dict. of dict.) 
                    idem for the GOterms from "molecular_function" namespace
            dict_GOp (dict. of dict.) 
                    idem for the GOterms from "biological_process" namespace

        Example:
            >>> go = GO.GO()
            >>> go.dict_GOc['GO:0000015']
            {'Name': 'phosphopyruvate hydratase complex',
             'Definition': 'A multimeric enzyme complex that catalyzes the 
             conversion of 2-phospho-D-glycerate to phosphoenolpyruvate and 
             water. [GOC:jl, ISBN:0198506732]'}
        """
        if obo_reload:
            os.system(f"wget http://purl.obolibrary.org/obo/go/go-basic.obo \
                -P {basedir}data/")

        self.dict_GO = {}
        with open(path2dir + filename, "r") as f:
            for line in f:
                line = line.strip("\n").replace("\"", "").split(": ")
                if "id" in line:
                    if "GO" in line[1]:
                        ID = line[1]
                    else:
                        ID = "other"
                    self.dict_GO[ID] = {}
                elif "name" in line:
                    self.dict_GO[ID]["Name"] = line[1]
                    # print(self.dict_GO)
                elif "namespace" in line:
                    self.dict_GO[ID]["Namespace"] = line[1]
                elif "def" in line:
                    self.dict_GO[ID]["Definition"] = line[1]

        if 'other' in self.dict_GO.keys():
            self.dict_GO.pop('other')

        self.dict_GOc, self.dict_GOp, self.dict_GOf = {}, {}, {}
        for GO in self.dict_GO.keys():
            if self.dict_GO[GO]["Namespace"] == "biological_process":
                self.dict_GOp[GO] = {"Name": self.dict_GO[GO]["Name"],
                                     "Definition": self.dict_GO[GO]["Definition"]}
            elif self.dict_GO[GO]["Namespace"] == "molecular_function":
                self.dict_GOf[GO] = {"Name": self.dict_GO[GO]["Name"],
                                     "Definition": self.dict_GO[GO]["Definition"]}
            elif self.dict_GO[GO]["Namespace"] == "cellular_component":
                self.dict_GOc[GO] = {"Name": self.dict_GO[GO]["Name"],
                                     "Definition": self.dict_GO[GO]["Definition"]}


def GO_enrichment_test(dict_GOterms,
                       dict_features,
                       targ_features,
                       all_features="all",
                       targ_cats="all",
                       all_cats="all",
                       min_nb_elements=4,
                       thresh_padj=0.1,
                       output_dir=f"{resdir}enrichment_test/",
                       output_file=f"enrich_test{datetime.now()}",
                       obo_reload=False,
                       GO_filename="go-basic.obo",
                       GO_path2dir=f"{basedir}data/"):
    '''
    Computes enrichment tests (hypergeometric test) of a target in a
    category.
    N.B.: performs a p-value correction for false discovery rate (See
          statsmodels.stats.multitest.fdrcorrection for more details)

    Required args:
        dict_GOterms (dict.): classification of each element in a dictionary 
                of shape {category:[elements]}. N.B.: the same element can be 
                associated to multiple features
                Example: {"GOterm1":["GeneA","GeneB"],
                          "GOterm2:["GeneA","GeneC"]}
        dict_features (dict.): feature corresponding to each element in a
                dictionary of shape {feature:[elements]}
                Example: {"act": ["GeneA","GeneC","GeneD"],"rep":["GeneB"]}
        targ_features (list or string): targeted feature(s)

    Optional args:
        all_features (list or string): features that will be used to calculate
                the proportion, including targ_features (default: all keys of 
                the dict_features dictionary)
        targ_cats (list): list of categories. The enrichment test is 
                performed for each category. (default: all keys of dict_cats)
        all_cats (list): list of categories used to compute the global
                proportion and the expected number in the selection, 
                including targ_cats.(default: all keys of dict_cats)
        min_nb_elements: Number of elements used as thresholds for the feature
                selection. If there is strictly less than min_nb_elements
                corresponding to a feature, the result corresponding to this 
                feature is not relevant and is therefore neither returned 
                nor reported in the output file.
                By default, min_nb_elements is set to 4.
        thresh_padj (float.) : threshold on the adjusted p-value, used to 
                create the file containing only the significantly enriched
                GO terms
        output_dir (str.): output directory
        output_file (str.): output filename
        GO_filename (str.): name of the Gene Ontology .obo file
        GO_path2dir (str.): path to the directory containing the .obo file
        obo_reload (Bool.) : if true, the Gene Ontology basic.obo file is 
                reloaded. This file allows each GO term to be associated with 
                a category (GOf, GOp or GOc), a name and a description
    Returns:
        3 files (one for GOf, one for GOp and one for GOc) containing the
        following columns:
            0. 'Term'(str.): GO term
            1. 'Selected_gene_nb'(int.): Nb of elements corresponding to this
                    feature in the selection
            2. 'Expected_selected_nb' (int.): Expected nb of elements 
                    corresponding to this feature in the selection
            3. 'Total_gene_nb'(int.): Nb of elements corresponding to this
                    feature in the dict_features
            4. 'Proportion'(float): ratio between Selected_gene_nb and 
                    Total_gene_nb
            5. 'Prop_conf_int' (np.array): 95% confidence interval with equal
                    areas around the proportion.
            6. 'p-value' (float): p-value obtained with the enrichment test
            7. 'Adj p-value (FDR)' (float): p-value corrected for false
                    discovery rate
            8. 'Global_proportion' (float): ratio between nb of elements in 
                    the selection and nb of elements in dict_features
            The data are ordered according to the adjusted pvalues.
        1 common file containing the significantly enriched GO terms

    Example:
        >>> g = Genome.Genome("ecoli")
        >>> g.load_GO()
        >>> g.load_state_from_FC()
        >>> GO.GO_enrichment_test(g.GO['GO'],g.statesFC['Blot'],'act')
    '''
    # Formats the objects given as arguments as lists containing no duplicates.
    dGO = GO(filename=GO_filename, path2dir=GO_path2dir, obo_reload=obo_reload)

    list_GO = ["GOc", "GOp", "GOf"]

    dGOp, dGOc, dGOf = {}, {}, {}

    unknown_terms = []
    for t in dict_GOterms.keys():
        try:
            sp = dGO.dict_GO[t]['Namespace']
            if sp == "biological_process":
                dGOp[t] = dict_GOterms[t]
            elif sp == "molecular_function":
                dGOf[t] = dict_GOterms[t]
            elif sp == "cellular_component":
                dGOc[t] = dict_GOterms[t]

        except BaseException:
            unknown_terms.append(t)
    if len(unknown_terms) != 0:
        print(f"{len(unknown_terms)} GO terms are not in the GO annotation")

    for n, dict_cats in enumerate([dGOc, dGOp, dGOf]):
        o_file = f"{output_file}_{list_GO[n]}"
        df_res = pd.DataFrame(stat_analysis.enrichment_test(dict_cats=dict_cats,
                                                            dict_features=dict_features,
                                                            targ_features=targ_features,
                                                            all_features=all_features,
                                                            targ_cats=targ_cats,
                                                            all_cats=all_cats,
                                                            min_nb_elements=min_nb_elements,
                                                            output_dir=output_dir,
                                                            output_file=o_file))
        terms = df_res["Category"]
        names = []
        defs = []
        for t in terms:
            names.append(dGO.dict_GO[t]['Name'])
            defs.append(dGO.dict_GO[t]['Definition'])
        df2 = df_res.assign(name=names, definition=defs)

        df2[['p-value',
             'Adj p-value (FDR)']] = df2[['p-value',
                                          'Adj p-value (FDR)']].applymap(lambda x: '{:.2E}'.format(x))
        df2[['Proportion',
             'Global_proportion',
             'Expected_selected_nb']] = df2[['Proportion',
                                             'Global_proportion',
                                             'Expected_selected_nb']].applymap(lambda x: '{:.3f}'.format(x))
        df2csv = df2
        df2csv.to_csv(f"{output_dir}{o_file}.csv", sep='\t', index=False)

    f_out = open(f"{output_dir}{output_file}_padj{thresh_padj}.csv", 'w')
    f_out.write(f"GOtype\tCategory\tSelected_gene_nb\tTotal_gene_nb\tExpected_selected_nb\tp-value\tAdj p-value (FDR)\tName\tDefinition\n")

    for GOtype in ["GOc", "GOf", "GOp"]:
        empty_test = True
        with open(f"{output_dir}{output_file}_{GOtype}.csv", "r") as f:
            header = next(f)
            for line in f:
                line = line.strip('\n').split('\t')
                if float(line[6]) < thresh_padj:
                    f_out.write(f"{GOtype}\t{line[0]}\t{line[1]}\t{line[2]}\t{line[8]}\t{line[5]}\t{line[6]}\t{line[9]}\t{line[10]}\n")
                    empty_test = False
        f.close()
        if empty_test:
            f_out.write(f"{GOtype}\tNA\n")
    f_out.close()
    print(f"Results saved in {output_dir}{output_file}_significant.csv")


def gaf2annot(path2file, startline, regexp_locus):
    '''
    converts a .gaf file into an .annot file which can then be used to load
    GO information with the load_go method of the genome class.

    Args:
        path2file (str.): path to the file containing the data
        startline (int.): file start line
        regexp_locus (str.): Regular expression that matches only the
            locus in the gaf file
    Output:
        .annot file saved in the folder of the initial .gaf file.
    Example:
    >>> GO_analysis.gaf2annot("path2file/ecocyc.gaf",35,"b\\d\\d\\d\\d")
    '''
    if ".gaf" in path2file:
        path2file = path2file[:-4]

    f_out = open(path2file + ".annot", "w")
    with open(path2file + ".gaf", "r") as f_in:
        i = 0
        while i < startline:
            header = next(f_in)
            i += 1
        for line in f_in:
            line = line.strip('\n').split('\t')
            try:
                GOterm = line[4]
                locus_part = line[10].split('|')
                for l in locus_part:
                    if re.compile(regexp_locus).search(l):
                        f_out.write(f"{l}\t{GOterm}\n")
            except Exception as e:
                print(e)
    f_in.close()
    f_out.close()
