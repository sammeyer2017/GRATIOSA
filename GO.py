#! /usr/bin/env python
# -*- coding: utf-8 -*-
from globvar import *
import wget
from datetime import datetime
import stat_analysis
#==============================================================================#

class GO:

    def __init__(self,
                filename="go-basic.obo",
                path2dir =f"{basedir}topo_database/",
                obo_reload=False):
        if obo_reload :
            wget.download("http://purl.obolibrary.org/obo/go/go-basic.obo",
                          out=f"{basedir}topo_database/")

        self.dict_GO = {}
        with open(path2dir+filename, "r") as f:
            for line in f:
                line = line.strip("\n").replace("\"","").split(": ")
                if "id" in line:
                    if "GO" in line[1]:
                        ID = line[1]
                    else :
                        ID ="other"
                    self.dict_GO[ID]={}
                elif "name" in line:
                    self.dict_GO[ID]["name"] = line[1]
                elif "namespace" in line:
                    self.dict_GO[ID]["namespace"] = line[1]
                elif "def" in line:
                    self.dict_GO[ID]["definition"] = line[1]
        self.dict_GO.pop('other')

        self.dict_GOc, self.dict_GOp, self.dict_GOf = {}, {}, {}
        for GO in self.dict_GO.keys() :
            if self.dict_GO[GO]["namespace"] == "biological_process":
                self.dict_GOp[GO]={"name":self.dict_GO[GO]["name"],
                              "definition":self.dict_GO[GO]["definition"]}
            elif self.dict_GO[GO]["namespace"] == "molecular_function":
                self.dict_GOf[GO]={"name":self.dict_GO[GO]["name"],
                              "definition":self.dict_GO[GO]["definition"]}
            elif self.dict_GO[GO]["namespace"] == "cellular_component":
                self.dict_GOc[GO]={"name":self.dict_GO[GO]["name"],
                              "definition":self.dict_GO[GO]["definition"]}

def GO_enrichment_test(
                    dict_GOterms,
                    dict_features,
                    targ_features,
                    path2dir =f"{basedir}topo_database/",
                    GO_filename ="go-basic.obo",
                    all_features="all",
                    targ_cats="all",
                    all_cats="all",
                    min_nb_elements=4,
                    output_dir=f"{resdir}enrichment_test/",
                    output_file=f"enrich_test{datetime.now()}"):
    ''' Computes enrichment tests (hypergeometric test) of a target in a 
    category.
    N.B.: performs a p-value correction for false discovery rate (See
          statsmodels.stats.multitest.fdrcorrection for more details)

    Required args:
        dict_GOterms (dict.): classification of each element in a dictionary of
                           shape {category:[elements]}
                           N.B.: the same element can be associated to multiple
                           features
                           Example: {"GOterm1":["GeneA","GeneB"],
                                      "GOterm2:["GeneA","GeneC"]}
        dict_features (dict.): feature corresponding to each element in a
                               dictionary of shape {feature:[elements]}
                               Example: {"act": ["GeneA","GeneC","GeneD"],
                                          "rep":["GeneB"]}
        targ_features (list or string): targeted feature(s)

    Optional args:
        all_features (list or string): features that will be used to calculate
                                  the proportion, including targ_features
                                  (default: all keys of the dict_features
                                  dictionary)
        targ_cats (list): list of categories. The enrichment test is performed
                      for each category. (default: all keys of dict_cats)
        all_cats (list): list of categories used to compute the global
                             proportion and the expected number in the selection, 
                             including targ_cats.
                             (default: all keys of dict_cats)
        min_nb_elements: Number of elements used as thresholds for the feature
                         selection. If there is strictly less than min_nb_elements
                         corresponding to a feature, the result corresponding to
                         this feature is not relevant and is therefore neither
                         returned nor reported in the output file.
                         By default, min_nb_elements is set to 4.
        output_dir (str.): output directory
        output_file (str.): output filename

    Returns:
        dataframe containing the following columns:
            0. 'Term'(str.): GO term
            1. 'Selected_gene_nb'(int.): Nb of elements corresponding to this
                                      feature in the selection
            2. 'Expected_selected_number' (int.): Expected nb of elements corresponding
                                        to this feature in the selection
            3. 'Total_gene_nb'(int.): Nb of elements corresponding to this
                                   feature in the dict_features
            4. 'Proportion'(float): ratio between Selected_gene_nb and Total_gene_nb
            5. 'Prop_conf_int' (np.array): 95% confidence interval with equal
                                            areas around the proportion.
            6. 'p-value' (float): p-value obtained with the enrichment test
            7. 'Adj p-value (FDR)' (float): p-value corrected for false
                                            discovery rate
            8. 'Genomic average' (float): ratio between nb of elements in the
                                          selection and nb of elements in
                                          dict_features
            N.B.:
                The DataFrame, ordered according to the adjusted pvalues,is
                reported in the output_file.

    Example:
        >>>  dict_features = {
                        "act":["B","D","E","H","I","M","P","Q","R","S","T","W"],
                        "rep":["C","F","G","U","X"],
                        "None":["A"],
                        "NA":["J"]}
        >>> dict_cats = {"GOterm1":["A","B","D","E","F","P","Q","R","S","T","U"],
                         "GOterm2":["C","E"],
                         "GOterm3":["A","B","F","G","H","I","M","U","V","W","X"],
                         "GOterm4":["C","F","G","J"]}
        >>> stat_analysis.enrichment_test(dict_cats,dict_features,
                                          targ_features=["act","None"],
                                          all_features=["act","None","rep","NA"],
                                          targ_cats=["GOterm1","GOterm2","GOterm3"],
                                          min_nb_elements=3,output_file="test")
          Term  Selected_gene_nb  Total_gene_nb  Proportion  Prop_conf_int  p-value  \
        0  GOterm1                9          11    0.818182  [0.4545, 1.0]  0.27206
        1  GOterm3                6          10    0.600000     [0.2, 1.0]  0.97059

          Adj p-value (FDR)  Global_proportion  Expected_selected_number
        0           0.54412        0.722222        7.944444
        1           0.97059        0.722222         7.22222
          Term Selected_gene_nb  Total_gene_nb  Proportion     Prop_conf_int  \
        0  GOterm1                9          11    0.818182  [0.6818, 0.9545]
        1  GOterm3                6          10    0.600000        [0.4, 0.8]
          p-value  Adj p-value (FDR)  Global_proportion  Expected_selected_number
        0 0.16563            0.33126        0.684211        7.526316
        1 0.90867            0.90867        0.684211        6.84210
        # GOterm2 was ignored because its nb of elements is less than 3.
        # GOterm4 was ignored because it was not selected in the "features"
    '''
    # Formats the objects given as arguments as lists containing no duplicates.
    

    dGO = GO(filename=GO_filename,path2dir=path2dir)

    list_GO = ["GOc","GOp","GOf"]

    dGOp, dGOc, dGOf = {},{},{}

    unknown_terms = [] 
    for t in dict_GOterms.keys():
        try :
            sp = dGO.dict_GO[t]['namespace']
            if sp == "biological_process" :
                dGOp[t] = dict_GOterms[t]
            elif sp == "molecular_function" :
                dGOf[t] = dict_GOterms[t]
            elif sp == "cellular_component":
                dGOc[t] = dict_GOterms[t]

        except :
            unknown_terms.append(t)
    print(f"{len(unknown_terms)} GO terms are not in the GO annotation")


    for n,dict_cats in enumerate([dGOc,dGOp,dGOf]):
        o_file = f"{output_file}_{list_GO[n]}"
        df_res = stat_analysis.enrichment_test(dict_cats = dict_cats, 
                                               dict_features = dict_features, 
                                               targ_features = targ_features, 
                                               all_features = all_features, 
                                               targ_cats = targ_cats,
                                               all_cats = all_cats,
                                               min_nb_elements = min_nb_elements,
                                               output_dir = output_dir,
                                               output_file = o_file)
        terms = df_res["Category"]
        names = []
        defs = []
        for t in terms :
            names.append(dGO.dict_GO[t]['name'])
            defs.append(dGO.dict_GO[t]['definition'])
        df2 = df_res.assign(name=names,definition=defs)
        df2csv = df2.round({'Proportion': 4, 'Global_proportion': 4})
        df2csv.to_csv(f"{output_dir}{output_file}{list_GO[n]}.csv", sep='\t', index=False)