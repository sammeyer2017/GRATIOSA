#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
This module allows the classification and the statistical analysis 
(enrichment test, proportion test, student test) of omics and spatial data 
loaded on the different objects created, in particular with the Genome, 
Transcriptome, Chipseq and HiC classes.
'''
import numpy as np
from statsmodels.stats import proportion, multitest
from scipy import stats
from GRATIOSA.globvar import *
from pathlib import Path
from datetime import datetime
import pandas as pd


def data_classification(data_x, data_y, class_nb, *args, **kwargs):
    '''
    Classification of data into fractions according to defined thresholds
    (using the thresholds argument), according to class sizes (using the
    class_sizes argument) or into classes of equal size (if neither the
    thresholds argument nor the class_sizes argument is specified).

    Args:
        data_x (list): list of elements such as positions or gene names
        data_y (list): list of data associated with each element such as 
                signal coverage or gene expression
        class_nb (int.): number of classes (fractions) to create
        class_names (Optionnal [list]): list of names to give to each class
                By default, each class will be named by a number between 0 
                and class_nb.
        thresholds (Optionnal [list of ints.]): Thresholds used for the classification.
                If not given, the classification will distribute the data in 
                fractions of equal sizes.
        class_sizes (Optionnal [list of ints.]): Number of elements to put in each class.
                The total must be equal to the number of  elements in data_x 
                and in data_y.

    Returns: 
        tuple: 2 dictionaries:
            * dict. of shape {class_name:list of elements}
            * dict. of shape {class_name: list of data associated with each element}

    Example:
        >>> stat_analysis.data_classification(["a","b","d","c","e"],[1,2,4,3,5], 
        ...                                   class_nb=3,thresholds = [1,5])
        ({0: ['a'], 1: ['b', 'c', 'd'], 2: []}, {0: [1], 1: [2, 3, 4], 2: []})
        >>> stat_analysis.data_classification(["a","b","d","c","e"],[1,2,4,3,5], 
        ...                                   class_nb=3,class_sizes=[1,1,3])
        ({0: ['a'], 1: ['b'], 2: ['c', 'd', 'e']}, {0: [1], 1: [2], 2: [3, 4, 5]})
        >>> stat_analysis.data_classification(["a","b","d","c","e","f"],[1,2,4,3,5,6], 
        ...                                   class_nb=3)
        ({0: ['a', 'b'], 1: ['c', 'd'], 2: ['e', 'f']}, {0: [1, 2], 1: [3, 4], 2: [5, 6]})
    '''

    class_names = kwargs.get('class_names', np.arange(class_nb))
    thresholds = kwargs.get('thresholds')
    nb_per_class = kwargs.get('class_sizes')
    data_y = list(map(float, data_y))

    # STEP 1: Checks the adequacy between arguments
    if len(class_names) != class_nb:
        class_nb = len(class_names)
        print("WARNING: the number of classes names isn't equal to the class_nb")
        print(f"WARNING: the number of classes was set to {len(class_names)}")
    if thresholds:
        if class_nb == len(thresholds) + 1:
            # If only internal thresholds were given,
            # adds the min and max to the thresholds.
            thresholds = \
                [min(data_y + thresholds)] + thresholds + [max(data_y)]
        elif class_nb != len(thresholds) - 1:
            sys.exit("Please adjust the nb_class to thresholds.")
        print("Performing the classification in the following classes:")
        for i in np.arange(class_nb):
            print(f"{class_names[i]}: {thresholds[i]} to {thresholds[i+1]}")
    if len(data_x) != len(data_y):
        sys.exit("data_x and data_y do not have the same length")
    if nb_per_class:
        if np.sum(nb_per_class) != len(data_x):
            sys.exit("The sum of nb_per_class is not equal to the length of data_x.")
        elif len(class_names) != len(nb_per_class):
            sys.exit("Please adjust the class_names to the number of classes described with class_sizes.")

    # STEP 2: Performs the classification
    # Sorts the data
    data_y, data_x = map(list, zip(*sorted(zip(data_y, data_x))))
    classif_y = {}
    classif_x = {}

    # Option 1: classification using the chosen thresholds
    if thresholds:
        t, i = 0, 0
        while t < class_nb:
            classif_y[class_names[t]] = []
            classif_x[class_names[t]] = []
            while data_y[i] <= thresholds[t + 1] and i + 1 < len(data_y):
                classif_y[class_names[t]].append(data_y[i])
                classif_x[class_names[t]].append(data_x[i])
                i += 1
            t += 1

    # Option 2: classification using the chosen class sizes
    elif nb_per_class:
        prev_class_sums = 0
        for c in np.arange(class_nb):
            classif_y[class_names[c]
                      ] = data_y[prev_class_sums:nb_per_class[c] + prev_class_sums]
            classif_x[class_names[c]
                      ] = data_x[prev_class_sums:nb_per_class[c] + prev_class_sums]
            prev_class_sums += nb_per_class[c]

    # Option 3: classification in equal-size fractions
    else:
        size_fraction = len(data_y) / class_nb
        print("Performing the classification in the following classes:")
        for c in np.arange(class_nb - 1):
            classif_y[class_names[c]] = \
                data_y[int(c * size_fraction):int((c + 1) * size_fraction)]
            classif_x[class_names[c]] = \
                data_x[int(c * size_fraction):int((c + 1) * size_fraction)]
        classif_y[class_names[c + 1]] = data_y[int((c + 1) * size_fraction):]
        classif_x[class_names[c + 1]] = data_x[int((c + 1) * size_fraction):]
        for k in classif_y.keys():
            print(f"{k}: from {min(classif_y[k])} to {max(classif_y[k])}")

    for k in classif_y.keys():
        print(f"{k}: {len(classif_x[k])} elements")

    return classif_x, classif_y


def proportion_test(dict_cats,
                    dict_features,
                    targ_features,
                    all_features="all",
                    cats="all",
                    alt_hyp="one-sided",
                    output_dir=f"{resdir}proportion_test/",
                    output_file=f"prop_test{datetime.now()}"):
    '''
    Compares the proportion of targ_features/all_features between categories.
    The proportion test is based on a normal test using the
    stats.proportion.proportions_ztest function from statsmodels.

    Args:
        dict_cats (dict.): classification of each element in a dictionary of
               shape {category:[elements]}
        dict_features (dict.): feature corresponding to each element in a
                dictionary of shape {feature:[elements]}
        targ_features (list or string): targeted feature(s)
        all_features (Optionnal [list or str.]): features that will be used to 
                calculate the proportion, including targ_features
                (default: all keys of the dict_features dictionary)
        cats (Optionnal [list]): list of categories to compare
                (default: all keys of dict_cats)
        alt_hyp (Optionnal ["two-sided" or "one-sided"]): alternative hypothesis.
                If "one-sided" is chosen, both one-sided tests will be 
                performed with the statsmodels.stats.proportions_ztest 
                function and the smaller p-value will be kept. 
                See statsmodels documentation for more details.
                (default: "one-sided")
        output_dir (Optionnal [str.]): output directory
        output_file (Optionnal [str.]): output filename

    Returns:
        Dictionary: dict. of shape {"categories":cats, "proportions":prop,
        "confidence intervals":(ci0,ci1), "p-values":{(cat1,cat2):pval}

    Note:
        Test results are saved in the output_file

    Example:
        >>> dict_cats = {"borders":["A","B","G","K","L"],
        ...              "loops":["C","D","F","H","J"],
        ...              "None":["E","I"]}
        >>> data = {"act":["A","B","C","F","H"],
        ...         "rep":["D","G","L","I"],
        ...         "None":["J","K"],"NA":["E"]}
        >>> res = stat_analysis.proportion_test(dict_cats,data,"act",
        ...                                     all_features=["act","rep"],
        ...                                     alt_hyp="two-sided")
        >>> res["categories"]
        ['borders', 'loops', 'None']
        >>> res["proportions"]
        [0.5, 0.75, 0.0]
        >>> res["confidence intervals"]
        [(0.010009, 0.98999),(0.32565, 1.0),(0.0, 0.0)]
        >>> res["p-values"]
        {'borders-loops': 0.4652088184521418,
         'borders-None': 0.3613104285261787,
         'loops-None': 0.17090352023079747}
    '''
    # Formats the objects given as arguments as lists containing no duplicates.
    if isinstance(targ_features, str):
        targ_features = [targ_features]
    else:
        targ_features = list(dict.fromkeys(targ_features))

    if isinstance(all_features, str):
        if all_features == "all":
            all_features = set(dict_features.keys())
    else:
        all_features = list(dict.fromkeys(all_features))

    if cats == "all":
        cats = list(dict_cats.keys())
    else:
        cats = list(dict.fromkeys(cats))

    # For each category, creates a dictionnary containing:
    #   one list of elements associated to at least one feature of all_features
    #   one list of elements associated to at least one feature of tar_features
    targ_feat_elements, all_feat_elements = [], []
    for f in dict_features.keys():
        if f in targ_features:
            targ_feat_elements.extend(dict_features[f])
        if f in all_features:
            all_feat_elements.extend(dict_features[f])
    list_per_cat = {}
    for c in cats:
        list_per_cat[c] = {}
        list_per_cat[c]["tar"] = set(dict_cats[c]) & set(targ_feat_elements)
        list_per_cat[c]["all"] = set(dict_cats[c]) & set(all_feat_elements)

    # For each category, computes the proportion between the number of elements
    # in targ_features and the number in all_features and its confidence interval.
    # Saves the results in a .txt file.
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    results = open(f'{output_dir}{output_file}.txt', 'w')
    results.write(f"Targeted_features: {targ_features}\n")
    results.write(f"All_features: {all_features}\n\n")
    results.write(f"Category\tTargeted_feature\tAll_features\tProportion\t"
                  "Confidence interval\n")
    prop, cim = [], []
    for c in cats:
        try:
            nb_A = len(list_per_cat[c]["tar"])
            nb_tot = len(list_per_cat[c]["all"])
            print(f"{c}: {nb_tot} elements")
            pexp = float(nb_A / nb_tot)
            std = np.array(stats.binom.std(nb_tot, pexp, loc=0)) / nb_tot
            ci = proportion.proportion_confint(nb_A, nb_tot, alpha=0.05,
                                               method='normal')
            ci0 = max(ci[0], 0)
            ci1 = min(ci[1], 1)
            cim.append((ci0, ci1))
            prop.append(pexp)
            results.write(f"{c}\t{nb_A}\t{nb_tot}\t{pexp}\t{(ci0,ci1)}\n")
        except Exception as ex:
            print(f"Error with {c}: {ex}")
    results.write("\n")

    # Tests pairwise if the proportions between categories are significantly
    # differents. Performs one-sided or two-sided proportions test depending 
    # on the "alt_hyp" argument.
    comp = []
    pvals = {}
    for i in np.arange(len(cats)):
        for j in np.arange(i + 1, len(cats)):
            comp.append([cats[i], cats[j]])
    for cat1, cat2 in comp:
        s1 = len(list_per_cat[cat1]["tar"])
        s2 = len(list_per_cat[cat2]["tar"])
        n1 = s1 + len(list_per_cat[cat1]["all"])
        n2 = s2 + len(list_per_cat[cat2]["all"])
        if alt_hyp == "one-sided":
            stat1, pval1 = proportion.proportions_ztest([s1, s2], [n1, n2],
                                                        alternative='larger')
            stat2, pval2 = proportion.proportions_ztest([s1, s2], [n1, n2],
                                                        alternative='smaller')
            pval = min(pval1, pval2)
            results.write(
                f"One-sided test\tp-value between {cat1} and {cat2}\n")
            results.write(f"Larger\t{str(pval1)}\n")
            results.write(f"Smaller\t{str(pval2)}\n")
            print(f"One-sided test\tp-value between {cat1} and {cat2}\n")
            print(f"Larger\t{str(pval1)}\n")
            print(f"Smaller\t{str(pval2)}\n")

        elif alt_hyp == "two-sided":
            stat, pval = proportion.proportions_ztest([s1, s2], [n1, n2],
                                                      alternative='two-sided')
            results.write(f"Two-sided test p-value between {cat1} and {cat2}\t"
                          f"{str(pval)}\n")
            print(f"Two-sided test p-value between {cat1} and {cat2}\t"
                  f"{str(pval)}\n")
        else:
            sys.exit("The alternative hypothesis \"alt_hyp\" has to be"
                     "\"one-sided\" or \"two-sided\"")
        pvals[(cat1, cat2)] = pval

    results.close()
    print(f"Results saved in {output_dir}{output_file}")
    return {"categories": cats, "proportions": prop,
            "confidence intervals": cim, "p-values": pvals}


def enrichment_test(dict_cats,
                    dict_features,
                    targ_features,
                    all_features="all",
                    targ_cats="all",
                    all_cats="all",
                    min_nb_elements=4,
                    output_dir=f"{resdir}enrichment_test/",
                    output_file=f"enrich_test{datetime.now()}"):
    '''
    Computes enrichment tests (hypergeometric test) of a target in a
    category.

    Args:
        dict_cats (dict.): classification of each element in a dictionary of
                shape {category:[elements]}
                N.B.: the same element can be associated to multiple features
        dict_features (dict.): feature corresponding to each element in a
                dictionary of shape {feature:[elements]}
        targ_features (list or string): targeted feature(s)
        all_features (Optional[list or str.]): features that will be used to 
                calculate the proportion, including targ_features
                (default: all keys of the dict_features dictionary)
        targ_cats (Optional[list]): list of categories. The enrichment test is 
                performed for each category. (default: all keys of dict_cats)
        all_cats (Optional[list]): list of categories used to compute the global
                proportion and the expected number in the selection,
                including targ_cats.
                (default: all keys of dict_cats)
        min_nb_elements (Optional[int.]): Number of elements used as thresholds 
                for the feature selection. If there is strictly less than 
                min_nb_elements corresponding to a feature, the result 
                corresponding to this feature is not relevant and is 
                therefore neither returned nor reported in the output file.
                By default, min_nb_elements is set to 4.
        output_dir (Optional[str.]): output directory
        output_file (Optional[str.]): output filename

    Returns: 
        DataFrame: DataFrame containing the following columns:
            * 'Category'(str.): category
            * 'Selected_gene_nb'(int.): Nb of elements corresponding to this
              feature in the selection
            * 'Expected_selected_nb' (int.): Expected nb of elements 
              corresponding to this feature in the selection
            * 'Total_gene_nb'(int.): Nb of elements corresponding to this
              feature in the dict_features
            * 'Proportion'(float): ratio between Selected_gene_nb and 
              Total_gene_nb
            * 'Prop_conf_int' (np.array): 95% confidence interval with equal
              areas around the proportion.
            * 'p-value' (float): p-value obtained with the enrichment test
            * 'Adj p-value (FDR)' (float): p-value corrected for false
              discovery rate
            * 'Global_proportion' (float): ratio between nb of elements in 
              the selection and nb of elements in dict_features

    Note: 
        This function performs a p-value correction for false discovery rate using
        statsmodels.stats.multitest.fdrcorrection 

    Note: 
        The created DataFrame, ordered according to the adjusted pvalues,
        is reported in the output_file.

    Example:
        >>> dict_features = {"act": ["B", "D", "E", "H", "I", "M", "P", "Q", "R", "S", "T", "W"], 
        ...                  "rep": ["C", "F", "G", "U", "X"], "None": ["A"], "NA": ["J"]}                
        >>> dict_cats = {"GOterm1": ["A", "B", "D", "E", "F", "P", "Q", "R", "S", "T", "U"], 
        ...              "GOterm2": ["C", "E"], 
        ...              "GOterm3": ["A", "B", "F", "G", "H", "I", "M", "U", "V", "W", "X"], 
        ...              "GOterm4": ["C", "F", "G", "J"]}               
        >>> stat_analysis.enrichment_test(dict_cats,
        ...                           dict_features,
        ...                           targ_features=["act", "None"],
        ...                           all_features=["act", "None", "rep", "NA"],
        ...                           targ_cats=["GOterm1", "GOterm2", "GOterm3"],
        ...                           min_nb_elements=3,
        ...                           output_file="test")
          Category  Selected_gene_nb  Total_gene_nb  Proportion  Prop_conf_int  p-value  
        0  GOterm1                9          11    0.818182  [0.4545, 1.0]  0.27206
        1  GOterm3                6          10    0.600000     [0.2, 1.0]  0.97059
          Adj p-value (FDR)  Global_proportion  Expected_selected_nb
        0           0.54412        0.722222        7.944444
        1           0.97059        0.722222         7.22222
          Category  Selected_gene_nb  Total_gene_nb  Proportion     Prop_conf_int 
        0  GOterm1                9          11    0.818182  [0.6818, 0.9545]
        1  GOterm3                6          10    0.600000        [0.4, 0.8]
          p-value  Adj p-value (FDR)  Global_proportion  Expected_selected_nb
        0 0.16563            0.33126        0.684211        7.526316
        1 0.90867            0.90867        0.684211        6.84210
        # GOterm2 was ignored because its nb of elements is less than 3.
        # GOterm4 was ignored because it was not selected in the "features"
    '''
    # Formats the objects given as arguments as lists containing no duplicates.
    if isinstance(targ_features, str):
        targ_features = [targ_features]
    else:
        targ_features = list(dict.fromkeys(targ_features))

    if isinstance(all_features, str):
        if all_features == "all":
            all_features = dict_features.keys()
    else:
        all_features = list(dict.fromkeys(all_features))

    if isinstance(targ_cats, str):
        if targ_cats == "all":
            targ_cats = list(dict_cats.keys())
        else:
            targ_cats = [targ_cats]
    else:
        targ_cats = list(dict.fromkeys(targ_cats))

    if isinstance(all_cats, str):
        if all_cats == "all":
            all_cats = list(dict_cats.keys())
    else:
        all_cats = list(dict.fromkeys(all_cats))

    if set(targ_cats) - set(all_cats):
        sys.exit("Every category listed in \"targ_cats\""
                 "has to be included in \"all_cats\".\n")

    # Keeps only valid elements ie elements that correspond:
    #   to a category listed in all_cats
    #   and to a feature listed in all_features
    cat_elements = set([el for c in all_cats for el in dict_cats[c]])
    all_targ_elements = set(
        [el for t in all_features for el in dict_features[t]])
    valid_elements = cat_elements & all_targ_elements
    M = len(valid_elements)    # M: total number of valid elements

    # Selects valid elements that correspond to a feature listed in
    # targ_features
    targ_elements = set(
        [el for tA in targ_features for el in dict_features[tA]])
    selection = valid_elements & targ_elements
    n = len(selection)    # n: Nb of valid elements in the selection
    global_average = n / M

    # Perfoms enrichment test for each category
    file = []
    for cat in targ_cats:
        # N: Nb of valid elements in the category
        N = len(set(dict_cats[cat]) & valid_elements)
        if len(targ_cats) < 20:
            print(f"{cat}: {N} valid elements")

        # x: Nb of valid elements associated to the category in the selection
        x = len(set(dict_cats[cat]) & selection)
        if N != 0:
            prop = x / N

            # exp_nb: expected number of elements associated to the category
            # in the selection
            exp_nb = global_average * N

            # Hypergeometric test (See scipy.stats.hypergeom documentation)
            pval = stats.hypergeom.sf(x - 1, M, n, N)

            # ci0 and ci1 are the 95% confidence interval with equal areas
            # around the proportion.
            ci = stats.hypergeom.interval(0.95, M, n, N, loc=0)
            diff_ci = (ci[1] - ci[0]) / float(N)
            ci0, ci1 = max(prop - diff_ci / 2, 0), min(prop + diff_ci / 2, 1)

            # Keeps relevant results i.e. categorys with more than
            # "min_nb_elements" elements
            if N >= min_nb_elements:
                file.append([cat, pval, x, N, prop,
                             np.round((ci0, ci1), 4), global_average, exp_nb])

    df = pd.DataFrame(data=file, columns=['Category', 'p-value', 'Selected_gene_nb',
                                          'Total_gene_nb', 'Proportion',
                                          'Prop_conf_int', 'Global_proportion',
                                          'Expected_selected_nb'])

    # pvalue correction for false discovery rate
    # (See statsmodels.stats.multitest.fdrcorrection documentation)
    df['Adj p-value (FDR)'] = multitest.fdrcorrection(df['p-value'])[1]

    # Returns and saves the results as a dataframe
    df = df[['Category', 'Selected_gene_nb', 'Total_gene_nb', 'Proportion', 'Prop_conf_int',
             'p-value', 'Adj p-value (FDR)', 'Global_proportion', 'Expected_selected_nb']]
    # Sorts the results according to the adjusted pvalues
    df = df.sort_values(by=['Adj p-value (FDR)'], ascending=True)

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    df2csv = df.round({'Proportion': 3,
                       'Global_proportion': 3,
                       'p-value': 3,
                       'Adj p-value (FDR)': 3,
                       'Expected_selected_nb': 3})
    df2csv.to_csv(f"{output_dir}{output_file}.csv", sep='\t', index=False)
    print(f"Results saved in {output_dir}{output_file}.csv")

    return df


def quantitative_data_student_test(dict_data, cats="all",
                                   alt_hyp="one-sided",
                                   output_dir=f"{resdir}student_test/",
                                   output_file=f"student_test{datetime.now()}"):
    '''
    Computes the T-test for the means of independants categories using the
    scipy.stats.ttest_ind.

    Args:
        dict_data (dict.): datapoints corresponding to each category in a
                dictionary of shape {category:list of datapoints}
        cats (Optional [list]): list of categories to compare
                (default: all keys of dict_data)
        alt_hyp (Optional ["two-sided" or "one-sided"]): alternative hypothesis.
                If "one-sided" is chosen, both one-sided tests will be 
                performed with the scipy.stats.ttest_ind function and the 
                smaller p-value will be kept. See scipy documentation for 
                more details. (default: "one-sided")
        output_dir (Optional[str.]): output directory
        output_file (Optional[str.]): output filename

    Returns:
        Dictionary: dict. of shape 
        {"categories":cats, "means":means,
        "confidence intervals":(ci0,ci1), "p-values":{(cat1,cat2):pval}
    Note:
        Test results are also reported in the output_file.

    Example:
        >>> dict_data = {'a':[1,2,5,6,19], 'b':[10,24,4,15]}
        >>> stat_analysis.quantitative_data_student_test(dict_data)
        {'categories': ['a', 'b'],
         'means': [6.6, 13.25],
         'confidence intervals': [(0.2610995224360444, 12.938900477563955),
                                  (4.958672795504613, 21.54132720449539)],
         'p-values': {('a', 'b'): 0.12169488735182109}}
    '''
    if cats == "all":
        cats = list(dict_data.keys())
    elif isinstance(cats, str):
        cats = [cats]

    # Saves the results in a .txt file.
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    results = open(f'{output_dir}{output_file}.txt', 'w')
    results.write(f"Category\tMean\tConfidence interval\n")
    means, cim = [], []
    for c in cats:
        m = np.mean(dict_data[c])
        means.append(m)
        ci = stats.norm.interval(alpha=0.95, loc=m,
                                 scale=stats.sem(dict_data[c]))
        cim.append((ci[0], ci[1]))
        results.write(f"{c}\t{m}\t{(ci[0],ci[1])}\n")
    results.write("\n")

    # Tests pairwise if the means between categories are significantly
    # differents. Performs one-sided or two-sided proportions test depending 
    # on the "alt_hyp" argument.
    comp = []
    pvals = {}
    for i in np.arange(len(cats)):
        for j in np.arange(i + 1, len(cats)):
            comp.append([cats[i], cats[j]])
    for cat1, cat2 in comp:
        stats.ttest_ind(dict_data[cat1], dict_data[cat2])
        if alt_hyp == "one-sided":
            stat1, pval1 = stats.ttest_ind(dict_data[cat1], dict_data[cat2],
                                           alternative='greater')
            stat2, pval2 = stats.ttest_ind(dict_data[cat1], dict_data[cat2],
                                           alternative='less')
            pval = min(pval1, pval2)
            results.write(
                f"One-sided test\tp-value between {cat1} and {cat2}\n")
            results.write(f"Larger\t{str(pval1)}\n")
            results.write(f"Smaller\t{str(pval2)}\n")
        elif alt_hyp == "two-sided":
            stat, pval = stats.ttest_ind(dict_data[cat1], dict_data[cat2],
                                         alternative='two-sided')
            results.write(f"Two-sided test p-value between {cat1} and {cat2}\t"
                          f"{str(pval)}\n")
        else:
            sys.exit("The alternative hypothesis \"alt_hyp\" has to be"
                     "\"one-sided\" or \"two-sided\"")
        pvals[(cat1, cat2)] = pval
    results.close()
    print(f"Results saved in {output_dir}{output_file}")
    return {"categories": cats, "means": means,
            "confidence intervals": cim, "p-values": pvals}