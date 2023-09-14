#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module allows to graphically represent (barplots) the statistical
analysis performed with the module stat_analysis.
"""

import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
from GRATIOSA.globvar import *
from GRATIOSA import stat_analysis

plt.rcParams.update({'font.size': 9,
                     'font.family': "Arial",
                     'xtick.labelsize': 8,
                     'ytick.labelsize': 8,
                     'axes.labelsize': 9,
                     'axes.titlesize': 10})


def significance(pval):
    """
    Converts p-values in stars annotation.
    """
    if pval <= 0.001:
        s = '***'
    elif pval <= 0.01:
        s = '**'
    elif pval <= 0.05:
        s = '*'
    else:
        s = 'ns'
    return s


def barplot_annotate_brackets(categories, y_up, dict_pval, *args, **kwargs):
    """
    Annotates barplot with p-values, using the significance function to 
    convert p-values in stars annotation.

    Args:
        categories (list): list of the plotted categories, in the order of
                their position on the plot
        y_up (list): list of the maximal y position of each bar (taking into
                account the confidence intervals)
        dict_pval (dict.): dictionnary of shape {"cat1-cat2":pvalue} with 
                cat1 and cat2 contained in the categories list given as
                argument
        linewidth (Optional [float.]): Linewidth of the brackets. (default: 1.5)
    """
    # gets the position of each category on the barplot
    pos = {}
    for p, cat in enumerate(categories):
        pos[cat] = p

    # creates a dictionnary containing, for each significant pvalue:
    # - its pvalue converted as "*" to represent it on the plot
    # - start and end position of the bracket (positions of the 2 categories)
    # - the minimal height of the bracket
    dict_br = {}
    for cats in list(dict_pval.keys()):
        if dict_pval[cats] <= 0.05:
            s = significance(dict_pval[cats])
            c0, c1 = cats
            x0, x1 = np.sort([pos[c0], pos[c1]])
            hmax = 0
            for x in np.arange(x0, x1 + 1):
                hmax = max(hmax, y_up[x])
            dict_br[cats] = {"x0": x0, "x1": x1, "hmax": hmax, "s": s}

    # gets the ylim to defines the necessary difference between the heights
    # of 2 brackets
    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh = (ax_y1 - ax_y0) / 10

    # plots the brackets without overlay
    h = 0
    linewidth = kwargs.get("linewidth", 1.5)
    while dict_br:
        h += dh
        pre_x1 = -1
        plotted_b = []
        for b in dict_br.keys():
            b_x0 = dict_br[b]["x0"]
            b_x1 = dict_br[b]["x1"]
            if b_x0 > pre_x1 and dict_br[b]["hmax"] < h:
                barx = [b_x0, b_x0, b_x1, b_x1]
                bary = [h, h + dh / 3, h + dh / 3, h]
                plt.plot(barx, bary, c='black', linewidth=linewidth)
                plt.text(float(b_x0 + b_x1) / 2, h + dh / 3, dict_br[b]["s"],
                         ha='center')
                pre_x1 = b_x1
                plotted_b.append(b)
        for b in plotted_b:
            del dict_br[b]


def plot_proportion_test(dict_cats,
                         dict_features,
                         targ_features,
                         all_features="all",
                         cats="all",
                         alt_hyp="one-sided",
                         output_dir=f"{resdir}proportion_test/",
                         output_file=f"prop_test{datetime.now()}",
                         file_extension=".pdf",
                         xlabel="",
                         ylabel="Proportion",
                         title="",
                         annot_brackets=True,
                         *args, **kwargs):
    '''
    Barplots of the proportion test: targ_features/all_features between
    categories. The proportion test, based on normal test, is computed with
    stat_analysis.proportion_test (see its documentation for more details)

    Args:
        dict_cats (dict.): classification of each elements in a dictionary of
                shape {category:[elements]}
                Example: {"border":["GeneA","GeneB"],"None":["GeneC","GeneD"]}
        dict_features (dict.): feature corresponding to each element in a
                dictionary of shape {feature:[elements]}
                Example: {"act": ["GeneA","GeneC","GeneD"], "rep":["GeneB"]}
        targ_features (list or string): targeted feature(s)
        all_features (Optional [list or string]): features that will be used to calculate
                the proportion, including targ_features default: all keys of 
                the dict_features dictionary)
        cats (Optional [list]): list of categories to compare (default: all keys of 
                dict_cats)
        alt_hyp (Optional ["two-sided" or "one-sided"]): alternative hypothesis.
                If "one-sided" is chosen, both one sided tests will be 
                perform with the statsmodels.stats.proportions_ztest 
                function and the smaller pvalue will be kept. 
                See statsmodels documentation for more details. 
                (default: "one-sided")
        output_dir (Optional [str.]): output directory
        output_file (Optional [str.]): output filename for the proportion test data 
                (.txt) and the plot
        file_extension (Optional [str.]): Graphic file extension type (.pdf by default)
        xlabel (Optional [str.]): label for the x-axis (default: empty)
        ylabel (Optional [str.]): label for the y-axis (default: "Proportion")
        title (Optional [str.]): general title for the figure
        annot_brackets (Optional [bool.]): if True, the barplot will be annotated according
                to the p-values using stars annotaion (default: True)
        brackets_linewidth (Optional [float.]): Linewidth of the brackets. (default: 1.5)
        ymin (Optional [float]): y-axis bottom limit
        ymax (Optional [float]): y-axis top limit
        figsize (Optional [(float,float)]): width and height in inches (by default: (w,2.2)
                with w dependent on the number of categories)
        xticks_rotation (Optional [int.]): x-ticks labels rotation in degrees
        xticks_labels (Optional [list.]): x-ticks labels (by default: cats)
        err_capsize (Optional [float.]): Length of the error bar caps in points
        bar_linewidth (Optional [float.]): Width of the bars edge. (default: 1.5)
        bar_width (Optional [float.]): Width of the bars. (default dependent on the number 
                of categories. If less than 5 cats: 0.7)


    Example:
        >>> import numpy as np
        >>> dict_cat = {"cat1":np.arange(100,154),"cat2":np.arange(1,100),
        ...             "cat3":np.arange(154,180),"cat4":np.arange(180,230)}
        >>> dict_features ={"act":list(np.arange(1,90))+list(np.arange(100,120))
        ...                +list(np.arange(154,158))+list(np.arange(180,220)),
        ...        "rep":list(np.arange(90,94))+list(np.arange(120,150))
        ...                +list(np.arange(158,177))+list(np.arange(220,225)),
        ...        "None":list(np.arange(94,100))+list(np.arange(150,154))
        ...                +list(np.arange(177,180))+list(np.arange(225,230))}
        >>> plot_stat_analysis.plot_proportion_test(dict_cat,dict_features,"act",
        ...                all_features=["act","rep"],alt_hyp="two-sided",output_file="test")
    '''
    res = stat_analysis.proportion_test(
        dict_cats=dict_cats, dict_features=dict_features, targ_features=targ_features, all_features=all_features,
        cats=cats, alt_hyp=alt_hyp, output_dir=output_dir, output_file=output_file)

    cats = res["categories"]
    prop = res["proportions"]

    ci = np.zeros((2, len(cats)))
    for c in np.arange(len(cats)):
        ci[0, c] = prop[c] - res["confidence intervals"][c][0]
        ci[1, c] = res["confidence intervals"][c][1] - prop[c]

    if len(cats) < 5:
        wb = 0.7
    elif len(cats) < 15:
        wb = 0.5
    else:
        wb = 0.4
    figsize = kwargs.get("figsize", (len(cats) * wb/1.5+0.2, 2.2))
    plt.figure(figsize=figsize)

    bar_width = kwargs.get("bar_width", wb)
    bar_linewidth = kwargs.get("bar_linewidth", 1.5)
    err_capsize = kwargs.get("err_capsize", bar_width*5)
    plt.bar(cats, prop, yerr=ci,
            color='white', edgecolor='black', ecolor='black',
            width=bar_width, linewidth=bar_linewidth, capsize=err_capsize)

    error_max = np.amax(res['confidence intervals'], axis=1)

    if annot_brackets:
        brackets_linewidth = kwargs.get("brackets_linewidth", 1.5)
        barplot_annotate_brackets(
            cats,
            error_max,
            res["p-values"],
            linewidth=brackets_linewidth)

    plt.ylabel(ylabel)
    if xlabel != "":
        plt.xlabel(xlabel)
    ymin = kwargs.get('ymin', plt.ylim()[0])
    ymax = kwargs.get('ymax', plt.ylim()[1])
    xticks_rotation = kwargs.get("xticks_rotation", 0)
    xticks_labels = kwargs.get("xticks_labels", cats)
    plt.xticks(ticks=np.arange(len(cats)),
               labels=xticks_labels,
               rotation=xticks_rotation)
    plt.ylim(ymin, ymax)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(
        f"{output_dir}{output_file}{file_extension}",
        transparent=False,
        facecolor='white')
    plt.show()
    plt.close()


def plot_enrichment_test(dict_cats,
                         dict_features,
                         targ_features,
                         all_features="all",
                         targ_cats="all",
                         all_cats="all",
                         min_nb_elements=4,
                         output_dir=f"{resdir}enrichment_test/",
                         output_file=f"enrich_test{datetime.now()}",
                         file_extension=".pdf",
                         xlabel="",
                         ylabel="Proportion",
                         title="",
                         legend_text="Global\nproportion",
                         legend_loc="best",
                         annot_star=True,
                         *args, **kwargs):
    '''
    Barplots of enrichment tests (hypergeometric test) of features in a 
    sublist. The test is performed with stat_analysis.enrichment_test (see 
    its documentation for more details)

    Args:
        dict_cats (dict.): classification of each elements in a dictionary 
                of shape {category:[elements]}
                N.B.: the same element can be associated to multiple features
                fcExample: {"GOterm1":["GeneA","GeneB"],"GOterm2:["GeneA","GeneC"]}
        dict_features (dict.): feature corresponding to each element in a
                dictionary of shape {feature:[elements]}
                Example: {"act": ["GeneA","GeneC","GeneD"], "rep":["GeneB"]}
        targ_features (list or string): targeted feature(s)
        all_features (Optional [list or str.]): features that will be used to calculate
                the proportion, including targ_features
                (Default: all keys of the dict_features dictionary)
        targ_cats (Optional [list]): list of categories. The enrichment test is 
                performed for each catergory. (default: all keys of dict_cats)
        all_cats (Optional [list]): list of categories used to compute the global
                proportion and the expected number in the selection. 
                All_cats includes targ_cats. (default: all keys of dict_cats)
        min_nb_elements (Optional [int.]): Number of elements used as thresholds for the 
                feature selection. If there is stricly less than 
                min_nb_elements corresponding to a feature, the result 
                corresponding to this feature is not relevant and is 
                therefore neither returned nor reported in the output file.
                By default, min_nb_elements is set to 4.
        output_dir (Optional [str.]): output directory
        output_file (Optional [str.]): output filename for the proportion test data 
                (.txt) and the plot
        file_extension (Optional [str.]): Graphic file extension type (.pdf by default)
        xlabel (Optional [str.]): label for the x-axis (default: empty)
        ylabel (Optional [str.]): label for the y-axis (default: "Proportion")
        title (Optional [str.]): general title for the figure
        annot_star (Optional [bool.]): if True, the barplot will be annotated according
                to the p-values using stars annotaion (default: True)
        ymin (Optional [float]): y-axis bottom limit
        ymax (Optional [float]): y-axis top limit
        legend_text (Optional [str.]): Legend text (default: "Global proportion")
                If set to None, no legend will be plotted.
        legend_loc (Optional [str.]): Location of the legend such as 'upper right',
                'lower right', 'lower left', 'lower left' and 'best' 
                (default: 'best'). See matplotlib.pyplot.legend for more 
                options
        figsize (Optional [(float,float)]): width and height in inches (by default: (w,2.2)
                with w dependent on the number of categories)
        xticks_rotation (Optional [int.]): x-ticks labels rotation in degrees
        xticks_labels (Optional [list.]): x-ticks labels (by default: targ_cats)
        err_capsize (Optional [float]): Length of the error bar caps in points
        bar_linewidth (Optional [float]): Width of the bars edge. (default: 1.5)
        bar_width (Optional [float]): Width of the bars. (default dependent on the 
                number of categories. If less than 5 cats: 0.7)

    Example:
        >>>  dataX = {"GOterm1":["A","B","D","E","F"],
        ...           "GOterm2":["C","E"],
        ...           "GOterm3":["A","B","F","G","H","I","M"],
        ...           "GOterm4":["C","F","G","J"]}
        >>>  = ["A","E","I","F","G","H","J"]
        >>> plot_stat_analysis.plot_enrichment_test(
        ...                         dict_cats,dict_features,
        ...                         targ_feature=["act","None"],
        ...                         all_features=["act","None","rep","NA"],
        ...                         targ_cats=["GOterm1","GOterm2","GOterm3"],
        ...                         min_nb_elements=3,output_file="test1")
    '''

    df_res = stat_analysis.enrichment_test(dict_cats=dict_cats,
                                           dict_features=dict_features,
                                           targ_features=targ_features,
                                           all_features=all_features,
                                           targ_cats=targ_cats,
                                           all_cats=all_cats,
                                           min_nb_elements=min_nb_elements,
                                           output_dir=output_dir,
                                           output_file=output_file)
    print(df_res[["Category",
                  "Selected_gene_nb",
                  "Expected_selected_nb",
                  "Adj p-value (FDR)"]])
    print(df_res)
    res = df_res.set_index("Category").to_dict('index')

    prop = []
    ci = np.zeros((2, len(targ_cats)))
    for n, t in enumerate(targ_cats):
        prop.append(res[t]['Proportion'])
        ci[0, n] = res[t]['Proportion'] - res[t]["Prop_conf_int"][0]
        ci[1, n] = res[t]["Prop_conf_int"][1] - res[t]['Proportion']

    if len(targ_cats) < 5:
        wb = 0.7
    elif len(targ_cats) < 15:
        wb = 0.5
    else:
        wb = 0.4
    figsize = kwargs.get("figsize", (len(targ_cats) * wb/1.5+0.2, 2.2))
    plt.figure(figsize=figsize)

    bar_width = kwargs.get("bar_width", wb)
    bar_linewidth = kwargs.get("bar_linewidth", 1.5)
    err_capsize = kwargs.get("err_capsize", bar_width*5)
    plt.bar(targ_cats, prop, yerr=ci,
            color='white', edgecolor='black', ecolor='black',
            width=bar_width, linewidth=bar_linewidth, capsize=err_capsize)
    
    # gets the ylim to defines the necessary difference between the
    # errorbar and the annotation in order to be readable
    ax_y0, ax_y1 = plt.gca().get_ylim()
    dh = (ax_y1 - ax_y0) / 40

    max_ci = 0
    for n, t in enumerate(targ_cats):
        max_ci = max(max_ci, ci[1, n] + prop[n] + 4 * dh)
        if res[t]["Adj p-value (FDR)"] < 0.05 and annot_star:
            s = significance(res[t]["Adj p-value (FDR)"])
            plt.text(n, ci[1, n] + prop[n] + dh /
                     2, s, ha='center', fontsize=14)

    ymin = kwargs.get('ymin', plt.ylim()[0])
    ymax = kwargs.get('ymax', max_ci*1.1)
    plt.ylim(ymin, ymax)
    plt.ylabel(ylabel)
    if xlabel != "":
        plt.xlabel(xlabel)
    xticks_rotation = kwargs.get("xticks_rotation", 0)
    xticks_labels = kwargs.get("xticks_labels", targ_cats)
    plt.xticks(ticks=np.arange(len(targ_cats)),
               labels=xticks_labels,
               rotation=xticks_rotation)

    # plots an horizontal line corresponing to the Global_proportion
    plt.axhline(df_res["Global_proportion"][0], c="blue", ls='--',
                label=legend_text)
    if legend_text not in [None,"None"] :
        plt.legend(loc=legend_loc)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(
        f"{output_dir}{output_file}{file_extension}",
        transparent=False,
        facecolor='white')
    plt.show()
    plt.close()


def plot_student_test(dict_data, cats="all",
                      alt_hyp="one-sided",
                      output_dir=f"{resdir}student_test/",
                      output_file=f"student_test{datetime.now()}",
                      file_extension=".pdf",
                      xlabel="",
                      ylabel="Mean(data)",
                      title="",
                      annot_brackets=True,
                      *args, **kwargs):
    '''
    Barplots of the student test computed with stat_analysis.quantitative_
    data_student_test (see its documentation for more details)

    Args:
        dict_data (dict.): feature corresponding to each element in a
                dictionary of shape {category:list of datapoints}
        cats (Optional [list]): list of categories to compare
                (default: all keys of dict_data)
        alt_hyp (Optional ["two-sided" or "one-sided"]): alternative hypothesis.
                If "one-sided" is chosen, both one-sided tests will be 
                performed with the scipy.stats.ttest_ind function and the 
                smaller p-value will be kept. See scipy documentation for 
                more details. (default: "one-sided")
        output_dir (Optional [str.]): output directory
        output_file (Optional [str.]): output filename for the student test data (.txt)
                and the plot
        file_extension (Optional [str.]): Graphic file extension type (.pdf by default)
        xlabel (Optional [str.]): label for the x-axis (default: empty)
        ylabel (Optional [str.]): label for the y-axis (default: "Mean(data)")
        title (Optional [str.]): general title for the figure
        annot_brackets (Optional [bool.]): if True, the barplot will be annotated 
                according to the p-values using stars annotaion
                (default: True)
        brackets_linewidth (Optional [float.]): Linewidth of the brackets.(default: 1.5)
        ymin (Optional [float]): y-axis bottom limit
        ymax (Optional [float]): y-axis top limit
        figsize (Optional [(float,float)]): width and height in inches (by default: (w,2.2)
                with w dependent on the number of categories)
        xticks_rotation (Optional [int.]): x-ticks labels rotation in degrees
        xticks_labels (Optional [list]): x-ticks labels (by default: cats)
        err_capsize (Optional [float]): Length of the error bar caps in points
        bar_linewidth (Optional [float]): Width of the bars edge. (default: 2.0)
        bar_width (Optional [float]): Width of the bars. (default: 0.6) (default 
                dependent on the number of categories. If less than 5 cats: 0.7)

    Example:
        >>> dict_data = {'a':[1,2,5,6,19], 'b':[10,24,4,15]}
        >>> plot_stat_analysis.plot_student_test(dict_data)
    '''
    res = stat_analysis.quantitative_data_student_test(dict_data, cats=cats,
                                                       alt_hyp=alt_hyp,
                                                       output_dir=output_dir,
                                                       output_file=output_file)

    cats = res["categories"]
    means = res["means"]

    ci = np.zeros((2, len(cats)))
    for c in np.arange(len(cats)):
        ci[0, c] = means[c] - res["confidence intervals"][c][0]
        ci[1, c] = res["confidence intervals"][c][1] - means[c]

    if len(cats) < 5:
        wb = 0.7
    elif len(cats) < 15:
        wb = 0.5
    else:
        wb = 0.4
    figsize = kwargs.get("figsize", (len(cats) * wb/1.5+0.2, 2.2))
    plt.figure(figsize=figsize)

    bar_width = kwargs.get("bar_width", wb)
    bar_linewidth = kwargs.get("bar_linewidth", 1.5)
    err_capsize = kwargs.get("err_capsize", bar_width*5)
    plt.bar([str(c) for c in cats], means, yerr=ci,
            color='white', edgecolor='black', ecolor='black',
            width=bar_width, linewidth=bar_linewidth, capsize=err_capsize)
    error_max = np.amax(res['confidence intervals'], axis=1)

    if annot_brackets:
        brackets_linewidth = kwargs.get("brackets_linewidth", 1.5)
        barplot_annotate_brackets(
            cats,
            error_max,
            res["p-values"],
            linewidth=brackets_linewidth)

    plt.ylabel(ylabel)
    if xlabel != "":
        plt.xlabel(xlabel)
    ymin = kwargs.get('ymin', plt.ylim()[0])
    ymax = kwargs.get('ymax', plt.ylim()[1])
    xticks_rotation = kwargs.get("xticks_rotation", 0)
    xticks_labels = kwargs.get("xticks_labels", cats)
    plt.xticks(ticks=np.arange(len(cats)),
               labels=xticks_labels,
               rotation=xticks_rotation)
    plt.xticks(rotation=xticks_rotation)
    plt.ylim(ymin, ymax)
    plt.axhline(0, color="black")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(
        f"{output_dir}{output_file}{file_extension}",
        transparent=False,
        facecolor='white')
    plt.show()
    plt.close()