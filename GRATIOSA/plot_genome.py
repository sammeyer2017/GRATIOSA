import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from datetime import datetime
from pathlib import Path
from GRATIOSA import Genome
from GRATIOSA.globvar import *

'''
This module plots experimental information (one subplot per experiment) and 
gene annotations for a region. The genes annotation subplot is created with
the subplot_genes function. RNASeq subplots and Chipseq subplots are created 
with the subplot_rnaseq_coverage and subplot_signal functions respectively.
'''


plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'legend.fontsize': 10})
plt.rcParams.update({'font.family': "Arial"})


def plot_region(gen, beg, end,
                RNASeq_cond=[], 
                signals_cond=[], 
                gene_names=True,
                output_dir=f"{resdir}Genome_plot/",
                output_file=f"Genome_{datetime.now()}",
                file_extension=".pdf",
                *args, **kwargs):
    """
    Plots experimental information (one subplot per experiment) and gene
    annotations for a region. The genes annotation subplot is created with
    the subplot_genes function. RNASeq subplots and Chipseq subplots
    are created with the subplot_rnaseq_coverage and subplot_signal functions
    respectively.

    Args:
        gen: Genome instance
        beg (int.): beginning of the region to be plotted
        end (int.): end of the region to be plotted
        RNASeq_cond (list of str.): list of the RNASeq conditions names to
            represent. Left empty by default ie no cover is plotted.
        tr_object (Transcriptome instance): Transcriptome instance with 
            loaded RNASeq coverages
        signals_cond (list of str.): list of the Chipseq conditions names to
            represent. Left empty by default ie no signal is plotted.
        ch_object (Chipseq instance): Chipseq instance with loaded signals.
        gene_names (Optional [Bool.]): If true, the name of the genes is noted on their
            representation. (True by default).
        output_dir (Optional [str.]): output directory
        output_file (Optional [str.]): output filename for plot
        file_extension (Optional [str.]): Graphic file extension type (.pdf by default)
        R_ylabels (Optional [list of str.]): y-axis label for RNASeq data subplots.
            By default, conditions names given in input with the RNASeq_cond
            argument.
        S_ylabels (Optional [list of str.]): y-axis label for RNASeq data subplots.
            By default, conditions names given in input with the signals_cond
            argument.
        vlines (Optional [dict.]): Dictionnary of shape {position: annotation text}.
            One green vertical line and its annotation are plotted for each
            item of this dictionnary.
        figsize (Optional [(float,float)]): width and height in inches
            (by default: (3.5,n) with n the number of suplots)
        hratios (Optional [array-like]): Relative heights of the subplots. If not given,
            all rows will have the same height except the genes subplot.
        hspace (Optional [float]): Height space between subplots, expressed as a 
            fraction of the average axis height. By default: 0.05
        figpos (Optional [array-like of length 4]): Position of the subplots as a 
            fraction of figure width or height. The array contains the 
            positions in the following order: bottom, top, left and right.

    
    Warning:
        This method needs a genomic annotation. If no annotation is
        loaded to the Genome instance, the load_annotation method of with
        the default "sequence.gff3" file is computed. To use another 
        annotation, please load an annotation to your Genome instance before 
        using this method.

    Warning:
        To plot a ChIPseq coverage, all the selected conditions must first be loaded 
        on the Chipseq instance using one of the ChipSeq method (for example
        load_signal or load_signals_average). This Chipseq instance must
        then be given as input with the argument "ch_object".

    Warning:
        To plot a RNAseq coverage, all the selected conditions must first be loaded 
        on the Transcriptome instance using tr.load_rnaseq_cov() and this
        instance must be given as input with the argument "tr_object".

    Example:
        >>> g = Genome.Genome("ecoli")
        >>> g.load_annotation()
        >>> tr = Transcriptome.Transcriptome("ecoli")
        >>> tr.load_rnaseq_cov()
        >>> ch = Chipseq.Chipseq("ecoli")
        >>> ch.load_signals_average(list_cond=["cond1","cond2"],
        ...                         average_name="cond12")
        >>> plot_genome.plot_region(g,2101000,2106000,
                                    RNASeq_cond = ["WT"],
                                    signals_cond=["cond12"],
                                    R_ylabels=["WT"],
                                    S_ylabels=["cond12"],
                                    tr_object=tr,
                                    ch_object=ch,
                                    hratios = [1,2,2,1],
                                    figsize=(6,5),
                                    figpos=[0.1,0.98,0.15,0.95])
    """

    # loads RNASeq data
    if len(RNASeq_cond) != 0:
        try:
            tr = kwargs.get('tr_object')
        except BaseException:
            sys.exit("Please give a Transcriptome instance with loaded RNASeq coverages in input with the 'tr_object' argument.")
        error = set(RNASeq_cond) - set(tr.rnaseq_cov_pos.keys())
        if len(error) != 0:
            sys.exit(
                f"{error} not found in the RNASeq coverages. Available coverages are {list(tr.rnaseq_cov_pos.keys())}")
        if isinstance(RNASeq_cond, str):
            RNASeq_cond = [RNASeq_cond]
        R_ylabels = kwargs.get('R_ylabels', RNASeq_cond)
        if isinstance(R_ylabels, str):
            R_ylabels = [R_ylabels]

    # loads ChipSeq data
    if len(signals_cond) != 0:
        try:
            ch = kwargs.get('ch_object')
            ch.get_all_signals()
        except BaseException:
            sys.exit("Please give a Chipseq instance with loaded signal in input with the 'ch_object' argument.")
        error = set(signals_cond) - set(ch.all_signals.keys())
        if len(error) != 0:
            sys.exit(f"{error} not found in the loaded signals. Available signals are {list(ch.all_signals.keys())}")
        if isinstance(signals_cond, str):
            signals_cond = [signals_cond]
        S_ylabels = kwargs.get('S_ylabels', signals_cond)
        if isinstance(S_ylabels, str):
            S_ylabels = [S_ylabels]

    # figure settings
    l = len(RNASeq_cond) + len(signals_cond)
    if l < 4:
        height_r = [1] + [2] * l
    else:
        height_r = [1] * (l + 1)
    hratios = kwargs.get("hratios", height_r)
    hspace = kwargs.get("hspace", 0.05)
    figsize = kwargs.get("figsize", (3.5, l + 1))
    fig = plt.figure(figsize=figsize)
    figpos = kwargs.get("figpos", [0.15, 0.98, 0.25, 0.9])
    gs = gridspec.GridSpec(l + 1, 1, height_ratios=hratios, bottom=figpos[0],
                           top=figpos[1], left=figpos[2], right=figpos[3],
                           hspace=hspace)
    # xticks
    dist = end - beg
    tickda = np.array([100, 200, 500, 1000, 2000, 5000, 10000,
                      20000, 50000, 100000, 200000, 500000, 1000000])
    tickdist = tickda[np.argmin(np.abs(1 - 5 * tickda / (end - beg)))]
    xticks = np.arange(
        (beg / tickdist) * tickdist,
        (end / tickdist + 1) * tickdist,
        tickdist)
    xticks = [None] * l + [xticks]

    # genes annotation subplot
    ax = plt.subplot(gs[0])
    subplot_genes(ax, gen, beg, end, xticks=xticks[0], gene_names=gene_names)
    vlines = kwargs.get('vlines', None)
    if vlines is not None:
        for v in vlines.keys():
            ax.axvline(x=v - beg, color="green", lw=1, ls="-")
            ax.annotate(vlines[v], xy=(v - beg, 0), color="green", va='center')

    n = 1
    # RNASeq subplots
    for nr, Rcond in enumerate(RNASeq_cond):
        ax = plt.subplot(gs[n])
        subplot_rnaseq_coverage(ax, tr, Rcond, beg, end,
            ylabel=R_ylabels[nr],
            xticks=xticks[n])
        n += 1

    # ChipSeq subplots
    for ns, Scond in enumerate(signals_cond):
        ax = plt.subplot(gs[n])
        subplot_signal(ax,ch,Scond,beg,end,
            ylabel=S_ylabels[ns],
            xticks=xticks[n])
        n += 1

    fig.align_ylabels()
    ax.set_xlabel("Genomic position")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    pathf = output_dir + output_file + file_extension
    plt.savefig(pathf,transparent=False,facecolor='white')
    print(f"Saved as {pathf}")
    plt.show()
    plt.close()


def subplot_genes(ax, gen, beg, end, xticks=None, gene_names=True):
    """
    Plots a subplot with the annoted oriented genes that are located in
    the selected region.

    Args:
        ax (matplotlib Axes instance)
        gen (Genome instance)
        beg (int.): beginning of the region to be plotted
        end (int.): end of the region to be plotted
        xticks (list of float): xaxis' tick locations.
            None by default ie. the xaxis is left empty.
        gene_names (Bool.): If true, the name of the genes is noted on their
            representation. (True by default).

    Warnings:
        This method needs a genomic annotation. If no annotation is
        loaded to the Genome instance, the load_annotation method of with
        the default "sequence.gff3" file is computed. To use another 
        annotation, please load an annotation to your Genome instance before using this method.
    """
    if not hasattr(gen, "genes"):
        gen.load_annotation()

    if xticks is not None:
        ax.set_xticks(xticks)
    else:
        ax.set_xticklabels([])
        ax.set_xticks([])

    length = end - beg
    text_up = False

   # take only genes in region
    gene_list = [g for g in gen.genes.values() if g.right >
                 beg and g.left < end]

    # sort list
    gene_start = np.argsort([(g.left + g.right) / 2 for g in gene_list])
    gene_list = [gene_list[i] for i in gene_start]

    for i, g in enumerate(gene_list):

        if g.strand:
            # when arrow is too little, a rectangle is drawn instead
            if g.end - g.start < length / 20:
                hlenght = 0
            else:
                hlenght = 1
            ax.quiver(g.start, 1, g.end - g.start, 0, angles='xy',
                scale_units='xy', scale=1, width=0.05, headwidth=1, 
                headlength=hlenght, headaxislength=1, minlength=0)

            if gene_names and beg < (g.end + g.start) / 2 < end:
                if g.end - g.start > length / 10:
                    ax.text(
                        (g.start + g.end) / 2,
                        1,
                        "$\\it{%s}$" % (g.name),
                        size=8,
                        horizontalalignment='center',
                        verticalalignment='center',
                        color="white")
                else:
                    text_up = True  # if a gene name is written above the arrow, ylim will be modified
                    ax.text(
                        (g.start + g.end) / 2,
                        2.5,
                        "$\\it{%s}$" % (g.name),
                        size=8,
                        horizontalalignment='center',
                        verticalalignment='center',
                        color="black")

        else:
            if g.start - g.end < length / 20:
                hlenght = 0
            else:
                hlenght = 1
            ax.quiver(
                g.start,
                -2,
                g.end - g.start,
                0,
                angles='xy',
                scale_units='xy',
                scale=1,
                width=0.05,
                headwidth=1,
                headlength=hlenght,
                headaxislength=1,
                minlength=0)
            if gene_names and beg < (g.end + g.start) / 2 < end:
                if g.start - g.end > length / 10:
                    ax.text((g.start + g.end) / 2, -2, "$\\it{%s}$" %
                            (g.name), size=8, horizontalalignment='center', 
                             verticalalignment='center', color="w")
                else:
                    ax.text((g.start + g.end) / 2, -0.5, "$\\it{%s}$" %
                            (g.name), size=8, horizontalalignment='center', 
                            verticalalignment='center', color="black")

    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xlim(beg, end)
    if text_up:
        ax.set_ylim(-3.5, 3.5)
    else:
        ax.set_ylim(-3.5, 2.5)


def subplot_rnaseq_coverage(ax, tr, cond, beg, end, ylabel=None, xticks=None):
    """
    Plots a subplot with the with the experimental data (RNASeq coverages on
    both DNA strands) of the selected region.

    Args:
        ax (matplotlib Axes instance)
        tr (Transcriptome instance): Transcriptome instance with loaded
            RNASeq coverages
        cond (str.): RNASeq condition name to represent. 
            WARNING: This condition must first be loaded on the Transcriptome 
            instance using tr.load_rnaseq_cov().
        beg (int.): beginning of the region to be plotted
        end (int.): end of the region to be plotted
        ylabel (str.): y-axis label. By default, the condition name given
            in input with the cond argument.
        xticks (list of float): xaxis' tick locations.
            None by default ie. the xaxis is left empty.
    """
    cov_pos = tr.rnaseq_cov_pos[cond]
    cov_neg = tr.rnaseq_cov_neg[cond]

    x = np.arange(beg, end)
    ax.fill_between(x, cov_pos[beg:end], color="blue")
    ax.fill_between(x, -cov_neg[beg:end], color="red")
    ax.axhline(y=0, color="black", lw=1)

    if xticks is not None:
        ax.set_xticks(xticks)
    else:
        ax.set_xticklabels([])
        ax.set_xticks([])
    ax.set_xlim(beg, end)
    if ylabel is None:
        ylabel = cond
    ax.set_ylabel(ylabel)


def subplot_signal(ax, ch, cond, beg, end, ylabel=None, xticks=None):
    """
    Plots a subplot with the with the experimental data (Chipseq signal) of 
    the selected region.

    Args:
        ax (matplotlib Axes instance)
        ch (Chipseq instance): Chipseq instance with loaded signals.
        cond (str.): Chipseq condition name to represent. WARNING: This
            condition must first be loaded on the Chipseq instance using one 
            of the ChipSeq method (for example load_signal or 
            load_signals_average).
        beg (int.): beginning of the region to be plotted
        end (int.): end of the region to be plotted
        ylabel (str.): y-axis label. By default, the condition name given
            in input with the cond argument.
        xticks (list of float): xaxis' tick locations.
            None by default ie. the xaxis is left empty.
    """
    if not hasattr(ch, "all_signals"):
        ch.get_all_signals()
    coverage = ch.all_signals[cond]

    y = coverage[beg:end]
    x = np.arange(beg, end)

    test_neg = [i <= 0 for i in y]
    test_pos = [i > 0 for i in y]
    ax.fill_between(x, y, where=test_pos, color='b')
    ax.fill_between(x, y, where=test_neg, color='r')
    ax.axhline(y=0, color="black", lw=1)

    if xticks is not None:
        ax.set_xticks(xticks)
    else:
        ax.set_xticklabels([])
        ax.set_xticks([])
    ax.set_xlim(beg, end)
    if ylabel is None:
        ylabel = cond
    ax.set_ylabel(ylabel)
