from globvar import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.font_manager as font_manager
from matplotlib import gridspec
import matplotlib.patches as pat
import numpy as np
from scipy import stats
import pandas as pd

plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'legend.fontsize': 10})
#plt.rcParams.update({'mathtex.fontset': "cm"})
plt.rcParams.update({'font.family': "Arial"})
#plt.rcParams.update({'usetex': True})


def subplot_genes(ax, gen, beg, end,buf=1000, xticks=None, gene_names=False):
    """ 
    plots a small subplot with oriented genes
    """
    if not hasattr(gen,"genes") :  
        print("loading annotation")
        gen.load_annotation()

    if xticks is not None:
        ax.set_xticks(xticks)
    else : 
        ax.set_xticklabels([])
        ax.set_xticks([])

    length = end-beg
    text_up = False

    # take only genes in region
    gene_list=[g for g in gen.genes.values() if (g.left+g.right)/2>(beg-buf) and (g.left+g.right)/2<(end+buf)]

    # sort list
    gene_start=np.argsort([(g.left+g.right)/2 for g in gene_list])
    gene_list=[gene_list[i] for i in gene_start]
     
    for i,g in enumerate(gene_list):

        if g.strand :
            #when arrow is too little, a rectangle is drawn instead
            if g.end-g.start < length/20:
                hlenght = 0
            else :
                hlenght = 1
            ax.quiver(g.start,1,g.end-g.start,0,angles='xy', scale_units='xy', scale=1,width=0.05,headwidth=1,headlength=hlenght,headaxislength=1,minlength=0)      
           
            if gene_names and beg < (g.end+g.start)/2 < end :
                if g.end-g.start > length/10 :
                    ax.text((g.start+g.end)/2,1,"$\it{%s}$"%(g.name),size=8,horizontalalignment='center',verticalalignment='center',color = "white")
                else : 
                    text_up = True #if a gene name is written above the arrow, ylim will be modified 
                    ax.text((g.start+g.end)/2,2.5,"$\it{%s}$"%(g.name),size=8,horizontalalignment='center',verticalalignment='center',color = "black")

        else :
            if g.start-g.end < length/20 :
                hlenght = 0
            else :
                hlenght = 1            
            ax.quiver(g.start,-2,g.end-g.start,0,angles='xy', scale_units='xy', scale=1,width=0.05,headwidth=1,headlength=hlenght,headaxislength=1,minlength=0) 
            if gene_names and beg < (g.end+g.start)/2 < end :
                if g.start-g.end > length/10 :
                    ax.text((g.start+g.end)/2,-2,"$\it{%s}$"%(g.name),size=8,horizontalalignment='center',verticalalignment='center',color = "w")
                else : 
                    ax.text((g.start+g.end)/2,-0.5,"$\it{%s}$"%(g.name),size=8,horizontalalignment='center',verticalalignment='center',color = "black")
  
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xlim(beg,end)
    if text_up :
        ax.set_ylim(-3.5,3.5)
    else : 
        ax.set_ylim(-3.5,2.5)


def subplot_rnaseq_coverage(ax,gen,cond,beg,end,ylabel,xticks=None):
    """ 
    plots experimental information for region
    gen= genome with loaded coverage
    """

    if not hasattr(gen,"cov_rnaseq_pos") : 
        gen.load_genomic_RNASeq_cov(cond)
    if cond in gen.cov_rnaseq_pos.keys() : 
        cov_pos = gen.cov_rnaseq_pos[cond]
        cov_neg = gen.cov_rnaseq_neg[cond]
    else :
        print(f"{cond} not found in the RNASeq coverages")

    x=np.arange(beg,end)
    ax.fill_between(x,cov_pos[beg:end],color="blue")
    ax.fill_between(x,-cov_neg[beg:end],color="red")
    ax.axhline(y=0,color="black",lw=1)
  
    if xticks is not None:
        ax.set_xticks(ticks)
    else : 
        ax.set_xticklabels([])
        ax.set_xticks([])
    ax.set_xlim(beg,end)
    ax.set_ylabel(ylabel)

def subplot_signal(ax,gen,cond,beg,end,ylabel,xticks=None):
    """ 
    plots experimental information for region
    gen= genome with loaded coverage
    """

    if not hasattr(gen,"all_signals") : 
        gen.load_genomic_signal()
    if cond in gen.all_signals.keys() : 
        coverage = gen.all_signals[cond]
    else :
        print(f"{cond} not found in the loaded signals")
    
    y = coverage[beg:end]
    x=np.arange(beg,end)

    test_neg = [i<=0 for i in y] 
    test_pos = [i> 0 for i in y] 
    ax.fill_between(x,y,where = test_pos, color = 'b')
    ax.fill_between(x,y,where = test_neg , color = 'r')
    ax.axhline(y=0,color="black",lw=1)

    if xticks is not None:
        ax.set_xticks(xticks) 
    else : 
        ax.set_xticklabels([])
        ax.set_xticks([]) 
    ax.set_xlim(beg,end)
    ax.set_ylabel(ylabel)


def plot_region(gen, beg,end,outname,RNASeq_cond=[],signals_cond=[],gene_names=True,*args, **kwargs):
    """ 
    plots experimental information and gene annotations for a region
    gen= genome with loaded coverage
    """  
  
    outdir = kwargs.get('outdir', resdir)
    R_ylabels = kwargs.get('R_ylabels', RNASeq_cond)
    S_ylabels = kwargs.get('S_ylabels', signals_cond)
    xlabel = kwargs.get('xlabel', "Genomic position")
    if type(R_ylabels) == str: R_ylabels = [R_ylabels]
    if type(S_ylabels) == str: S_ylabels = [S_ylabels]
    vlines = kwargs.get('vlines',None)     #{3925859:"OriC"}
    if type(RNASeq_cond) == str: RNASeq_cond = [RNASeq_cond]
    if type(signals_cond) == str: signals_cond = [signals_cond]
    
    l = len(RNASeq_cond)+len(signals_cond)
    if l < 4 :
        height_r = [1] + [2]*l
    else : 
        height_r = [1]*(l+1)
    hratios = kwargs.get("hratios",height_r)
    hspace = kwargs.get("hspace",dh)
    figsize = kwargs.get("figsize",(w,l+1))
    fig=plt.figure(figsize=figsize)
    figpos=kwargs.get("figpos",[0.15,0.98,0.25,0.9])
    gs = gridspec.GridSpec(l+1, 1, height_ratios=hratios, bottom=figpos[0], 
                            top=figpos[1], left=figpos[2], right=figpos[3], hspace=hspace)

    # xticks
    dist = end - beg
    tickda=np.array(tickd)
    tickdist=tickda[np.argmin(np.abs(1-(ticknb+1)*tickda/(end-beg)))]
    xticks=np.arange((beg/tickdist)*tickdist,(end/tickdist+1)*tickdist,tickdist)
    xticks = [None]*l+[xticks]

    # annotation
    ax=plt.subplot(gs[0])
    subplot_genes(ax, gen, beg, end, buf=2000, xticks=xticks[0], gene_names=gene_names)
    if vlines is not None :
        for v in vlines.keys() : 
            ax.axvline(x=v-beg,color="green",lw=1,ls="-")
            ax.annotate(vlines[v], xy=(v-beg, 0),color="green",va='center')

    n = 1
    for nr, Rcond in enumerate(RNASeq_cond) : 
        ax=plt.subplot(gs[n])
        subplot_rnaseq_coverage(ax,gen,Rcond,beg,end,ylabel=R_ylabels[nr],xticks=xticks[n])
        n += 1

    for ns,Scond in enumerate(signals_cond) : 
        ax=plt.subplot(gs[n])
        subplot_signal(ax,gen,Scond,beg,end,ylabel=S_ylabels[ns],xticks=xticks[n])
        n += 1
    fig.align_ylabels()
    ax.set_xlabel(xlabel)
    plt.savefig(outdir+outname+".pdf")
    print(f"Saved as {outdir}{outname}.pdf")
    plt.show()
    plt.close()
