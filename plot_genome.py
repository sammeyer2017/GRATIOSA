from globvar import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.font_manager as font_manager
from matplotlib import gridspec
import matplotlib.patches as pat
import numpy as np

plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'legend.fontsize': 10})
#plt.rcParams.update({'mathtex.fontset': "cm"})
plt.rcParams.update({'font.family': "Arial"})
#plt.rcParams.update({'usetex': True})


def subplot_annotation_list(ax,gene_list,beg,end, ticks=None, plot_gene_names=False):
    """ 
    plots a small subplot with oriented genes, with different colors corresponding to operons
    """
    if ticks is not None:
        ax.set_xticks(ticks)
    # define colors corresponding to operons
    operonvals=list(set([x.operon for x in gene_list]))
    cols=dict([(x, cm.rainbow(float(i)/len(operonvals))) for i,x in enumerate(operonvals)])
    ax.axhline(0, color='black')
    gw=0.2
    for i,g in enumerate(gene_list):
        print(g.name, g.start, g.end)
        if g.orientation:
            #ax.arrow(g.start,(1+i%2)*gw,g.end,(1+i%2)*gw,width=.1,fc=cols[g.operon],ec=cols[g.operon])
            # old version with ax.plot: rectangles are not well defined
            # ax.plot([g.start,g.end], [(1+i%2)*gw,(1+i%2)*gw],linewidth=2,color=cols[g.operon])
            rect=pat.Rectangle((g.start,(.75+i%2)*gw),g.end-g.start,gw,color=cols[g.operon])
            ax.add_artist(rect)
            if plot_gene_names:
                ax.text((g.start+g.end)/2,3*gw,g.name,size=6,horizontalalignment='center', color=cols[g.operon])
        else:
            rect=pat.Rectangle((g.start,-(1.75+i%2)*gw),g.end-g.start,gw,color=cols[g.operon])
            ax.add_artist(rect)
#            ax.plot([g.start,g.end], [-(1+i%2)*gw,-(1+i%2)*gw],linewidth=2,color=cols[g.operon])
            if plot_gene_names:
                ax.text((g.start+g.end)/2,-3*gw-0.3,g.name,size=6,horizontalalignment='center', color=cols[g.operon])
    ax.set_xticklabels([])
    ax.set_yticks([])
    #ax.set_yticklabels([])
    ax.set_xlim(beg,end)
    ax.set_ylim(-1.5,1.5)
    return 0
#float(end-beg)/20

# ,head_width=.1,head_length=.1
# ,width=.5
def subplot_annotation(ax, gen, beg, end, buf=1000, ticks=None, plot_gene_names=False):
    if ticks is not None:
        ax.set_xticks(ticks)
    # take only genes in region
    gene_list=[g for g in gen.genes.values() if (g.left+g.right)/2>(beg-buf) and (g.left+g.right)/2<(end+buf)]
    # sort list
    gene_start=np.argsort([(g.left+g.right)/2 for g in gene_list])
    gene_list=[gene_list[i] for i in gene_start]
    subplot_annotation_list(ax,gene_list,beg,end,ticks=ticks, plot_gene_names=plot_gene_names)
    return 0    

def subplot_exp(ax,gen,cond,beg,end, ticks=None):
    """ 
    plots experimental information for region
    gen= genome with loaded data: 
     - coverage
     - TSS
    """
    if not hasattr(gen,"cov_plus"):
        print("ERROR: plotting requires coverage")
        return 1
    x=np.arange(beg,end)
    ax.fill_between(x,gen.cov_plus[cond][beg:end],color="blue")
    ax.fill_between(x,-gen.cov_minus[cond][beg:end],color="blue")
    ax.axhline(y=0,color="black",lw=2)
    if hasattr(gen,"TSS_complete"):
        for t in gen.TSS_complete["+"]:
            if t>beg and t<end:
                ax.axvline(x=t, ymin=0, color="red", ls="-", alpha=.5)
        for t in gen.TSS_complete["-"]:
            if t>beg and t<end:
                ax.axvline(x=t, ymax=0, color="red", ls="-", alpha=.5)
    if ticks is not None:
        ax.set_xticks(ticks)
    ax.set_xticklabels([])
    ax.set_xlim(beg,end)
    ax.set_ylabel("coverage (exp)")
    return 0

def subplot_model(ax, gen, cond, beg, end, ticks=None):
    x=np.arange(beg,end)
    ax.plot(x, np.zeros(len(x)), color="black")
    # ticks: around 5
    if ticks is not None:
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks/1000)
    ax.set_xlabel("position (kb)")
    ax.set_ylabel("coverage (model)")
    ax.set_xlim(beg,end)
    return 1


def plot_region(gen, cond, beg,end, name=None, plot_gene_names=False):
    if name is None:
        name=basedir+"plots/%d_%d_%s"%(beg,end,cond)
    fig=plt.figure(figsize=(w,sum(h)+2*dh))
    gs = gridspec.GridSpec(3, 1, height_ratios=h, bottom=.15, top=.98, left= .21, right=.98, hspace=dh)
    # ticks
    tickda=np.array(tickd)
    tickdist=tickda[np.argmin(np.abs(1-(ticknb+1)*tickda/(end-beg)))]
    ticks=np.arange((beg/tickdist)*tickdist,(end/tickdist+1)*tickdist,tickdist)
    # annotation
    ax=plt.subplot(gs[0])
    subplot_annotation(ax, gen, beg, end, buf=2000, ticks=ticks, plot_gene_names=plot_gene_names)
    # exp
    ax=plt.subplot(gs[1])
    subplot_exp(ax,gen,cond,beg,end, ticks=ticks)
    # model
    ax=plt.subplot(gs[2])
    subplot_model(ax,gen,cond,beg,end, ticks=ticks)
    #plt.tight_layout()
    for ext in exts:
        print(name+ext)
        plt.savefig(name+ext)
    plt.close()
    return 0

def plot_TU(gen, cond, TU_index):
    plot_region(gen, cond, (gen.TU[TU_index].left-1000,gen.TU[TU_index].right+1000), name=basedir+"calibration/plots/TU_%d_%s"%(TU_index,cond))
    return 0


