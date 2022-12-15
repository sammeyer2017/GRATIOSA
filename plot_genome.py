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


def subplot_genes(ax, gen, beg, end,buf=1000, ticks=None, plot_gene_names=False):
    """ 
    plots a small subplot with oriented genes
    """
    if ticks is not None:
        ax.set_xticks(ticks)

    length = end-beg
    text_up = False

    # take only genes in region
    gene_list=[g for g in gen.genes.values() if (g.left+g.right)/2>(beg-buf) and (g.left+g.right)/2<(end+buf)]

    # sort list
    gene_start=np.argsort([(g.left+g.right)/2 for g in gene_list])
    gene_list=[gene_list[i] for i in gene_start]
     
    for i,g in enumerate(gene_list):
        #print(g.name, g.start, g.end)
        if g.orientation == 1 :
            #when arrow is too little, a rectangle is drawn instead
            if g.end-g.start < length/30 :
                hlenght = 0
            else :
                hlenght = 1
            ax.quiver(g.start,1,g.end-g.start,0,angles='xy', scale_units='xy', scale=1,width=0.03,headwidth=1,headlength=hlenght,headaxislength=1,minlength=0)      
           
            if plot_gene_names and beg < (g.end+g.start)/2 < end :
                if g.end-g.start > length/10 :
                    ax.text((g.start+g.end)/2,1,"$\it{%s}$"%(g.name),size=5,horizontalalignment='center',verticalalignment='center',color = "w")
                else : 
                    text_up = True #if a gene name is written above the arrow, ylim will be modified 
                    ax.text((g.start+g.end)/2,2,"$\it{%s}$"%(g.name),size=5,horizontalalignment='center',verticalalignment='center',color = "black")

        elif g.orientation == -1 :
            if g.start-g.end < length/30 :
                hlenght = 0
            else :
                hlenght = 1            
            ax.quiver(g.start,-1,g.end-g.start,0,angles='xy', scale_units='xy', scale=1,width=0.03,headwidth=1,headlength=hlenght,headaxislength=1,minlength=0) 
            if plot_gene_names and beg < (g.end+g.start)/2 < end :
                if g.start-g.end > length/10 :
                    ax.text((g.start+g.end)/2,-1,"$\it{%s}$"%(g.name),size=5,horizontalalignment='center',verticalalignment='center',color = "w")
                else : 
                    ax.text((g.start+g.end)/2,0,"$\it{%s}$"%(g.name),size=5,horizontalalignment='center',verticalalignment='center',color = "black")
  
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xlim(beg,end)
    if text_up :
        ax.set_ylim(-2.5,3)
    else : 
        ax.set_ylim(-2.5,2.5)


def subplot_rnaseq_coverage(ax,gen,cond,beg,end,binsize,ticks=None):
    """ 
    plots experimental information for region
    gen= genome with loaded coverage
    ##INTEGRER LES DONNEES SANS TRAITEMENET ET SMOOTH
    """
    f_name = cond + "_" + str(binsize) + "b"

    if not hasattr(gen,"cov_rnaseq_pos_bin") : 
        gen.load_cov_rnaseq_bin(binsize,cond)

    elif f_name not in gen.cov_rnaseq_pos_bin.keys() :
        gen.load_cov_rnaseq_bin(binsize,cond)

    x=np.arange(beg,end,binsize)
    ax.fill_between(x,gen.cov_rnaseq_pos_bin[f_name][int(beg/binsize):int(end/binsize)],color="blue")
    ax.fill_between(x,-gen.cov_rnaseq_neg_bin[f_name][int(beg/binsize):int(end/binsize)],color="red")
    ax.axhline(y=0,color="black",lw=1)
  
    if ticks is not None:
        ax.set_xticks(ticks)
    #ax.set_xticklabels([])
    ax.set_xlim(beg,end)
    ax.set_ylabel("coverage (exp)")

def subplot_chipseq_coverage(ax,gen,cond,beg,end,data_treatment,replicates = True,ticks=None,*args, **kwargs):
    """ 
    plots experimental information for region
    gen= genome with loaded coverage
    """
        
    if data_treatment == "bin" :
        size = kwargs.get('binsize')
        c_name = cond+"_"+data_treatment+str(size)+"b"
        if size == None :
            sys.exit("please give a binsize")

    elif data_treatment == "smooth" :
        size = kwargs.get('window')
        c_name = cond+"_"+data_treatment+str(size)+"b"
        if size == None : 
            sys.exit("please give a window size")

    else :
        print("no smoothing nor binning")
        c_name = cond
        size = 1

    #to plot a single experimental coverage
    if replicates == False :
        if data_treatment == "bin" :
            if not hasattr(gen,"cov_chipseq_bin") :
                gen.load_cov_chipseq_bin(size,cond)
            elif c_name not in gen.cov_chipseq_bin.keys() :
                gen.load_cov_chipseq_bin(binsize,cond)
            y=gen.cov_chipseq_bin[c_name][int(beg/size):int((end+size)/size)]
            y=np.repeat(y,size)
       
        elif data_treatment == "smooth" :
            if not hasattr(gen,"cov_chipseq_smooth") :
                gen.load_cov_chipseq_smooth(size,cond)
            elif c_name not in gen.cov_chipseq_smooth.keys() :
                gen.load_cov_chipseq_smooth(size,cond)
            y=gen.cov_chipseq_smooth[c_name][beg:end]
        else :
            if not hasattr(gen,"cov_chipseq") :
                gen.load_cov_chipseq(cond)
            elif c_name not in gen.cov_chipseq.keys() :
                gen.load_cov_chipseq(cond)
            y=gen.cov_chipseq[c_name][beg:end+1]

    #to plot the mean coverage of replicates 
    else : 
        if not hasattr(gen,"cov_chipseq_average") :
            gen.load_cov_chipseq_average(cond,data_treatment=data_treatment,binsize = size,window = size)
        elif c_name not in gen.cov_chipseq_average.keys() :
            gen.load_cov_chipseq_average(cond,data_treatment=data_treatment,binsize = size,window = size)

        if data_treatment == "bin" :
            y=np.repeat(gen.cov_chipseq_average[c_name][int(beg/size):int((end+size)/size)],size)
        else : 
            y=gen.cov_chipseq_average[c_name][beg:end]

    x=np.arange(beg,len(y)+beg,1)

    test_neg = [i<=0 for i in y] 
    test_pos = [i> 0 for i in y] 
    ax.fill_between(x,y,where = test_pos, color = 'b')
    ax.fill_between(x,y,where = test_neg , color = 'r')
    ax.axhline(y=0,color="black",lw=1)

    if ticks is not None:
        ax.set_xticks(ticks) 
    #ax.set_xticklabels([])
    ax.set_xlim(beg,end)
    ax.set_ylabel("log2FC(pdown/mock)")


def plot_region(gen, cond, beg,end,exp,data_treatment=None,replicates=False, name=None, plot_gene_names=True,*args, **kwargs):
    """ 
    plots experimental information and gene annotations for a region
    gen= genome with loaded coverage
    exp can be "rnaseq" (with cov_rnaseq_pos and cov_rnaseq_neg) or "chipseq" (only one cov)
    """

    if not hasattr(gen,"genes") :  ## A METTRE DANS SUBPLOT_GENES OU CHANGER SCHEMA
        print("loading annotation")
        gen.load_annotation()
        
  
    fig=plt.figure(figsize=(w,sum(h)+2*dh))
    
    if type(cond) == str: 
        cond = [cond]
    l = len(cond)

    if l < 4 :
        hr = [1] + [2]*l
    else : 
        hr = [1]*(l+1)

    gs = gridspec.GridSpec(l+1, 1, height_ratios=hr, bottom=.15, top=.98, left= .15, right=.95, hspace=dh)

    # ticks
    dist = end - beg
    tickda=np.array(tickd)
    tickdist=tickda[np.argmin(np.abs(1-(ticknb+1)*tickda/(end-beg)))]
    ticks=np.arange((beg/tickdist)*tickdist,(end/tickdist+1)*tickdist,tickdist)
    #if dist >= 5000 : 
     #   ticks = [t/1000 for t in ticks]

    # annotation
    ax=plt.subplot(gs[0])
    #subplot_genes(ax, gen, beg, end, buf=2000, ticks=ticks, plot_gene_names=plot_gene_names)
    #coli_ori = 3925859
    #ax.axvline(x=coli_ori,color="green",lw=1,ls="-")
    #ax.annotate('OriC', xy=(coli_ori, 0),color="green",va='center')

    if data_treatment == "bin" :
        size = kwargs.get('binsize')
        if size == None :
            sys.exit("please give a binsize")

    elif data_treatment == "smooth" :
        size = kwargs.get('window')
        if size == None : 
            sys.exit("please give a window size")

    else :
        print("no smoothing nor binning")
        c_name = cond
        size = 1

    if name is None:
        list_cond = []
        for co in cond :
            if len(cond) > 1 :
                list_cond = list_cond + ["_"] + [co]
            else : 
                list_cond = co
        if data_treatment in ["smooth","bin"] : 
            name=basedir+"plots/%s_region_%d_%d_%s_%s_%db"%(gen.name,beg,end,list_cond,data_treatment,size)
        else : 
            name=basedir+"plots/%s_region_%d_%d_%s"%(gen.name,beg,end,list_cond)
      
    # exp
    if exp == "rnaseq" : 
        for n,c in enumerate(cond):
            ax=plt.subplot(gs[n+1])
            subplot_rnaseq_coverage(ax,gen,c,beg,end,binsize=binsize,ticks=ticks)
            
    elif exp == "chipseq":
        for n,c in enumerate(cond):
            ax=plt.subplot(gs[n+1])
            subplot_chipseq_coverage(ax,gen,c,beg,end,data_treatment = data_treatment,binsize=size,window=size,ticks=ticks,replicates=replicates)
    
    # model
    #ax=plt.subplot(gs[2])
    #subplot_model(ax,gen,cond,beg,end, ticks=ticks)

    #for ext in exts:
    #    print(name+ext)
    plt.savefig(name+".pdf")
    plt.close()
    return 0


############### SAM ######################

def old_subplot_exp(ax,gen,cond,beg,end, ticks=None):
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

def old_subplot_model(ax, gen, cond, beg, end, ticks=None):
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


def old_plot_region(gen, cond, beg,end, name=None, plot_gene_names=False):
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
    subplot_genes(ax, gen, beg, end, buf=2000, ticks=ticks, plot_gene_names=plot_gene_names)
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

def old_plot_TU(gen, cond, TU_index):
    plot_region(gen, cond, (gen.TU[TU_index].left-1000,gen.TU[TU_index].right+1000), name=basedir+"calibration/plots/TU_%d_%s"%(TU_index,cond))
    return 0


def old_subplot_annotation_list(ax,gene_list,beg,end, ticks=None, plot_gene_names=False):
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
 
def old_subplot_annotation(ax, gen, beg, end, buf=1000, ticks=None, plot_gene_names=False):
    if ticks is not None:
        ax.set_xticks(ticks)
    # take only genes in region
    gene_list=[g for g in gen.genes.values() if (g.left+g.right)/2>(beg-buf) and (g.left+g.right)/2<(end+buf)]
    # sort list
    gene_start=np.argsort([(g.left+g.right)/2 for g in gene_list])
    gene_list=[gene_list[i] for i in gene_start]
    subplot_annotation_list(ax,gene_list,beg,end,ticks=ticks, plot_gene_names=plot_gene_names)
    return 0 


def old_subplot_chipseq_coverage(ax,gen,cond,beg,end,coverage,ticks=None):
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