#! /usr/bin/env python
# -*- coding: utf-8 -*-
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import GC
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from matplotlib import patches
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
from globvar import *
from datetime import datetime
import seaborn as sns
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats
import itertools as it
from useful_functions import *

params = {
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
   'axes.labelsize': 11,
   'font.size': 11,
   'legend.fontsize': 9,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'axes.linewidth':1.5, #0.8
   'axes.titlesize':11,
   'axes.spines.top':True,
   'axes.spines.right':True,
   'font.family': "Arial"
   }
plt.rcParams.update(params)

def load_domains(gen):
    '''
    Assigns each gene to its corresponding domain
    '''      
    if not hasattr(gen,"genes"):
        gen.load_annotation()

    d = []
    # Load domain names and coordinates
    with open(basedir+"data/"+gen.name+"/domains.txt","r") as f:
        header = next(f)
        for line in f:
            line=line.strip()
            line=line.split('\t')
            d.append((line[0], int(line[1]),int(line[2])))
            # shape (name, left, right)
    for gene in gen.genes.keys():
        g = gen.genes[gene]
        for name,start,end in d:
            if start < end:
                if g.start > start and g.start < end:
                    g.domain = name
            else: # domain overlapping on start / end
                if g.start > start or g.start < end:
                    g.domain = name

    # export file with each gene and its corresponding domain

    # pathDB = basedir+"data/"+gen.name+"/GO_analysis/domains.txt"                
    # file = open(pathDB,"w")
    # for gene in gen.genes.keys():
    #     file.write(gene+"\td"+gen.genes[gene].domain+"\n")
    # file.close()
    

def Shannon_entropy(seq):
    '''
    Computes Shannon entropy which quantifies randomness of codons utilization in a sequence
    Unit: bit / codon
    '''      
    H = 0  
    codons = [''.join(i) for i in it.product('ATCG', repeat = 3)] # all codons
    res = dict.fromkeys(codons, 0) 
    for el in [seq[i:i+3] for i in range(0,len(seq),1)]:
    # counts the number of overlapping codons in windows
        try:    
            res[el] += 1
        except: # remaining seq of length != 3 or unknown codon if e.g. sequencing incertitude
            pass
    tot = float(sum(res.values())) # nb of total codons in windows
    for k in res.keys(): # for each codon
        px = res[k]/tot # probability of having this particular codon
        if px != 0:
            H += (px*np.log2(px))
    return -H

def Gibbs_entropy(seq):
    '''
    Computes Gibbs entropy which quantifies thermodynamic stability of a sequence
    '''      
    # Unified nearest-neighbor free energy parameters for di-nucleotides using thermodynamic stability parameters of Watson-Crick bp in 1M NaCl at 37°C. Units = DG(Kcal/mol)
    dinucleotides = {'AA':-1,'TT':-1,'AT':-0.88,'TA':-0.58,'CA':-1.45,'TG':-1.45,'GT':-1.44,'AC':-1.44,'CT':-1.28,'AG':-1.28,'GA':-1.30, 'TC':-1.30,'CG':-2.17, 'GC':-2.24, 'GG':-1.84, 'CC':-1.84}
    codons = [''.join(i) for i in it.product('ATCG', repeat = 3)] # all codons
    kb1 = 1.38*10**-23 #  kb Boltzmann constant in J/K
    kb2 = 0.0019872041 # kb Boltzmann constant in kcal / (mol.K)
    T = 310.15 # 37°C in K

    energies = dict.fromkeys(codons, 0) # associates each codon a free energy depending on its dinucleotides e.g. ATG = AT + TG
    for k in energies.keys():
        energies[k] = (dinucleotides[k[0:2]] + dinucleotides[k[1:3]])

    # each codon is a microstate (64 microstates) and its probability of existence depends on its dinucleotides thermodynamic stabilities
    res = dict.fromkeys(codons, 0)
    for el in [seq[i:i+3] for i in range(0,len(seq),1)]:
    # counts the number of overlapping codons in windows            
        try:    
            res[el] += 1
        except: # remaining seq of length != 3
            pass

    for codon in res.keys(): # converts the nb of codons in probability depending on thermodynamic stability and occurence
        res[codon] = res[codon] * np.exp(-(energies[codon]/(kb2*T)))

    tot = sum(res.values()) # all proba
    G = 0
    for codon in res.keys():
        PGi = res[codon]/tot
        if PGi != 0:
            G += (PGi*np.log(PGi))
    return (-kb1*G)      

def Melting_energy(seq):
    '''
    Computes melting energy of a sequence
    '''      
    # Unified nearest-neighbor free energy parameters for di-nucleotides using thermodynamic stability parameters of Watson-Crick bp in 1M NaCl at 37°C. Units = DG(Kcal/mol)
    dinucleotides = {'AA':-1.,'TT':-1.,'AT':-0.88,'TA':-0.58,'CA':-1.45,'TG':-1.45,'GT':-1.44,'AC':-1.44,'CT':-1.28,'AG':-1.28,'GA':-1.30, 'TC':-1.30,'CG':-2.17, 'GC':-2.24, 'GG':-1.84, 'CC':-1.84}

    # melting energy 6 bp-averaged
    G = (np.sum([dinucleotides[k]*seq.count(k) for k in dinucleotides.keys()])/len(seq))*6
    return G

def computes_ME_genome(gen):
    '''
    Computes melting energy of all genes separately
    '''
    if not hasattr(gen,'genes'):
        gen.load_annotation()
    if not hasattr(gen,'seq'):
        gen.load_seq()

    # Unified nearest-neighbor free energy parameters for di-nucleotides using thermodynamic stability parameters of Watson-Crick bp in 1M NaCl at 37°C. Units = DG(Kcal/mol)
    dinucleotides = {'AA':-1.,'TT':-1.,'AT':-0.88,'TA':-0.58,'CA':-1.45,'TG':-1.45,'GT':-1.44,'AC':-1.44,'CT':-1.28,'AG':-1.28,'GA':-1.30, 'TC':-1.30,'CG':-2.17, 'GC':-2.24, 'GG':-1.84, 'CC':-1.84}
    for g in gen.genes.keys():
        seq = gen.seq[gen.genes[g].start:gen.genes[g].end] if gen.genes[g].strand else gen.seqcompl[gen.genes[g].end:gen.genes[g].start]
        gen.genes[g].ME = (np.sum([dinucleotides[k]*seq.count(k) for k in dinucleotides.keys()])/len(seq))

    # mean melting energy of genome
    seq = gen.seq    
    ME1 = np.sum([dinucleotides[k]*seq.count(k) for k in dinucleotides.keys()])/len(seq)
    seq = gen.seqcompl
    ME2 = np.sum([dinucleotides[k]*seq.count(k) for k in dinucleotides.keys()])/len(seq)
    gen.ME = np.mean([ME1,ME2])

def thermodynamic_properties(gen,windows=250000, increment=4000):
    '''
    Computes DNA thermodynamic properties (melting energy / GC content / Shannon / Gibbs entropies) all over the genome on sliding windows
    '''
    if not hasattr(gen,'seq'):
        gen.load_seq()

    bins = [] # bins = windows of the genome : [start coordinate,end coordinate]
    for i in range(1,len(gen.seq),increment): # create bins depending on windows size and increment value
        if (i+windows) <= len(gen.seq): # enough length to create a bin
            bins.append([i,i+windows])
        else: # i + windows > genome size, overlap with the beginning (circular chromosome)
            bins.append([i, windows - (len(gen.seq)-i)])
    unknown
    bins = np.array(bins) # convert to .npy

    new_bins = [] # bins with DNA properties
    for start,end in bins:
        ov = False if start < end else True # overlapping or not
        s = gen.seq[start-1:enunknownd] if not ov else gen.seq[start-1:] + gen.seq[0:end] # sequence of windows
        ME = Melting_energy(s) # melting energy of the windows
        H = Shannon_entropy(s) # shannon entropy of the windows
        GCcont = GC(s) # GC content of the windows
        new_bins.append([start,end,H,GCcont,ME])
        # G = Gibbs_entropy(s)
        # new_bins.append([start,end,ME,H,G])

    new_bins = np.array(new_bins) # convert to .npy
    res = (new_bins, windows, increment)
    return res

def genomic_properties(gen,windows=500000, increment=4000):
    '''
    Computes DNA genomic properties (orientation usage, leading strand usage, melting energy usage)
    all over the genome on sliding windows depending on condition
    '''
    if not hasattr(gen,'seq'):
        gen.load_seq()

    if not hasattr(gen,"genes_valid"):
        gen.load_fc_pval()
    if not hasattr(gen,"ME"): # melting energy required
        computes_ME_genome(gen) 
    if not hasattr(gen,"orientation"): # orientation required
        gen.load_gene_orientation()
    if not hasattr(gen,"annot_detailed"): # leading strand usage required
        load_annot_detailed(gen)

    bins = [] # bins = windows of the genome : [start coordinate,end coordinate,nb of activated genes,nb of repressed genes,nb of genes not affected]
    bins_overlap = [] # bins where end coordinate < start coordinate (overlap circular chromosome)
    for i in range(1,len(gen.seq),increment): # create bins depending on windows size and increment value
        if (i+windows) <= len(gen.seq): # enough length to create a bin
            start = i ; end = i+windows
            genes = [g for g in gen.genes.keys() if gen.genes[g].start >= start and gen.genes[g].start <= end]

        else: # i + windows > genome size, overlap with the beginning (circular chromosome)
            start = i ; end = windows - (len(gen.seq)-i)
            genes = [g for g in gen.genes.keys() if gen.genes[g].start >= start or gen.genes[g].start <= end]
        # genes = genes of windows
        bins.append([start,end,genes]) 

    def features(glist):
        '''
        Starting from a list of genes, return the melting energy, leading strand usage, and orientation of the genes
        '''
        orient = {"convergent":0, "divergent":0, "tandem":0, "isolated":0}
        nb = len(glist) ; MEs = [] ; LSs = [] 
        for g in glist:
            try:
                MEs.append(gen.genes[g].ME) # melting energy = float
                LSs.append(gen.genes[g].leading_strand) # leading strand usage = 0 (lagging) or 1 (leading)
                orient[gen.genes[g].orientation] += 1
            except:
                pass
        # nb of genes in list, melting energy, leading strand, orientation        
        return nb,MEs,LSs,orient

    # all over the genome
    allnb, allME, allLS, allorient = features(gen.genes.keys())
    pexpLS = np.mean(allLS) # proportion of leading strand among leading and lagging
    # proportion of convergent genes among convergent and divergent
    pexporient = float(allorient["convergent"]) / (allorient["convergent"] + allorient["divergent"])

    new_bins = []
    # for all windows
    for start,end,genes in bins:
        nbwind, MEwind, LSwind, orientwind = features(genes)
        # zscore quantifies difference compared to genome
        # see in useful functions
        zscore_LS = compute_zscore_binomial(len(LSwind), LSwind.count(1), pexpLS)
        zscore_orient = compute_zscore_binomial(orientwind["convergent"] + orientwind["divergent"], orientwind["convergent"], pexporient)
        zscore_ME = compute_zscore_2means(np.mean(MEwind), np.mean(allME), 0, np.std(MEwind), np.std(allME), len(MEwind), len(allME))

        new_bins.append([start,end,nbwind,zscore_LS,zscore_orient,zscore_ME])

    res = (np.array(new_bins), windows, increment)
    return res

def draw_thermo_graphes(gen, bins, *args, **kwargs): 
    '''
    Starting from thermodynamic_properties results (bins containing windows coordinates and features),
    draws either GC content / melting energy / Shannon / Gibbs entropies
    Possibility to draw domain borders to evaluate correlation with those features
    '''    
    new_bins = [] 
    for start,end,H,GCcont,ME in bins:
        if start < end: 
            coor = (float(start) + float(end)) / 2 # x coordinate of windows = middle
        else: # end of genome (circular)
            size = (float(len(gen.seq)) - float(start) + float(end))/2 
            coor = start + size if (start + size) < len(gen.seq) else end - size
        
        new_bins.append([coor,H,GCcont,ME])


    new_bins = np.array(new_bins)
    new_bins[:,0] = (new_bins[:,0])/1000000 # converts coordinates from pb to Mb for clearness

    data = np.split(new_bins, np.where(np.diff(new_bins[:,0]) < 0)[0]+1) # [10,20,1,2] -> [[10,20],[1,2]] to plot properly
   
    dom = True  # whether or not we represent domains borders stored in domains.txt

    width = 5 ; height = 3
    linew = 1.25 
    fig, ax1 = plt.subplots()

    # color = "red"
    # ax1.set_xlabel('Position (Mb)')
    # ax1.set_ylabel('Melting energy (Kcal / mol)', color=color)
    # #ax1.set_ylim(5.93,5.95)
    # #ax1.set_xlim(0,5)    
    # for d in data:
    #     ax1.plot(d[:,0], d[:,3], color=color, linewidth= linew)
    # ax1.axhline(y = np.mean(new_bins[:,3]), linewidth = 0.5, color=color, linestyle="dashed")
    # ax1.tick_params(axis='y', labelcolor=color)

    color = "red"
    ax1.set_xlabel('Position (Mb)')
    ax1.set_ylabel('GC content (%)', color=color)
    #ax1.set_ylim(5.93,5.95)
    #ax1.set_xlim(0,5)    
    for d in data:
        ax1.plot(d[:,0], d[:,2], color=color, linewidth= linew)
    ax1.axhline(y = np.mean(new_bins[:,2]), linewidth = 0.5, color=color, linestyle="dashed")
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis = coordinate
    color = "blue"
    ax2.set_xlabel('Position (Mb)')
    ax2.set_ylabel('Sequence randomness (bit/codon)', color=color)
    #ax1.set_ylim(5.93,5.95)
    # ax1.set_xlim(0,5)    
    for d in data:
        ax2.plot(d[:,0], d[:,1], color=color, linewidth= linew)

    ax2.tick_params(axis='y', labelcolor=color)
    # mean shannon entropy
    ax2.axhline(y = np.mean(new_bins[:,1]), linewidth = 0.5, color=color, linestyle="dashed")

    if dom: # draws domain borders
        d = [] # shape [(name, domain border 1, domain border 2)]
        with open(basedir+"data/"+gen.name+"/domains.txt","r") as f:
            header = next(f)
            for line in f:
                line=line.strip()
                line=line.split('\t')
                d.append((line[0], int(line[1]),int(line[2])))

        for name,front1,front2 in d:
            # converts into Mb
            front1 = front1/float(1000000)
            front2 = front2/float(1000000)
            mid = (front1 + front2) /2 if front2 > front1 else (len(gen.seq)/float(1000000) + front1)/2 
            ax1.axvline(x=front1, linestyle="dashed", linewidth = 0.5, color='black')        
            ax1.axvline(x=front2, linestyle="dashed", linewidth = 0.5, color='black')
            ax1.text(x=mid,y = 52, s = name, fontsize=8, horizontalalignment='center', fontweight="bold")

    ax1.set_title(gen.name)
    fig.set_size_inches(width, height)
    fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
    plt.tight_layout()
    plt.savefig(basedir+"data/"+gen.name+"/annotation/thermodynamic.png", dpi=500, transparent=False)        
    plt.savefig(basedir+"data/"+gen.name+"/annotation/thermodynamic.svg")        
    plt.close()


def draw_thermo_circles(gen, bins, wind, incr, *args, **kwargs):
    '''
    Starting from thermodynamic_properties results (bins containing windows coordinates and features),
    draws either GC content / melting energy / Shannon / Gibbs entropies in the shape of rings
    opt arg : colormap, vmin, vmax
    '''
    colormap= kwargs.get('colormap','jet') # default value
    try:
        cScale_fc = plt.get_cmap(colormap)
    except:
        print 'Incorrect colormap, please check https://matplotlib.org/users/colormaps.html'
        print 'Loading default'
        cScale_fc = plt.get_cmap('jet')

    typ = kwargs.get('type','thermodynamic')

    # if we start from thermodynamic_properties results
    if typ == "thermodynamic":
        ME = bins[0:,4]
        H = bins[0:,2]
        GCcont = bins[0:,3]
        #Gibbs = bins[0:,3]
        lgd = ["ME","H","GC"]
        lab = ["Melting energy (Kcal/mol)","Shannon entropy (bit/codon)","GC content (%)"]
        params = [ME,H,GCcont]

    # if we start from genomic_properties results
    elif typ == "genomic":
        LS = bins[0:,3]
        orient = bins[0:,4]
        ME = bins[0:,5]
        lgd = ["LS","conv","ME"]
        lab = ["Leading strand usage (red leading blue lagging)","Orientation usage (red conv blue div)","Melting energy usage (red more AT rich)"]
        params = [LS,orient,ME]

    name = gen.name.split("_")[1] if len(gen.name.split("_")) > 1 else gen.name.split("_")[0]

    z= 0 # idx for legend etc.
    for p in params:
        vmin = kwargs.get('vmin', min(p))
        vmax = kwargs.get('vmax', max(p))          
        # normalisation of colors
        #vmin = -4 ; vmax = 4
        cNorm_fc  = colors.Normalize(vmin=vmin, vmax=vmax) 
        # map which assigns a colour depending on value between vmin and vmax
        cMap_fc = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_fc) 
        # config, see globvar for more
        # init plots
        width = 3.5 ; height = width/1
        fig, ax = plt.subplots()
        #plt.axis([0, width, 0, height]) ; 
        ax.set_axis_off() 

        angle = 360.0/len(p) # angle between two fragments
        # display colour in the middle of the windows
        start_angle = angle * (wind/2 - incr/2)  / incr

        i=0
        for value in p:
            # edgecolor = assign a colour depending on value using cMap
            # draw arc on circle
            arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_fc.to_rgba(value),lw=7)
            ax.add_patch(arc)
            i+= angle

        cMap_fc._A = [] # fake array to print color map
        cbar = fig.colorbar(cMap_fc,fraction=0.025, pad=0.04)#,shrink=0.3)
        cbar.set_label(lab[z])
        tick_locator = ticker.MaxNLocator(nbins=4)
        cbar.locator = tick_locator
        cbar.update_ticks() 
        cbar.ax.tick_params(direction='out', labelleft=False, labelright=True,right=True, left=False)
        plt.annotate('OriC', xy=(center_x,radius), xytext=(center_x,radius+0.1),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
        plt.annotate('Ter', xy=(center_x,0), xytext=(center_x,-0.3),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
        plt.annotate(gen.name, xy=(center_x,radius), xytext=(center_x,center_y),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
        plt.axis("equal")
        fig.set_size_inches(width, height)
        plt.tight_layout()
        plt.savefig(basedir+"data/"+gen.name+"/annotation/thermodynamic{}.png".format(lgd[z]), dpi=500, transparent=True)        
        plt.savefig(basedir+"data/"+gen.name+"/annotation/thermodynamic{}.svg".format(lgd[z]))        
        plt.close()
        z += 1

def draw_expression_circles(gen, *arg, **kwargs):
    '''
    generate density circles based on FC and pvalues
    opt arguments : colormap = see matplotlib doc, vmin vmax (color normalisation), windows, increment
    meth (= act/repr/actrepr/overlay), format (= png/pdf/svg)
    '''
    windows= kwargs.get('windows', 500000)
    increment= kwargs.get('increment', 4000)

    path = basedir+"data/"+gen.name+"/fold_changes/circles-"+str(datetime.now())
    os.makedirs(path)

    if not hasattr(gen, 'genes_valid'): # if no fc loaded 
        try:
            print 'Trying to load FC...'
            gen.load_fc_pval()
            print 'FC loaded'
        except:
            print 'Unable to load fc'
            sys.exit()

    if not hasattr(gen, 'seq'): # if no seq loaded
        try:
            print 'Trying to load seq...'
            gen.load_seq()
            print 'seq loaded'
        except:
            print'Unable to load seq'
            sys.exit()

    # if we need to have precise ring size for each condition, for easier ring assembly
    sizes = {"WT_expo_vs_ihfA_expo":2.3, 
            "WT_stat_vs_ihfA_stat":2.7,
            "WT_nov_expo_vs_ihfA_nov_expo":3.1,
            "WT_nov_stat_vs_ihfA_nov_stat":3.5,
            "WT_nov_expo_vs_WT_expo":2.3,
            "WT_nov_stat_vs_WT_stat":2.7,
            "ihfA_nov_expo_vs_ihfA_expo":3.1,
            "ihfA_nov_stat_vs_ihfA_stat":3.5
            }   

    # sizes = {"WT_PGA_stat_vs_WT_stat":3.1, 
    #         "ihfA_PGA_stat_vs_ihfA_stat":3.5,
    #         "WT_PGA_nov_stat_vs_WT_PGA_stat":3.1,
    #         "ihfA_PGA_nov_stat_vs_ihfA_PGA_stat":3.5,
    #         "WT_PGA_stat_vs_ihfA_PGA_stat":3.1,
    #         "WT_PGA_nov_stat_vs_ihfA_PGA_nov_stat":3.5,
    #         }   

    for cond in gen.genes_valid.keys():
        print 'Computing condition',cond
        #gen_states = compute_state_genes_operon_correction(gen,cond)
        gen_states = compute_state_genes(gen,cond) # np array of shape [[gene start, state in condition]]
        # 0 = non affected, 1 = activated, -1 = repressed
        bins = count_genes_in_windows(gen.seq,cond, gen_states, windows, increment)
        print 'Windows on genome :\n',bins
        try:

            tot_act = len(gen_states[np.where(gen_states[:,1] == 1)])
            print 'Total activated genes on genome :',tot_act
            tot_repr = len(gen_states[np.where(gen_states[:,1] == -1)])
            print 'Total repressed genes on genome :',tot_repr
            tot_non = len(gen_states[np.where(gen_states[:,1] == 0)])
            print 'Total non affected genes on genome :',tot_non

            meth = kwargs.get('meth', 'actrepr')

            if meth == "act":
                zscores = compute_zscores_act(tot_act,tot_repr,tot_non,bins)
            elif meth == "repr":
                zscores = compute_zscores_repr(tot_act,tot_repr,tot_non,bins)
            elif meth == "actrepr":
                zscores = compute_zscores_actrepr(tot_act,tot_repr,tot_non,bins)
            elif meth == "overlay":
                zscores_act = compute_zscores_act(tot_act,tot_repr,tot_non,bins)
                zscores_repr = compute_zscores_repr(tot_act,tot_repr,tot_non,bins)
            elif meth == "genes":
                zscores_act = compute_act(tot_act,tot_repr,tot_non,bins)
                zscores_repr = compute_repr(tot_act,tot_repr,tot_non,bins)                
            elif meth == "genes_all":
                zscores = compute_all(tot_act,tot_repr,tot_non,bins)
            else:
                print "Unknown method, computing default..."
                zscores = compute_zscores_actrepr(tot_act,tot_repr,tot_non,bins)

        except Exception as e:
            print e
            print 'Invalid data (e.g. no genes repressed nor activated)'
            sys.exit()

        # init plots
        try:
            width = sizes[cond] ; height = width/1
        except:
            width = 3.5 ; height = width/1

        fig, ax = plt.subplots()
        ax.set_axis_off() 
       
        vmin = kwargs.get('vmin', -4)
        vmax = kwargs.get('vmax', +4)

        if meth == "overlay":
            # normalisation of colours
            try:
                cNorm_fc = colors.Normalize(vmin=vmin, vmax=vmax)
            except:
                print 'Incorrect values for vmin / vmax'
                print 'Loading default'
                cNorm_fc = colors.Normalize(vmin=-4, vmax=4)

            cScale_act = plt.get_cmap("Reds")
            cScale_repr = plt.get_cmap("Blues")
            cMap_act = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_act) 
            cMap_repr = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_repr) 
            angle = 360.0/len(zscores_act) # angle between two fragments
            # display colour in the middle of the windows
            start_angle = angle * (windows/2 - increment/2)  / increment 
            i=0
            for value in zscores_act:
                # edgecolor = assign a colour depending on value using cMap
                # draw arc on circle
                arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_act.to_rgba(value), lw=7, alpha = 0.7)
                ax.add_patch(arc)
                i+= angle
            angle = 360.0/len(zscores_repr) 
            start_angle = angle * (windows/2 - increment/2)  / increment
            i=0
            for value in zscores_repr:
                # edgecolor = assign a colour depending on value using cMap
                # draw arc on circle
                arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_repr.to_rgba(value), lw=7, alpha = 0.7)
                ax.add_patch(arc)
                i+= angle

            cMap_act._A = [] # fake array to print map
            cMap_repr._A = [] # fake array to print map
            cbar = plt.colorbar(cMap_act,shrink=0.3)
            cbar.set_label("Z-score act",size=9)
            cbar.ax.tick_params(labelsize=7)
            
            cbar = plt.colorbar(cMap_repr,shrink=0.3)
            cbar.set_label("Z-score repr",size=9)
            cbar.ax.tick_params(labelsize=7)
        
        elif meth == "genes":
            # normalisation of colours
            vmin = 0 ; vmax= 50
            cNorm_fc = colors.Normalize(vmin=min(min(zscores_act),min(zscores_repr)), vmax=max(max(zscores_act),max(zscores_repr)))
            cScale_act = plt.get_cmap("Reds")
            cScale_repr = plt.get_cmap("Blues")
            cMap_act = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_act) 
            cMap_repr = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_repr) 
            angle = 360.0/len(zscores_act) # angle between two fragments
            # display colour in the middle of the windows
            start_angle = angle * (windows/2 - increment/2)  / increment 
            i=0
            for value in zscores_act:
                # edgecolor = assign a colour depending on value using cMap
                # draw arc on circle
                arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_act.to_rgba(value), lw=7, alpha = 0.7)
                ax.add_patch(arc)
                i+= angle
            angle = 360.0/len(zscores_repr) 
            start_angle = angle * (windows/2 - increment/2)  / increment
            i=0
            for value in zscores_repr:
                # edgecolor = assign a colour depending on value using cMap
                # draw arc on circle
                arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_repr.to_rgba(value), lw=7, alpha = 0.7)
                ax.add_patch(arc)
                i+= angle

            cMap_act._A = [] # fake array to print map
            cMap_repr._A = [] # fake array to print map
            cbar = plt.colorbar(cMap_act,shrink=0.3)
            cbar.set_label("Nb act",size=9)
            cbar.ax.tick_params(labelsize=7)
            
            cbar = plt.colorbar(cMap_repr,shrink=0.3)
            cbar.set_label("Nb repr",size=9)
            cbar.ax.tick_params(labelsize=7)

            plt.annotate('OriC', xy=(center_x,radius), xytext=(center_x,radius+0.1),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
            plt.annotate('Ter', xy=(center_x,0), xytext=(center_x,-0.3),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
            plt.annotate(cond, xy=(center_x,radius), xytext=(center_x,center_y),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)

        elif meth == "genes_all":
            # Colormap for Zscore
            colormap= kwargs.get('colormap','Reds') # default value
            colormap= kwargs.get('colormap','jet') # default value
            cScale_fc = plt.get_cmap(colormap)
            # normalisation of colours
            cNorm_fc = colors.Normalize(vmin=min(zscores), vmax=max(zscores))
            cNorm_fc = colors.Normalize(vmin=20, vmax=150)
            # map which assigns a colour depending on value between vmin and vmax
            cMap_fc = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_fc) 
            angle = 360.0/len(zscores) # angle between two fragments
            # display colour in the middle of the Windows
            start_angle = angle * (windows/2 - increment/2)  / increment
            i=0
            for value in zscores:
                # edgecolor = assign a colour depending on value using cMap
                # draw arc on circle
                arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_fc.to_rgba(value), lw=7)
                ax.add_patch(arc)
                i+= angle

            cmap=True
            if cmap:           
                cMap_fc._A = [] # fake array to print map
                cbar = fig.colorbar(cMap_fc,fraction=0.025, pad=0.04)#,shrink=0.3)
                cbar.set_label("Nb DE genes")
                tick_locator = ticker.MaxNLocator(nbins=4)
                cbar.locator = tick_locator
                cbar.update_ticks()    
                cbar.ax.tick_params(direction='out', labelleft=False, labelright=True,right=True, left=False)

                plt.annotate('OriC', xy=(center_x,radius), xytext=(center_x,radius+0.1),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
                plt.annotate('Ter', xy=(center_x,0), xytext=(center_x,-0.3),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
                plt.annotate(cond, xy=(center_x,radius), xytext=(center_x,center_y),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)        


        else:        
            vmin = kwargs.get('vmin', -4)
            vmax = kwargs.get('vmax', +4)

            # Colormap for Zscore
            colormap= kwargs.get('colormap','jet') # default value
            try:
                cScale_fc = plt.get_cmap(colormap)
            except:
                print 'Incorrect colormap, please check https://matplotlib.org/users/colormaps.html'
                print 'Loading default'
                cScale_fc = plt.get_cmap('jet')
            # normalisation of colours
            try:
                cNorm_fc = colors.Normalize(vmin=vmin, vmax=vmax)
            except:
                print 'Incorrect values for vmin / vmax'
                print 'Loading default'
                cNorm_fc = colors.Normalize(vmin=-4, vmax=4)

            # map which assigns a colour depending on value between vmin and vmax
            cMap_fc = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_fc) 

            angle = 360.0/len(zscores) # angle between two fragments
            # display colour in the middle of the Windows
            start_angle = angle * (windows/2 - increment/2)  / increment
            i=0
            for value in zscores:
                # edgecolor = assign a colour depending on value using cMap
                # draw arc on circle
                arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_fc.to_rgba(value), lw=7)
                ax.add_patch(arc)
                i+= angle

            cmap=False
            if cmap:           
                cMap_fc._A = [] # fake array to print map
                cbar = fig.colorbar(cMap_fc,fraction=0.025, pad=0.04)#,shrink=0.3)
                cbar.set_label("Z-score")
                tick_locator = ticker.MaxNLocator(nbins=4)
                cbar.locator = tick_locator
                cbar.update_ticks()    
                cbar.ax.tick_params(direction='out', labelleft=False, labelright=True,right=True, left=False)

                plt.annotate('OriC', xy=(center_x,radius), xytext=(center_x,radius+0.1),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
                plt.annotate('Ter', xy=(center_x,0), xytext=(center_x,-0.3),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
                plt.annotate(cond, xy=(center_x,radius), xytext=(center_x,center_y),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
        
        plt.axis("equal")
        fig.set_size_inches(width, height)
        plt.tight_layout()
        form = kwargs.get('format','png')
        if form == "png":
            plt.savefig(path+"/circle-"+cond.replace("/","-")+".png", dpi=500, transparent=True)
        elif form == "svg":
            plt.savefig(path+"/circle-"+cond.replace("/","-")+".svg", transparent=True) 
        elif form == "pdf":
            plt.savefig(path+"/circle-"+cond.replace("/","-")+".pdf", transparent=False)
        else:
            print 'Unknown format, computing default...'
            plt.savefig(path+"/circle-"+cond.replace("/","-")+".svg",  transparent=False)
        plt.close()

def compute_state_genes(gen, cond): 
    '''
    returns a npy where a row is a gene caracterised by a start pos and a gene state
    gene is considered activated above a given fc, repressed below a given fc
    '''
    gen_states = []
    for gene in gen.genes_valid[cond]:
        # if activated
        if gen.genes[gene].fc_pval[cond][0] >= fc_treshold_pos and gen.genes[gene].fc_pval[cond][1] <= pval_treshold:
            gen_states.append([gen.genes[gene].start,1])
        # if repressed
        elif gen.genes[gene].fc_pval[cond][0] <= fc_treshold_neg and gen.genes[gene].fc_pval[cond][1] <= pval_treshold:
            gen_states.append([gen.genes[gene].start,-1])
        # if not affected
        else:
            gen_states.append([gen.genes[gene].start,0])
    
    gen_states = np.array(gen_states)
    return gen_states

def compute_state_genes_operon_correction(gen, cond): 
    '''
    returns a npy where a row is a gene caracterised by a start pos and a gene state
    gene is considered activated above a given fc, repressed below a given fc
    '''
    gen.load_operon()
    not2do = []
    nb_operon = 0
    gen_states = []

    for ope in gen.operon_complete.keys():
        operon = gen.operon_complete[ope]
        expr = []
        expr_none = []
        try:
            for g in operon.genes:
                gene = gen.genes[g]
                try:
                    if gene.fc_pval[cond][1] <= pval_treshold:
                        expr.append(gene.fc_pval[cond][0])
                        not2do.append(g)
                    else:
                        expr_none.append(gene.fc_pval[cond][0])
                        not2do.append(g)

                except:
                    pass

            if expr != []:
                if np.mean(expr) <= fc_treshold_neg:
                    gen_states.append([operon.start,operon.start,-1])                
                elif np.mean(expr) >= fc_treshold_pos:
                    gen_states.append([operon.start,operon.start,1])                
                else:
                    if expr_none == []:
                        gen_states.append([operon.start,operon.start,0])
                nb_operon += 1

            elif expr_none != []:
                gen_states.append([operon.start,operon,operon.start,0])
                nb_operon += 1
        except:
            pass


    for gene in gen.genes_valid[cond]:
        if gene not in not2do:
            # if activated
            if gen.genes[gene].fc_pval[cond][0] >= fc_treshold_pos and gen.genes[gene].fc_pval[cond][1] <= pval_treshold:
                gen_states.append([gene,gen.genes[gene].start,1])
            # if repressed
            elif gen.genes[gene].fc_pval[cond][0] <= fc_treshold_neg and gen.genes[gene].fc_pval[cond][1] <= pval_treshold:
                gen_states.append([gene,gen.genes[gene].start,-1])
            # if not affected
            else:
                gen_states.append([gene,gen.genes[gene].start,0])
    
    gen_states = np.array(gen_states)

    print "NB OPERON :",nb_operon," CORRESPONDING TO GENES :",len(not2do)

    return gen_states



def count_genes_in_windows(seq, cond, gen_states, windows, increment):
    '''
    compute bins on the genome depending on windows size and increment, and
    calculate the nb of activated / repressed / non affected genes in each bin for further zscore
    '''
    bins = [] # bins = windows of the genome : [start coordinate,end coordinate,nb of activated genes,nb of repressed genes,nb of genes not affected]
    bins_overlap = [] # bins where end coordinate < start coordinate (overlap circular chromosome)
    for i in range(1,len(seq),increment): # create bins depending on windows size and increment value
        if (i+windows) <= len(seq): # enough length to create a bin
            bins.append([i,i+windows,0,0,0])
        else: # i + windows > genome size, overlap with the beginning (circular chromosome)
            bins_overlap.append([i, windows - (len(seq)-i),0,0,0])
    bins_overlap = np.array(bins_overlap)
    bins = np.array(bins) # convert to .npy
    
    for start,state in gen_states: # reminder gene state : a row = beginning of gene, state (activated, repressed or not affected)
        if state == 1: # activated
        # test to which bins the gene belongs to, and add one to the nb of activated genes of these bins
            bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),2] += 1
            bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),2] += 1
        elif state == -1: # repressed gene
            bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),3] += 1
            bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),3] += 1
        elif state == 0: # not affected gene
            bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),4] += 1
            bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),4] += 1
    
    bins = np.concatenate((bins,bins_overlap))
    return bins


def compute_zscores_actrepr(tot_act, tot_repr, tot_non, bins):
    '''
    p_exp = nb of total activated genes / nb of total activated + repressed genes on the whole genome
    '''
    zscores = []
    p_exp = float(tot_act) / float((tot_act + tot_repr))
    print 'g+ / (g+ + g-) on genome :',p_exp
    # compute zscore for each bin
    for start,end,nb_act,nb_repr,nb_null in bins:
        nb_act = float(nb_act) ; nb_repr = float(nb_repr)
        try:
            zscore =(nb_act - (nb_act+nb_repr)*p_exp) / (np.sqrt((nb_act+nb_repr)*p_exp*(1-p_exp)))
        except: # division by zero if no genes activated nor repressed
            zscore = 0
        zscores.append(zscore)
    return zscores

def compute_zscores_act(tot_act, tot_repr, tot_non, bins):
    '''
    p_exp = nb of total activated genes / nb of total genes on the whole genome
    '''
    zscores = []
    p_exp = float(tot_act) / float((tot_act + tot_repr + tot_non))
    print 'g+ / gtot on genome :',p_exp
    # compute zscore for each bin
    for start,end,nb_act,nb_repr,nb_null in bins:
        nb_act = float(nb_act) ; nb_tot = float((nb_act+nb_repr+nb_null))
        try:
            zscore =(nb_act - nb_tot*p_exp) / (np.sqrt(nb_tot*p_exp*(1-p_exp)))
        except: # division by zero if no genes activated nor repressed
            zscore = 0
        zscores.append(zscore)
    return zscores

def compute_zscores_repr(tot_act, tot_repr, tot_non, bins):
    '''
    p_exp = nb of total repressed genes / nb of total genes on the whole genome
    '''
    zscores = []
    p_exp = float(tot_repr) / float((tot_act + tot_repr + tot_non))
    print 'g- / gtot on genome :',p_exp
    # compute zscore for each bin
    for start,end,nb_act,nb_repr,nb_null in bins:
        nb_repr = float(nb_repr) ; nb_tot = float((nb_act+nb_repr+nb_null))
        try:
            zscore =(nb_repr - nb_tot*p_exp) / (np.sqrt(nb_tot*p_exp*(1-p_exp)))
        except: # division by zero if no genes activated nor repressed
            zscore = 0
        zscores.append(zscore)
    return zscores

def compute_act(tot_act, tot_repr, tot_non, bins):
    zscores = []
    for start,end,nb_act,nb_repr,nb_null in bins:
        zscores.append(nb_act)
    return zscores

def compute_repr(tot_act, tot_repr, tot_non, bins):
    zscores = []
    for start,end,nb_act,nb_repr,nb_null in bins:
        zscores.append(nb_repr)
    return zscores

def compute_all(tot_act, tot_repr, tot_non, bins):
    zscores = []
    for start,end,nb_act,nb_repr,nb_null in bins:
        zscores.append(nb_act+nb_repr)
    return zscores



def draw_all_circles(gen, *arg, **kwargs):
    '''
    Rings for ihfA paper
    For each expression condition, generates 3 rings: orientation usage, leading strand usage, melting energy usage
    '''
    windows= kwargs.get('windows', 500000)
    increment= kwargs.get('increment', 4000)

    path = basedir+"data/"+gen.name+"/fold_changes/circles-"+str(datetime.now())
    os.makedirs(path)

    if not hasattr(gen, 'genes_valid'): # if no fc loaded 
        gen.load_fc_pval()

    meth = kwargs.get("meth","genome")
    for cond in gen.genes_valid.keys():
        print 'Computing condition',cond
        #gen_states = compute_state_genes_operon_correction(gen,cond)
        bins = properties_windows(gen, cond, windows, increment, meth=meth)
        sizes = {}
        # init plots
        try:
            width = sizes[cond] ; height = width/1
        except:
            width = 3.5 ; height = width/1

        lgd = ["Leading strand usage (red = leading, blue = lagging)", "Orientation usage (red = conv, blue = div)", "Melting energy usage (red = low, blue = high)",
        "Leading strand usage (red = leading, blue = lagging)", "Orientation usage (red = conv, blue = div)", "Melting energy usage (red = low, blue = high)",
         "Leading strand usage (red = leading, blue = lagging)", "Orientation usage (red = conv, blue = div)", "Melting energy usage (red = low, blue = high)"]
        names = [cond.replace("/","-") +"-"+ meth + "-LS-actonly",cond.replace("/","-") +"-"+ meth + "-conv-actonly",cond.replace("/","-") +"-"+ meth + "-ME-actonly", cond.replace("/","-") +"-"+ meth + "-LS-reponly",cond.replace("/","-") +"-"+ meth + "-conv-reponly",cond.replace("/","-") +"-"+ meth + "-ME-reponly", cond.replace("/","-") +"-"+ meth + "-LS-actrep",cond.replace("/","-") +"-"+ meth + "-conv-actrep",cond.replace("/","-") +"-"+ meth + "-ME-actrep"]

        sizes_cond = {"WT_expo_vs_ihfA_expo":2.3, 
                "WT_stat_vs_ihfA_stat":2.7,
                "WT_nov_expo_vs_ihfA_nov_expo":3.1,
                "WT_nov_stat_vs_ihfA_nov_stat":3.5,
                "WT_nov_expo_vs_WT_expo":2.3,
                "WT_nov_stat_vs_WT_stat":2.7,
                "ihfA_nov_expo_vs_ihfA_expo":3.1,
                "ihfA_nov_stat_vs_ihfA_stat":3.5
                }   
        z = 0 ; sizes = [sizes_cond[cond]]*9
        for zscores in [bins[:,5],bins[:,6],bins[:,7],bins[:,8],bins[:,9],bins[:,10],bins[:,11],bins[:,12],bins[:,13]]:
            fig, ax = plt.subplots()
            ax.set_axis_off() 
           
            vmin = kwargs.get('vmin', -4)
            vmax = kwargs.get('vmax', +4)

            # Colormap for Zscore
            colormap= kwargs.get('colormap','jet') # default value
            cScale_fc = plt.get_cmap(colormap)
            # normalisation of colours
            cNorm_fc = colors.Normalize(vmin=vmin, vmax=vmax)
            # map which assigns a colour depending on value between vmin and vmax
            cMap_fc = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_fc) 

            angle = 360.0/len(zscores) # angle between two fragments
            # display colour in the middle of the Windows
            start_angle = angle * (windows/2 - increment/2)  / increment
            i=0
            for value in zscores:
                # edgecolor = assign a colour depending on value using cMap
                # draw arc on circle
                arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_fc.to_rgba(value), lw=7)
                ax.add_patch(arc)
                i+= angle

            cmap=False
            if cmap:           
                cMap_fc._A = [] # fake array to print map
                cbar = fig.colorbar(cMap_fc,fraction=0.025, pad=0.04)#,shrink=0.3)
                cbar.set_label("Z-score")
                tick_locator = ticker.MaxNLocator(nbins=4)
                cbar.locator = tick_locator
                cbar.update_ticks()    
                cbar.ax.tick_params(direction='out', labelleft=False, labelright=True,right=True, left=False)

                plt.annotate('OriC', xy=(center_x,radius), xytext=(center_x,radius+0.1),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
                plt.annotate('Ter', xy=(center_x,0), xytext=(center_x,-0.3),verticalalignment = 'center', horizontalalignment = 'center', wrap=True)
                plt.annotate(lgd[z], xy=(center_x,radius), xytext=(center_x,center_y),verticalalignment = 'center', horizontalalignment = 'center', wrap=True, fontsize=6)
        
            plt.axis("equal")
            width = sizes[z] ; height = sizes[z]
            fig.set_size_inches(width, height)
            plt.tight_layout()
            form = kwargs.get('format','png')
            if form == "png":
                plt.savefig(path+"/"+names[z]+".png", dpi=500, transparent=True)
            elif form == "svg":
                plt.savefig(path+"/"+names[z]+".svg", transparent=True) 
            elif form == "pdf":
                plt.savefig(path+"/"+names[z]+".pdf", transparent=False)
            else:
                print 'Unknown format, computing default...'
                plt.savefig(path+"/"+names[z]+".svg",  transparent=False)
            plt.close()

            z += 1


def properties_windows(gen, cond, windows, increment,meth="genome"):
    '''
    compute bins on the genome depending on windows size and increment, and
    calculate the nb of activated / repressed / non affected genes in each bin for further zscore
    computes leading strand usage / melting energy usage / orientation usage
    '''
    if not hasattr(gen,"genes_valid"):
        gen.load_fc_pval()
    if not hasattr(gen,"ME"):
        computes_ME_genome(gen)
    if not hasattr(gen,"orientation"):
        gen.load_gene_orientation()
    if not hasattr(gen,"annot_detailed"):
        load_annot_detailed(gen)

    bins = [] # bins = windows of the genome : [start coordinate,end coordinate,nb of activated genes,nb of repressed genes,nb of genes not affected]
    bins_overlap = [] # bins where end coordinate < start coordinate (overlap circular chromosome)
    for i in range(1,len(gen.seq),increment): # create bins depending on windows size and increment value
        if (i+windows) <= len(gen.seq): # enough length to create a bin
        # act = activated genes in windows, rep = repressed, non = non affected
            start = i ; end = i+windows
            act = [g for g in gen.genes_valid[cond] if gen.genes[g].start >= start and gen.genes[g].start <= end and gen.genes[g].fc_pval[cond][0] >= fc_treshold_pos and gen.genes[g].fc_pval[cond][1] <= pval_treshold]
            rep = [g for g in gen.genes_valid[cond] if gen.genes[g].start >= start and gen.genes[g].start <= end and gen.genes[g].fc_pval[cond][0] <= fc_treshold_neg and gen.genes[g].fc_pval[cond][1] <= pval_treshold]
            non = [g for g in gen.genes_valid[cond] if gen.genes[g].start >= start and gen.genes[g].start <= end and gen.genes[g].fc_pval[cond][1] > pval_treshold] + [g for g in gen.genes_valid[cond] if gen.genes[g].start >= start and gen.genes[g].start <= end and gen.genes[g].fc_pval[cond][1] <= pval_treshold and gen.genes[g].fc_pval[cond][0] < fc_treshold_pos and gen.genes[g].fc_pval[cond][0] > fc_treshold_neg]

        else: # i + windows > genome size, overlap with the beginning (circular chromosome)
            start = i ; end = windows - (len(gen.seq)-i)
            act = [g for g in gen.genes_valid[cond] if gen.genes[g].start >= start or gen.genes[g].start <= end and gen.genes[g].fc_pval[cond][0] >= fc_treshold_pos and gen.genes[g].fc_pval[cond][1] <= pval_treshold]
            rep = [g for g in gen.genes_valid[cond] if gen.genes[g].start >= start or gen.genes[g].start <= end and gen.genes[g].fc_pval[cond][0] <= fc_treshold_neg and gen.genes[g].fc_pval[cond][1] <= pval_treshold]
            non = [g for g in gen.genes_valid[cond] if gen.genes[g].start >= start or gen.genes[g].start <= end and gen.genes[g].fc_pval[cond][1] > pval_treshold] + [g for g in gen.genes_valid[cond] if gen.genes[g].start >= start or gen.genes[g].start <= end and gen.genes[g].fc_pval[cond][1] <= pval_treshold and gen.genes[g].fc_pval[cond][0] < fc_treshold_pos and gen.genes[g].fc_pval[cond][0] > fc_treshold_neg]

        bins.append([start,end,act,rep,non]) 

    def features(glist):
        '''
        Starting from a list of genes, return the melting energy, leading strand usage, and orientation of the genes
        '''
        orient = {"convergent":0, "divergent":0, "tandem":0, "isolated":0}
        nb = len(glist) ; MEs = [] ; LSs = [] 
        for g in glist:
            try:
                MEs.append(gen.genes[g].ME)
                LSs.append(gen.genes[g].leading_strand)
                orient[gen.genes[g].orientation] += 1
            except:
                pass
        return nb,MEs,LSs,orient

    # for whole condition
    allact = [g for g in gen.genes_valid[cond] if gen.genes[g].fc_pval[cond][0] >= fc_treshold_pos and gen.genes[g].fc_pval[cond][1] <= pval_treshold]
    allrep = [g for g in gen.genes_valid[cond] if gen.genes[g].fc_pval[cond][0] <= fc_treshold_neg and gen.genes[g].fc_pval[cond][1] <= pval_treshold]
    allnon = [g for g in gen.genes_valid[cond] if gen.genes[g].fc_pval[cond][1] > pval_treshold] + [g for g in gen.genes_valid[cond] if gen.genes[g].fc_pval[cond][0] > fc_treshold_neg and gen.genes[g].fc_pval[cond][0] < fc_treshold_pos]

    allnbact, allMEact, allLSact, allorientact = features(allact)
    allnbrep, allMErep, allLSrep, allorientrep = features(allrep)
    allnbnon, allMEnon, allLSnon, allorientnon = features(allnon)

    # Meth = whether we compare windows features to genome static features or to 
    # features of all genes which respond in the condition for zscore generation
    if meth == "genome":
        allnb, allME, allLS, allorient = features(gen.genes.keys())
        pexpLS = np.mean(allLS)
        pexporient = float(allorient["convergent"]) / (allorient["convergent"] + allorient["divergent"])

    elif meth =="cond":
        pexpLS = np.mean(allLSact + allLSrep)
        pexporient = (float(allorientact["convergent"]) + allorientrep["convergent"]) / (allorientact["convergent"] + allorientact["divergent"] + allorientrep["convergent"] + allorientrep["divergent"])
        allME = allMEact + allMErep

    new_bins = []
    for start,end,act,rep,non in bins:
        nbact, MEact, LSact, orientact = features(act)
        nbrep, MErep, LSrep, orientrep = features(rep)
        nbnon, MEnon, LSnon, orientnon = features(non)

        # if we consider only activated genes
        zscore_LSact = compute_zscore_binomial(len(LSact), LSact.count(1), pexpLS)
        zscore_orientact = compute_zscore_binomial(orientact["convergent"] + orientact["divergent"], orientact["convergent"], pexporient)
        zscore_MEact = compute_zscore_2means(np.mean(MEact), np.mean(allME), 0, np.std(MEact), np.std(allME), len(MEact), len(allME))
        # if we consider only repressed genes
        zscore_LSrep = compute_zscore_binomial(len(LSrep), LSrep.count(1), pexpLS)
        zscore_orientrep = compute_zscore_binomial(orientrep["convergent"] + orientrep["divergent"], orientrep["convergent"], pexporient)
        zscore_MErep = compute_zscore_2means(np.mean(MErep), np.mean(allME), 0, np.std(MErep), np.std(allME), len(MErep), len(allME))
        # if we consider both activated and repressed genes
        zscore_LSactrep = compute_zscore_binomial(len(LSact) + len(LSrep), LSact.count(1) + LSrep.count(1), pexpLS)
        zscore_orientactrep = compute_zscore_binomial(orientact["convergent"] + orientact["divergent"] + orientrep["convergent"] + orientrep["divergent"], orientact["convergent"] + orientrep["convergent"], pexporient)
        zscore_MEactrep = compute_zscore_2means(np.mean(MEact+MErep), np.mean(allME), 0, np.std(MEact+MErep), np.std(allME), len(MEact) + len(MErep), len(allME))


        new_bins.append([start,end,nbact,nbrep,nbnon, zscore_LSact,zscore_orientact,zscore_MEact,zscore_LSrep,zscore_orientrep,zscore_MErep, zscore_LSactrep,zscore_orientactrep,zscore_MEactrep])

    return np.array(new_bins)


