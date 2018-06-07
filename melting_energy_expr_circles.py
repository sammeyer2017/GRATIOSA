#! /usr/bin/env python
# -*- coding: utf-8 -*-
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import numpy as np
from matplotlib import patches
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
from globvar import *
from datetime import datetime

def compute_melting_energy(seqgen,windows=500000, increment=4000):
    ''' 
    Compute melting energy on genome windows with a specific increment
    '''

    melting_energy = []

    bins = [] # bins = windows of the genome : [start coordinate,end coordinate]
    bins_overlap = [] # bins where end coordinate < start coordinate (overlap circular chromosome)

    for i in range(1,len(seqgen),increment): # create bins depending on windows size and increment value
        if (i+windows) <= len(seqgen): # enough length to create a bin
            bins.append([i,i+windows])
        else: # i + windows > genome size, overlap with the beginning (circular chromosome)
            bins_overlap.append([i, windows - (len(seqgen)-i)])
    
    bins = np.array(bins) # convert to .npy
    bins_overlap = np.array(bins_overlap)

    # compute melting energy on bins
    for start,end in bins:
        seq = Seq(seqgen[start-1:end])
        melting_energy.append(mt.Tm_GC(seq, strict=False))
        # if len(seq) < 1000:
        #     melting_energy = mt.Tm_NN(seq,strict=False)
        # else:
        #     melting_energy = mt.Tm_GC(seq,strict=False)

    for start,end in bins_overlap:
        seq = Seq(seqgen[start-1:] + seqgen[0:end])
        melting_energy.append(mt.Tm_GC(seq, strict=False))

    d = {'melting_energy':melting_energy,'windows':windows,'increment':increment}
    return d


def draw_melting_energy_circle(melting_energy,name, *args, **kwargs):

    '''
    Draw melting energy circle from melting energy
    opt arg : colormap, vmin, vmax
    '''

    colormap= kwargs.get('colormap','jet') # default value
    try:
        cScale_fc = plt.get_cmap(colormap)
    except:
        print 'Incorrect colormap, please check https://matplotlib.org/users/colormaps.html'
        print 'Loading default'
        cScale_fc = plt.get_cmap('jet')

    m_e = melting_energy['melting_energy']
    wind = melting_energy['windows']
    incr = melting_energy['increment']

    vmin = kwargs.get('vmin', min(m_e))
    vmax = kwargs.get('vmax', max(m_e))          
    # normalisation of colors
    cNorm_fc  = colors.Normalize(vmin=vmin, vmax=vmax) 
    # map which assigns a colour depending on value between vmin and vmax
    cMap_fc = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_fc) 
    # config, see globvar for more
    # init plots
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.subplot(1,1,1)
    plt.axis([0, fig_width, 0, fig_height]) ; ax.set_axis_off() 

    angle = 360.0/len(m_e) # angle between two fragments
    # display colour in the middle of the windows
    start_angle = angle * (wind/2 - incr/2)  / incr

    i=0
    for value in m_e:
        # edgecolor = assign a colour depending on value using cMap
        # draw arc on circle
        arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_fc.to_rgba(value), lw=7)
        ax.add_patch(arc)
        i+= angle

    cMap_fc._A = [] # fake array to print map
    cbar = plt.colorbar(cMap_fc,shrink=0.3)
    cbar.set_label("Melting temp °C",size=9)
    cbar.ax.tick_params(labelsize=7)
    plt.axis("equal")
    plt.annotate('OriC', xy=(center_x,radius), xytext=(center_x,radius+0.4),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=8)
    plt.annotate('Ter', xy=(center_x,0), xytext=(center_x,0),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=8)
    plt.annotate('Melting energy '+name, xy=(center_x,radius), xytext=(center_x,center_y),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=7)
    plt.savefig(basedir+"data/"+name+"/annotation/melting_energy.pdf", transparent=False)        
    plt.close()


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

    for cond in gen.genes_valid.keys():
        print 'Computing condition',cond
        gen_states = compute_state_genes(gen,cond)
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
            else:
                print "Unknown method, computing default..."
                zscores = compute_zscores_actrepr(tot_act,tot_repr,tot_non,bins)

        except:
            print 'Invalid data (e.g. no genes repressed nor activated)'
            sys.exit()

        # init plots
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.subplot(1,1,1)
        plt.axis([0, fig_width, 0, fig_height]) ; ax.set_axis_off() 
       
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

        else:
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
            # display colour in the middle of the windows
            start_angle = angle * (windows/2 - increment/2)  / increment
            i=0
            for value in zscores:
                # edgecolor = assign a colour depending on value using cMap
                # draw arc on circle
                arc = patches.Arc((center_x,center_y), radius, radius, angle=90,theta1=-(start_angle+i+angle), theta2=-(start_angle+i), edgecolor=cMap_fc.to_rgba(value), lw=7)
                ax.add_patch(arc)
                i+= angle

            cMap_fc._A = [] # fake array to print map
            cbar = plt.colorbar(cMap_fc,shrink=0.3)
            cbar.set_label("Z-score",size=9)
            cbar.ax.tick_params(labelsize=7)

        plt.axis("equal")

        plt.annotate('OriC', xy=(center_x,radius), xytext=(center_x,radius+0.4),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=8)
        plt.annotate('Ter', xy=(center_x,0), xytext=(center_x,0),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=8)
        plt.annotate(cond, xy=(center_x,radius), xytext=(center_x,center_y),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=7)
        form = kwargs.get('format','svg')
        if form == "png":
            plt.savefig(path+"/circle-"+cond+".png", dpi=400, transparent=False)
        elif form == "svg":
            plt.savefig(path+"/circle-"+cond+".svg", transparent=False) 
        elif form == "pdf":
            plt.savefig(path+"/circle-"+cond+".pdf", transparent=False)
        else:
            print 'Unknown format, computing default...'
            plt.savefig(path+"/circle-"+cond+".svg",  transparent=False)
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
            zscore =(nb_act - (nb_act+nb_repr)*p_exp) / (sqrt((nb_act+nb_repr)*p_exp*(1-p_exp)))
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
            zscore =(nb_act - (nb_tot*p_exp)) / sqrt(nb_tot*p_exp*(1-p_exp))
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
            zscore =(nb_repr - (nb_tot*p_exp)) / sqrt(nb_tot*p_exp*(1-p_exp))
        except: # division by zero if no genes activated nor repressed
            zscore = 0
        zscores.append(zscore)
    return zscores
