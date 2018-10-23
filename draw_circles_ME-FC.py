#! /usr/bin/env python
# -*- coding: utf-8 -*-
from Bio.SeqUtils import MeltingTemp as mt
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

from scipy import stats

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
    cbar.set_label("Melting temp C",size=9)
    cbar.ax.tick_params(labelsize=7)
    plt.axis("equal")
    plt.annotate('OriC', xy=(center_x,radius), xytext=(center_x,radius+0.1),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=7)
    plt.annotate('Ter', xy=(center_x,0), xytext=(center_x,-0.3),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=7)
    plt.annotate('Melting energy '+name, xy=(center_x,radius), xytext=(center_x,center_y),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=7)
    plt.savefig(basedir+"data/"+name+"/annotation/melting_energy.svg", transparent=False)        
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
        gen_states = compute_state_genes_operon_correction(gen,cond)
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
        ax = plt.subplot(1,1,1)

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

            cMap_fc._A = [] # fake array to print map
            cbar = plt.colorbar(cMap_fc,shrink=0.3)
            cbar.set_label("Z-score",size=7)
            cbar.ax.tick_params(labelsize=5)

        plt.axis("equal")

        plt.annotate('OriC', xy=(center_x,radius), xytext=(center_x,radius+0.1),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=7)
        plt.annotate('Ter', xy=(center_x,0), xytext=(center_x,-0.3),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=7)
        plt.annotate(cond, xy=(center_x,radius), xytext=(center_x,center_y),verticalalignment = 'center', horizontalalignment = 'center',fontstyle='italic', wrap=True, fontsize=7)
        form = kwargs.get('format','svg')
        if form == "png":
            plt.savefig(path+"/circle-"+cond.replace("/","-")+".png", dpi=400, transparent=False)
        elif form == "svg":
            plt.savefig(path+"/circle-"+cond.replace("/","-")+".svg", transparent=False) 
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
                    gen_states.append([operon.start,-1])                
                elif np.mean(expr) >= fc_treshold_pos:
                    gen_states.append([operon.start,1])                
                else:
                    if expr_none == []:
                        gen_states.append([operon.start,0])
                nb_operon += 1

            elif expr_none != []:
                gen_states.append([operon.start,0])
                nb_operon += 1
        except:
            pass


    for gene in gen.genes_valid[cond]:
        if gene not in not2do:
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
    return zscores_repr

def compute_domains_scores(gen,*arg, **kwargs):
    '''
    Starting from domain coordinates on genome in domains.txt and FC values of genes for several conditions,
    computes zscore of act / rep genes in each domain for each condition and draws a summary table.
    '''
    if not hasattr(gen, 'genes_valid'): # if no fc loaded 
        try:
            print 'Trying to load FC...'
            gen.load_fc_pval()
            print 'FC loaded'
        except:
            print 'Unable to load FC'
            sys.exit()
    conds = []
    domains = {}
    try:
        with open(basedir+"data/"+gen.name+"/domains.txt","r") as f:
            header = next(f)
            for line in f:
                line=line.strip()
                line=line.split('\t')
                domains[line[0]] = {}
                domains[line[0]]['coord'] = [int(line[1]),int(line[2])]
            f.close()
    except:
        print 'Error please load domains.txt'
        sys.exit()

    dnames = sorted(domains, key=lambda k: domains[k]['coord'])

    for cond in gen.genes_valid.keys():
        print 'Computing condition',cond
        conds.append(cond)

        bins = [] # bins = domains of the genome : [start coordinate,end coordinate,nb of activated genes,nb of repressed genes,nb of genes not affected]
        bins_overlap = [] # bins where end coordinate < start coordinate (overlap circular chromosome)
        for coor in sorted([domains[d]['coord'] for d in domains]):
            if coor[0] < coor[1]:
                bins.append(coor+[0,0,0])
            else:
                bins_overlap.append(coor+[0,0,0])
        bins = np.array(bins) ; bins_overlap = np.array(bins_overlap)

        gen_states = compute_state_genes_operon_correction(gen,cond)  
        for start,state in gen_states: # reminder gene state : a row = beginning of gene, state (activated, repressed or not affected)
            if state == 1: # activated
            # test to which bins the gene belongs to, and add one to the nb of activated genes of these bins
                bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),2] += 1
                if bins_overlap != []:
                    bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),2] += 1
            elif state == -1: # repressed gene
                bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),3] += 1
                if bins_overlap != []:
                    bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),3] += 1
            elif state == 0: # not affected gene
                bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),4] += 1
                if bins_overlap != []:
                    bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),4] += 1
        
        if bins_overlap != []:
            bins = np.concatenate((bins,bins_overlap))

        print 'Domains on genome :\n',bins
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
        
        print zscores


        i=0
        for start,end,nb_act,nb_repr,nb_null in bins:
            for d in domains.keys():
                if domains[d]['coord'] == [start,end]:
                    domains[d][cond] = float(np.round(zscores[i],3))
                    i += 1


    vmin = kwargs.get('vmin', -4)
    vmax = kwargs.get('vmax', +4)
    # Colormap for Zscore
    colormap= kwargs.get('colormap','jet') # default value


    df = pd.DataFrame.from_dict(domains,dtype=float)
    df.drop('coord',0,inplace=True) ; df.sort_index(axis=1,inplace=True) 
    df.fillna(0.0,inplace=True)

    #col =  cMap_fc.to_rgba(np.array(df.values,dtype=float))

    fig_width = 13.5 ; fig_height = 5
    fig = plt.figure(figsize=(fig_width,fig_height))
    #plt.tick_params(labelsize=16)
    plt.tight_layout()
    ax = sns.heatmap(df,linewidths=.5,annot=True, vmin=vmin, vmax=vmax, cmap=colormap)

    plt.savefig("/home/raphael/Documents/test.svg", transparent=False)
    

def correlation_matrix(gen,*arg, **kwargs):
    if not hasattr(gen, 'genes_valid'): # if no fc loaded 
        try:
            print 'Trying to load FC...'
            gen.load_fc_pval()
            print 'FC loaded'
        except:
            print 'Unable to load FC'
            sys.exit()

    gene_dict = {}
    for gene in gen.genes.keys():
        try:
            gene_dict[gene] = gen.genes[gene].fc_pval

        except:
            pass

    df = pd.DataFrame.from_dict(gene_dict,orient='index', dtype=float)
    df.dropna(inplace=True) ; df.sort_index(axis=1,inplace=True)
    df = df.applymap(lambda x: x[0])

    new_index = {}
    for g in list(df.index):
        new_index[g] = gen.genes[g].start/1000

    df.rename(index=new_index, inplace=True) ; df.sort_index(inplace=True)

    size= kwargs.get('size',250) # default value

    i = 0
    dfs = [] ;
    for i in range(1,df.shape[0],size): # create bins depending on windows size and increment value
        if (i+size) <= df.shape[0]: # enough length to create a bin
            dfs.append(df[i:i+size])
        else: # i + windows > genome size, overlap with the beginning (circular chromosome)
            dfs.append(df[i:])

    colormap= kwargs.get('colormap','jet') # default value

    i=0
    for d in dfs:
        name = "corr_{}_{}.png".format(str(size),str(i))
        fig = plt.figure()
        ax = sns.heatmap(d.T.corr(),linewidths=0,annot=False, vmin=-1, vmax=+1, cmap=colormap)
        plt.savefig("/home/raphael/Documents/"+name, transparent=False)
        plt.close("all")
        i+=1



def domains_common_genes(gen,*arg, **kwargs):
    '''
    '''
    if not hasattr(gen, 'genes_valid'): # Load FC
        try:
            print 'Trying to load FC...'
            gen.load_fc_pval()
            print 'FC loaded'
        except:
            print 'Unable to load FC'
            sys.exit()

    conds = [] # FC conds, useful for later
    domains = {} # domains ; shape : {d1 : {'coordinates':[start,end], 'genes':[g1,g2,g3,g4,g5,g6,g7], 'cond1':{'act':[g1,g2],'rep':[g3,g4],'non affected':[g5,g6,g7]}, 'cond2':{...}}, d2:{...}}
    
    try:
        with open(basedir+"data/"+gen.name+"/domains.txt","r") as f:
            header = next(f)
            for line in f:
                line=line.strip()
                line=line.split('\t')
                domains[line[0]] = {} # FC conditions, act rep and non genes per cond
                domains[line[0]]['coord'] = [int(line[1]),int(line[2])] # coordinates
                domains[line[0]]['genes'] = [] # gene names in this domain
            f.close()
    except:
        print 'Error please load domains.txt'
        sys.exit()

    conv = {} # dict of shape {gene_start:'gene_name'}
    # Load genes in domains
    for gene in gen.genes.keys():
        g = gen.genes[gene]
        try:
            conv[g.start] = gene
        except:
            pass
            
        for d in domains.keys():
            if domains[d]['coord'][0] < domains[d]['coord'][1]:
                if domains[d]['coord'][0] <= g.start and domains[d]['coord'][1] >= g.start:
                    domains[d]['genes'].append(gene)
            else:
                if domains[d]['coord'][0] <= g.start or domains[d]['coord'][1] >= g.start:
                    domains[d]['genes'].append(gene)

    # For each condition, load act, rep and non genes in each domain
    for cond in gen.genes_valid.keys():
        print 'Computing condition',cond
        conds.append(cond)
        
        for d in domains.keys():
            domains[d][cond] = {}
            domains[d][cond]['act'] = []
            domains[d][cond]['rep'] = []
            domains[d][cond]['non'] = []

        gen_states = compute_state_genes(gen,cond)

        for start,state in gen_states:
            for d in domains.keys():
                if domains[d]['coord'][0] < domains[d]['coord'][1]:
                    if domains[d]['coord'][0] <= start and domains[d]['coord'][1] >= start:
                        if state == 1:
                            domains[d][cond]['act'].append(conv[start])
                        elif state == -1:
                            domains[d][cond]['rep'].append(conv[start])
                        elif state ==0:
                            domains[d][cond]['non'].append(conv[start])

                else:
                    if start <= st or end >= st:
                        if state == 1:
                            domains[d][cond]['act'].append(conv[start])
                        elif state == -1:
                            domains[d][cond]['rep'].append(conv[start])
                        elif state ==0:
                            domains[d][cond]['non'].append(conv[start])
        

        bins = [] # bins = domains of the genome : [start coordinate,end coordinate,nb of activated genes,nb of repressed genes,nb of genes not affected]
        bins_overlap = [] # bins where end coordinate < start coordinate (overlap circular chromosome)
        for coor in sorted([domains[d]['coord'] for d in domains]):
            if coor[0] < coor[1]:
                bins.append(coor+[0,0,0])
            else:
                bins_overlap.append(coor+[0,0,0])
        bins = np.array(bins) ; bins_overlap = np.array(bins_overlap)

        gen_states = compute_state_genes_operon_correction(gen,cond)  
        for start,state in gen_states: # reminder gene state : a row = beginning of gene, state (activated, repressed or not affected)
            if state == 1: # activated
            # test to which bins the gene belongs to, and add one to the nb of activated genes of these bins
                bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),2] += 1
                if bins_overlap != []:
                    bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),2] += 1
            elif state == -1: # repressed gene
                bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),3] += 1
                if bins_overlap != []:
                    bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),3] += 1
            elif state == 0: # not affected gene
                bins[np.where((bins[:,0] <= start) & (bins[:,1] > start)),4] += 1
                if bins_overlap != []:
                    bins_overlap[np.where((bins_overlap[:,0] <= start) | (bins_overlap[:,1] > start)),4] += 1
        
        if bins_overlap != []:
            bins = np.concatenate((bins,bins_overlap))

        print 'Domains on genome :\n',bins

        tot_act = len(gen_states[np.where(gen_states[:,1] == 1)])
        print 'Total activated genes on genome :',tot_act
        tot_repr = len(gen_states[np.where(gen_states[:,1] == -1)])
        print 'Total repressed genes on genome :',tot_repr
        tot_non = len(gen_states[np.where(gen_states[:,1] == 0)])
        print 'Total non affected genes on genome :',tot_non

        zscores = compute_zscores_actrepr(tot_act,tot_repr,tot_non,bins)        

        i=0
        for start,end,nb_act,nb_repr,nb_null in bins:
            for d in domains.keys():
                if domains[d]['coord'] == [start,end]:
                    domains[d][cond]['Z'] = float(np.round(zscores[i],3))
                    i += 1
    
    res = {}
    labels = {}
    Zthresh = kwargs.get('Z',2)
    meth = kwargs.get('meth','all')


    if meth == 'all':       
        for d in domains.keys(): 
            res[d] = {}    
            labels[d] = {}
            for cond1 in conds:
                if abs(domains[d][cond1]['Z']) >= Zthresh:
                    res[d][cond1] = {}
                    labels[d][cond1] = {}
                    for cond2 in conds:
                        if abs(domains[d][cond2]['Z']) >= Zthresh:
                            l1 = domains[d][cond1]['act'] + domains[d][cond1]['rep']
                            l2 = domains[d][cond2]['act'] + domains[d][cond2]['rep']
                            tot = min(len(l1),len(l2))
                            intersect = len(list(set(l1) & set(l2)))
                            res[d][cond1][cond2] = (float(intersect) / float(tot))*100

                            N = len(domains[d]['genes']) ; n = min(len(l1),len(l2)) ; pA = max(len(l1),len(l2)) ; k = intersect - 1
                            rv = stats.hypergeom(N,pA,n) ; pval = 1 - rv.cdf(k)
                            print 'N',N,'n',n,'pA',pA,'k',k,'pval',pval

                            if pval <= 0.001:
                                s = '***'
                            elif pval <= 0.01:
                                s = '**' 
                            elif pval <= 0.05:
                                s = '*' 
                            else:
                                s = ''
                            
                            labels[d][cond1][cond2] = '{}={}_{}-{}{}'.format(str(int(res[d][cond1][cond2])),str(intersect),str(len(l1)),str(len(l2)),s)

    elif meth == 'actrep':       
        for d in domains.keys(): 
            res[d] = {}    
            labels[d] = {}
            for cond1 in conds:
                if abs(domains[d][cond1]['Z']) >= Zthresh:
                    res[d][cond1] = {}
                    labels[d][cond1] = {}
                    s1 = np.sign(domains[d][cond1]['Z'])
                    for cond2 in conds:
                        if abs(domains[d][cond2]['Z']) >= Zthresh:
                            s2 = np.sign(domains[d][cond2]['Z'])
                            if s1 < 0:
                                l1 = domains[d][cond1]['rep']
                            else:
                                l1 = domains[d][cond1]['act']
                            if s2 < 0:
                                l2 = domains[d][cond2]['rep']
                            else:
                                l2 = domains[d][cond2]['act']

                            tot = min(len(l1),len(l2))
                            intersect = len(list(set(l1) & set(l2)))
                            res[d][cond1][cond2] = (float(intersect) / float(tot))*100

                            N = len(domains[d]['genes']) ; n = min(len(l1),len(l2)) ; pA = max(len(l1),len(l2)) ; k = intersect - 1
                            rv = stats.hypergeom(N,pA,n) ; pval = 1 - rv.cdf(k)
                            print 'N',N,'n',n,'pA',pA,'k',k,'pval',pval

                            if pval <= 0.001:
                                s = '***'
                            elif pval <= 0.01:
                                s = '**' 
                            elif pval <= 0.05:
                                s = '*' 
                            else:
                                s = ''
                            
                            labels[d][cond1][cond2] = '{}={}_{}-{}{}'.format(str(int(res[d][cond1][cond2])),str(intersect),str(len(l1)),str(len(l2)),s)


    for d in domains.keys():
        df = pd.DataFrame.from_dict(res[d],dtype=float)
        lab = pd.DataFrame.from_dict(labels[d])
        vmin = kwargs.get('vmin', 0)
        vmax = kwargs.get('vmax', 100)
        # Colormap for Zscore
        colormap= kwargs.get('colormap','binary') # default value

        fig_width = 16 ; fig_height = 12
        fig = plt.figure(figsize=(fig_width,fig_height))

        mask = np.zeros_like(df)
        mask[np.triu_indices_from(mask)] = True
        with sns.axes_style("white"):
            ax = sns.heatmap(df,linewidths=.5,annot=lab, fmt="", vmin=vmin, vmax=vmax, cmap=colormap, mask=mask)
        #ax = sns.heatmap(df,linewidths=.5,annot=lab, fmt="", vmin=vmin, vmax=vmax, cmap=colormap)
        #plt.tick_params(labelsize=16)
        plt.tight_layout()
        plt.savefig("/home/raphael/Documents/{}.svg".format(str(d)), transparent=False)
    


def table_genes(gen):

    res = {}
    if not hasattr(gen, 'genes_valid'): # if no fc loaded 
        gen.load_fc_pval()

    for cond in gen.genes_valid.keys():
        for gene in gen.genes_valid[cond]:
            g = gen.genes[gene]
            if g.start not in res.keys():
                res[g.start] = {}

            # if activated
            if gen.genes[gene].fc_pval[cond][0] >= fc_treshold_pos and gen.genes[gene].fc_pval[cond][1] <= pval_treshold:
                res[g.start][cond] = 1
            # if repressed
            elif gen.genes[gene].fc_pval[cond][0] <= fc_treshold_neg and gen.genes[gene].fc_pval[cond][1] <= pval_treshold:
                res[g.start][cond] = -1
            # if not affected
            else:
                res[g.start][cond] = 0
    
    
    df = pd.DataFrame.from_dict(res,orient='index')
    df.dropna(inplace=True) ; df.sort_index(axis=1,inplace=True)

    new_index = {}
    for g in list(df.index):
        new_index[g] = g/1000

    df.rename(index=new_index, inplace=True) ; df.sort_index(axis=1,inplace=True)
    df.to_csv('/home/raphael/Documents/test.csv')


def genes_domains(gen):
    gen.load_annotation()
    gen.load_gene_orientation() # compute gene orientation
    domains = {}
    with open(basedir+"data/"+gen.name+"/domains.txt","r") as f:
        header = next(f)
        for line in f:
            line=line.strip()
            line=line.split('\t')
            domains[line[0]] = {}
            domains[line[0]]['coord'] = [int(line[1]),int(line[2])]
            domains[line[0]]['genes'] = []
            domains[line[0]]['orientation'] = {'convergent':0,'divergent':0,'tandem':0,'isolated':0,}
            domains[line[0]]['leading strand'] = {'forward':0,'reverse':0}
            domains[line[0]]['melting energy'] = []

        f.close()

    for gene in gen.genes.keys():
        g = gen.genes[gene]
        belong = False
        for d in domains.keys():
            if domains[d]['coord'][0] < domains[d]['coord'][1]:
                if domains[d]['coord'][0] <= g.start and domains[d]['coord'][1] >= g.start:
                    domains[d]['genes'].append(gene)
                    belong = True
            else:
                if start <= st or end >= st:
                    domains[d]['genes'].append(gene)
                    belong = True
            if belong:
                domains[d]['orientation'][g.orientation] += 1
                #domains[d]['melting energy'].append()

                if g.strand:        
                    domains[d]['leading strand']['forward'] += 1
                else:
                    domains[d]['leading strand']['reverse'] += 1
                

    df = pd.DataFrame.from_dict(domains,orient='index')
    df.to_csv('/home/raphael/Documents/test.csv')
