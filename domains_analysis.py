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
from draw_circles_ME_FC import *
from statsmodels.stats.proportion import proportions_ztest
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
    pathdb = '{}data/{}/domains'.format(basedir,gen.name)
    if not os.path.exists(pathdb):
      os.makedirs(pathdb)

    fig_width = 13.5 ; fig_height = 5
    fig = plt.figure(figsize=(fig_width,fig_height))
    #plt.tick_params(labelsize=16)
    sns.set(font_scale=1.2)
    ax = sns.heatmap(df,linewidths=.5,annot=True, vmin=vmin, vmax=vmax, cmap=colormap)
    plt.yticks(rotation=0)
    ax.tick_params(labelsize=20)
    plt.tight_layout()
    plt.savefig(pathdb+"/domains_zscores.png",dpi=200, transparent=False)
    

def correlation_matrix(gen,*arg, **kwargs):
    '''
    Draws correlation matrix for all genes starting from FC (/! not expression)
    '''
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
    dfs = [] # In order to generate comprehensive images of correlation matrix, divides it into smaller matrix
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

        # load genes in domains    
        for d in domains.keys(): # if non overlapping domain
            if domains[d]['coord'][0] < domains[d]['coord'][1]: # if gene in domain
                if domains[d]['coord'][0] <= g.start and domains[d]['coord'][1] >= g.start:
                    domains[d]['genes'].append(gene)
            else: # if overlapping domain
                if domains[d]['coord'][0] <= g.start or domains[d]['coord'][1] >= g.start:
                    domains[d]['genes'].append(gene)

    # for each condition, load act, rep and non genes in each domain
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
            for d in domains.keys(): # if non overlapping domain
                if domains[d]['coord'][0] < domains[d]['coord'][1]: # if gene in domain
                    if domains[d]['coord'][0] <= start and domains[d]['coord'][1] >= start:
                        if state == 1: # if gene act
                            domains[d][cond]['act'].append(conv[start])
                        elif state == -1: # if gene rep
                            domains[d][cond]['rep'].append(conv[start])
                        elif state ==0: # if non
                            domains[d][cond]['non'].append(conv[start])

                else: # if overlapping domain
                    if start <= st or end >= st:
                        if domains[d]['coord'][0] <= start or domains[d]['coord'][1] >= start:
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
    pvals = {}
    Zthresh = kwargs.get('Z',1.9)
    meth = kwargs.get('meth','actrep')

    if meth == 'all':       
        for d in domains.keys(): 
            res[d] = {}    
            labels[d] = {}
            pvals[d] = {}
            print d
            for cond1 in conds:
                if abs(domains[d][cond1]['Z']) >= Zthresh:
                    res[d][cond1] = {}
                    labels[d][cond1] = {}
                    pvals[d][cond1] = {}
                    for cond2 in conds:
                        if abs(domains[d][cond2]['Z']) >= Zthresh:
                            l1 = domains[d][cond1]['act'] + domains[d][cond1]['rep']
                            l2 = domains[d][cond2]['act'] + domains[d][cond2]['rep']
                            tot = min(len(l1),len(l2))
                            intersect = len(list(set(l1) & set(l2)))
                            res[d][cond1][cond2] = int((float(intersect) / float(tot))*100)
                            N = len(domains[d]['genes']) ; n = min(len(l1),len(l2)) ; pA = max(len(l1),len(l2)) ; k = intersect - 1
                            rv = stats.hypergeom(N,pA,n) ; pval = 1 - rv.cdf(k)
                            pvals[d][cond1][cond2] = pval

                            #print 'N',N,'n',n,'pA',pA,'k',k,'pval',pval
                            s = ''
                            if pval == 0:
                                s = '***'
                            if pval <= 0.001:
                                s = '***'
                            elif pval <= 0.01:
                                s = '**' 
                            elif pval <= 0.05:
                                s = '*' 
                            else:
                                s = 'ns'
                            
                            labels[d][cond1][cond2] = '{}% {}'.format(str(res[d][cond1][cond2]),s)

    elif meth == 'actrep':       
        for d in domains.keys(): 
            res[d] = {}    
            labels[d] = {}
            pvals[d] = {}
            print d
            for cond1 in conds:
                if abs(domains[d][cond1]['Z']) >= Zthresh:
                    res[d][cond1] = {}
                    labels[d][cond1] = {}
                    pvals[d][cond1] = {}
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
                            res[d][cond1][cond2] = int((float(intersect) / float(tot))*100)

                            N = len(domains[d]['genes']) ; n = min(len(l1),len(l2)) ; pA = max(len(l1),len(l2)) ; k = intersect - 1
                            rv = stats.hypergeom(N,pA,n) ; pval = 1 - rv.cdf(k)
                            pvals[d][cond1][cond2] = pval
                            #print 'N',N,'n',n,'pA',pA,'k',k,'pval',pval
                            s = ''
                            if pval <= 0.001:
                                s = '***'
                            elif pval <= 0.01:
                                s = '**' 
                            elif pval <= 0.05:
                                s = '*' 
                            else:
                                s = 'ns'
                            
                            #labels[d][cond1][cond2] = '{}={}_{}-{}{}'.format(str(int(res[d][cond1][cond2])),str(intersect),str(len(l1)),str(len(l2)),s)
                            labels[d][cond1][cond2] = '{}% {}'.format(str(res[d][cond1][cond2]),s)


    for d in domains.keys():
        df = pd.DataFrame.from_dict(pvals[d],dtype=float)
        lab = pd.DataFrame.from_dict(labels[d])
        vmin = kwargs.get('vmin', 0)
        vmax = kwargs.get('vmax', 0.1)
        # Colormap for Zscore
        colormap= kwargs.get('colormap','Greens_r') # default value

        pathdb = '{}data/{}/domains'.format(basedir,gen.name)
        if not os.path.exists(pathdb):
          os.makedirs(pathdb)

        fig_width = 12 ; fig_height = fig_width/1.25
        fig = plt.figure(figsize=(fig_width,fig_height))

        mask = np.zeros_like(df)
        mask[np.triu_indices_from(mask)] = True
        with sns.axes_style("dark"):
            ax = sns.heatmap(df,linewidths=.5,annot=lab, fmt="", vmin=vmin, vmax=vmax, cmap=colormap, mask=mask, xticklabels=False)
        sns.set(font_scale=1.4)

        plt.yticks(rotation=0)
        plt.xticks(rotation=90)
        ax.tick_params(labelsize=20)
        #ax = sns.heatmap(df,linewidths=.5,annot=lab, fmt="", vmin=vmin, vmax=vmax, cmap=colormap)
        #plt.tick_params(labelsize=16)
        plt.tight_layout()
        plt.savefig(pathdb+"/{}_genes.png".format(str(d)), dpi=200, transparent=False)
        plt.close('all')
    

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


def domains_features(gen,*arg, **kwargs):
    gen.load_gene_orientation() # compute gene orientation
    gen.load_seq()
    prom = kwargs.get('prom',200)
    domains = {}
    domains['genome'] = {'genes':[], 'orientation':{'convergent':0,'divergent':0,'tandem':0,'isolated':0}, 'leading strand':{'forward':0,'reverse':0}, 'melting energy':{'prom':[],'all':[]}}
    with open(basedir+"data/"+gen.name+"/domains.txt","r") as f:
        header = next(f)
        for line in f:
            line=line.strip()
            line=line.split('\t')
            domains[line[0]] = {}
            domains[line[0]]['coord'] = [int(line[1]),int(line[2])]
            domains[line[0]]['genes'] = []
            domains[line[0]]['orientation'] = {'convergent':0,'divergent':0,'tandem':0,'isolated':0}
            domains[line[0]]['leading strand'] = {'forward':0,'reverse':0}
            domains[line[0]]['melting energy'] = {'prom':[],'all':[]}
        f.close()

    for gene in gen.genes.keys():
        try:            
            g = gen.genes[gene]
            if g.strand:        
                domains['genome']['leading strand']['forward'] += 1
                seq_prom = gen.seq[g.start-prom:g.start]
                seq = gen.seq[g.start:g.end]
            elif not g.strand:
                domains['genome']['leading strand']['reverse'] += 1
                seq_prom = gen.seqcompl[g.start:g.start+prom]
                seq = gen.seq[g.end:g.start]

            me_prom = mt.Tm_GC(Seq(seq_prom), strict=False)
            me = mt.Tm_GC(Seq(seq), strict=False)
            domains['genome']['melting energy']['prom'].append(me_prom)
            domains['genome']['melting energy']['all'].append(me)
            domains['genome']['orientation'][g.orientation] += 1
            belong = False
            for d in domains.keys():
                if d != 'genome':
                    if domains[d]['coord'][0] < domains[d]['coord'][1]:
                        if domains[d]['coord'][0] <= g.start and domains[d]['coord'][1] >= g.start:
                            domains[d]['genes'].append(gene)
                            belong = True
                    else:
                        if domains[d]['coord'][0] <= g.start or domains[d]['coord'][1] >= g.start:
                            domains[d]['genes'].append(gene)
                            belong = True
                    if belong:
                        if g.strand:        
                            domains[d]['leading strand']['forward'] += 1
                        elif not g.strand:
                            domains[d]['leading strand']['reverse'] += 1
                        domains[d]['melting energy']['prom'].append(me_prom)
                        domains[d]['melting energy']['all'].append(me)
                        domains[d]['orientation'][g.orientation] += 1
        except Exception as e:
            print e
            pass
    res = {}
    for d in domains.keys():
        res[d] = {}
        res[d]['genes'] = domains[d]['genes']
        sobs = domains[d]['leading strand']['forward']
        sexp = domains['genome']['leading strand']['forward']
        nobs = domains[d]['leading strand']['forward'] + domains[d]['leading strand']['reverse']
        nexp = domains['genome']['leading strand']['forward'] + domains['genome']['leading strand']['reverse']
        pobs = int((float(sobs)/nobs)*100)
        stat, pval = proportions_ztest([sobs,sexp], [nobs,nexp])
        res[d]['p+'] = '{}% {}'.format(str(pobs),significance(pval))

        sobs = domains[d]['orientation']['convergent']
        sexp = domains['genome']['orientation']['convergent']
        nobs = domains[d]['orientation']['convergent'] + domains[d]['orientation']['divergent'] + domains[d]['orientation']['tandem']
        nexp = domains['genome']['orientation']['convergent'] + domains['genome']['orientation']['divergent'] + domains['genome']['orientation']['tandem']
        pobs = int((float(sobs)/nobs)*100)
        stat, pval = proportions_ztest([sobs,sexp], [nobs,nexp])
        res[d]['pconv'] = '{}% {}'.format(str(pobs),significance(pval))

        sobs = domains[d]['orientation']['divergent']
        sexp = domains['genome']['orientation']['divergent']
        nobs = domains[d]['orientation']['convergent'] + domains[d]['orientation']['divergent'] + domains[d]['orientation']['tandem']
        nexp = domains['genome']['orientation']['convergent'] + domains['genome']['orientation']['divergent'] + domains['genome']['orientation']['tandem']
        pobs = int((float(sobs)/nobs)*100)
        stat, pval = proportions_ztest([sobs,sexp], [nobs,nexp])
        res[d]['pdiv'] = '{}% {}'.format(str(pobs),significance(pval))

        pval = stats.ttest_ind(domains[d]['melting energy']['all'],domains['genome']['melting energy']['all'],equal_var=False)[1]
        res[d]['p_me_all'] = significance(pval)

        pval = stats.ttest_ind(domains[d]['melting energy']['prom'],domains['genome']['melting energy']['prom'],equal_var=False)[1]
        res[d]['p_me_prom'] = significance(pval)
            # res[d] = {}
            # res[d]['p+'] = float(domains[d]['leading strand']['forward']) / (float(domains[d]['leading strand']['forward']+domains[d]['leading strand']['reverse']))
            # res[d]['pconv'] = float(domains[d]['orientation']['convergent']) / (float(domains[d]['orientation']['convergent']+domains[d]['orientation']['divergent']+domains[d]['orientation']['tandem']))
            # res[d]['pdiv'] = float(domains[d]['orientation']['divergent']) / (float(domains[d]['orientation']['convergent']+domains[d]['orientation']['divergent']+domains[d]['orientation']['tandem']))
            # res[d]['me_all'] = np.mean(domains[d]['melting energy']['all'])
            # res[d]['me_prom'] = np.mean(domains[d]['melting energy']['prom'])

    # res['genome'] = {}
    # res['genome']['p+'] = float(domains['genome']['leading strand']['forward']) / (float(domains['genome']['leading strand']['forward']+domains['genome']['leading strand']['reverse']))
    # res['genome']['pconv'] = float(domains['genome']['orientation']['convergent']) / (float(domains['genome']['orientation']['convergent']+domains['genome']['orientation']['divergent']+domains['genome']['orientation']['tandem']))
    # res['genome']['pdiv'] = float(domains['genome']['orientation']['divergent']) / (float(domains['genome']['orientation']['convergent']+domains['genome']['orientation']['divergent']+domains['genome']['orientation']['tandem']))
    # res['genome']['me_all'] = np.mean(domains['genome']['melting energy']['all'])
    # res['genome']['me_prom'] = np.mean(domains['genome']['melting energy']['prom'])

    pathdb = '{}data/{}/domains'.format(basedir,gen.name)
    if not os.path.exists(pathdb):
      os.makedirs(pathdb)

    df = pd.DataFrame.from_dict(res,orient='index')
    df.to_csv(pathdb+'/domains_features.csv')

