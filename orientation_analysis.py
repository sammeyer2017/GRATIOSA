from globvar import *
from scipy import stats
#from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats import proportion
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from useful_functions import *
plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})
plt.rcParams.update({'font.size': 11})
plt.rcParams.update({'legend.fontsize': 9})
plt.rcParams.update({'font.family': "Arial"})

def compute_orientation_proportion(gen,cond_fc,*args,**kwargs):
    '''
    Compute gene orientation and state (act, non, rep) depending on FC data, 
    to check whether or not act / rep / non genes are more convergent / divergent / tandem.
    '''

    bound = kwargs.get('bound',5000) # maximal distance for seeking neighbour, either left or right
    couple = kwargs.get('couple',3)
    thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
    thresh_fc = kwargs.get('thresh_fc', 0) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated
    filtering = kwargs.get('filtering', False)
    name_filt = kwargs.get('name_filt', None)

    pathdb = '{}data/{}/orientation'.format(basedir,gen.name)
    if not os.path.exists(pathdb):
      os.makedirs(pathdb)

    titl = pathdb+'/{}_{}_FC{}_PVAL{}'.format(gen.name,cond_fc.replace('/',''),thresh_fc,thresh_pval)

    gen.compute_state_from_fc(thresh_pval=thresh_pval,thresh_fc=thresh_fc) 
    gen.load_gene_orientation(bound=bound,couple=couple) # compute gene orientation
    # number of convergent, divergent and tandem genes for each state 
    states = ['act','non','rep']
    orientations = ['convergent','tandem','divergent']
    res = {}
    for state in states:
        res[state] = {}
        for orientation in orientations:
            res[state][orientation] = []

    for gene in gen.genes:
        try:               
            g = gen.genes[gene]
            if filtering:
                if g.expression['baseMean']/g.length > 18:
                    res[g.state[cond_fc]][g.orientation].append(gene)
                  
                # if g.length > 750 and g.name[0:len(name_filt)] != name_filt:
                #     res[g.state[cond_fc]][g.orientation].append(gene)

            else:
                res[g.state[cond_fc]][g.orientation].append(gene)
        except:
            pass
    
    for state in states:
        for orientation in orientations:
            try:
                print '{}: {}'.format(orientation, str(len(res[state][orientation])))
            except:
                pass

    means = [] ; stds = [] ; cim = []
    meth = kwargs.get('meth','sam')
    # meth = method to use to display results, e.g. (convergent) / (convergent + divergent) for each state
    if couple == 2: 
        meth = 'couple'
            
    if meth == 'couple':
        for state in ['act','non','rep']:
            try:
                pexp = float(len(res[state]['divergent'])) / (len(res[state]['divergent'])+len(res[state]['tandem']))
                tot = len(res[state]['divergent'])+len(res[state]['tandem'])
                std = np.array(stats.binom.std(tot, pexp, loc=0))/tot
                means.append(pexp)
                stds.append(std)
            except:
                pass

    if meth == 'sam':
        cis = []
        results = open(titl+'.txt','w')
        results.write('Orientation\tAct+Rep\tAct\tRep\tNon\n')
        for orientation in orientations:
            try:
                pexp = float(len(res['act'][orientation])) / (len(res['act'][orientation])+len(res['rep'][orientation]))
                tot = len(res['act'][orientation])+len(res['rep'][orientation])
                std = np.array(stats.binom.std(tot, pexp, loc=0))/tot
                ci = proportion.proportion_confint(len(res['act'][orientation]), tot, alpha=0.05, method='normal')
                print "CI",ci
                cim.append((ci[1] - ci[0])/2)

                means.append(pexp)
                stds.append(std)
                results.write('{}\t{}\t{}\t{}\t{}\n'.format(orientation,str(len(res['act'][orientation])+len(res['rep'][orientation])), str(len(res['act'][orientation])), str(len(res['rep'][orientation])), str(len(res['non'][orientation]))))
            except Exception as e:
                print e
                pass
        
        s1 = len(res['act']['convergent'])
        s2 = len(res['act']['divergent'])
        n1 = len(res['act']['convergent']) + len(res['rep']['convergent'])
        n2 = len(res['act']['divergent']) + len(res['rep']['divergent'])

        stat1, pval1 = proportion.proportions_ztest([s1,s2], [n1,n2], alternative='larger') 
        stat2, pval2 = proportion.proportions_ztest([s1,s2], [n1,n2], alternative='smaller')
        pval = min(pval1,pval2)
        s = significance(pval)

        s1 = significance(pval1)
        s2 = significance(pval2)
        results.write('One-sided test\tp-value\n')        
        results.write('Larger\t{} {}\n'.format(s1,str(round(pval1,5))))
        results.write('Smaller\t{} {}'.format(s2,str(round(pval2,5))))        
        results.close()


    draw = kwargs.get("draw",True) # whether or not graphes have to be created
    org = kwargs.get('org', gen.name) # organism name to use for plot title

    if draw:
        fig_width = 2 ; fig_height = 2.2
        fig = plt.figure(figsize=(fig_width,fig_height))
        if meth == 'sam':
            xval = [1,2,3] ; labs = ['conv','tand','div']
            yval = means
            if s == 'ns':
                fs = 13
                dt = 0.01
            else:
                fs = 18
                dt = 0

            plt.bar(xval, yval, yerr = cim, width = 0.6, color = 'white',edgecolor = 'black', linewidth = 2, ecolor = 'black', capsize = 3.5)
            plt.xticks(xval,labs,rotation = 45)
            plt.xlim(0.5, 3.5)
            ylim1 = min(yval) - max(cim) - 0.01
            ylim2 = max(yval) + max(cim) + 0.15
            plt.ylim(ylim1,ylim2)
            barplot_annotate_brackets(0, 2, s,xval,yval,yerr=cim,fs=fs, dt=dt)#,dh=0.05, barh=.05, fs=14, dt=0)
            fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
            plt.ylabel('proportion of \nactivated genes')
            plt.xticks(xval,labs)
            #plt.title('{}, {}'.format(org,cond_fc),fontsize=11)
            #plt.title('{}, {}'.format(gen.name,cond_fc),fontsize=11)

            plt.tight_layout()
            plt.savefig(titl+'.svg',transparent=False)
            plt.close('all')

        else:

            try:
                xval = [0,1,2] ; labs = ['act','non','rep']
                yval = means
                plt.plot(xval, yval, 'rD', markersize=7, linestyle='None')
                plt.errorbar(xval, yval,yerr=stds,mec='black', capsize=10, elinewidth=1,mew=1,linestyle='None', color='black')        
            except:
                xval = [0,1] ; labs = ['act','rep']
                yval = means
                plt.plot(xval, yval, 'rD', markersize=7, linestyle='None')
                plt.errorbar(xval, yval,yerr=stds,mec='black', capsize=10, elinewidth=1,mew=1,linestyle='None', color='black')        
            
            fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
            if meth == 'couple':
                ax.set_ylabel('diverg / (tandem + diverg)',fontweight='bold')
            else:
                ax.set_ylabel('conv / (conv + diverg)',fontweight='bold')

            ax.set_xlim(xval[0]-1,xval[-1]+1)
            #ax.set_ylim(yval[0]-0.1,yval[-1]+0.05)
            plt.xticks(xval,labs,fontweight='bold')
            ax.set_title('{}, {}'.format(org,cond_fc),fontweight='bold')
            fig.set_size_inches(width, height)
            plt.tight_layout()
            plt.savefig(titl+'.svg',transparent=False)
            plt.close('all')


    write_genes = kwargs.get('write_genes',False)
    if write_genes:
        f = open(titl+'.txt',"w")
        f.write('convergent\tdivergent\ttandem\n')
        for state in ['act','non','rep']:
            f.write(state+'\n')
            f.write(",".join(str(x) for x in res[state]['convergent'])+'\t'+",".join(str(x) for x in res[state]['divergent'])+'\t'+",".join(str(x) for x in res[state]['tandem'])+'\n')
        f.close()

    write = kwargs.get('write',False)
    if write:
        results = open(titl+'.txt','w')
        results.write('Gene\tOrientation\tState\tLeft Neighbour state\tRight neighbour state\n')
        for gene in gen.genes:
            try:
                g = gen.genes[gene]
                lg = gen.genes[g.left_neighbour]
                rg = gen.genes[g.right_neighbour]
                results.write('{}-{}\t{}\t{}\t{}\t{}\n'.format(gene,g.name,g.orientation,g.state[cond_fc],lg.state[cond_fc],rg.state[cond_fc]))
            except:
                pass               
        results.close()


def compare_gene_orientation(self,c1,c2,*args,**kwargs):
    self.load_fc_pval()
    self.load_gene_orientation()
    thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
    thresh_fc = kwargs.get('thresh_fc', 0) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated

    res = {'convergent':{}, 'divergent':{}, 'tandem':{}}
    states = ['+*','+','-*','-']
    for orientation in res.keys():        
        for st1 in states:
            res[orientation][st1] = {}
            for st2 in states:
                res[orientation][st1][st2] = 0

    for gene in self.genes.keys():
        try:
            g = self.genes[gene]
            if g.fc_pval[c1][1] <= thresh_pval:
                stat1 = '*'
            else:
                stat1 = ''
            if g.fc_pval[c1][0] > thresh_fc:
                sign1 = '+'
            else:
                sign1 = '-'


            if g.fc_pval[c2][1] <= thresh_pval:
                stat2 = '*'
            else:
                stat2 = ''
            if g.fc_pval[c2][0] > thresh_fc:
                sign2 = '+'
            else:
                sign2 = '-'

            st1 = sign1+stat1 ; st2 = sign2+stat2
            res[g.orientation][st1][st2] += 1
        except:
            pass

    for orientation in res.keys():        
        df = pd.DataFrame.from_dict(res[orientation],orient='index')
        #df.drop(['non'],axis=1,inplace=True)
        df.sort_index(axis=1,inplace=True) ; df.sort_index(axis=0,inplace=True)
        df.to_csv('/home/raphael/Documents/test_{}.csv'.format(orientation))

def compute_orientation_dickeya(self,*args,**kwargs):
    bound = kwargs.get('bound',5000) # maximal distance for seeking neighbour, either left or right
    self.load_annotation()  
    self.load_gene_orientation(bound=bound)
    #files = ['act_puce_expo.txt','act_puce.txt','act_expo.txt','rep_puce_expo.txt','rep_puce.txt','rep_expo.txt']
    files = ['act_common.txt','rep_common.txt']
    #files = ['act_expo.txt','rep_expo.txt']
    #files = ['act_puce_expo.txt','rep_puce_expo.txt']
    res = {'act':{'tandem':0,'convergent':0,'divergent':0},'rep':{'tandem':0,'convergent':0,'divergent':0}}
    for name in files:
        file = open('{}{}'.format(basedir,name),'r')
        for line in file: 
            line = line.strip()
            try:
                res[name[0:3]][self.genes[line].orientation] += 1
            except:
                pass
        
    print res
    means = [] ; stds = []
    for state in ['act','rep']:
        try:
            pexp = float(res[state]['convergent']) / (res[state]['convergent']+res[state]['divergent'])
            tot = res[state]['convergent']+res[state]['divergent']
            std = np.array(stats.binom.std(tot, pexp, loc=0))/tot
            means.append(pexp)
            stds.append(std)
        except:
            pass

    titl = '/home/raphael/Documents/orientation'
    width = 5 ; height = width / 1.618
    fig, ax = plt.subplots()
    xval = [0,1] ; labs = ['rel','hyp']
    yval = means
    plt.plot(xval, yval, 'rD', markersize=7, linestyle='None')
    plt.errorbar(xval, yval,yerr=stds,mec='black', capsize=10, elinewidth=1,mew=1,linestyle='None', color='black')        
    
    fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
    ax.set_ylabel('conv / (conv + diverg)',fontweight='bold')
    ax.set_xlim(xval[0]-1,xval[-1]+1)
    #ax.set_ylim(0.4,0.6)
    plt.xticks(xval,labs,fontweight='bold')
    #ax.set_title('{}, {} et al.'.format(org,cond_fc),fontweight='bold')
    fig.set_size_inches(width, height)
    plt.tight_layout()
    plt.savefig(titl+'.svg',transparent=False)
    plt.close('all')

