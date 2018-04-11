from scipy import stats
import statsmodels.api as sm
from globvar import *
import numpy as np
import matplotlib.pyplot as plt

def compute_spacer_response(gen,cond_fc,cond_tss,*arg,**kwargs):
    '''For all TSS in TSSs[cond_tss], compute spacer length for each sigma factor from promoter site coordinates.
    Then, associate expression data available in cond_fc for genes associated to TSS to spacer lengths, which allows
    the analysis of the link between FC and spacer length. kwargs = thresh_pval = only consider genes with pval below, default 0.05. Returns a dict oh shape {[sigma_factor]:{[sp_length]:{[genes]:[g1,g2],[expr]:[1,2]}}}
    '''
    if not hasattr(gen, 'genes_valid'): # if no fc loaded 
        gen.load_fc_pval()
    if not hasattr(gen, 'TSSs'): # if no TSS loaded
        gen.load_TSS() 
    if not hasattr(gen, 'seq'): # if no seq loaded
        gen.load_seq() 

    thresh_pval = kwargs.get('thresh_pval', 1)
    try:
        thresh_pval = float(thresh_pval)
    except:
        print 'Invalid thresh_pval, setting default 1...'
        thresh_pval = 1

    spacer_sigma = {} # dict of shape {[sigma_factor]:{[sp_length]:{[genes]:[g1,g2],[expr]:[1,2]}}}

    try :
        for TSSpos in gen.TSSs[cond_tss].keys(): # for all TSS in cond_tss
            TSSu = gen.TSSs[cond_tss][TSSpos] # single TS
            for sig in TSSu.promoter.keys(): # for all sigma factors binding to this TSS
                try:
                    if sig not in spacer_sigma.keys(): # add sigma factor to sigma dict
                        spacer_sigma[sig] = {'genes_valid':[],'genes_tot':[],'expr':[]}
                    try: # if site coordinates
                        TSSu.compute_magic_prom(gen.seq,gen.seqcompl)
                        spacer = len(TSSu.promoter[sig]['spacer'])
                    except: # if spacer length instead of sites coordinates
                        spacer = TSSu.promoter[sig]['sites'][0]

                    # valid gene : gene with pval < thresh_pval
                    expr = [] # expression values for valid genes in that TSS 
                    valid_genes = [] # valid genes in that TSS
                    tot_genes = [] # genes that appear in FC + TSS data but with pval > thresh_pval

                    for gene in TSSu.genes:
                        try:
                            if gene in gen.genes_valid[cond_fc]:
                                tot_genes.append(gene) # gene in FC data + TSS data

                                if gen.genes[gene].fc_pval[cond_fc][1] <= thresh_pval:
                                # gene in FC data + TSS data and valid (pval < thresh)
                                    valid_genes.append(gene) 
                                    expr.append(gen.genes[gene].fc_pval[cond_fc][0])
                        except:
                            pass

                    if tot_genes != []:
                        if spacer not in spacer_sigma[sig].keys(): # init spacer in dict
                            spacer_sigma[sig][spacer] = {'genes_valid':[],'genes_tot':[],'expr':[]} # shape 15:{[expr]:[expr1,expr2],[genes_valid]:[g1,g2],[genes_tot]:[g1,g2,g3]} 
                        # genes_valid and genes_tot implemented in order to keep in memory all
                        # genes with that spacer length
                        spacer_sigma[sig]['genes_tot'].append(tot_genes) # gene in FC data + TSS data
                        spacer_sigma[sig][spacer]['genes_tot'].append(tot_genes)

                        if expr != [] and valid_genes != []:
                            spacer_sigma[sig]['expr'].append(expr)
                            spacer_sigma[sig]['genes_valid'].append(valid_genes)                            
                            spacer_sigma[sig][spacer]['expr'].append(expr)                               
                            spacer_sigma[sig][spacer]['genes_valid'].append(valid_genes)

                        print 'Success, FC',expr,'for genes',valid_genes,'added to sigma',sig,'with',spacer,'pb spacer'
                except:
                    pass

    except Exception as e:
        print 'Error',e   

    if spacer_sigma == {}:
        print 'Unable to compute spacer_sigma'
    else:
        return spacer_sigma

def draw_results(spacer_sigma,*arg,**kwargs):
    '''Starting from spacer_sigma, the dict associating expression values to each spacer
    length for each sigma factor, compute results and graphes. kwargs : thresh_fc = for log2FC above,
    gene considered activated, for FC below, repressed (default 0). Aggregation = True or False. If True, for each TSS, associate mean of expressions of genes, 
    otherwise associate all single expression data. l_genes = write in results the list of genes associated 
    to each sigma factor and spacer length, default False.
    '''
    thresh_fc = kwargs.get('thresh_fc', 0.0)
    try:
        thresh_fc = float(thresh_fc)
    except:
        print 'Invalid thresh_fc, setting default 0...'
        thresh_fc = 0.0

    aggregation = kwargs.get('aggregation', True)
    if type(aggregation) != bool:
        print 'Invalid aggregation (bool needed), setting default True...'
        aggregation = True

    l_genes = kwargs.get('l_genes', False)
    if type(l_genes) != bool:
        print 'Invalid l_genes (bool needed), setting default False...'
        l_genes = False

    # compute and write results
    titl = 'spacer-fc-'+str(thresh_fc)+'-agg-'+str(aggregation)+'-genes-'+str(l_genes)

    results = open(basedir+"data/"+titl+'.txt','w')

    results.write('SigmaFactor|Spacerlength\tTotal TSS associated\tValid TSS associated\tTotal genes\t')
    results.write('Valid genes\tExpr mean\tActivated TSS\tRepressed TSS\tList total genes\tList valid genes\n')
    
    for sigma in spacer_sigma.keys():
        sigtot_tss = str(len(spacer_sigma[sigma]['genes_tot']))
        sigvalid_tss = str(len(spacer_sigma[sigma]['genes_valid']))
        sigtot_genes = str(np.sum([len(l) for l in spacer_sigma[sigma]['genes_tot']]))
        sigvalid_genes = str(np.sum([len(l) for l in spacer_sigma[sigma]['genes_valid']]))

        if aggregation:
            sig_val = np.array([np.mean(l) for l in spacer_sigma[sigma]['expr']])
        else:
            sig_val = np.array([item for sublist in spacer_sigma[sigma]['expr'] for item in sublist])

        sigact_tss = str(sig_val[sig_val > thresh_fc].shape[0])
        sigrepr_tss = str(sig_val[sig_val < thresh_fc].shape[0])
        sigexpr_mean = str(np.mean(sig_val))

        results.write(str(sigma)+'\t'+sigtot_tss+'\t'+sigvalid_tss+'\t'+sigtot_genes+'\t'+sigvalid_genes+'\t')
        results.write(sigexpr_mean+'\t'+sigact_tss+'\t'+sigrepr_tss)

        if l_genes:
            results.write('\t'+str(spacer_sigma[sigma]['genes_tot'])[1:-1].replace('\'',''))
            results.write('\t'+str(spacer_sigma[sigma]['genes_valid'])[1:-1].replace('\'','')+'\n')
        else:
            results.write('\n')

        for sp_len in spacer_sigma[sigma].keys():
            if type(sp_len) == int: # for keys that are not genes valid, genes tot or expr
                sptot_tss = str(len(spacer_sigma[sigma][sp_len]['genes_tot']))
                sptot_genes = str(np.sum([len(l) for l in spacer_sigma[sigma][sp_len]['genes_tot']]))

                if spacer_sigma[sigma][sp_len]['expr'] == []: 
                # spacer length with genes in FC + TSS data but without valid FC
                    results.write(str(sp_len)+'\t'+sptot_tss+'\t'+'NaN'+'\t'+sptot_genes+'\t'+'NaN'+'\t')
                    results.write('NaN'+'\t'+'NaN'+'\t'+'NaN')
                    if l_genes:
                        results.write('\t'+str(spacer_sigma[sigma][sp_len]['genes_tot'])[1:-1].replace('\'',''))
                        results.write('\tNaN\n')
                    else:
                        results.write('\n')

                else:    
                    spvalid_tss = str(len(spacer_sigma[sigma][sp_len]['genes_valid']))
                    spvalid_genes = str(np.sum([len(l) for l in spacer_sigma[sigma][sp_len]['genes_valid']]))

                    if aggregation:
                        sp_val = np.array([np.mean(l) for l in spacer_sigma[sigma][sp_len]['expr']])
                    else:
                        sp_val = np.array([item for sublist in spacer_sigma[sigma][sp_len]['expr'] for item in sublist])

                    spact_tss = str(sp_val[sp_val > thresh_fc].shape[0])
                    sprepr_tss = str(sp_val[sp_val < thresh_fc].shape[0])
                    spexpr_mean = str(np.mean(sp_val))

                    results.write(str(sp_len)+'\t'+sptot_tss+'\t'+spvalid_tss+'\t'+sptot_genes+'\t'+spvalid_genes+'\t')
                    results.write(spexpr_mean+'\t'+spact_tss+'\t'+sprepr_tss)

                    if l_genes:
                        results.write('\t'+str(spacer_sigma[sigma][sp_len]['genes_tot'])[1:-1].replace('\'',''))
                        results.write('\t'+str(spacer_sigma[sigma][sp_len]['genes_valid'])[1:-1].replace('\'','')+'\n')
                    else:
                        results.write('\n')

    results.close()

    # stats and graphes
    if aggregation:
        spacer_sigma[sigfactor]['expr'] = np.array([np.mean(l) for l in spacer_sigma[sigfactor]['expr']])
    else:
        spacer_sigma[sigfactor]['expr'] = np.array([item for sublist in spacer_sigma[sigfactor]['expr'] for item in sublist])

    for el in spacers:
        if aggregation:
            spacer_sigma[sigfactor][el]['expr'] = np.array([np.mean(l) for l in spacer_sigma[sigfactor][el]['expr']])
        else:
            spacer_sigma[sigfactor][el]['expr'] = np.array([item for sublist in spacer_sigma[sigfactor][el]['expr'] for item in sublist])

    fig = plt.figure(figsize=(6,6))
    plt.subplot(1,1,1)
    d = spacer_sigma[sigfactor]
    plt.boxplot([d[x]['expr'] for x in spacers],labels=spacers)
    plt.ylabel('log2(FC)')
    plt.xlabel('spacer length')
    plt.title(titl,fontsize = 9.5)
    plt.savefig(basedir+"data/"+titl+"-boxplot.pdf",transparent=False)
    plt.close('all')

    rows = ['expr',15,16,17,18,19]
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(1,1,1)
    i = 1
    for row in rows:
        fig.add_subplot(3,3,i)
        if row == 'expr':
            x = spacer_sigma[sigfactor][row]
            plt.title('All expression values')
        else:
            x = spacer_sigma[sigfactor][row]['expr']
            plt.title('Spacer '+str(row)+' pb')

        plt.hist(x,bins=10, range=(-2.5,2.5),color = 'yellow',edgecolor = 'red')
        i+=1
    
    ax.axis('off')
    ax.text(0.5, 0.25, '(x) Log2(FC), (y) Nb of values',fontsize = 18, rotation = 0, horizontalalignment = 'center', transform = ax.transAxes)

    plt.savefig(basedir+"data/"+titl+"-distributions.pdf",transparent=False)
    plt.close('all')

    xval = []
    yval = []
    meanstd = []
    for row in spacers:
        val = spacer_sigma[sigfactor][row]['expr']
        act = val[val > thresh_fc].shape[0]
        rep = val[val < thresh_fc].shape[0]
        tot = float(val.shape[0])
        pexp = float(act) / tot

        ci = np.array(stats.binom.interval(0.95, tot, pexp, loc=0))/tot
        meanval = np.array(stats.binom.mean(tot, pexp, loc=0))/tot
        std = np.array(stats.binom.std(tot, pexp, loc=0))/tot
        yval.append(meanval)
        meanstd.append(std)
        xval.append(row)

    #xvalreg = [xval[0]]+xval[2:] ; yvalreg = [yval[0]]+yval[2:] ; erreg = [meanstd[0]]+meanstd[2:]
    xvalreg = xval ; yvalreg = yval ; erreg = meanstd   
    # weights : 1 / var
    X = sm.add_constant(xvalreg)
    wls_model = sm.WLS(yvalreg, X, weights=1/np.power(np.array(erreg),2))
    results = wls_model.fit()
    print results.summary()
    slope = results.params[1]
    OR = results.params[0]
    pval = results.pvalues[1]

    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(111)

    plt.plot(xvalreg, yvalreg, 'rD', markersize=9, linestyle='None')
    #plt.plot(xval[1],yval[1], marker='D', color='lightgray', markersize=9, linestyle='None')
    plt.errorbar(xvalreg, yvalreg,yerr=erreg,mec='black', capsize=12.5, elinewidth=1.25,mew=1.5,linestyle='None', color='black')
    #plt.errorbar(xval[1], yval[1],yerr=meanstd[1],mec='lightgray', capsize=12.5, elinewidth=1.25,mew=1.5,linestyle='None', color='gray')

    xvalreg = [14] + xvalreg + [20]
    plt.plot(xvalreg,np.array(xvalreg)*slope + OR, linestyle='dashed', color='black')
    plt.text(0.5,0.75, r'Lin reg : slope '+str(round(slope,3))+', pval '+str(round(pval,3)),fontsize = 12,fontweight = 'bold',horizontalalignment = 'center',transform = ax.transAxes)

    plt.ylabel('Activated genes proportion',fontweight = 'bold', fontsize = 15)
    plt.xlabel('Spacer length',fontweight = 'bold', fontsize = 15)
    plt.xlim(14,20)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.savefig(basedir+"data/"+titl+"-binom.svg", transparent=False)
    plt.close('all')

    allact = []
    allrepr = []
    rows = [15,16,17,18,19]
    for row in rows:
        val = spacer_sigma[sigfactor][row]['expr']
        act = val[val > thresh_fc].shape[0]
        rep = val[val < thresh_fc].shape[0]
        tot = float(val.shape[0])

        allact.extend([row]*act)
        allrepr.extend([row]*rep)


    fig = plt.figure(figsize=(10,7))
    barWidth = 0.35
    y1 = [round(100*allact.count(x)/float(len(allact)),3) for x in rows]
    y2 = [round(100*allrepr.count(x)/float(len(allrepr)),3) for x in rows]
    r1 = range(len(y1))
    r2 = [x + barWidth for x in r1]
    b1 = plt.bar(r1, y1, width = barWidth, color = ['red' for i in y1],edgecolor = ['black' for i in y1], linewidth = 2)
    b2 = plt.bar(r2, y2, width = barWidth, color = ['blue' for i in y2],edgecolor = ['black' for i in y2], linewidth = 2)

    plt.xticks([r + barWidth for r in range(len(y1))], ['15', '16', '17', '18', '19'], fontsize=14)
    plt.yticks(fontsize=14)

    plt.xlabel("Spacer length",fontweight = 'bold', fontsize = 15)
    plt.ylabel("% Genes",fontweight = 'bold', fontsize = 15)

    plt.legend([b1, b2], [str(len(allact))+' Activated', str(len(allrepr))+' Repressed'], title= 'Genes', ncol = 1, fontsize = 15)
    # plt.text(r1[0]+barWidth/2,30, 't-test comparison : spacer length mean between activated and repressed genes signicantly different'+str(round(float(stats.mannwhitneyu(allact,allrepr)[1]),3)),fontsize = 15, style = 'italic')
    print stats.mannwhitneyu(allact,allrepr)
    # plt.annotate("", xy=(r2[2]+barWidth/2,36), xycoords='data',xytext=(r1[4]-barWidth/2,36), textcoords='data',arrowprops=dict(arrowstyle="-", ec='#aaaaaa',connectionstyle="bar,fraction=0.2"))
    # plt.text(r2[3],35, round(float(stats.mannwhitneyu(allact,allrepr)[1]),3), ha = 'center', va = 'center')
    plt.savefig(basedir+"data/"+titl+"-histfin.svg", transparent=False)
    plt.close('all')
    print stats.ttest_ind(allact,allrepr)