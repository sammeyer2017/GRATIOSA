from scipy import stats
import statsmodels.api as sm
from globvar import *
import numpy as np
import matplotlib.pyplot as plt


def compute_spacer_response(gen,cond_fc,cond_tss,*arg,**kwargs):
    '''For all TSS in TSSs[cond_tss], compute spacer length for each sigma factor from promoter site coordinates.
    Then, associate expression data available in cond_fc for genes associated to TSS to spacer lengths, which allows
    the analysis of the link between FC and spacer length. kwargs = thresh_pval = only consider genes with pval below, 
    default 0.05. Returns a dict oh shape {[sigma_factor]:{[sp_length]:{[genes]:[g1,g2],[expr]:[1,2]}}}
    '''
    if not hasattr(gen, 'genes_valid'): # if no fc loaded 
        gen.load_fc_pval()
    if not hasattr(gen, 'TSSs'): # if no TSS loaded
        gen.load_TSS() 
    if not hasattr(gen, 'seq'): # if no seq loaded
        gen.load_seq() 

    thresh_pval = kwargs.get('thresh_pval', 0.05)
    spacer_sigma = {'org':gen.name} # dict of shape {[sigma_factor]:{[sp_length]:{[genes]:[g1,g2],[expr]:[1,2]}}}
    try :
        for TSSpos in gen.TSSs[cond_tss].keys(): # for all TSS in cond_tss
            TSSu = gen.TSSs[cond_tss][TSSpos] # single TS
            # valid gene : gene with pval < thresh_pval
            valid_expr = [] # expression values for valid genes in that TSS
            non_expr = [] # expression values for non valid genes in that TSS
            valid_genes = [] # valid genes in that TSS
            non_genes = [] # genes that appear in FC + TSS data but with pval > thresh_pval
            for gene in TSSu.genes:
                try:
                    if gen.genes[gene].fc_pval[cond_fc][1] <= thresh_pval:
                        valid_genes.append(gene) 
                        valid_expr.append(gen.genes[gene].fc_pval[cond_fc][0])
                    else:
                        non_genes.append(gene) 
                        non_expr.append(gen.genes[gene].fc_pval[cond_fc][0])
                except:
                    pass

            for sig in TSSu.promoter.keys(): # for all sigma factors binding to this TSS
                try:
                    if sig not in spacer_sigma.keys(): # add sigma factor to sigma dict
                        spacer_sigma[sig] = {'genes_valid':[],'genes_non':[],'expr_valid':[],'expr_non':[]}

                    try: # if site coordinates
                        TSSu.compute_magic_prom(gen.seq,gen.seqcompl)
                        spacer = len(TSSu.promoter[sig]['spacer'])
                    except: # if spacer length instead of sites coordinates
                        spacer = TSSu.promoter[sig]['sites'][0]

                    if spacer not in spacer_sigma[sig].keys(): # init spacer in dict
                        spacer_sigma[sig][spacer] = {'genes_valid':[],'genes_non':[],'expr_valid':[],'expr_non':[]} 
                    # shape 15:{[expr]:[expr1,expr2],[genes_valid]:[g1,g2],[genes_tot]:[g1,g2,g3]} 
                    # genes_valid and genes_tot implemented in order to keep in memory all
                    # genes with that spacer length
                    if non_genes != []:
                        spacer_sigma[sig]['genes_non'].append(non_genes) # gene in FC data + TSS data
                        spacer_sigma[sig]['expr_non'].append(non_expr) # gene in FC data + TSS data
                        spacer_sigma[sig][spacer]['genes_non'].append(non_genes)                               
                        spacer_sigma[sig][spacer]['expr_non'].append(non_expr)

                    if valid_genes != []:
                        spacer_sigma[sig]['genes_valid'].append(valid_genes) # gene in FC data + TSS data
                        spacer_sigma[sig]['expr_valid'].append(valid_expr) # gene in FC data + TSS data
                        spacer_sigma[sig][spacer]['genes_valid'].append(valid_genes)                               
                        spacer_sigma[sig][spacer]['expr_valid'].append(valid_expr)
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
    length for each sigma factor, compute txt results and graphes. kwargs : thresh_fc = for log2FC above,
    gene considered activated, for FC below, repressed (default 0). Aggregation = True or False. If True, 
    for each TSS, associate mean of expressions of genes, otherwise associate all single expression data. 
    l_genes = write in results the list of genes associated to each sigma factor and spacer length, default False.
    '''
    thresh_fc = kwargs.get('thresh_fc', 0.0)
    aggregation = kwargs.get('aggregation', True)
    l_genes = kwargs.get('l_genes', False)
    miscell = kwargs.get('miscell', False)
    # compute and write results
    # title for results
    titl = '{}-fc{}-agg{}-genes{}'.format(spacer_sigma['org'],str(thresh_fc),str(aggregation),str(l_genes))
    # path to database
    pathdb = '{}data/{}/spacer'.format(basedir,spacer_sigma['org'])
    if not os.path.exists(pathdb):
        os.makedirs(pathdb)
    # txt results
    results = open(pathdb+'/'+titl,'w')
    results.write('SigmaFactor|Spacerlength\tNon valid TSS associated\tValid TSS associated\tNon valid genes\t')
    results.write('Valid genes\tExpr mean\tActivated TSS\tRepressed TSS\tList non valid genes\tList valid genes\n')
    for sigma in spacer_sigma.keys():
        if sigma != 'org':
            if aggregation:                
                spacer_sigma[sigma]['expr_valid'] = np.array([np.mean(l) for l in spacer_sigma[sigma]['expr_valid']])
                spacer_sigma[sigma]['expr_non'] = np.array([np.mean(l) for l in spacer_sigma[sigma]['expr_non']])
            else:
                spacer_sigma[sigma]['expr_valid'] = np.array([item for sublist in spacer_sigma[sigma]['expr_valid'] for item in sublist])
                spacer_sigma[sigma]['expr_non'] = np.array([item for sublist in spacer_sigma[sigma]['expr_non'] for item in sublist])

            signon_tss = str(len(spacer_sigma[sigma]['genes_non']))
            sigvalid_tss = str(len(spacer_sigma[sigma]['genes_valid']))
            signon_genes = str(np.sum([len(l) for l in spacer_sigma[sigma]['genes_non']]))
            sigvalid_genes = str(np.sum([len(l) for l in spacer_sigma[sigma]['genes_valid']]))

            sig_val = spacer_sigma[sigma]['expr_valid']
            sigact_tss = str(sig_val[sig_val > thresh_fc].shape[0])
            sigrepr_tss = str(sig_val[sig_val < thresh_fc].shape[0])
            sigexpr_mean = str(np.mean(sig_val))

            results.write(str(sigma)+'\t'+signon_tss+'\t'+sigvalid_tss+'\t'+signon_genes+'\t'+sigvalid_genes+'\t')
            results.write(sigexpr_mean+'\t'+sigact_tss+'\t'+sigrepr_tss)

            if l_genes:
                results.write('\t'+str(spacer_sigma[sigma]['genes_non'])[1:-1].replace('\'',''))
                results.write('\t'+str(spacer_sigma[sigma]['genes_valid'])[1:-1].replace('\'','')+'\n')
            else:
                results.write('\n')

            for sp_len in spacer_sigma[sigma].keys():
                if type(sp_len) == int: # for keys that are not e.g. list of genes and expr
                    spnon_tss = str(len(spacer_sigma[sigma][sp_len]['genes_non']))
                    spnon_genes = str(np.sum([len(l) for l in spacer_sigma[sigma][sp_len]['genes_non']]))

                    if spacer_sigma[sigma][sp_len]['expr_valid'] == [] and spacer_sigma[sigma][sp_len]['expr_non'] != []: 
                    # spacer length with genes in FC + TSS data but without valid FC
                        results.write(str(sp_len)+'\t'+spnon_tss+'\t'+'NaN'+'\t'+spnon_genes+'\t'+'NaN'+'\t')
                        results.write('NaN'+'\t'+'NaN'+'\t'+'NaN')
                        if l_genes:
                            results.write('\t'+str(spacer_sigma[sigma][sp_len]['genes_non'])[1:-1].replace('\'',''))
                            results.write('\tNaN\n')
                        else:
                            results.write('\n')

                    else:
                        spvalid_tss = str(len(spacer_sigma[sigma][sp_len]['genes_valid']))
                        spvalid_genes = str(np.sum([len(l) for l in spacer_sigma[sigma][sp_len]['genes_valid']]))

                        if aggregation:
                            spacer_sigma[sigma][sp_len]['expr_valid'] = np.array([np.mean(l) for l in spacer_sigma[sigma][sp_len]['expr_valid']])
                            spacer_sigma[sigma][sp_len]['expr_non'] = np.array([np.mean(l) for l in spacer_sigma[sigma][sp_len]['expr_non']])
                        
                        else:
                            spacer_sigma[sigma][sp_len]['expr_valid'] = np.array([item for sublist in spacer_sigma[sigma][sp_len]['expr_valid'] for item in sublist])
                            spacer_sigma[sigma][sp_len]['expr_non'] = np.array([item for sublist in spacer_sigma[sigma][sp_len]['expr_non'] for item in sublist])

                        sp_val = spacer_sigma[sigma][sp_len]['expr_valid']
                        spact_tss = str(sp_val[sp_val > thresh_fc].shape[0])
                        sprepr_tss = str(sp_val[sp_val < thresh_fc].shape[0])
                        spexpr_mean = str(np.mean(sp_val))

                        results.write(str(sp_len)+'\t'+spnon_tss+'\t'+spvalid_tss+'\t'+spnon_genes+'\t'+spvalid_genes+'\t')
                        results.write(spexpr_mean+'\t'+spact_tss+'\t'+sprepr_tss)

                        if l_genes:
                            results.write('\t'+str(spacer_sigma[sigma][sp_len]['genes_non'])[1:-1].replace('\'',''))
                            results.write('\t'+str(spacer_sigma[sigma][sp_len]['genes_valid'])[1:-1].replace('\'','')+'\n')
                        else:
                            results.write('\n')
        
    results.close()

    # stats and figures

    # sigfactor in globvar : sigfactor of interest
    # spacers in globvar : spacers of interests
    def save_single_fig(graphtype):
        plt.savefig('{}/{}-{}.svg'.format(pathdb,titl,graphtype),transparent=False)
        plt.close('all')
    def add_labels(*arg,**kwargs):
        xlab = kwargs.get('xlab','Log2(FC)')
        ylab = kwargs.get('ylab','Spacer length')
        tit =  kwargs.get('tit',titl)
        plt.xlabel(xlab, fontweight = 'bold')
        plt.ylabel(ylab, fontweight = 'bold')
        plt.title(tit, fontweight='bold')

    def good_shape():
        plt.rcParams.update({'pdf.fonttype': 42})
        plt.rcParams.update({'ps.fonttype': 42})
        plt.rcParams.update({'font.size': 11})
        plt.rcParams.update({'legend.fontsize': 11})
        plt.rcParams.update({'font.family': "Arial"})
        # plt.rcParams.update('xtick', labelsize=11)
        # plt.rcParams.update('ytick', labelsize=11)
        # plt.rcParams.update('axes', labelsize=11)

    if miscell == True:
        fig = plt.figure(figsize=(6,6))
        d = spacer_sigma[sigfactor]
        plt.boxplot([d[x]['expr'] for x in spacers],labels=spacers)
        add_labels()
        save_single_fig('boxplot')

        rows = ['expr'] ; rows.extend(spacers)
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
        save_single_fig('hists')


    xval = []
    yval = []
    meanstd = []

    for sp in spacers:
        val = spacer_sigma[sigfactor][sp]['expr_valid']
        act = val[val > thresh_fc].shape[0]
        rep = val[val < thresh_fc].shape[0]
        tot = float(val.shape[0])
        pexp = float(act) / tot

        ci = np.array(stats.binom.interval(0.95, tot, pexp, loc=0))/tot
        meanval = np.array(stats.binom.mean(tot, pexp, loc=0))/tot
        std = np.array(stats.binom.std(tot, pexp, loc=0))/tot
        
        yval.append(meanval)
        meanstd.append(std)
        xval.append(sp)


    X = sm.add_constant(xval)
    wls_model = sm.WLS(yval, X, weights=1/np.power(np.array(meanstd),2))
    results = wls_model.fit()
    slope = results.params[1]
    OR = results.params[0]
    pval = results.pvalues[1]

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['legend.fontsize'] = 11
    plt.rcParams['font.size'] = 11
    plt.rcParams['font.family'] = 'Arial'

    width = 6
    height = width / 1.618

    fig, ax = plt.subplots()
    fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)

    ax.set_ylabel('Activated promoters (%)')
    ax.set_xlabel('Spacer length (nt)')
    ax.set_xlim(14,20)
    fig.set_size_inches(width, height)
    
    plt.plot(xval, yval, 'rD', markersize=9, linestyle='None')
    plt.errorbar(xval, yval,yerr=meanstd,mec='black', capsize=12.5, elinewidth=1.25,mew=1.5,linestyle='None', color='black')
    xval = [14] + xval + [20] # delimit plot
    plt.plot(xval,np.array(xval)*slope + OR, linestyle='dashed', color='black')
    plt.text(0.5,0.75, r'Lin reg : slope '+str(round(slope,3))+', pval '+str(round(pval,3)),fontsize = 8,fontweight = 'bold',horizontalalignment = 'center',transform = ax.transAxes)

    good_shape()
    save_single_fig('spacers')

    # compute spacer lengths for activated and repressed genes
    allact = []
    allrep = []
    allnone = []
    for sp in spacers:
        val = spacer_sigma[sigfactor][sp]['expr_valid']
        act = val[val > thresh_fc].shape[0]
        rep = val[val < thresh_fc].shape[0]
        non = spacer_sigma[sigfactor][sp]['expr_non'].shape[0]

        allact.extend([sp]*act)
        allrep.extend([sp]*rep)
        allnone.extend([sp]*non)

    mact = np.mean(allact)
    mrep = np.mean(allrep)

    eact = np.std(allact)/np.sqrt(len(allact))
    erep = np.std(allrep)/np.sqrt(len(allrep))
    
    if allnone != []:
        mnone = np.mean(allnone)
        enone = np.std(allnone)/np.sqrt(len(allnone))
        xval = [1,2,3]
        yval = [mact,mnone,mrep]
        stdval = [eact,enone,erep]
        col = ['red','black','blue']
        labs = ['Activated','Non affected','Repressed']
    else:
        xval = [1,2]
        yval = [mact,mrep]
        stdval = [eact,erep]
        col = ['red','blue']
        labs = ['Activated','Repressed']

    plt.plot(xval,yval, markersize=9,marker='D',linestyle='None')#,color=col)
    plt.errorbar(xval, yval,yerr=stdval,linestyle='None', capsize=12.5, ecolor='black', elinewidth=1.25, mew=1.5)
    
    plt.xticks(xval,labs)
    plt.xlim(xval[0]-1,xval[-1]+1)
    
    plt.show()
