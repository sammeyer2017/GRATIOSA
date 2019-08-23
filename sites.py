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
# Sites: NAP binding, gyrase cleavage...
def load_sites(gen):
    '''
    Load sites from sites.info: left and right binding coordinates on genome, sequence of site, score, strand
    '''
    gen.sites = {}
    if os.path.exists(basedir+"data/"+gen.name+"/sites/sites.info"):
        with open(basedir+"data/"+gen.name+"/sites/sites.info","r") as f:
            skiphead = next(f) # skip head
            for header in f:
                header=header.strip()
                header=header.split('\t')
                gen.sites[header[0]] = load_sites_cond(basedir+"data/"+gen.name+"/sites/"+header[1], header[0], header[2], int(header[3]), int(header[4]), int(header[5]), int(header[6]), int(header[7]), int(header[8]))
        f.close()

def load_sites_cond(filename, condition, separator, start_line, left, right, strand, seq, score, *args, **kwargs):
    gen_sites = {}
    with open(filename, 'r') as f:
        i=1
        while i < start_line:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n')
            if separator == '\\t':
                line = line.split('\t')
            else:
                line=line.split(separator)
            try:
                gen_sites[(int(line[left]) + int(line[right]))/2] = {"left":int(line[left]),"right":int(line[right]),"score":float(line[score]), "strand":line[strand], "seq":line[seq]}
            except:
                pass
    f.close()
    return gen_sites


def characterize_sites(gen,cond,*args, **kwargs):
    '''
    For each site, find the 2 nearest neighbours on genome whatever the strand is and associates FC data
    '''
    if not hasattr(gen, 'sites'):
        load_sites(gen)

    pathdb = '{}data/{}/sites/results'.format(basedir,gen.name)
    if not os.path.exists(pathdb):
      os.makedirs(pathdb)

    # genome FC values
    gen.load_fc_pval()
    # export results
    export_csv = kwargs.get('export',True)
    # FC data if present
    cond_FC = kwargs.get('cond_FC',"No expr")

    genes = {}
    # create dic with genes and position (res[position] = gene)
    for i in gen.genes:
        genes[int(gen.genes[i].start)]=i
    # sort by position
    genes_pos=[x[0] for x in sorted(list(genes.items()), key=operator.itemgetter(0))]

    for s in gen.sites[cond].keys():
        try:
            site = gen.sites[cond][s]
            # sort genes depending on distance to site
            l = sorted([(abs(x-s),x) for x in genes_pos])
            # g1 = closest to site, g2 = 2nd closest
            g1 = gen.genes[genes[l[0][1]]] ; g2 = gen.genes[genes[l[1][1]]]
            # if no FC data
            if cond_FC == "No expr":
                expr1 = ("NA","NA") ; expr2 = ("NA","NA")
            else: # try to find corresponding FC data, shape (FC,pvalue)
                try:
                    expr1 = g1.fc_pval[cond_FC]
                except:
                    expr1 = ("NA","NA")
                try:
                    expr2 = g2.fc_pval[cond_FC]
                except:
                    expr2 = ("NA","NA")

            # extract relevant informations from g1 and g2 for site: locus tag, name, distance to site, strand, FC data
            site["g1"] = [g1.gene_id, g1.name, abs(g1.start - s), g1.strand, expr1]
            site["g2"] = [g2.gene_id, g2.name, abs(g2.start - s), g2.strand, expr2]

        except Exception as e:
            print e
            pass

    if export_csv:
    	# export data
        titl = "{}".format(cond)
        res = open(pathdb+"/{}.csv".format(titl),'w')
        res.write('Left\tRight\tStrand\tSequence\tScore\tGene 1 locus tag\tGene 1 name\tDistance to ATG\tStrand\t{} log2(FC)\t{} adj pvalue\tGene 2 locus tag\tGene 2 name\tDistance to ATG\tStrand\t{} log2(FC)\t{} adj pvalue\n'.format(cond_FC,cond_FC,cond_FC,cond_FC))
        
        sites = gen.sites[cond].keys() ; sites.sort()
        for s in sites:
            site = gen.sites[cond][s]
            try:
                res.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(str(site["left"]),str(site["right"]), site["strand"], site["seq"], site["score"], site["g1"][0], site["g1"][1], str(site["g1"][2]), str(site["g1"][3]), str(site["g1"][4][0]), str(site["g1"][4][1]), site["g2"][0], site["g2"][1], str(site["g2"][2]), str(site["g2"][3]), str(site["g2"][4][0]), str(site["g2"][4][1])))
            except Exception as e:
                print e
                pass

        res.close()


# THIS PART WAS USED TO CHARACTERIZE GYRASE CLEAVAGE SITES, AN ANALYSIS IMPOSED BY NAR REVIEWERS FOR THE 
# PAPER Bacterial genome architecture shapes global transcriptional regulation by DNA supercoiling

# def load_sites(gen):
#     '''
#     Load sites from file and for each site, find nearest neighbours (left and right) on genome and computes orientation
#     '''
#     gen.sites = {}
#     if os.path.exists(basedir+"data/"+gen.name+"/GCS/GCS.info"):
#         with open(basedir+"data/"+gen.name+"/GCS/GCS.info","r") as f:
#             skiphead = next(f) # skip head
#             for header in f:
#                 header=header.strip()
#                 header=header.split('\t')
#                 gen.sites[header[0]] = load_sites_cond(basedir+"data/"+gen.name+"/GCS/"+header[1], header[0], int(header[4]), int(header[5]), int(header[6]), header[3], int(header[2]))
#         f.close()

# def load_sites_cond(filename, condition, coor_col, N_col, score_col, separator, start_line, *args, **kwargs):
#     gen_sites = {}
#     with open(filename, 'r') as f:
#         i=1
#         while i < start_line:
#             header=next(f)
#             i+=1
#         for line in f:
#             line = line.strip('\n')
#             if separator == '\\t':
#                 line = line.split('\t')
#             else:
#                 line=line.split(separator)
#             try:
#                 gen_sites[int(line[coor_col])] = {"pos":int(line[coor_col]),"N3E":float(line[N_col]), "score":float(line[score_col])}
#             except:
#                 pass
#     f.close()
#     return gen_sites


# def compare_sites_orientation(gen,*args, **kwargs):
#     if not hasattr(gen, 'sites'):
#         load_sites(gen)

#     pathdb = '{}data/{}/orientation'.format(basedir,gen.name)
#     if not os.path.exists(pathdb):
#       os.makedirs(pathdb)

#     # genome stats
#     gen.load_gene_orientation()

#     res = {}
#     # create dic with genes and position (res[position] = gene)
#     for i in gen.genes:
#         res[int(gen.genes[i].start)]=i
#     # sort dic per position of genes
#     l_res=sorted(list(res.items()), key=operator.itemgetter(0))

#     intergenics = []
#     for i in range(1,len(l_res)):
#         g1 = gen.genes[l_res[i-1][1]]
#         g2 = gen.genes[l_res[i][1]]

#         if g1.strand and g2.strand:
#             intergenics.append([g1.end, g2.start,"tandem"])
#         elif g1.strand and not g2.strand:
#             intergenics.append([g1.end, g2.end,"convergent"])
#         elif not g1.strand and g2.strand:
#             intergenics.append([g1.start, g2.start,"divergent"])
#         elif not g2.strand and not g2.strand:
#             intergenics.append([g1.start, g2.end,"tandem"])

#     res_genome = {"intergenic":{"tandem":0,"convergent":0,"divergent":0}, "intragenic":{"tandem":0,"convergent":0,"divergent":0}, "intra-inter":{"tandem":0,"convergent":0,"divergent":0}}
#     for gene in gen.genes.keys():
#         try:
#             g = gen.genes[gene]
#             res_genome["intragenic"][g.orientation] += g.length
#             res_genome["intra-inter"][g.orientation] += g.length
           
#         except:
#             pass

#     for region in intergenics:
#         if region[1] - region[0] > 0: # if non overlapping
#             res_genome["intergenic"][region[2]] += (region[1] - region[0])
#             res_genome["intra-inter"][region[2]] += (region[1] - region[0])

#     # number of convergent, divergent and tandem genes for each state 
#     conds = gen.sites.keys()
#     orientations = ['convergent','tandem','divergent']
#     res_sites = {}
#     for cond in conds:
#         res_sites[cond] = {"intergenic":{},"intragenic":{},"intra-inter":{}}
#         for orientation in orientations:
#             res_sites[cond]["intergenic"][orientation] = []
#             res_sites[cond]["intragenic"][orientation] = []
#             res_sites[cond]["intra-inter"][orientation] = []

#     for cond in conds:
#         for site in gen.sites[cond].keys():
#             try:               
#                 s = gen.sites[cond][site]
#                 intergenic = False
#                 for region in intergenics:
#                     if region[1] - region[0] > 0 and region[0] <= site <= region[1]:# if non overlapping
#                         res_sites[cond]["intergenic"][region[2]].append(site)
#                         res_sites[cond]["intra-inter"][region[2]].append(site)
#                         intergenic = True

#                 if not intergenic:
#                     g = res[min([x for x in res.keys()], key=lambda x:abs(x-site))]
#                     res_sites[cond]["intragenic"][gen.genes[g].orientation].append(site)
#                     res_sites[cond]["intra-inter"][gen.genes[g].orientation].append(site)

#             except:
#                 pass
                   
#     typ = kwargs.get("typ","intra-inter")
#     print 'Type:',typ

#     gen_tand = res_genome[typ]["tandem"]
#     gen_conv = res_genome[typ]["convergent"]
#     gen_div = res_genome[typ]["divergent"]
#     gen_tot = float(gen_tand + gen_conv + gen_div)

#     gen_ptand = gen_tand / gen_tot
#     gen_pconv = gen_conv / gen_tot
#     gen_pdiv = gen_div / gen_tot
#     print "Genome",gen_tot,"conv:",gen_conv,"div:",gen_div,"tand:",gen_tand
#     print "pconv",gen_pconv,"ptand",gen_ptand,"pdiv",gen_pdiv

#     pconv = [] ; ptand = [] ; pdiv = []
#     nbconv = [] ; nbtand = [] ; nbdiv = []
#     scores_conv = [] ; scores_tand = [] ; scores_div = []
#     seqs_conv = [] ; seqs_tand = [] ; seqs_div = []

#     ciconv = [] ; citand = [] ; cidiv = []

#     varconv = [] ; vartand = []; vardiv = []
#     pvals = []


#     for cond in ["Micro","Oxo","Cfx","RifCfx"]:
#         conv = float(len(res_sites[cond][typ]['convergent']))
#         div = float(len(res_sites[cond][typ]['divergent']))
#         tand = float(len(res_sites[cond][typ]['tandem']))

#         nbconv.append(conv) ; nbdiv.append(div) ; nbtand.append(tand)

#         tot = conv + div + tand
#         print cond, tot, "conv",conv,"div",div,"tand",tand

#         ci = proportion.proportion_confint(conv, tot, alpha=0.05, method='normal')
#         pconv.append(conv/tot) ; ciconv.append(tot*(ci[1] - ci[0])/2)
#         varconv.append(stats.binom.var(tot, conv/tot, loc=0))

#         ci = proportion.proportion_confint(tand, tot, alpha=0.05, method='normal')
#         ptand.append(tand/tot) ; citand.append(tot*(ci[1] - ci[0])/2)
#         vartand.append(stats.binom.var(tot, tand/tot, loc=0))

#         ci = proportion.proportion_confint(div, tot, alpha=0.05, method='normal')
#         pdiv.append(div/tot) ; cidiv.append(tot*(ci[1] - ci[0])/2)
#         vardiv.append(stats.binom.var(tot, div/tot, loc=0))
#         print "VAR",stats.binom.var(tot, div/tot, loc=0)

#         stat, pval = proportion.proportions_ztest([conv,gen_conv], [conv+div,gen_conv+gen_div], alternative = "larger")
#         pvals.append(pval)
#         print "pval",pval

#         score_conv = [gen.sites[cond][x]['N3E'] for x in res_sites[cond][typ]["convergent"]]
#         score_div = [gen.sites[cond][x]['N3E'] for x in res_sites[cond][typ]["divergent"]] 
#         score_tand = [gen.sites[cond][x]['N3E'] for x in res_sites[cond][typ]["tandem"]] 

#         scores_conv.append(np.sum(score_conv)) ; scores_tand.append(np.sum(score_tand)) ; scores_div.append(np.sum(score_div))

#         seq_conv = [gen.sites[cond][x]['score'] for x in res_sites[cond][typ]["convergent"]]
#         seq_div = [gen.sites[cond][x]['score'] for x in res_sites[cond][typ]["divergent"]] 
#         seq_tand = [gen.sites[cond][x]['score'] for x in res_sites[cond][typ]["tandem"]] 

#         seqs_conv.append(np.sum(seq_conv)) ; seqs_tand.append(np.sum(seq_tand)) ; seqs_div.append(np.sum(seq_div))

#     norm = kwargs.get("norm",True)
#     final = kwargs.get("final",1)
#     fact = 1000

#     if norm:
#         fconv = fact / float(res_genome[typ]["convergent"]) ; fdiv = fact / float(res_genome[typ]["divergent"]) ; ftand = fact / float(res_genome[typ]["tandem"])
       
#         nbconv = np.asarray(nbconv) * fconv ; nbtand = np.asarray(nbtand) * ftand ; nbdiv = np.asarray(nbdiv) * fdiv
#         ciconv = np.asarray(ciconv) * fconv ; citand = np.asarray(citand) * ftand ; cidiv = np.asarray(cidiv) * fdiv
#         scores_conv = np.asarray(scores_conv) * fconv ; scores_tand = np.asarray(scores_tand) * ftand ; scores_div = np.asarray(scores_div) * fdiv
#         seqs_conv = np.asarray(seqs_conv) * fconv ; seqs_tand = np.asarray(seqs_tand) * ftand ; seqs_div = np.asarray(seqs_div) * fdiv

#         varconv = fconv**2*np.asarray(varconv) ; vartand = ftand**2*np.asarray(vartand) ; vardiv = fdiv**2*np.asarray(vardiv)        

#         if final == 1:
#             titl1 = 'global_{}_norm_final1.png'.format(typ)
#             titl2 = 'global_{}_norm_final2.png'.format(typ)
           
#             width = 4 ; height = 2
#             fig, ax = plt.subplots()

#             barWidth = 0.2
#             y1 = nbconv ; y2 =nbtand ; y3 = nbdiv ; ci1 = ciconv ; ci2 = citand ; ci3 = cidiv
#             scal = kwargs.get("scal",True)
#             if scal: # scalability
#                 y1[0] = y1[0]*5 ; y2[0] = y2[0]*5 ; y3[0] = y3[0]*5 
#                 ci1[0] = ci1[0]*5 ; ci2[0] = ci2[0]*5 ; ci3[0] = ci3[0]*5

#             r1 = range(len(y1))
#             r2 = [x + barWidth for x in r1]
#             r3 = [x + barWidth for x in r2]
#             plt.bar(r1, y1, width = barWidth, linewidth = 2, label="conv", yerr = ci1)
#             plt.bar(r2, y2, width = barWidth, linewidth = 2, label="tand", yerr = ci2)
#             plt.bar(r3, y3, width = barWidth, linewidth = 2, label="div", yerr = ci3)

#             r = list([x[0] for x in [r1,r2,r3]] + [x[1] for x in [r1,r2,r3]] + [x[2] for x in [r1,r2,r3]] + [x[3] for x in [r1,r2,r3]])
#             y = list([x[0] for x in [y1,y2,y3]] + [x[1] for x in [y1,y2,y3]] + [x[2] for x in [y1,y2,y3]] + [x[3] for x in [y1,y2,y3]])
#             yerr = list([x[0] for x in [ci1,ci2,ci3]] + [x[1] for x in [ci1,ci2,ci3]] + [x[2] for x in [ci1,ci2,ci3]] + [x[3] for x in [ci1,ci2,ci3]])

#             i = 0
#             for pval in pvals:
#                 s = significance(pval)
#                 barplot_annotate_brackets(i,i+2,s,r,y,dh=0.1,barh=0.02)#,dh=0.05, barh=.05, fs=14, dt=0)
#                 i += 3

#             plt.ylim(top=max(y)+0.5)
#             plt.xticks([r + barWidth / 1.5 for r in range(len(y1))], ["Micro","Oxo","Cfx","RifCfx"])        
#             plt.ylabel('binding sites density\n(kb-1)')
#             plt.legend()

#             fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
#             fig.set_size_inches(width, height)
#             plt.tight_layout()
#             plt.savefig(pathdb+"/"+titl1)
#             plt.close('all')

#             width = 1.75 ; height = 2
#             fig, ax = plt.subplots()

#             # computes ratio of variances, ratio of conv / div and p-value
#             ratiovar1 =(nbconv[-2]**2/nbdiv[-2]**2)*((varconv[-2]/nbconv[-2]**2) + (vardiv[-2]/nbdiv[-2]**2))
#             ratiovar2 = (nbconv[-1]**2/nbdiv[-1]**2)*((varconv[-1]/nbconv[-1]**2) + (vardiv[-1]/nbdiv[-1]**2))

#             ratio1 = nbconv[-2] / nbdiv[-2]
#             ratio2 = nbconv[-1] / nbdiv[-1]

#             Zscore = abs((ratio1 - ratio2)/np.sqrt(ratiovar1 + ratiovar2))
#             pval = 1 - stats.norm.cdf(Zscore)

#             print "Z-score",Zscore,"Pval",pval

#             v = [2*np.sqrt(ratiovar1),2*np.sqrt(ratiovar2)]
#             y = [ratio1,ratio2]
#             r = [0,0.5]
#             plt.bar(r[0], y[0], width = barWidth, linewidth = 2, yerr = v[0], color = 'white',edgecolor = 'black')
#             plt.bar(r[1], y[1], width = barWidth, linewidth = 2, yerr = v[1], color = 'white',edgecolor = 'black')
#             #plt.bar(xval, yval, yerr = cim, width = 0.6, color = 'white',edgecolor = 'black', linewidth = 2, ecolor = 'black', capsize = 5)

#             plt.xticks([0,0.5],["Cfx","RifCfx"])
#             plt.ylabel('conv/div\nenrichment')
#             plt.ylim(2.5,top=max(y)+1.5)
#             plt.xlim((-0.25,0.75 ))

#             s = significance(pval)
#             barplot_annotate_brackets(0,1,s,r,y,yerr=v,dh=0.1,barh=0.02)#,dh=0.05, barh=.05, fs=14, dt=0)

#             fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
#             fig.set_size_inches(width, height)
#             plt.tight_layout()
#             plt.savefig(pathdb+"/"+titl2)
#             plt.close('all')


#         if final == 2:
#             titl = 'global_{}_norm.png'.format(typ)
#             width = 6 ; height = 4.5
#             fig, ax = plt.subplots()
            
#             plt.subplot(311)
#             barWidth = 0.2

#             y1 = nbconv ; y2 =nbtand ; y3 = nbdiv ; ci1 = cisconv ; ci2 = cistand ; ci3 = cisdiv
#             print "Sites:",y1,y2,y3
#             y1[0] = y1[0]*5 ; y2[0] = y2[0]*5 ; y3[0] = y3[0]*5 # scalability
#             ci1[0] = ci1[0]*5 ; ci2[0] = ci2[0]*5 ; ci3[0] = ci3[0]*5 # scalability

#             r1 = range(len(y1))
#             r2 = [x + barWidth for x in r1]
#             r3 = [x + barWidth for x in r2]

#             y_r1 = [y1[i] - ci1[i][0] for i in range(len(ci1))]
#             y_r2 = [y2[i] - ci2[i][0] for i in range(len(ci2))]
#             y_r3 = [y3[i] - ci3[i][0] for i in range(len(ci3))]
#             plt.bar(r1, y1, width = barWidth, linewidth = 2, label="conv", yerr = y_r1)
#             plt.bar(r2, y2, width = barWidth, linewidth = 2, label="tand", yerr = y_r2)
#             plt.bar(r3, y3, width = barWidth, linewidth = 2, label="div", yerr = y_r3)
#             plt.xticks([r + barWidth / 1.5 for r in range(len(y1))], ["Micro","Oxo","Cfx","RifCfx"])        
#             plt.ylabel('Binding sites\ndensity (kb-1)',fontweight='bold')

#             r = [x[0] for x in [r1,r2,r3]] + [x[1] for x in [r1,r2,r3]] + [x[2] for x in [r1,r2,r3]] + [x[3] for x in [r1,r2,r3]]
#             y = [x[0] for x in [y1,y2,y3]] + [x[1] for x in [y1,y2,y3]] + [x[2] for x in [y1,y2,y3]] + [x[3] for x in [y1,y2,y3]]
#             i = 0
#             for pval in pvals:
#                 s = significance(pval)
#                 barplot_annotate_brackets(i,i+2,s,r,y,dh=0.1)#,dh=0.05, barh=.05, fs=14, dt=0)
#                 i += 3

#             plt.ylim(top=max(y)+6)

#             plt.subplot(312)
#             y1 = scores_conv ; y2 = scores_tand ; y3 = scores_div

#             print "Scores:",y1,y2,y3

#             r1 = range(len(y1))
#             r2 = [x + barWidth for x in r1]
#             r3 = [x + barWidth for x in r2]

#             plt.bar(r1, y1, width = barWidth, linewidth = 2, label="conv")#, yerr = sdsconv)
#             plt.bar(r2, y2, width = barWidth, linewidth = 2, label="tand")#, yerr = sdstand)
#             plt.bar(r3, y3, width = barWidth, linewidth = 2, label="div")#, yerr = sdsdiv)
#             plt.legend(loc="upper left")
#             plt.xticks([r + barWidth / 1.5 for r in range(len(y1))], ["Micro","Oxo","Cfx","RifCfx"])        
#             plt.ylabel('Binding signal\ndensity (kb-1)',fontweight='bold')

#             plt.subplot(313)
#             y1 = seqs_conv ; y2 = seqs_tand ; y3 = seqs_div
#             print "Sequence:",np.asarray(y1)/np.asarray(y3)
#             r1 = range(len(y1))
#             r2 = [x + barWidth for x in r1]
#             r3 = [x + barWidth for x in r2]

#             plt.bar(r1, y1, width = barWidth, linewidth = 2, label="conv")
#             plt.bar(r2, y2, width = barWidth, linewidth = 2, label="tand")
#             plt.bar(r3, y3, width = barWidth, linewidth = 2, label="div")
#             plt.xticks([r + barWidth / 1.5 for r in range(len(y1))], ["Micro","Oxo","Cfx","RifCfx"])        
#             plt.ylabel('Sequence score\ndensity (kb-1)',fontweight='bold')

#             fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
#             fig.set_size_inches(width, height)
#             plt.tight_layout()
#             plt.savefig(pathdb+"/"+titl)
#             plt.close('all')
    
#     if not norm:
#         titl = 'global_{}.png'.format(typ)
#         width = 6 ; height = 4.5
#         fig, ax = plt.subplots()
        
#         plt.subplot(311)
#         barWidth = 0.2
#         y1 = mconv ; y2 = mtand ; y3 = mdiv
#         r1 = range(len(y1))
#         r2 = [x + barWidth for x in r1]
#         r3 = [x + barWidth for x in r2]

#         plt.bar(r1, y1, width = barWidth, linewidth = 2, label="conv", yerr = sdconv)
#         plt.bar(r2, y2, width = barWidth, linewidth = 2, label="tand", yerr = sdtand)
#         plt.bar(r3, y3, width = barWidth, linewidth = 2, label="div", yerr = sddiv)
#         plt.xticks([r + barWidth / 1.5 for r in range(len(y1))], ["Micro","Oxo","Cfx","RifCfx"])        
#         plt.ylabel('sites\n / (total sites)',fontweight='bold')

#         r = [x[0] for x in [r1,r2,r3]] + [x[1] for x in [r1,r2,r3]] + [x[2] for x in [r1,r2,r3]] + [x[3] for x in [r1,r2,r3]]
#         y = [x[0] for x in [y1,y2,y3]] + [x[1] for x in [y1,y2,y3]] + [x[2] for x in [y1,y2,y3]] + [x[3] for x in [y1,y2,y3]]

#         i = 0
#         for pval in pvals:
#             s = significance(pval)
#             #plt.annotate(i,i+2,s,r,y)
#             barplot_annotate_brackets(i,i+2,s,r,y,dh=0.55)#,dh=0.05, barh=.05, fs=14, dt=0)
#             i += 3
        
#         plt.ylim(top=max(y)+0.4)
#         plt.subplot(312)
#         y1 = scores_conv ; y2 = scores_tand ; y3 = scores_div
#         r1 = range(len(y1))
#         r2 = [x + barWidth for x in r1]
#         r3 = [x + barWidth for x in r2]

#         plt.bar(r1, y1, width = barWidth, linewidth = 2, label="conv", yerr = sdsconv)
#         plt.bar(r2, y2, width = barWidth, linewidth = 2, label="tand", yerr = sdstand)
#         plt.bar(r3, y3, width = barWidth, linewidth = 2, label="div", yerr = sdsdiv)
#         plt.legend(loc="upper left")
#         plt.xticks([r + barWidth / 1.5 for r in range(len(y1))], ["Micro","Oxo","Cfx","RifCfx"])        
#         plt.ylabel('N3E',fontweight='bold')

#         plt.subplot(313)
#         y1 = seqs_conv ; y2 = seqs_tand ; y3 = seqs_div
#         r1 = range(len(y1))
#         r2 = [x + barWidth for x in r1]
#         r3 = [x + barWidth for x in r2]

#         plt.bar(r1, y1, width = barWidth, linewidth = 2, label="conv")
#         plt.bar(r2, y2, width = barWidth, linewidth = 2, label="tand")
#         plt.bar(r3, y3, width = barWidth, linewidth = 2, label="div")
#         plt.xticks([r + barWidth / 1.5 for r in range(len(y1))], ["Micro","Oxo","Cfx","RifCfx"])        
#         plt.ylabel('Score (sequence)',fontweight='bold')

#         fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
#         fig.set_size_inches(width, height)
#         plt.tight_layout()
#         plt.savefig(pathdb+"/"+titl)
#         plt.close('all')
