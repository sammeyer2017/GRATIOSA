from Bio.SeqUtils import GC
from scipy import stats
import statsmodels.api as sm
from globvar import *
import numpy as np
import matplotlib.pyplot as plt

def compute_at_content(gen,cond_fc,cond_tss,*arg,**kwargs):

	gen.compute_magic_prom()
	if not hasattr(gen, 'genes_valid'): # if no fc loaded 
		gen.load_fc_pval()
	
	thresh_pval = kwargs.get('thresh_pval', 1.0)
	# wind = kwargs.get('windows', 6)
	# incr = kwargs.get('increment', 1)
	meth = kwargs.get('meth', 'at')
	sp = kwargs.get('sp', False)
	at = {}
	for TSSpos in gen.TSSs[cond_tss].keys(): # for all TSS in cond_tss
		TSSu = gen.TSSs[cond_tss][TSSpos] # single TS
		for sig in TSSu.promoter.keys(): # for all sigma factors binding to this TSS
			try:
				if sig not in at.keys(): # add sigma factor to sigma dict
					at[sig] = {'act':[],'rep':[]}

				at_val = [] # AT content for : -35, spacer, -10, discriminator
				if meth == 'at':
					if sp:
						if len(TSSu.promoter[sig]['spacer']) in spacers:
							at_val.append(100 - GC(TSSu.promoter[sig]['minus35']))		
							at_val.append(100 - GC(TSSu.promoter[sig]['spacer']))		
							at_val.append(100 - GC(TSSu.promoter[sig]['minus10']))		
							at_val.append(100 - GC(TSSu.promoter[sig]['discriminator']))
					else:
						at_val.append(100 - GC(TSSu.promoter[sig]['minus35']))		
						at_val.append(100 - GC(TSSu.promoter[sig]['spacer']))		
						at_val.append(100 - GC(TSSu.promoter[sig]['minus10']))		
						at_val.append(100 - GC(TSSu.promoter[sig]['discriminator']))

				if meth == 'len':
					if sp:
						if len(TSSu.promoter[sig]['spacer']) in spacers:
							at_val.append(len(TSSu.promoter[sig]['minus35']))		
							at_val.append(len(TSSu.promoter[sig]['spacer']))		
							at_val.append(len(TSSu.promoter[sig]['minus10']))		
							at_val.append(len(TSSu.promoter[sig]['discriminator']))
					else:
						at_val.append(len(TSSu.promoter[sig]['minus35']))		
						at_val.append(len(TSSu.promoter[sig]['spacer']))		
						at_val.append(len(TSSu.promoter[sig]['minus10']))		
						at_val.append(len(TSSu.promoter[sig]['discriminator']))

				expr = []
				for gene in TSSu.genes:
					try:
						if gen.genes[gene].fc_pval[cond_fc][1] <= thresh_pval:
							expr.append(gen.genes[gene].fc_pval[cond_fc][0])
					except:
						pass

				if expr != [] and at_val != []:
					if np.mean(expr) < 0:
						at[sig]['rep'].append(at_val)
					elif np.mean(expr) > 0:
						at[sig]['act'].append(at_val)
			except:
				pass
	return at


def draw_results(at):

	title = kwargs.get('title', 'AT content')

	lab = ['minus35','spacer','minus10','discriminator']
	nb_act = len(at[sigfactor]['act'])
	nb_rep = len(at[sigfactor]['rep'])

	at_act = zip(*at[sigfactor]['act'])
	at_rep = zip(*at[sigfactor]['rep'])

	print 'Number of activated promoters :',nb_act,'VS repressed :',nb_rep
	for l,a,r in zip(lab,at_act,at_rep):
		try:
			print l
			print 'Mean',title,'activated :',np.mean(a),'VS repressed :',np.mean(r)
			print 'Student :',stats.ttest_ind(a,r)
			print 'MW != :',stats.mannwhitneyu(a,r,alternative='two-sided')
			print 'MW > :',stats.mannwhitneyu(a,r,alternative='greater')
			print 'MW < :',stats.mannwhitneyu(a,r,alternative='less')
			plt.hist([a,r], label=['Act','Rep'],normed=True)
			plt.legend()
			plt.xlabel(title,fontweight='bold')
			plt.ylabel('Frequency',fontweight='bold')
			plt.title(l+' '+title,fontweight='bold')
			plt.show()
		except:
			pass


    # plt.hist([d[0],d[1]], label=['Biocyc','bTSS'], range = (4,9), bins = 6,normed=True)

	# rep = plt.boxplot([x for x in zip(*at[sigfactor]['rep'])], labels = lab, patch_artist = True,positions=[2,5,8,11])
	# for patch in rep['boxes']:
	# 	patch.set_color('skyblue')
	# 	patch.set_edgecolor('black')
	# for st in at[sigfactor].keys():
	# 	vals = [x for x in zip(*at[sigfactor][st])]
	# 	for val in vals:
	# 		plt.hist(val)
	# 		plt.title(st)
	# plt.xlim(0,12)
	# plt.legend((act, rep),('Act','Rep'))
	# plt.xlabel('Sequence')
	# plt.ylabel('AT content')
	# plt.title('Discriminator')    
    
	# act = plt.boxplot([x for x in zip(*at[sigfactor]['act'])], labels = lab, patch_artist = True,positions=[1,4,7,10])
	# for patch in act['boxes']:
	# 	patch.set_color('lightcoral')
	# 	patch.set_edgecolor('black')
