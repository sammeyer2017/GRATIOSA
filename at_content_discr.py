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


def compute_at_windows(gen,cond_fc,cond_tss,*arg,**kwargs):

	gen.compute_magic_prom()
	if not hasattr(gen, 'genes_valid'): # if no fc loaded 
		gen.load_fc_pval()
	if not hasattr(gen, 'seq'): # if no fc loaded 
		gen.load_seq()

	thresh_pval = kwargs.get('thresh_pval', 1.0)
	wind2 = kwargs.get('windows',6)
	wind = wind2/2
	align = kwargs.get('align',-10)
	before = kwargs.get('bef',50)
	after = kwargs.get('aft',20)
	at = {}
	for TSSpos in gen.TSSs[cond_tss].keys(): # for all TSS in cond_tss
		TSSu = gen.TSSs[cond_tss][TSSpos] # single TS
		for sig in TSSu.promoter.keys(): # for all sigma factors binding to this TSS
			try:
				if sig not in at.keys(): # add sigma factor to sigma dict
					at[sig] = {'act':[],'rep':[],'none':[],'lact':[],'lrep':[], 'lnone':[]}

				at_val = [] # AT content for each position
				# AT content for positions -before to +after
				if align == -10:
					ref = (TSSu.promoter[sig]['sites'][0] + TSSu.promoter[sig]['sites'][1])/2
				elif align == -35:
					ref = (TSSu.promoter[sig]['sites'][2] + TSSu.promoter[sig]['sites'][3])/2
				elif align == 1:
					ref = TSSu.pos
				
				if TSSu.strand == True:
					refb = ref - (before + align)
					refe = ref + (after - align)
					for i in range(refb,refe+1):
						at_val.append(100 - GC(gen.seq[i-1-wind:i+wind]))

				elif TSSu.strand == False:
					refb = ref - (after - align)
					refe = ref + (before + align)
					for i in range(refb,refe+1):
						at_val.append(100 - GC(gen.seqcompl[i-1-wind:i+wind]))
					at_val = at_val[::-1]

				expr = [] ; expr_none = []
				for gene in TSSu.genes:
					try:
						if gen.genes[gene].fc_pval[cond_fc][1] <= thresh_pval:
							expr.append(gen.genes[gene].fc_pval[cond_fc][0])
						else:
							expr_none.append(gen.genes[gene].fc_pval[cond_fc][0])
					except:
						pass
				
				if expr != [] and at_val != []:
					if np.mean(expr) < 0:
						at[sig]['rep'].append(at_val)
						at[sig]['lrep'].append(TSSu.pos)
					elif np.mean(expr) > 0:
						at[sig]['act'].append(at_val)
						at[sig]['lact'].append(TSSu.pos)
				
				elif expr_none != [] and at_val != []:
					at[sig]['none'].append(at_val)
					at[sig]['lnone'].append(TSSu.pos)

			except Exception as e:
				print e
				pass

	titl = '{}, {} et al., p-value = {}'.format(cond_tss,cond_fc,thresh_pval)

	sigma = kwargs.get('sigma',sigfactor)
	nb_act = len(at[sigma]['act'])
	nb_rep = len(at[sigma]['rep'])
	nb_none = len(at[sigma]['none'])

	at_act = zip(*at[sigma]['act'])
	at_rep = zip(*at[sigma]['rep'])
	at_none = zip(*at[sigma]['none'])

	mact = [np.mean(x) for x in at_act]
	mrep = [np.mean(x) for x in at_rep]
	mnone = [np.mean(x) for x in at_none]

	pos = [x for x in range(-before,after+1)]
	print nb_none
	draw = kwargs.get('draw', True)	
	if draw == True:
		act = plt.plot(pos,mact,'rD',linestyle='solid', color='red', markersize=3,label=str(nb_act)+' Activated')
		rep = plt.plot(pos,mrep,'bD',linestyle='solid', color='blue',markersize=3,label=str(nb_rep)+' Repressed')
		if nb_none != 0:
			plt.plot(pos,mnone,'kD',linestyle='dashed', color='black',markersize=3,label=str(nb_none)+' Non')
		plt.xlabel('position',fontweight='bold')
		plt.ylabel('AT content',fontweight='bold')
		plt.title(titl)
		plt.legend(title= 'Promoters', ncol = 1, fontsize='medium')
		plt.arrow(align,0,0, max(max(mact),max(mrep)),linestyle='dashed',color='gray')

		pval = []
		for a,r in zip(at_act,at_rep):
			pval.append(stats.ttest_ind(a,r)[1])
		i = - before
		for val in pval:
			if val <= 0.001:
				s = '***'
			elif val <= 0.01:
				s = '**' 
			elif val <= 0.05:
				s = '*' 
			else:
				s = ''
			if s != '':
				plt.text(i,min(min(mact),min(mrep))+1,s,fontweight='bold')
			i += 1

    	plt.savefig("{}res/jet3/{}-{}_{}_{}.svg".format(basedir,titl,before,after,sigma),transparent=False)
    	plt.close('all')
	# spc = kwargs.get('spc',False)
	# if spc == True:
	# 	spact = [item for sublist in at_act[20:41] for item in sublist]
	# 	sprep = [item for sublist in at_rep[20:41] for item in sublist]
	# 	disact = [item for sublist in at_act[40:51] for item in sublist]
	# 	disrep = [item for sublist in at_rep[40:51] for item in sublist]
	# 	print 'Spacer :',stats.ttest_ind(spact,sprep)[1]
	# 	print 'Discriminator :',stats.ttest_ind(disact,disrep)[1]

def corr_at_content(gen,cond_fc,cond_tss,*arg,**kwargs):

	gen.compute_magic_prom()
	if not hasattr(gen, 'genes_valid'): # if no fc loaded 
		gen.load_fc_pval()
	
	thresh_pval = kwargs.get('thresh_pval', 1.0)
	# wind = kwargs.get('windows', 6)
	# incr = kwargs.get('increment', 1)
	at = {}
	for TSSpos in gen.TSSs[cond_tss].keys(): # for all TSS in cond_tss
		TSSu = gen.TSSs[cond_tss][TSSpos] # single TS
		for sig in TSSu.promoter.keys(): # for all sigma factors binding to this TSS
			try:
				if sig not in at.keys(): # add sigma factor to sigma dict
					at[sig] = {'act':[],'rep':[]}

				at_val = [] # AT content for : -35, spacer, -10, discriminator
				len_val = []
				at_val.append(100 - GC(TSSu.promoter[sig]['minus35']))
				at_val.append(100 - GC(TSSu.promoter[sig]['spacer']))		
				at_val.append(100 - GC(TSSu.promoter[sig]['minus10']))	
				at_val.append(100 - GC(TSSu.promoter[sig]['discriminator']))		

				len_val.append(len(TSSu.promoter[sig]['minus35']))	
				len_val.append(len(TSSu.promoter[sig]['spacer']))	
				len_val.append(len(TSSu.promoter[sig]['minus10']))
				len_val.append(len(TSSu.promoter[sig]['discriminator']))
				
				expr = []
				for gene in TSSu.genes:
					try:
						if gen.genes[gene].fc_pval[cond_fc][1] <= thresh_pval:
							expr.append(gen.genes[gene].fc_pval[cond_fc][0])
					except:
						pass

				if expr != [] and at_val != []:
					if np.mean(expr) < 0:
						at[sig]['rep'].append([at_val,len_val])
					elif np.mean(expr) > 0:
						at[sig]['act'].append([at_val,len_val])
			except:
				pass
	return at

def corr_results(at,*arg,**kwargs):

	titl = kwargs.get('title','')

	nb_act = len(at[sigfactor]['act'])
	nb_rep = len(at[sigfactor]['rep'])

	at_act = zip(*[x[0] for x in at[sigfactor]['act']])
	at_rep = zip(*[x[0] for x in at[sigfactor]['rep']])
	len_act = zip(*[x[1] for x in at[sigfactor]['act']])
	len_rep = zip(*[x[1] for x in at[sigfactor]['rep']])

	print stats.pearsonr(at_act[1],at_act[3])
	plt.plot(at_act[1],at_act[3], 'rD',markersize=2,linestyle='None')
	plt.title(stats.pearsonr(at_act[1],at_act[3]))
	plt.xlabel('AT spacer')
	plt.ylabel('AT discr')
	# plt.show()
	print stats.pearsonr(at_rep[1],at_rep[3])
	plt.plot(at_rep[1],at_rep[3], 'bD',markersize=2,linestyle='None')
	plt.title(stats.pearsonr(at_rep[1],at_rep[3]))
	plt.xlabel('AT spacer')
	plt.ylabel('AT discr')
	# plt.show()

	print stats.pearsonr(len_act[1],len_act[3])
	plt.plot(len_act[1],len_act[3], 'rD',markersize=2,linestyle='None')
	plt.title(stats.pearsonr(len_act[1],len_act[3]))
	plt.xlabel('len spacer')
	plt.ylabel('len discr')
	# plt.show()

	print stats.pearsonr(len_rep[1],len_rep[3])
	plt.plot(len_rep[1],len_rep[3], 'bD',markersize=2,linestyle='None')
	plt.title(stats.pearsonr(len_rep[1],len_rep[3]))
	plt.xlabel('len spacer')
	plt.ylabel('len discr')
	# plt.show()

	# for a,r in zip(at_act,at_rep):
	# 	print stats.pearsonr(a,r)
	# 	plt.plot(a,r, linestyle='None')
	# 	plt.show()
	# 	plt.plot()
	# if draw == True:
	# 	act = plt.plot(pos,mact,'rD',linestyle='solid', color='red', markersize=3,label=str(nb_act)+' Activated')
	# 	rep = plt.plot(pos,mrep,'bD',linestyle='solid', color='blue',markersize=3,label=str(nb_rep)+' Repressed')
	# 	plt.xlabel('position',fontweight='bold')
	# 	plt.ylabel('AT content',fontweight='bold')
	# 	plt.title(titl)
	# 	plt.legend(title= 'Promoters', ncol = 1, fontsize='medium')
	# 	plt.arrow(align,0,0, max(max(mact),max(mrep)),linestyle='dashed',color='gray')

	# 	if stat == True:
	# 		i = -50
	# 		for val in pval:
	# 			if val <= 0.001:
	# 				s = '***'
	# 			elif val <= 0.01:
	# 				s = '**' 
	# 			elif val <= 0.05:
	# 				s = '*' 
	# 			else:
	# 				s = ''

	# 			if s != '':# and i > -40 and i<1:
	# 				plt.text(i,min(min(mact),min(mrep))+1,s,fontweight='bold')

	# 			i += 1


