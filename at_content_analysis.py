from Bio.SeqUtils import GC
from scipy import stats
#import statsmodels.api as sm
from globvar import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

params = {
	'pdf.fonttype': 42,
	'ps.fonttype': 42,
   'axes.labelsize': 22,
   'font.size': 11,
   'font.family':'Arial',
   'legend.fontsize': 9,
   'xtick.labelsize': 11,
   'ytick.labelsize': 11,
   'text.usetex': False,
   'axes.linewidth':1.5, #0.8
   'axes.titlesize':11,
   'axes.spines.top':True,
   'axes.spines.right':True,
   }
plt.rcParams.update(params)

def compute_at_windows(gen,cond_fc,cond_tss,*arg,**kwargs):
	'''
	Compute AT content on each position of promoter on x bp centered windows. Promoter is then
	classified activated, repressed or non regulated based on FC / pvalues data and AT contents 
	between groups are compared. gen : genome object, cond_fc : FC condition, cond_tss : list of TSS used
	'''
	# Requirements
	gen.load_fc_pval()
	gen.load_TSS() 
	gen.load_seq() 
	# kwargs
	thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
	thresh_fc = kwargs.get('thresh_fc', 0) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated
	align = kwargs.get('align',-10) # position where promoters are aligned for comparison
	before = kwargs.get('bef',35) # number of bp to look below TSS position
	after = kwargs.get('aft',10) # number of bp to look above TSS position
	methstat = kwargs.get('methstat','actvsrep') # stat method to compute pvalue between groups. Actvsrep : ttest(activated > repressed). Actvsnone : ttest(act > none)*ttest(rep < none).
	statw = kwargs.get('statw',False) # True : write pval results for each position in a text file 
	draw = kwargs.get('draw', 'CI') # std : draw AT curves with annotated pvalues. If draw == 'CI', draw CI
	org = kwargs.get('org', gen.name) # organism name to use for plot title
	wind2 = kwargs.get('windows',6) # length of windows to compute AT contents
	wind = wind2/2
	shift = kwargs.get('shift',1) # shift of windows, per default one windows computed every nt 

	at = {'all':{'act':[],'rep':[],'none':[]}} # all : TSS aligned on +1 regardless of sigma factor
	# general shape of dict is at[sigma]:{[at contents for activated promoters],[at contents for repressed promoters], [at contents for non regulated promoters]}
	# general shape of [at contents for promoters] : list of lists where each list corresponds to the values of a promoter for each position
	for TSSpos in gen.TSSs[cond_tss].keys(): # for all TSS in cond_tss
		TSSu = gen.TSSs[cond_tss][TSSpos] # single TSS
		at_val = [] # at values of the promoter
		# alignment on TSS
		if TSSu.strand == True:
			for i in range(TSSpos-before,TSSpos+after+1,shift):
				at_val.append(100 - GC(gen.seq[i-1-wind:i+wind]))
		elif TSSu.strand == False:
			for i in range(TSSpos-after,TSSpos+before+1,shift):
				at_val.append(100 - GC(gen.seqcompl[i-1-wind:i+wind]))
			at_val = at_val[::-1]
		# decides if the promoter is activated, repressed or non regulated
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
			if np.mean(expr) < 0 - thresh_fc:
				at['all']['rep'].append(at_val)
			elif np.mean(expr) > 0 + thresh_fc:
				at['all']['act'].append(at_val)
			else:
				if expr_none == []:
					at['all']['none'].append(at_val)

		elif expr_none != [] and at_val != []:
			at['all']['none'].append(at_val)
		
		for sig in TSSu.promoter.keys(): # for all sigma factors binding to this TSS
			try:
				if sig not in at.keys(): # add sigma factor to sigma dict
					at[sig] = {'act':[],'rep':[],'none':[]}

				at_val = [] # AT content for each position
				# AT content for positions -before to +after
				if align == -10:
					ref = (TSSu.promoter[sig]['sites'][0] + TSSu.promoter[sig]['sites'][1])//2
				elif align == -35:
					ref = (TSSu.promoter[sig]['sites'][2] + TSSu.promoter[sig]['sites'][3])//2
				elif align == 0:
					ref = TSSu.pos
				
				if TSSu.strand == True:
					refb = ref - (before + align)
					refe = ref + (after - align)
					for i in range(refb,refe+1,shift):
						at_val.append(100 - GC(gen.seq[i-1-wind:i+wind]))

				elif TSSu.strand == False:
					refb = ref - (after - align)
					refe = ref + (before + align)
					for i in range(refb,refe+1,shift):
						at_val.append(100 - GC(gen.seqcompl[i-1-wind:i+wind]))
					at_val = at_val[::-1]
				
				if expr != [] and at_val != []:
					if np.mean(expr) < 0 - thresh_fc:
						at[sig]['rep'].append(at_val)
					elif np.mean(expr) > 0 + thresh_fc:
						at[sig]['act'].append(at_val)
					else:
						if expr_none == []:
							at[sig]['none'].append(at_val)

				elif expr_none != [] and at_val != []:
					at[sig]['none'].append(at_val)

			except:
				pass

	sigma = kwargs.get('sigma',sigfactor)
	if sigma == 'all':
		align = 0

	titl = '{}-{}-{}-{}-FC{}-PVAL{}'.format(cond_tss,cond_fc,sigma,methstat,thresh_fc,thresh_pval)

	nb_act = len(at[sigma]['act'])
	nb_rep = len(at[sigma]['rep'])
	nb_none = len(at[sigma]['none'])
	# lists of AT contents per position instead of lists of AT contents per promoter
	at_act = zip(*at[sigma]['act'])
	at_rep = zip(*at[sigma]['rep'])
	at_none = zip(*at[sigma]['none'])
	# mean of AT contents per position for act, rep and non regulated promoters
	mact = [np.mean(x) for x in at_act]
	mrep = [np.mean(x) for x in at_rep]
	mnone = [np.mean(x) for x in at_none]
	pvals = []
	if methstat == 'actvsrep':
		for a,r in zip(at_act,at_rep): # H0 : ATact <= ATrep, H1 : ATact > ATrep, t-test(act,rep), pval / 2
			val = stats.ttest_ind(a,r,equal_var=False)[1] / 2
			pvals.append(val)

	elif methstat == 'actvsnone': # H0 : ATact <= ATnone and ATrep >= ATnone, H1 : ATact > ATnone and ATrep < ATnone
		for a,r,n in zip(at_act,at_rep,at_none):
			val = (stats.ttest_ind(a,n,equal_var=False)[1] / 2) * (stats.ttest_ind(n,r,equal_var=False)[1] / 2)
			pvals.append(val)

	if statw == True: # write pvalues in file for each position
		res = open("{}res/jet4/{}-bef{}-aft{}.txt".format(basedir,titl,before,after),'w')
		for val in pvals:
			res.write('{}\n'.format(str(val)))
		res.close()

	pos = [x for x in range(-before,after+1,shift)]

	if draw == 'std':
		fig, ax = plt.subplots()
		act = plt.plot(pos,mact,'rD',linestyle='solid', color='red', markersize=3,label=str(nb_act)+' Activated')
		rep = plt.plot(pos,mrep,'bD',linestyle='solid', color='blue',markersize=3,label=str(nb_rep)+' Repressed')
		try:
			if nb_none > 1:
				plt.plot(pos,mnone,'kD',linestyle='solid', color='black',markersize=3,label=str(nb_none)+' Non affected')
		except:
			pass
		plt.arrow(align,0,0, max(max(mact),max(mrep)),linestyle='dashed',color='gray')

		i = - before
		for val in pvals: 
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
			i += shift

		width = 6
		height = width / 1.618
		fig.set_size_inches(width, height)

		ax.set_ylabel('AT content (%)',fontweight='bold')
		ax.set_xlabel('Position (nt)',fontweight='bold')
		ax.set_xlim(-before,after)
		ax.set_title(titl,fontweight='bold')
		fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
		plt.tight_layout()

		if sigma == 'all':
			leg = 'All Promoters'
		else:
			leg = '$\sigma {}$ Promoters'.format(sigma.replace('sigma',''))
		
		ax.legend(title=leg,ncol=1,loc='upper left')
		plt.savefig("/home/raphael/Documents/topo/results/stage/{}-bef{}-aft{}.svg".format(titl,before,after),transparent=False)
		plt.close('all')

	if draw == 'CI':
		fig, ax = plt.subplots()		
		ciact = [(np.mean(x)-np.std(x)/np.sqrt(len(x)),np.mean(x)+np.std(x)/np.sqrt(len(x))) for x in at_act]
		cirep = [(np.mean(x)-np.std(x)/np.sqrt(len(x)),np.mean(x)+np.std(x)/np.sqrt(len(x))) for x in at_rep]
		cinone = [(np.mean(x)-np.std(x)/np.sqrt(len(x)),np.mean(x)+np.std(x)/np.sqrt(len(x))) for x in at_none]

		plt.fill_between(pos,[x[0] for x in ciact], [x[1] for x in ciact], facecolor = 'red', alpha = 0.5, label=str(nb_act)+' Activated',linewidth=0)
		plt.fill_between(pos,[x[0] for x in cirep], [x[1] for x in cirep], facecolor = 'blue', alpha = 0.5, label=str(nb_rep)+' Repressed',linewidth=0)

		plt.plot(pos,mact, linestyle = 'solid', color = 'red', linewidth = 1.5, alpha=0.8)
		plt.plot(pos,mrep, linestyle = 'solid', color = 'blue', linewidth = 1.5, alpha=0.8)
		try:
			if nb_none > 1:
				plt.plot(pos,mnone, linestyle = 'solid', color = 'black', label=str(nb_none)+' Non affected', linewidth = 1.5)
		except:
			pass

		width = 10
		height = width / 1.618
		fig.set_size_inches(width, height)

		ax.set_ylabel('AT content (%)',fontweight='bold')
		ax.set_xlabel('Position (nt)',fontweight='bold')
		ax.set_xlim(-before,after)
		ax.set_title(org,fontweight='bold')
		fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
		plt.tight_layout()

		if sigma == 'all':
			leg = 'All Promoters'
		else:
			leg = '$\sigma {}$ Promoters'.format(sigma.replace('sigma',''))
		
		ax.legend(title=leg,ncol=1,loc='upper left')
		plt.arrow(align,0,0, max(max(mact),max(mrep)),linestyle='dashed',color='gray')
		plt.savefig("/home/raphael/Documents/topo/results/stage/{}-bef{}-aft{}.svg".format(titl,before,after),transparent=False)
		plt.close('all')

def compute_at_spacer(gen,cond_tss,*arg,**kwargs):

	gen.compute_magic_prom()
	sig = kwargs.get('sig', 'sigma70') # below, gene considered valid, above, gene considered non regulated
	wind2 = kwargs.get('windows',6) # length of windows to compute AT contents
	wind = wind2//2
	align = kwargs.get('align',-10) # position where promoters are aligned for comparison
	before = kwargs.get('bef',50) # number of bp to look below TSS position
	after = kwargs.get('aft',20) # number of bp to look above TSS position

	at = {}
	for sp in spacers:
		at[sp] = []
	
	for TSSpos in gen.TSSs[cond_tss].keys(): # for all TSS in cond_tss
		TSSu = gen.TSSs[cond_tss][TSSpos] # single TS
		try:
			at_val = [] # AT content for each position
			# AT content for positions -before to +after
			if align == -10:
				ref = (TSSu.promoter[sig]['sites'][0] + TSSu.promoter[sig]['sites'][1])//2
			elif align == -35:
				ref = (TSSu.promoter[sig]['sites'][2] + TSSu.promoter[sig]['sites'][3])//2
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
			
			sp = len(TSSu.promoter[sig]['spacer'])
			at[sp].append(at_val)

		except Exception as e:
			pass	

	titl = '{},{}'.format(cond_tss,align)

	pos = [x for x in range(-before,after+1)]

	fig, ax = plt.subplots()

	cols = ['cyan','purple','blue','green','red']
	i = 0
	for sp in spacers:
		vals = zip(*at[sp])
		mvals = [np.mean(x) for x in vals]
		plt.plot(pos,mvals,'kD',linestyle='solid', color=cols[i], markersize=3,label='{} bp ({} prom.)'.format(sp,len(at[sp])))
		i += 1

	plt.arrow(align,0,0, 75,linestyle='dashed',color='gray')
	width = 6
	height = width / 1.618
	fig.set_size_inches(width, height)

	ax.legend()
	ax.set_ylabel('AT content (%)',fontweight='bold')
	ax.set_xlabel('Position (nt)',fontweight='bold')
	ax.set_xlim(-35,+10)
	ax.set_title(cond_tss,fontweight='bold')
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
	plt.tight_layout()
	plt.savefig("{}res/{}-{}_{}.svg".format(basedir,titl,before,after),transparent=False)
	plt.close('all')



def compute_at_discr(gen,*arg,**kwargs):
	'''
	'''

	gen.load_fc_pval()
	gen.load_TSS() 
	gen.load_seq() 

	thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
	thresh_fc = kwargs.get('thresh_fc', 0) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated
	wind = kwargs.get('wind', 6)
	cond_tss = kwargs.get('cond_tss', 'biocyc')
	cond_fc = kwargs.get('cond_fc', 'Blot')

	ATactref = kwargs.get('actref', 0.08)
	ATrepref = kwargs.get('repref', 0.02)
	ATref = 100 - GC(gen.seq)
	ATact = []
	ATrep = []

	for TSS in gen.TSSs[cond_tss].keys():
		try:
			TSSu = gen.TSSs[cond_tss][TSS]
			if TSSu.strand == True:
				AT = 100 - GC(gen.seq[TSSu.promoter[sigfactor]['sites'][1]:TSSu.promoter[sigfactor]['sites'][1]+wind])
			elif TSSu.strand == False:
				AT = 100 - GC(gen.seq[TSSu.promoter[sigfactor]['sites'][0]-wind-1:TSSu.promoter[sigfactor]['sites'][0]-1])

			expr = []
			for gene in TSSu.genes:
				try:
					if gen.genes[gene].fc_pval[cond_fc][1] <= thresh_pval:
						expr.append(gen.genes[gene].fc_pval[cond_fc][0])
				except:
					pass

			if expr != []:
				if np.mean(expr) < 0 - thresh_fc:
					ATact.append(AT)
				elif np.mean(expr) > 0 + thresh_fc:
					ATrep.append(AT)
		except:
			pass

	ATs = ATact + ATrep
	titl = '{}_{}_FC{}_PVAL{}'.format(cond_tss,cond_fc,str(thresh_fc),str(thresh_pval))

	fig = plt.figure()
	fig.set_size_inches(8,10)

	plt.subplot(321)
	plt.hist(ATs,normed=True,align='left',edgecolor='black',color='lightgray')
	plt.title('All promoters', fontweight='bold')
	plt.xlim(0,100)
	
	plt.subplot(323)
	plt.hist(ATact,normed=True,align='left',edgecolor='black',color='red')
	plt.title('Activated promoters', fontweight='bold')
	plt.xlim(0,100)

	plt.subplot(325)
	plt.hist(ATrep,normed=True,align='left',edgecolor='black',color='blue')
	plt.ylabel('Proportion', fontweight='bold')
	plt.xlabel('AT %', fontweight='bold')
	plt.title('Repressed promoters', fontweight='bold')
	plt.xlim(0,100)

#####

	plt.subplot(322)
	plt.hist(np.random.binomial(len(ATs),ATref/100),normed=True,align='left',edgecolor='black',color='lightgray')
	plt.title('All promoters', fontweight='bold')
	
	plt.subplot(324)
	plt.hist(np.random.binomial(len(ATact),ATactref),normed=True,align='left',edgecolor='black',color='red')
	plt.title('Activated promoters', fontweight='bold')

	plt.subplot(326)
	plt.hist(np.random.binomial(len(ATrep),ATrepref),normed=True,align='left',edgecolor='black',color='blue')
	plt.ylabel('Proportion', fontweight='bold')
	plt.xlabel('AT %', fontweight='bold')
	plt.title('Repressed promoters', fontweight='bold')
	
	plt.savefig("{}res/jet4/{}.svg".format(basedir,titl),transparent=False)
	plt.close('all')

def write_genes(gen,cond_fc):
	act = open(basedir+'data/genes_act.txt','w')
	rep = open(basedir+'data/genes_rep.txt','w')
	for gene in gen.genes.keys():
		try:
			if gen.genes[gene].fc_pval[cond_fc][1] <= 0.05:
				if gen.genes[gene].fc_pval[cond_fc][0] > 0:
					act.write(gene+'\n')
				else:
					rep.write(gene+'\n')
		except:
			pass
	
	act.close()
	rep.close()
