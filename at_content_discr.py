from Bio.SeqUtils import GC
from scipy import stats
import statsmodels.api as sm
from globvar import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

params = {
   'axes.labelsize': 11,
   'font.size': 11,
   'legend.fontsize': 11,
   'xtick.labelsize': 11,
   'ytick.labelsize': 11,
   'text.usetex': False,
   'axes.linewidth':1.5, #0.8
   'axes.spines.top':True,
   'axes.spines.right':True,
   }

plt.rcParams.update(params)

def compute_at_windows(gen,cond_fc,cond_tss,*arg,**kwargs):
	'''
	Compute AT content on each position of promoter on x bp centered windows. Promoter is then classified activated, repressed or non regulated based on FC data and AT contents between groups are compared.
	'''
	gen.load_fc_pval()
	gen.load_TSS() 
	gen.load_seq() 

	thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
	thresh_fc = kwargs.get('thresh_fc', 0) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated
	align = kwargs.get('align',-10) # position where promoters are aligned for comparison
	before = kwargs.get('bef',50) # number of bp to look below TSS position
	after = kwargs.get('aft',20) # number of bp to look above TSS position
	methstat = kwargs.get('methstat','actvsrep') # stat method to compute pvalue between groups. Actvsrep : ttest(activated > repressed). Actvsnone : ttest(act > none)*ttest(rep < none).
	statw = kwargs.get('statw',True) # write pval results 
	draw = kwargs.get('draw', True) # draw AT curves. If draw == 'CI', draw CI
	org = kwargs.get('org', gen.name) # real organism name for plot title

	wind2 = kwargs.get('windows',6) # length of windows to compute AT contents
	wind = wind2/2

	at = {'all':{'act':[],'rep':[],'none':[]}} # all TSS aligned on +1 regardless of sigma factor
	# general shape of dict is at[sigma]:{[at contents for activated promoters],[at contents for repressed promoters], [at contents for non regulated promoters]}
	for TSSpos in gen.TSSs[cond_tss].keys(): # for all TSS in cond_tss
		TSSu = gen.TSSs[cond_tss][TSSpos] # single TS
		at_val = []
		if TSSu.strand == True:
			for i in range(TSSpos-before-1,TSSpos+after-1+1):
				at_val.append(100 - GC(gen.seq[i-1-wind:i+wind]))
		elif TSSu.strand == False:
			for i in range(TSSpos-after+1,TSSpos+before+1+1):
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
		align = +1

	titl = '{}, {} et al., {}, {}'.format(cond_tss,cond_fc,sigma,methstat)

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
		res = open("{}res/jet4/{}-{}_{}_{}.txt".format(basedir,titl,before,after,sigma),'w')
		for val in pvals:
			res.write('{}\n'.format(str(val)))
		res.close()

	pos = [x for x in range(-before,after+1)]

	if draw == True:
		act = plt.plot(pos,mact,'rD',linestyle='solid', color='red', markersize=3,label=str(nb_act)+' Activated')
		rep = plt.plot(pos,mrep,'bD',linestyle='solid', color='blue',markersize=3,label=str(nb_rep)+' Repressed')
		try:
			plt.plot(pos,mnone,'kD',linestyle='solid', color='black',markersize=3,label=str(nb_none)+' Non affected')
		except:
			pass
		plt.title(titl, fontweight='bold')
		plt.legend(title='Promoters', ncol = 1, fontsize='medium')
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
			i += 1

    	plt.savefig("{}res/jet4/{}-{}_{}_{}.svg".format(basedir,titl,before,after,sigma),transparent=False)
    	plt.close('all')

	if draw == 'CI':
		
		ciact = [(np.mean(x)-np.std(x)/np.sqrt(len(x)),np.mean(x)+np.std(x)/np.sqrt(len(x))) for x in at_act]
		cirep = [(np.mean(x)-np.std(x)/np.sqrt(len(x)),np.mean(x)+np.std(x)/np.sqrt(len(x))) for x in at_rep]
		cinone = [(np.mean(x)-np.std(x)/np.sqrt(len(x)),np.mean(x)+np.std(x)/np.sqrt(len(x))) for x in at_none]

		width = 5.5
		height = width / 1.618
		fig, ax = plt.subplots()

		plt.fill_between(pos,[x[0] for x in ciact], [x[1] for x in ciact], facecolor = 'red', alpha = 0.5, label=str(nb_act)+' Activated',linewidth=0)
		plt.fill_between(pos,[x[0] for x in cirep], [x[1] for x in cirep], facecolor = 'blue', alpha = 0.5, label=str(nb_rep)+' Repressed',linewidth=0)

		plt.plot(pos,mact, linestyle = 'solid', color = 'red', linewidth = 1.5, alpha=0.8)
		plt.plot(pos,mrep, linestyle = 'solid', color = 'blue', linewidth = 1.5, alpha=0.8)
		try:
			plt.plot(pos,mnone, linestyle = 'solid', color = 'black', label=str(nb_none)+' Non affected', linewidth = 1.5)
		except:
			pass

		ax.set_ylabel('AT content (%)',fontweight='bold')
		ax.set_xlabel('Position (nt)',fontweight='bold')
		ax.set_xlim(-35,10)
		ax.set_title(org,fontweight='bold')
		fig.set_size_inches(width, height)
		fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)

		if sigma == 'all':
			leg = 'All Promoters'
		else:
			leg = '$\sigma{}$ Promoters'.format(sigma.replace('sigma',''))

		plt.legend(title=leg, ncol = 1, loc='upper left')
		plt.arrow(align,0,0, max(max(mact),max(mrep)),linestyle='dashed',color='gray')
		plt.savefig("{}res/jet4/{}-{}_{}_{}.svg".format(basedir,titl,before,after,sigma),transparent=False)
		plt.close('all')

		# plt.arrow(-4,0,0, max(max(mact),max(mrep)),linestyle='solid',color='black',linewidth =2)
		# plt.arrow(1,0,0, max(max(mact),max(mrep)),linestyle='solid',color='black',linewidth =2)
		# rec = patches.Rectangle((-2,65),6,3.5,fill=True, facecolor='lightgray', edgecolor='black',hatch='\\')
		# ax.add_patch(rec)

