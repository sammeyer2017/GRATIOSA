from Bio.SeqUtils import GC
from scipy import stats
import statsmodels.api as sm
from globvar import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from useful_functions import *
from collections import Counter
from statsmodels.stats import proportion

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

def compute_at_windows(gen,cond_fc,cond_tss,*arg,**kwargs):
	'''
	Compute AT content on each position of promoter on x bp centered windows. Promoter is then
	classified activated, repressed or non regulated based on FC / pvalues data and AT contents 
	between groups are compared. gen : genome object, cond_fc : FC condition, cond_tss : list of TSS used
	'''
	# requirements
	#gen.load_TSS()
	gen.load_seq()
	if not hasattr(gen,'genes_valid'):
            gen.load_fc_pval()
	if not hasattr(gen,'TSSs'):
            gen.load_TSS()	#gen.load_fc_pval()

	# kwargs
	thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
	thresh_fc = kwargs.get('thresh_fc', 0) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated
	align = kwargs.get('align',-10) # position where promoters are aligned for comparison ; either -10 element, -35 or TSS
	before = kwargs.get('bef',35) # number of bp to look below TSS position
	after = kwargs.get('aft',10) # number of bp to look above TSS position
	org = kwargs.get('org', gen.name) # organism name to use for plot title

	wind2 = kwargs.get('windows',4) # length of windows to compute AT contents
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
		# Classifies promoter act, rep or non affected depending on gene expresion values
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
	# if sigma = all, aligns on TSS
	if sigma == 'all':
		align = 0

	pathdb = '{}data/{}/discriminator'.format(basedir,gen.name)
	if not os.path.exists(pathdb):
		os.makedirs(pathdb)

	titl = '{}-{}-{}-FC{}-PVAL{}'.format(cond_tss,cond_fc.replace('/',''),sigma,thresh_fc,thresh_pval)

	# Number of act rep and non affected promoters 
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
	pos = [x for x in range(-before,after+1,shift)]

	methstat = kwargs.get('methstat','actvsrep') # stat method to compute pvalue between groups. Actvsrep : ttest(activated > repressed). Actvsnone : ttest(act > none)*ttest(rep < none).
	statw = kwargs.get('statw',True) # True : write pval results for each position in a text file 
	if statw == True: # write pvalues in file for each position
		if methstat == 'actvsrep':
			for a,r in zip(at_act,at_rep): # H0 : ATact <= ATrep, H1 : ATact > ATrep, t-test(act,rep), pval / 2
				val = stats.ttest_ind(a,r,equal_var=False)[1] / 2
				pvals.append(val)

		# elif methstat == 'actvsnone': # H0 : ATact <= ATnone and ATrep >= ATnone, H1 : ATact > ATnone and ATrep < ATnone
		# 	for a,r,n in zip(at_act,at_rep,at_none):
		# 		val = (stats.ttest_ind(a,n,equal_var=False)[1] / 2) * (stats.ttest_ind(n,r,equal_var=False)[1] / 2)
		# 		pvals.append(val)

		res = open(pathdb+"/{}.txt".format(titl),'w')
		res.write('position\tp-value unilateral t-test\n')
		for p,val in zip(pos,pvals):
			res.write('{}\t{}\n'.format(str(p),str(val)))
		res.close()

	draw = kwargs.get('draw', 'CI') # std : draw AT curves with annotated pvalues. If draw == 'CI', draw CI
	if draw == 'CI':
		fig, ax = plt.subplots()		
		fact = 1 
		ciact = [(np.mean(x)-fact*(np.std(x)/np.sqrt(len(x))),np.mean(x)+fact*(np.std(x)/np.sqrt(len(x)))) for x in at_act]
		cirep = [(np.mean(x)-fact*(np.std(x)/np.sqrt(len(x))),np.mean(x)+fact*(np.std(x)/np.sqrt(len(x)))) for x in at_rep]
		cinone = [(np.mean(x)-fact*(np.std(x)/np.sqrt(len(x))),np.mean(x)+fact*(np.std(x)/np.sqrt(len(x)))) for x in at_none]

		plt.fill_between(pos,[x[0] for x in ciact], [x[1] for x in ciact], facecolor = 'red', alpha = 0.5, label=str(nb_act)+' Activated',linewidth=0)
		plt.fill_between(pos,[x[0] for x in cirep], [x[1] for x in cirep], facecolor = 'blue', alpha = 0.5, label=str(nb_rep)+' Repressed',linewidth=0)

		plt.plot(pos,mact, linestyle = 'solid', color = 'red', linewidth = 1.5, alpha=0.8)
		plt.plot(pos,mrep, linestyle = 'solid', color = 'blue', linewidth = 1.5, alpha=0.8)
		try:
			if nb_none > 20:
				plt.plot(pos,mnone, linestyle = 'solid', color = 'black', label=str(nb_none)+' Non affected', linewidth = 1.5)
		except:
			pass

		ax.set_ylabel('AT content (%)',fontweight='bold')
		ax.set_xlabel('Position (nt)',fontweight='bold')
		ax.set_xlim(-before,after)
		ax.set_title(org+'_{}'.format(cond_fc),fontweight='bold')
		fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
		plt.tight_layout()

		if sigma == 'all':
			leg = 'All Promoters'
		else:
			leg = '$\sigma {}$ Promoters'.format(sigma.replace('sigma',''))
		
		ax.legend(title=leg,ncol=1,loc='upper left')
		#plt.arrow(align,0,0, max(max(mact),max(mrep)),linestyle='dashed',color='gray')

		width = 5 ; height = width / 1.618
		fig.set_size_inches(width, height)
		plt.tight_layout()
		plt.savefig(pathdb+"/{}.svg".format(titl),transparent=False)
		plt.close('all')

	fullplots = kwargs.get('fullplots',True)
	# Creates additional graphes: 
	#	% activated promoters for each possible AT content for given windows for given position, default -2 (typical )
	#	Comparison of mean AT contents between activated, non affected and repressed promoters for given position
	locs = kwargs.get('locs',[-2]) # positions where AT contents are compared
	if fullplots:
		for loc in locs: # for each position to be compared
			idx = before + loc # index of at_act list corresponding to loc position
			a = Counter(at_act[idx]) # dict of nb of occurences of each element
			a['idx'] = 'act'
			r = Counter(at_rep[idx])
			r['idx'] = 'rep'
			if nb_none > 20:
				n = Counter(at_none[idx])
				n['idx'] = 'non'
			else:
				n = {}
			res_loc = {}	
			for d in [a,r,n]:
				for k in d.keys():
					if k != 'idx':
						if k not in res_loc.keys():
							res_loc[k] = {'act':0,'rep':0,'non':0}
						res_loc[k][d['idx']] += d[k]

			x = [] ; y = [] ; stds = [] ; cim = []
			for k in res_loc.keys():
				tot = res_loc[k]['act']+res_loc[k]['rep']
				if tot > 7:
					pexp = float(res_loc[k]['act']) / float(tot)
					std = np.array(stats.binom.std(tot, pexp, loc=0))/tot
					ci = proportion.proportion_confint(res_loc[k]['act'], tot, alpha=0.05, method='normal')
					cim.append((ci[1] - ci[0])/2)

					x.append(k)
					y.append(pexp)
					stds.append(std)

			X = sm.add_constant(x)
			wls_model = sm.WLS(y, X, weights=1/np.power(np.array(stds),2)) # weights proportional to the inverse of std2
			results = wls_model.fit()
			slope = results.params[1]
			OR = results.params[0]
			pval = str(round(results.pvalues[1],3))
			spearman = str(round(stats.spearmanr(x,y)[1],3))
			pearson = str(round(stats.pearsonr(x,y)[1],3))
			width = 3 ; height = width / 1.4
			fig, ax = plt.subplots()
			plt.plot(x, y, marker='o', color='black',markersize=6, linestyle='None', fillstyle='none')
			plt.errorbar(x, y,yerr=cim,mec='black', capsize=5, elinewidth=1,mew=1,linestyle='None', color='black')
			plt.plot([-20]+x+[120],np.array([-20]+x+[120])*slope + OR, color='black',linewidth=0.5)   
			ax.set_ylabel("proportion of\nactivated promoters")
			ax.set_xlabel('AT %')
			ax.set_xlim(-10,110)
			ax.set_title("Reg: {}, Pear.: {}, Spear.: {}".format(pval,spearman,pearson),fontsize=9)
			fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
			fig.set_size_inches(width, height)
			plt.tight_layout()
			plt.savefig(pathdb+"/"+titl+"_pos"+str(loc)+"_reg.svg",transparent=False)
			plt.close('all')

			print x,y,cim
			fact = 1.96
			if nb_none > 20:
				labs = ["act","non","rep"]
				yval = [np.mean(at_act[idx]),np.mean(at_none[idx]),np.mean(at_rep[idx])]
				stdval = [fact*np.std(at_act[idx])/np.sqrt(len(at_act[idx])),fact*np.std(at_none[idx])/np.sqrt(len(at_none[idx])),fact*np.std(at_rep[idx])/np.sqrt(len(at_rep[idx]))]
				xval = [0,1,2]
				col = ['red','black','blue']
			else:
				labs = ["act","rep"]
				yval = [np.mean(at_act[idx]),np.mean(at_rep[idx])]
				stdval = [fact*np.std(at_act[idx])/np.sqrt(len(at_act[idx])),fact*np.std(at_rep[idx])/np.sqrt(len(at_rep[idx]))]
				xval = [0,1]
				col = ['red','blue']

			s = significance(pvals[idx])
			print pvals[idx],idx
			if s == 'ns':
				fs = 14
				dt = 0.01
			else:
				fs = 18
				dt = 0

			fig_width = 2 ; fig_height = 2.2
			fig = plt.figure(figsize=(fig_width,fig_height))
			plt.bar(xval, yval, yerr = stdval, width = 0.6, color = 'white',edgecolor = 'black', linewidth = 2, ecolor = 'black', capsize = 3.5)
			plt.xticks(xval,labs,rotation = 45)
			plt.xlim(xval[0]-0.5,xval[-1]+0.5)
			if len(xval) == 2:
				barplot_annotate_brackets(0,1, s, xval,yval,yerr=stdval,fs=fs, dt=dt,barh=0.01)
			else:
				barplot_annotate_brackets(0,2, s, xval,yval,yerr=stdval,fs=fs, dt=dt,barh=0.01)

			plt.ylabel("AT %")
			ylim1 = min(yval) - 5
			ylim2 = max(yval) + 10
			plt.ylim(ylim1,ylim2)

			fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
			#plt.title('{}'.format(cond_fc),fontsize=11)
			plt.tight_layout()
			plt.savefig(pathdb+"/"+titl+"_pos"+str(loc)+"_group.svg",transparent=False)
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


def compute_at_evolution2(gen,cond_fc,cond_tss,*arg,**kwargs):
	'''
	'''

	thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
	thresh_fc = kwargs.get('thresh_fc', 0) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated
	
	at = {'act':{'FC':[],'AT':[]},'rep':{'FC':[],'AT':[]},'non':{'FC':[],'AT':[]}}
	at_vals = {}
	for TSS in gen.TSSs[cond_tss].keys():
		TSSu = gen.TSSs[cond_tss][TSS]
		expr = [] ; expr_none = [] ; at_val = []
		try:
			at_val = 100 - GC(TSSu.promoter['region'])
			print at_val
			if at_val not in at_vals.keys():
				at_vals[at_val] = []
			for gene in TSSu.genes:
				try:
					if gen.genes[gene].fc_pval[cond_fc][1] <= thresh_pval:
						expr.append(gen.genes[gene].fc_pval[cond_fc][0])
					else:
						expr_none.append(gen.genes[gene].fc_pval[cond_fc][0])
				except:
					pass
			if expr != []:
				at_vals[at_val].append(np.mean(expr))
				if np.mean(expr) < 0 - thresh_fc:
					at['rep']['AT'].append(at_val)
					at['rep']['FC'].append(np.mean(expr))
				elif np.mean(expr) > 0 + thresh_fc:
					at['act']['AT'].append(at_val)
					at['act']['FC'].append(np.mean(expr))
				else:
					if expr_none == []:
						at['non']['AT'].append(at_val)
						at['non']['FC'].append(np.mean(expr))
			elif expr_none != []:
				at_vals[at_val].append(np.mean(expr_none))
				at['non']['AT'].append(at_val)
				at['non']['FC'].append(np.mean(expr_none))
		except:
			pass

	nb_act = len(at['act']['AT'])
	nb_rep = len(at['rep']['AT'])
	nb_none = len(at['non']['AT'])
	x = [] ; y = [] ; stds = []
	for at_val in at_vals.keys():
		try:
			val = np.array(at_vals[at_val])
			if val.shape[0] > 10:
				act = val[val > 0 + thresh_fc].shape[0]
				rep = val[val < 0 - thresh_fc].shape[0]
				tot = act + rep
				pexp = float(act) / float(tot)
				std = np.array(stats.binom.std(tot, pexp, loc=0))/tot
				x.append(at_val)
				y.append(pexp)
				stds.append(std)
		except:
			pass
	X = sm.add_constant(x)
	wls_model = sm.WLS(y, X, weights=1/np.power(np.array(stds),2)) # weights proportional to the inverse of std2
	results = wls_model.fit()
	slope = results.params[1]
	OR = results.params[0]
	pval = results.pvalues[1]

	# lists of AT contents per position instead of lists of AT contents per promoter
	# mean of AT contents per position for act, rep and non regulated promoters
	titl = '{}_{}_FC{}_PVAL{}'.format(cond_tss,cond_fc.replace('/',''),str(thresh_fc),str(thresh_pval))
	pathdb = '{}data/{}/discriminator'.format(basedir,gen.name)
	if not os.path.exists(pathdb):
		os.makedirs(pathdb)

	meth = kwargs.get('meth',2)
	if meth ==1:
		x = at['act']['AT'] + at['rep']['AT']# + at['non']['AT']
		y = at['act']['FC'] + at['rep']['FC']# + at['non']['FC']

	if not os.path.exists(pathdb+'/stats.txt'):
		res = open(pathdb+'/stats.txt','w')
		res.close()
	res = open(pathdb+'/stats.txt','a')
	spearman = 	stats.spearmanr(x,y)
	pearson = stats.pearsonr(x,y)
	s0 = str(round(spearman[0],4)) ; s1 = str(round(spearman[1],4))
	p0 = str(round(pearson[0],4)) ; p1 = str(round(pearson[1],4))
	res.write('{}\tPearson = {}|{}\tSpearman= {}|{}\tWLR = {}\n'.format(cond_fc.replace('/',''), s0, s1, p0, p1,str(round(pval,4))))
	res.close()


	if meth == 1:
		fig, ax = plt.subplots()
		plt.scatter(at['act']['AT'], at['act']['FC'], label=str(nb_act)+' Activated', color='red')		
		plt.scatter(at['rep']['AT'], at['rep']['FC'], label=str(nb_rep)+' Repressed', color='blue')		
		plt.scatter(at['non']['AT'], at['non']['FC'], label=str(nb_none)+' Non affected', color='lightgray')		

		ax.set_ylabel('log2FC')
		ax.set_xlabel('AT %')
		ax.set_title(cond_fc)
		fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
		leg = 'All Promoters'	
		ax.legend(title=leg,ncol=1,loc='upper left')
		width = 4 ; height = width / 1.618
		fig.set_size_inches(width, height)
		plt.tight_layout()
		plt.savefig(pathdb+"/"+titl+".svg",transparent=False)
		plt.close('all')
	elif meth == 2:
		fig, ax = plt.subplots()
		plt.plot(x, y, 'rD', markersize=7, linestyle='None')
		plt.errorbar(x, y,yerr=stds,mec='black', capsize=10, elinewidth=1,mew=1,linestyle='None', color='black')
		plt.plot(x,np.array(x)*slope + OR, linestyle='dashed', color='black')    
		ax.set_ylabel("Activated promoters proportion")
		ax.set_xlabel('AT %')
		ax.set_title(cond_fc)
		fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
		width = 4 ; height = width / 1.618
		fig.set_size_inches(width, height)
		plt.tight_layout()
		plt.savefig(pathdb+"/"+titl+".svg",transparent=False)
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
