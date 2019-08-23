from Bio.SeqUtils import GC
from scipy import stats
import statsmodels.api as sm
from globvar import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from useful_functions import *
from collections import Counter
import os
import subprocess

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


def compute_bubble_opening(gen,condTSS,*arg,**kwargs):
	'''
	Starting from transcription bubble sequence (default: length 21 bp, 12 bp upstream TSS pos and 8 bp downstream), 
	computes opening probability for given SC lvls using Twist DNA, initiation rate, normalized initiation rate
	'''
	SCs = kwargs.get('SC',[-0.03,-0.06])
	path2twist = kwargs.get('path','/home/raphael/Twist-DNA_1.1/')
	temp = kwargs.get('t',37)
	salt = kwargs.get('salt',0.1)
	bsizemin = kwargs.get('bmin',1)
	bsizemax = kwargs.get('bmax',0)
	bstep = kwargs.get('bstep',1)
	bthresh = kwargs.get('bthresh',-3)
	flank = 100 # number of poly G to add at each side

	pos = 11 # eq to -2 pos

	B = kwargs.get('B',1)
	m = kwargs.get('m',2)

	gen.compute_magic_prom()

	res = {'info':{'condTSS':condTSS,'SCs':SCs}} 
	kall = {}
	for SC in SCs:
		kall[SC] = 0

	for TSSpos in gen.TSSs[condTSS].keys():
		try:
			TSS = gen.TSSs[condTSS][TSSpos]
			fasta = open(path2twist+'bubble.fasta','w')
			fasta.write('>bubble\n')

			fasta.write('G'*flank)
			fasta.write(TSS.promoter[sigfactor]['discr_model'])
			fasta.write('G'*flank)
			fasta.close()

			for SC in SCs:
				param = open(path2twist+'input.dat','w')
				param.write(str(temp)+'\t::Temperature T in Celsius (default 37 C)\n')
				param.write(str(SC)+'\t::Superhelical density Sigma (default -0.06)\n')
				param.write(str(salt)+'\t::Salt concentration CNa (default 0.1 M)\n')
				param.write(str(bsizemin)+'\t::Smallest bubble size bsizemin\n')
				param.write(str(bsizemax)+'\t::Largest bubble size bsizemax\n')
				param.write(str(bstep)+'\t::Bubble size step bstep\n')
				param.write(str(bthresh)+'\t::Bubble storing threshold (default -3)')
				param.close()


				#twistDNA = '{}TwistDNA < {}bubble.fasta'.format(path2twist,path2twist)
				subprocess.Popen("./TwistDNA < bubble.fasta", cwd=path2twist, shell=True)
				
				f = open(path2twist+'outputs/openproba.bed', 'r')
				lines = f.readlines()
				val = float(lines[pos+flank].strip().split()[3])
				val = val / np.log10(np.exp(1))
				#val = np.power(10,val)
				f.close()

				k = np.exp(-B*(val+m))
				
				kall[SC] += k

				try:
					res[TSSpos][SC] = {'p':val,'k':k,'kn':None}
				except:
					res[TSSpos] = {}
					res[TSSpos][SC] = {'p':val,'k':k,'kn':None}

		except Exception as e:
			pass

	for TSSpos in gen.TSSs[condTSS].keys():
		try:
			TSS = gen.TSSs[condTSS][TSSpos]
			for SC in SCs:
				res[TSSpos][SC]['kn'] = res[TSSpos][SC]['k'] / kall[SC]
		except:
			pass

	return res


def predict_FC_from_opening(gen,res,condFC,*arg,**kwargs):
	'''
	Computes free energy (~ level of expression) for each promoter, based on total theta and stiffness.
	Args : NPd,ABCd from previous step. Kwargs : var (i.e. which parameter is a variable of sequence)
	= 'all','theta','kt','none'. sig = default level of sigma supercoiling level, dsig = delta sigma caused
	by relaxation. Returns a dict of shape {'NP':{splen:[predictedlog2FC1]}}.
	'''
	# results of shape {'NP':{17:[logFC1, logFC2]}}
	contrasts = kwargs.get('contrasts',['-0.03_-0.06'])
	thresh_pval = kwargs.get('thresh_pval', 0.05) # below, gene considered valid, above, gene considered non regulated
	thresh_fc = kwargs.get('thresh_fc', 0) # 0 +- thresh_fc : below, gene considered repressed, above, gene considered activated, between, gene considered non regulated

	FCs = {'all':[], 'valid':[], 'non':[]}
	condTSS = res['info']['condTSS']

	for TSSpos in res.keys():
		try:
			TSS = gen.TSSs[condTSS][TSSpos]

			valid_expr = [] # expression values for valid genes in that TSS i.e. pval < thresh
			non_expr = [] # expression values for non valid genes in that TSS i.e. pval > thresh
			all_expr = [] # expression values for all genes in that TSS
			for gene in TSS.genes:
				try:
					if gen.genes[gene].fc_pval[condFC][1] <= thresh_pval:
						valid_expr.append(gen.genes[gene].fc_pval[condFC][0])
					else:
						non_expr.append(gen.genes[gene].fc_pval[condFC][0])

					all_expr.append(gen.genes[gene].fc_pval[condFC][0])
				except:
					pass

			for contrast in contrasts:
				SC1 = float(contrast.split('_')[0])
				SC2 = float(contrast.split('_')[1])

			FCpred = np.log2(res[TSSpos][SC2]['kn']/res[TSSpos][SC1]['kn'])

			if valid_expr != []:
				FCs['valid'].append([np.mean(valid_expr),FCpred])

			if all_expr != []:
				FCs['all'].append([np.mean(all_expr),FCpred])

			if non_expr != []:
				FCs['non'].append([np.mean(non_expr),FCpred])

		except Exception as e:
			pass


	# FCreal = [a[0] for a in FCs['valid']]
	# FCpred = [a[1] for a in FCs['valid']]
	# x = FCpred
	# y = FCreal
	
	non = [a[1] for a in FCs['non']]
	act = [a[1] for a in FCs['valid'] if a[0] > 0 + thresh_fc]
	rep = [a[1] for a in FCs['valid'] if a[0] < 0 - thresh_fc]
	print len(non),len(act),len(rep)
	pathdb = '{}data/{}/discriminator/'.format(basedir,gen.name)
	if not os.path.exists(pathdb):
		os.makedirs(pathdb)
	titl = 'test'

	meth = kwargs.get('meth','raph')
	org = kwargs.get('org', gen.name)

	if meth =='raph':
		mact = np.mean(act)
		mrep = np.mean(rep) 
		mnon = np.mean(non)
		means = [mact,mnon,mrep]

		eact = np.std(act)/np.sqrt(len(act)) # 95 CI on mean
		erep = np.std(rep)/np.sqrt(len(rep)) # 95 CI on mean
		enon = np.std(non)/np.sqrt(len(non)) # 95 CI on mean
		stds = [eact,enon,erep]

		pval = stats.ttest_ind(rep,act,equal_var=False)[1] / 2
		s = significance(pval)

		fig_width = 2 ; fig_height = 2.2
		fig = plt.figure(figsize=(fig_width,fig_height))
		xval = [1,2,3] ; labs = ['act','non','rep']
		yval = means
		if s == 'ns':
			fs = 13
			dt = 0.01
		else:
			fs = 18
			dt = 0

		plt.bar(xval, yval, yerr = stds, width = 0.6, color = 'white',edgecolor = 'black', linewidth = 2, ecolor = 'black', capsize = 5)
		plt.xticks(xval,labs,rotation = 45)
		plt.xlim(0.5, 3.5)
		ylim1 = min(yval) - 0.1
		ylim2 = max(yval) + 0.1
		plt.ylim(ylim1,ylim2)
		barplot_annotate_brackets(0, 2, s,xval,yval,yerr=stds,fs=fs, dt=dt)#,dh=0.05, barh=.05, fs=14, dt=0)
		fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
		plt.ylabel('Predicted log2(FC)')
		plt.title('{}, {}'.format(gen.name,condFC),fontsize=11)
		plt.tight_layout()
		plt.savefig(pathdb+titl+'.svg',transparent=False)
		plt.close('all')


	if meth == 'all':
		# X = sm.add_constant(x)
		# wls_model = sm.WLS(y, X, weights=1/np.power(np.array(stds),2)) # weights proportional to the inverse of std2
		# results = wls_model.fit()
		# slope = results.params[1]
		# OR = results.params[0]
		# pval = str(round(results.pvalues[1],3))
		spearman = str(round(stats.spearmanr(x,y)[1],3))
		pearson = str(round(stats.pearsonr(x,y)[1],3))
		pval = None

		width = 3 ; height = width / 1.4
		fig, ax = plt.subplots()
		plt.plot(x, y, marker='o', color='black',markersize=7, linestyle='None', fillstyle='none')
		# plt.errorbar(x, y,yerr=stds,mec='black', capsize=10, elinewidth=1,mew=1,linestyle='None', color='black')
		# plt.plot(x,np.array(x)*slope + OR, color='black',linewidth=0.5)   
		ax.set_ylabel("FCreal")
		ax.set_xlabel('FCpred')
		# ax.set_xlim(min(x)-10,max(x)+10)
		ax.set_title(condFC+", Reg: {}, Pear.: {}, Spear.: {}".format(pval,spearman,pearson),fontsize=9)
		fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
		fig.set_size_inches(width, height)
		plt.tight_layout()
		plt.savefig(pathdb+titl+"_pred"+"_reg.svg",transparent=False)
		plt.close('all')


def export_opening_prob(gen,condTSS,*arg,**kwargs):
	SCs = np.linspace(0,-0.1,11)
	path2twist = kwargs.get('path','/home/raphael/Twist-DNA_1.1/')
	temp = kwargs.get('t',37)
	salt = kwargs.get('salt',0.1)
	bsizemin = kwargs.get('bmin',1)
	bsizemax = kwargs.get('bmax',0)
	bstep = kwargs.get('bstep',1)
	bthresh = kwargs.get('bthresh',-3)
	flank = 100 # number of poly G to add at each side

	pos = 11 # eq to -2 pos

	B = kwargs.get('B',1)
	m = kwargs.get('m',2)

	gen.compute_magic_prom()

	res = []
	i = 0
	for TSSpos in gen.TSSs[condTSS].keys():
		if i < 100 and TSSpos == 1287665:
			try:
				val = []
				TSS = gen.TSSs[condTSS][TSSpos]
				fasta = open(path2twist+'bubble.fasta','w')
				fasta.write('>bubble\n')
				fasta.write('G'*flank)
				fasta.write(TSS.promoter[sigfactor]['discr_model'])
				fasta.write('G'*flank)
				fasta.close()

				for SC in SCs:
					param = open(path2twist+'input.dat','w')
					param.write(str(temp)+'\t::Temperature T in Celsius (default 37 C)\n')
					param.write(str(SC)+'\t::Superhelical density Sigma (default -0.06)\n')
					param.write(str(salt)+'\t::Salt concentration CNa (default 0.1 M)\n')
					param.write(str(bsizemin)+'\t::Smallest bubble size bsizemin\n')
					param.write(str(bsizemax)+'\t::Largest bubble size bsizemax\n')
					param.write(str(bstep)+'\t::Bubble size step bstep\n')
					param.write(str(bthresh)+'\t::Bubble storing threshold (default -3)')
					param.close()

					subprocess.Popen("./TwistDNA < bubble.fasta", cwd=path2twist, shell=True)
					
					f = open(path2twist+'outputs/openproba.bed', 'r')
					lines = f.readlines()
					v = float(lines[pos+flank].strip().split()[3])
					v = v / np.log10(np.exp(1))
					val.append(v)
					f.close()
				res.append(val)
			except Exception as e:
				pass
	
	return SCs,res
