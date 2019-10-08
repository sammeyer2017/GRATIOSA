import numpy as np
import pandas as pd
from globvar import *
from TU import TU
import matplotlib.pyplot as plt
from collections import OrderedDict 
# Project of TU identification, first in dickeya

########################################################################
# FUNCTIONS FOR TU PREDICTION
########################################################################

def predict_TU_distance(gen,*arg, **kwargs):
	"""
	Defines potential TUs in genome object: successive isodirectional genes (same strand) with small intergenic distances
	"""
	d_thresh = kwargs.get("d_thresh", 500) # maximal intergenic distance for successive genes to belong to the same TU
	if not hasattr(gen, 'TUs'):
		gen.load_TU()

	gen.TUs["pred distance"] = {}

	lgenes = {}
	lgenes[0] = sorted([(gen.genes[g].start,gen.genes[g].end,g) for g in gen.genes.keys() if not gen.genes[g].strand])[::1]
	lgenes[1] = sorted([(gen.genes[g].start,gen.genes[g].end,g) for g in gen.genes.keys() if gen.genes[g].strand])

	for s in lgenes.keys():
		i = 0
		while i < len(lgenes[s]):
			j = i
			# while successive genes are close enough, extends TU
			if s:
				try:
					while lgenes[s][j+1][0] - lgenes[s][j][1] < d_thresh:
						j += 1
				except: # list index out of range: all genes processed
					pass
			elif not s:
				try:
					while lgenes[s][j][1] - lgenes[s][j+1][0] < d_thresh:
						j += 1
				except: # list index out of range: all genes processed
					pass
				

			j += 1
			genes = lgenes[s][i:j]
			names = [x[2] for x in genes]

			start = genes[0][0] ; stop = genes[-1][1]

			gen.TUs["pred distance"][start] = TU(start=start, stop=stop, orientation=s, genes = names)
			#gen.TUs["pred correlation"][start].add_correlation(pairs)
			i = j

def predict_TU_correlation(gen,*arg, **kwargs):
	"""
	Defines potential TUs in genome object: successive isodirectional genes (same strand) with high correlation coefficients
	"""
	corr_thresh = kwargs.get("corr_thresh", 0.75) # minimum corr coeff for successive genes to belong to the same TU
	if not hasattr(gen, 'TUs'):
		gen.load_TU()

	if not hasattr(gen, 'genes_valid_expr'):
		gen.load_expression()
	
	gen.TUs["pred correlation"] = {}

	gene_expression = {}
	dfs = {}
	for s in [0,1]:
		gene_expression[s] = {}
		dfs[s] = {}

	for gene in gen.genes.keys():
		try:
			g = gen.genes[gene]
			gene_expression[g.strand][(g.start,g.end,gene)] = g.expression
		except: # list index out of range: all genes processed
			pass

	for s in gene_expression.keys():
		dfs[s] = pd.DataFrame.from_dict(gene_expression[s], orient='index', dtype=float)
		dfs[s].sort_index(inplace=True) ; dfs[s] = dfs[s].T.corr()


	for s in dfs.keys():
		i = 0
		while i < (dfs[s].shape[0]):
			j = i
			# while successive genes have enough correlation values, extends TU
			try:
				while dfs[s].iloc[j,j+1] > corr_thresh and j < dfs[s].shape[0]-2:
					j += 1
			except:
				pass

			j += 1
			genes = dfs[s].index[i:j]			
			names = [x[2] for x in genes]

			if s:
				start = genes[0][0] ; stop = genes[-1][1]
			elif not s:
				start = genes[-1][0] ; stop = genes[0][1]
				names = names[::-1]

			gen.TUs["pred correlation"][start] = TU(start=start, stop=stop, orientation=s, genes = names)
			#gen.TUs["pred correlation"][start].add_correlation(pairs)

			i = j

def predict_TU_cov(gen,*arg, **kwargs):
	"""
	Defines potential TUs in genome object: successive isodirectional genes (same strand) with cov > 0 in intergenic regions
	"""
	if not hasattr(gen, 'cov_pos') or not hasattr(gen, 'cov_neg'):
		gen.load_cov()
	
	if not hasattr(gen, 'TUs'):
		gen.load_TU()

	UTR_thresh = kwargs.get("UTR", 50) # mean length of UTR
	totalcov_thresh = kwargs.get("totalcov", 100) # minimal coverage required between UTR for a condition to be taken into account
	successive_thresh = kwargs.get("successive", 10) # nb of successive positions where cov = 0 between intergenic regions required to end a TU
	nbconds_thresh = kwargs.get("nbcond", 5) # nb of conditions where successive coverage between intergenic is 0 necessary to end a TU
	
	gen.TUs["pred cov"] = {}

	cov = {} ; cov[0] = gen.cov_neg ; cov[1] = gen.cov_pos

	lgenes = {}
	lgenes[0] = sorted([(gen.genes[g].start,gen.genes[g].end,g) for g in gen.genes.keys() if not gen.genes[g].strand])[::1]
	lgenes[1] = sorted([(gen.genes[g].start,gen.genes[g].end,g) for g in gen.genes.keys() if gen.genes[g].strand])
	for s in lgenes.keys():
		i = 0
		while i < len(lgenes[s]):
			j = i
			extend = True
			while extend and j < len(lgenes[s]) -2 :
				extend = False
				g1 = gen.genes[lgenes[s][j][2]]
				g2 = gen.genes[lgenes[s][j+1][2]]

				start = g1.end if s else g1.start # end 3'UTR first gene
				stop = g2.start if s else g2.end  # start 5'UTR second gene

				# List of all intergenic cov (one for each RNAseq condition) if genes enough expressed
				ig_cov = [cov[s][c][start:stop] for c in cov[s].keys() if np.sum(cov[s][c][start-UTR_thresh:stop+UTR_thresh]) > totalcov_thresh]
				# List of all positions where cov = 0 for each RNAseq condition
				counts = [np.where(np.array(c) == 0)[0] for c in ig_cov]
				# for each RNAseq condition, split successive positions where cov = 0 into arrays 
				successive = [np.split(c, np.where(np.diff(c) != 1)[0]+1) for c in counts]

				tot = 0 # Number of RNAseq condition where coverage = 0 on successive positions
				for condRNAseq in successive:					
					nullcov = [np.where(np.shape(c)[0] >= successive_thresh) for c in condRNAseq]
					if np.sum([np.shape(z[0])[0] for z in nullcov]) >= 1:
						tot +=1

				if tot < nbconds_thresh:
					extend = True
					j += 1



				# counts = [np.shape(np.where(np.array(c) == 0))[1] for c in ig_cov]
				# if np.shape(np.where(np.array(counts) == 0))[1] > nbconds_thresh:
				# 	extend = True
				# 	j += 1

			j += 1

			genes = lgenes[s][i:j]
			names = [x[2] for x in genes] if s else [x[2] for x in genes][::-1]
			TUstart = genes[0][0] if s else genes[-1][0] 
			TUstop = genes[-1][1] if s else genes[0][1]   
			gen.TUs["pred cov"][TUstart] = TU(start=TUstart, stop=TUstop, orientation=s, genes = names)
			i = j

########################################################################
# DIVERSE FUNCTIONS 
########################################################################

def export_TUs(gen):
	if not hasattr(gen, 'TUs'):
		gen.load_TU()

	conds = gen.TUs.keys()		
	#conds = ['cov_ratio_corr_pred distance', 'cov_ratio_corr_operons','operons']

	for cond in conds:
		try:
			res = {}
			for idx in gen.TUs[cond].keys():
				gen.TUs[cond][idx].genes = [gen.genes[x].name for x in gen.TUs[cond][idx].genes]
				if len(gen.TUs[cond][idx].genes) > 1:
					res[idx] = gen.TUs[cond][idx].__dict__
			df = pd.DataFrame.from_dict(res,orient='index')
			df.sort_index(inplace=True)
			#df.to_csv(basedir+"data/"+gen.name+"/TU/res/"+cond+".csv",sep='\t',encoding='utf-8',columns=['start', 'stop', 'left', 'right','orientation','genes'], index=False)
			#df.to_csv(basedir+"data/"+gen.name+"/TU/res/"+cond+".csv",sep='\t',encoding='utf-8', index=False, columns=['start', 'stop','orientation','genes','TSS','TTS'],)
			df.to_csv(basedir+"data/"+gen.name+"/TU/res/"+cond+".csv",sep='\t',encoding='utf-8', index=False)
		except:
			pass

def compare_TUs(gen):
	if not hasattr(gen, 'TUs'):
		gen.load_TU()

	conds = gen.TUs.keys()
	res = {}
	for cond in conds:
		res[cond] = {} ; res[cond]["all"] = 0
		for idx in gen.TUs[cond].keys():
			tu = gen.TUs[cond][idx]
			if len(tu.genes) not in res[cond].keys():
				res[cond][len(tu.genes)] = 0

			res[cond][len(tu.genes)] += 1
			res[cond]["all"] += 1
		
		if 1 in res[cond].keys():
			print cond,res[cond][1],"TUs with one gene"

		print cond,np.sum([res[cond][k] for k in res[cond].keys() if k != 1]),"all TUs"

	df = pd.DataFrame.from_dict(res,orient='index')
	df.sort_index(axis=1, inplace=True)
	df.to_csv(basedir+"data/"+gen.name+"/TU/res/compare_lists.csv",sep='\t',encoding='utf-8')

	paired_cond = True
	if paired_cond:
		ctest = "pred cov"
		cref = "operons"
		res = {}
		for idx in gen.TUs[cref].keys():
			TUref = gen.TUs[cref][idx]
			try:
				TUcand = gen.TUs[ctest][idx]
				diff = len(TUcand.genes) - len(TUref.genes)
				if diff not in res.keys():
					res[diff] = 0
				res[diff] += 1
			except:
				pass
		res['all'] = np.sum([res[k] for k in res.keys()])
		df = pd.DataFrame.from_dict(res,orient='index')
		df.sort_index(axis=1, inplace=True)
		df.to_csv(basedir+"data/"+gen.name+"/TU/res/compare_conds.csv",sep='\t',encoding='utf-8')


		# pairs = [(conds[i],conds[j]) for i in range(len(conds)) for j in range(i+1, len(conds))]
		# for c1,c2 in pairs:
		# 	res[(c1,c2)] = {}
		# 	starts1 = [gen.TUs[c1][tu].start for tu in gen.TUs[c1].keys()]
		# 	starts2 = [gen.TUs[c2][tu].start for tu in gen.TUs[c2].keys()]
		# 	stops1 = [gen.TUs[c1][tu].stop for tu in gen.TUs[c1].keys()]
		# 	stops2 = [gen.TUs[c2][tu].stop for tu in gen.TUs[c2].keys()]

		# 	res[(c1,c2)]["TUs"] = (len(starts1), len(starts2))
		# 	res[(c1,c2)]["start"] = len(set(starts1).intersection(set(starts2)))
		# 	res[(c1,c2)]["stop"] = len(set(stops1).intersection(set(stops2)))  

		# 	if c1 in ["cov_ratio_corr_operons","operons"] and c2 in ["cov_ratio_corr_operons","operons"]:
		# 		print c1,c2,res[(c1,c2)]
		# 	if c1 in ["cov_ratio_corr_pred distance","operons"] and c2 in ["cov_ratio_corr_pred distance","operons"]:
		# 		print c1,c2,res[(c1,c2)]

def plot_cov(gen,*arg,**kwargs):
	"""
	Plot cov of a given genomic region
	"""
	if not hasattr(gen, 'cov_pos') or not hasattr(gen, 'cov_neg'):
		gen.load_cov()
	
	if not hasattr(gen, 'TUs'):
		gen.load_TU()

	gnames = {} # dict {gene_name:gene_id} to map names to ID for genome
	for g in gen.genes.keys():
		try:
			gnames[gen.genes[g].name] = g
		except:
			pass

	cov = {} ; cov[0] = gen.cov_neg ; cov[1] = gen.cov_pos

	lgenes = ["Dda3937_00879","Dda3937_00880","Dda3937_00881","Dda3937_00882","Dda3937_00883","Dda3937_00884","Dda3937_04310","Dda3937_04488","Dda3937_01884","Dda3937_01885","Dda3937_01886","Dda3937_01887","Dda3937_01888","Dda3937_01889","Dda3937_01890","Dda3937_01891","Dda3937_01892","Dda3937_01893","Dda3937_01894","Dda3937_01897","Dda3937_01898"]
	
	lgenes = ["queA", "tgt", "yajC", "secD", "secF"]
	lgenes = [gnames[g] for g in lgenes]

	lgenes = ["yajG", "ampG", "cyoA", "cyoB", "cyoC", "cyoD", "cyoE"]
	lgenes = [gnames[g] for g in lgenes]

	lgenes = ["Dda3937_02416", "Dda3937_02417", "Dda3937_02418", "Dda3937_02419", "Dda3937_02420"]

	s = gen.genes[lgenes[0]].strand
	s = 1 if s else 0
	pos = (gen.genes[lgenes[0]].start - 500, gen.genes[lgenes[-1]].end + 500) if s else (gen.genes[lgenes[-1]].end - 500, gen.genes[lgenes[0]].start + 500)

	colors = ["green"]*4 + ["blue"]*3
	colors = ["blue"]*5
	for c in cov[s].keys():
		width = 5 ; height = 3
		fig,ax = plt.subplots()
		plt.plot(np.arange(pos[0], pos[1]+1, 1), cov[s][c][pos[0]-1:pos[1]])
		plt.axhline(y=0, color = 'red')		
		wid = (max(cov[s][c][pos[0]-1:pos[1]]) - min(cov[s][c][pos[0]-1:pos[1]]))/15.0
		i = 0
		for g in lgenes:
			plt.axvline(x=gen.genes[g].start, color = 'green', linestyle='dashed',ymax=0.15, alpha = 0.75)
			plt.axvline(x=gen.genes[g].end, color = 'red', linestyle='dashed',ymax=0.15, alpha = 0.75)
			plt.arrow(gen.genes[g].start, 0, gen.genes[g].end - gen.genes[g].start, 0, width=wid, head_width=wid, head_length=200, alpha = 0.25, color=colors[i])
			i +=1
		plt.xlabel("Positions") ; plt.ylabel("Coverage") ; plt.title(c)
		fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97) ; fig.set_size_inches(width, height) ; plt.tight_layout()
		plt.savefig(basedir+"data/"+gen.name+"/TU/res/cov/{}{}.png".format(str(pos),c)) ; plt.close('all')



########################################################################
# FUNCTIONS TO ADD INFORMATIONS ON TU
########################################################################

def compute_TUs_corr(gen,*arg, **kwargs):
	if not hasattr(gen, 'TUs'):
		gen.load_TU()

	if not hasattr(gen, 'genes_valid_expr'):
		gen.load_expression()

	draw_distrib = kwargs.get("draw_distrib",True)
	notlog = kwargs.get("log",False)
	draw_general_distrib = kwargs.get("draw_general_distrib",False)

	def generate_hist(data,cond):
		weights = np.ones_like(data)/float(len(data))
		plt.hist(data, weights=weights,range = (0,1),bins=20,align = 'left')
		plt.xlabel("correlation coefficient") ; plt.ylabel("probability")
		plt.xticks(np.arange(0, 1.1, 0.1))
		plt.savefig(basedir+"data/"+gen.name+"/TU/res/"+cond+"_correlation.svg",bbox_inches="tight")
		plt.close('all')
	
	gene_expression = {}
	for g in gen.genes.keys():
		try:
			gene_expression[g] = gen.genes[g].expression
		except:
			pass

	# creates DF from dict
	df = pd.DataFrame.from_dict(gene_expression,orient='index', dtype=float)
	# Delete log or not
	if notlog:
		df = df.apply(lambda x : 2**x)
	# transforms DF in expression correlation matrix
	df = df.T.corr() 

	for cond in gen.TUs.keys():
		all_corr = [] # correlation values for all the TU of this condition
		condTU = gen.TUs[cond]
		for idx in gen.TUs[cond].keys():
			sTU = gen.TUs[cond][idx]
			g = sTU.genes
			# all possible pairs of genes among TU + their correlation value
			pairs = [(g[i],g[j],df[g[i]][g[j]]) for i in range(len(g)) for j in range(i+1, len(g))]
			gen.TUs[cond][idx].add_correlation(pairs) # add correlation data to TU

			if len(g) > 1:
				all_corr += [y[2] for y in pairs]

		if draw_distrib:
			generate_hist(all_corr,cond)
	
	if draw_general_distrib:
		# correlation among all genes
		# deletes diagonal i.e. gene corr with itself
		data =  df.values ; np.fill_diagonal(data, np.nan)
		data = data[~np.isnan(data)]
		generate_hist(data,"all")


def compute_TUs_cov(gen,*arg, **kwargs):
	if not hasattr(gen, 'TUs'):
		gen.load_TU()

	if not hasattr(gen, 'cov_pos') or not hasattr(gen, 'cov_neg'):
		gen.load_cov()
	min_cov = 100
	draw_distrib = kwargs.get("draw_distrib",True)
	draw_general_distrib = kwargs.get("draw_general_distrib",False)

	def generate_hist(data,cond):
		weights = np.ones_like(data)/float(len(data))
		plt.hist(data, weights=weights,range = (0,100),bins=20,align = 'left')
		plt.xlabel("coverage") ; plt.ylabel("probability")
		plt.savefig(basedir+"data/"+gen.name+"/TU/res/"+cond+"_cov.svg",bbox_inches="tight")
		plt.close('all')
	

	for cond in gen.TUs.keys():
		condTU = gen.TUs[cond]
		all_cov = []
		for idx in gen.TUs[cond].keys():
			sTU = gen.TUs[cond][idx]
			g = sTU.genes
			pairs = []
			for i in range(len(g)-1):
				g1 = gen.genes[g[i]]
				g2 = gen.genes[g[i+1]]
				if sTU.orientation:
					if g2.start - g1.end > 0:
						pairs.append((g[i],g[i+1],[gen.cov_pos[c][g1.end:g2.start] for c in gen.cov_pos.keys() if np.sum(gen.cov_pos[c][g1.end-50:g2.start+50] > min_cov)]))
					else:
						pairs.append((g[i],g[i+1],[gen.cov_pos[c][g2.start:g1.end] for c in gen.cov_pos.keys() if np.sum(gen.cov_pos[c][g2.start-50:g1.end+50] > min_cov)]))

				elif not sTU.orientation:
					if g1.end - g2.start > 0:
						pairs.append((g[i],g[i+1],[gen.cov_neg[c][g2.start:g1.end] for c in gen.cov_neg.keys() if np.sum(gen.cov_neg[c][g2.start-50:g1.end+50] > min_cov)]))
					else:
						pairs.append((g[i],g[i+1],[gen.cov_neg[c][g1.end:g2.start] for c in gen.cov_neg.keys() if np.sum(gen.cov_neg[c][g1.end-50:g2.start+50] > min_cov)]))
			
			gen.TUs[cond][idx].add_intergenic_cov(pairs)

			if len(g) > 1:
				for p in pairs:
					all_cov += [np.mean(c) for c in p[2]]

		if draw_distrib:
			generate_hist(all_cov,cond)


	if draw_general_distrib:
		gene_pos = {0:[],1:[]}
		all_cov = []

		for gene in gen.genes.keys():
			try:
				g = gen.genes[gene]
				gene_pos[g.strand].append((g.start,gene))
			except:
				pass

		for s in [0,1]:
			cov = gen.cov_pos if s else gen.cov_neg
			gene_pos[s] = sorted(gene_pos[s])
			for i in range(len(gene_pos[s])-1):
				try:
					g1 = gen.genes[gene_pos[s][i][1]]
					g2 = gen.genes[gene_pos[s][i+1][1]]

					if s:
						covs = [np.mean(x) for x in [cov[c][g1.end-50:g2.start+50] for c in cov.keys() if np.sum(cov[c][g1.end-50:g2.start+50]) > min_cov]]

					else:					
						covs = [np.mean(x) for x in [cov[c][g2.start-50:g1.end+50] for c in cov.keys() if np.sum(cov[c][g2.start-50:g1.end+50]) > min_cov]]
					
					for c in covs:
						all_cov += c

				except:
					pass
		generate_hist(all_cov,"all")
	

def compute_TUs_expression_ratio(gen,*arg, **kwargs):
	if not hasattr(gen, 'TUs'):
		gen.load_TU()

	if not hasattr(gen, 'genes_valid_expr'):
		gen.load_expression()

	draw_distrib = kwargs.get("draw_distrib",True)
	cond_RNAseq = kwargs.get("cond_RNAseq","log2rpkm_rnaseq_PE_grouped.csv")
	notlog = kwargs.get("log",True)
	draw_general_distrib = kwargs.get("draw_general_distrib",False)

	def generate_hist(data,cond):
		weights = np.ones_like(data)/float(len(data))
		plt.hist(data, weights=weights,range = (0,5),bins=100,align = 'left')
		plt.xlabel("ratio expression") ; plt.ylabel("probability")
		plt.xticks(np.arange(0, 5.1, 0.5))
		plt.savefig(basedir+"data/"+gen.name+"/TU/res/"+cond+"_expression_ratio.svg",bbox_inches="tight")
		plt.close('all')

	for cond in gen.TUs.keys():
		all_expr = []
		condTU = gen.TUs[cond]
		for idx in gen.TUs[cond].keys():
			try:
				sTU = gen.TUs[cond][idx]
				g = sTU.genes
				pairs = [(g[i],g[j]) for i in range(len(g)) for j in range(i+1, len(g))]
				pairs_expr = []
				for g1,g2 in pairs:
					if notlog:
						values = [[2**gen.genes[g1].expression[c],2**gen.genes[g2].expression[c]] for c in gen.genes_valid_expr[cond_RNAseq]["conditions"]]
					else:
						values = [[gen.genes[g1].expression[c],gen.genes[g2].expression[c]] for c in gen.genes_valid_expr[cond_RNAseq]["conditions"]]

					
					pairs_expr.append((g1,g2,[v[0]/v[-1] for v in values]))

				gen.TUs[cond][idx].add_expression_ratio(pairs_expr)
				if len(g) > 1:
					all_expr += [item for sublist in [y[2] for y in pairs_expr] for item in sublist]
			except:
				pass	

		if draw_distrib:
			generate_hist(all_expr,cond)

	if draw_general_distrib:
		g = gen.genes.keys()
		pairs = [(g[i],g[j]) for i in range(len(g)) for j in range(i+1, len(g))]
		pairs_expr = np.array([])
		# expression ratio among all genes
		for g1,g2 in pairs:
			try:
				if notlog:
					values = [[2**gen.genes[g1].expression[c],2**gen.genes[g2].expression[c]] for c in gen.genes_valid_expr[cond_RNAseq]["conditions"]]
				else:
					values = [[gen.genes[g1].expression[c],gen.genes[g2].expression[c]] for c in gen.genes_valid_expr[cond_RNAseq]["conditions"]]

				val = [v[0]/v[-1] for v in values]
				pairs_expr = np.append(pairs_expr,val)
			except:
				pass

		generate_hist(pairs_expr,"all")




def compute_TUs_idx_corr_ratio(gen,*arg, **kwargs):
	compute_TUs_expression_ratio(gen,draw_distrib = False)
	compute_TUs_corr(gen,draw_distrib = False)
	draw_distrib = kwargs.get("draw_distrib",True)
	def generate_hist(data,cond):
		weights = np.ones_like(data)/float(len(data))
		plt.hist(data, weights=weights,range = (0,1),bins=20,align = 'left')
		plt.xlabel("idx ratio*cor") ; plt.ylabel("probability")
		plt.xticks(np.arange(0, 1.1, 0.1))
		plt.savefig(basedir+"data/"+gen.name+"/TU/res/"+cond+"_ratio*cor.svg",bbox_inches="tight")
		plt.close('all')

	for cond in gen.TUs.keys():
		all_idxs = []
		condTU = gen.TUs[cond]
		for idx in gen.TUs[cond].keys():
			try:
				sTU = gen.TUs[cond][idx]
				g = sTU.genes
				pairs = [(g[i],g[j]) for i in range(len(g)) for j in range(i+1, len(g))]
				pairs_idx = []
				for g1,g2 in pairs:				
					idxpair = [ratio[2]*corr[2] for corr in sTU.correlation for ratio in sTU.expression_ratio if g1 in corr and g2 in corr and g1 in ratio and g2 in ratio]
					pairs_idx.append((g1,g2,idxpair[0]))

				gen.TUs[cond][idx].add_idx_corr_ratio(pairs_idx)
				if len(g) > 1:
					all_idxs += [y[2] for y in pairs_idx]
			except Exception as e:
				print e
				pass	

		if draw_distrib:
			generate_hist(all_idxs,cond)


	notlog = kwargs.get("log",False)
	draw_general_distrib = kwargs.get("draw_general_distrib",False)
	cond_RNAseq = kwargs.get("cond_RNAseq","log2rpkm_rnaseq_PE_grouped.csv")

	if draw_general_distrib:
		gene_expression = {}
		for g in gen.genes.keys():
			try:
				gene_expression[g] = gen.genes[g].expression
			except:
				pass
		# creates DF from dict
		df = pd.DataFrame.from_dict(gene_expression,orient='index', dtype=float)
		# Delete log or not
		if notlog:
			df = df.apply(lambda x : 2**x)
		# transforms DF in expression correlation matrix
		df = df.T.corr() 

		g = gen.genes.keys()
		pairs = [(g[i],g[j]) for i in range(len(g)) for j in range(i+1, len(g))]
		pairs_idxs = []
		# expression ratio among all genes
		for g1,g2 in pairs:
			try:
				print gen.genes[g1].expression,gen.genes[g2].expression

				if notlog:
					values = [[2**gen.genes[g1].expression[c],2**gen.genes[g2].expression[c]] for c in gen.genes_valid_expr[cond_RNAseq]["conditions"]]
				else:
					values = [[gen.genes[g1].expression[c],gen.genes[g2].expression[c]] for c in gen.genes_valid_expr[cond_RNAseq]["conditions"]]

				pairs_idxs.append((g1,g2,np.mean([min(v)/max(v) for v in values])*df[g1][g2]))	
			
			except Exception as e:
				pass

		generate_hist([y[2] for y in pairs_idxs],"all")




########################################################################
# FUNCTIONS FOR TU FILTERING
########################################################################

def filter_correlation_TU(gen,*arg, **kwargs):
	compute_TUs_corr(gen,draw_distrib = False)
	corr_thresh = kwargs.get("corr_thresh", 0.9) # minimum corr coeff for successive genes to belong to the same TU
	tolerance_thresh = kwargs.get("tolerance_thresh", 0.05) # minimum corr coeff for successive genes to belong to the same TU
	clustering = kwargs.get("clustering", 0)
	conds = kwargs.get("conds", gen.TUs.keys())
	for cond in conds:
		newcond = "clust"+str(clustering)+str(corr_thresh)+"corr_"+cond
		gen.TUs[newcond] = {}
		for idx in gen.TUs[cond].keys():
			try:
				sTU = gen.TUs[cond][idx]
				newTU = []
				i = 0 ; j = 0
				newTU.append([sTU.genes[i]])
				if len(sTU.genes) > 1:				
					while i < len(sTU.genes)-1:
						candidate = sTU.genes[i+1]
						correlations = [corr[2] for corr in sTU.correlation for g in newTU[j] if g in corr and candidate in corr]
						
						pairs = [(newTU[j][a],newTU[j][b]) for a in range(len(newTU[j])) for b in range(a+1, len(newTU[j]))]
						correlations_TU = [corr[2] for corr in sTU.correlation for a,b in pairs if a in corr and b in corr]

						med = np.median(correlations)

						if clustering:
							thresh = corr_thresh if correlations_TU == [] else np.median(correlations_TU) - tolerance_thresh
							if med >= thresh and med >= corr_thresh:
								newTU[j].append(candidate)
							else:
								newTU.append([candidate])	
								j += 1
						else:
							if med >= corr_thresh:
								newTU[j].append(candidate)
							else:
								newTU.append([candidate])			
								j += 1
						i += 1
				
				for x in newTU:
					gen.TUs[newcond][gen.genes[x[0]].start] = TU(start=gen.genes[x[0]].start, stop = gen.genes[x[-1]].end, orientation = sTU.orientation, genes = x)
			except:
				pass


def filter_expression_ratio_TU(gen,*arg, **kwargs):
	compute_TUs_expression_ratio(gen,draw_distrib = False)
	bound1 = 0.5 ; bound2 = 2 ; cond_thresh = 3 ; counts_thresh = 0
	ratio_thresh = kwargs.get("ratio_thresh", 0.5) # minimum expression ratio value for successive genes to belong to the same TU
	tolerance_thresh = kwargs.get("tolerance_thresh", 0.1) # minimum corr coeff for successive genes to belong to the same TU
	clustering = kwargs.get("clustering", False)
	counting = kwargs.get("counting", True)
	conds = kwargs.get("conds", gen.TUs.keys())
	for cond in conds:
		newcond = "clust_"+str(clustering)+"_ratio_"+cond
		gen.TUs[newcond] = {}
		for idx in gen.TUs[cond].keys():
			try:
				sTU = gen.TUs[cond][idx]
				newTU = []
				i = 0 ; j = 0
				newTU.append([sTU.genes[i]])
				if len(sTU.genes) > 1:				
					while i < len(sTU.genes)-1:
						candidate = sTU.genes[i+1]
						ratios = [ratio[2] for ratio in sTU.expression_ratio for g in newTU[j] if g in ratio and candidate in ratio]

						pairs = [(newTU[j][a],newTU[j][b]) for a in range(len(newTU[j])) for b in range(a+1, len(newTU[j]))]
						ratios_TU = [ratio[2] for ratio in sTU.expression_ratio for a,b in pairs if a in ratio and b in ratio]

						med = np.median(ratios)

						if counting:
							counts = [np.shape(r[(r > bound1) & (r < bound2)])[0] for r in [np.array(x) for x in ratios]]
							if np.shape(np.where(np.array(counts) > cond_thresh))[1] > counts_thresh:
								newTU[j].append(candidate)
							else:
								newTU.append([candidate])	
								j += 1

						elif clustering:
							thresh = ratio_thresh if ratios_TU == [] else np.median(ratios_TU) - tolerance_thresh
							if med >= thresh and med >= ratio_thresh:
								newTU[j].append(candidate)
							else:
								newTU.append([candidate])	
								j += 1
						else:
							if med >= ratio_thresh:
								newTU[j].append(candidate)
							else:
								newTU.append([candidate])			
								j += 1
						
						i += 1

				for x in newTU:
					gen.TUs[newcond][gen.genes[x[0]].start] = TU(start=gen.genes[x[0]].start, stop = gen.genes[x[-1]].end, orientation = sTU.orientation, genes = x)
			except Exception as e:
				pass

def filter_idx_corr_ratio(gen,*arg, **kwargs):
	compute_idx_corr_ratio(gen,draw_distrib = False)
	idx_thresh = kwargs.get("idx_thresh", 0.3) # minimum value for successive genes to belong to the same TU
	tolerance_thresh = kwargs.get("tolerance_thresh", 0.1) # minimum corr coeff for successive genes to belong to the same TU
	clustering = kwargs.get("clustering", False)

	conds = kwargs.get("conds", gen.TUs.keys())
	for cond in conds:
		newcond = "clust_"+str(clustering)+"_ratio*corr_"+cond
		gen.TUs[newcond] = {}
		for idx in gen.TUs[cond].keys():
			try:
				sTU = gen.TUs[cond][idx]
				newTU = []
				i = 0 ; j = 0
				newTU.append([sTU.genes[i]])
				if len(sTU.genes) > 1:				
					while i < len(sTU.genes)-1:
						candidate = sTU.genes[i+1]
						idxs = [idxpair[2] for idxpair in sTU.idx_corr_ratio for g in newTU[j] if g in ratio and candidate in ratio]

						pairs = [(newTU[j][a],newTU[j][b]) for a in range(len(newTU[j])) for b in range(a+1, len(newTU[j]))]
						idxs_TU = [idxpair[2] for idxpair in sTU.idx_corr_ratio for a,b in pairs if a in idxpair and b in idxpair]
					
						med = np.median(idxs)
						if clustering:
							if idxs_TU == []:
								if med > idx_thresh:
									newTU[j].append(candidate)
								else:
									newTU.append([candidate])			
									j += 1
								
							else:
								medTU = np.median(idxs_TU)
								if med >= (medTU - tolerance_thresh):
									newTU[j].append(candidate)
								else:
									newTU.append([candidate])			
									j += 1

						else:

							if med > idx_thresh:
								newTU[j].append(candidate)
							else:
								newTU.append([candidate])			
								j += 1
						
						i += 1


				for x in newTU:
					gen.TUs[newcond][gen.genes[x[0]].start] = TU(start=gen.genes[x[0]].start, stop = gen.genes[x[-1]].end, orientation = sTU.orientation, genes = x)
			except:
				pass



def filter_coverage_TU(gen,*arg, **kwargs):
	compute_TUs_cov(gen, draw_distrib=False)
	cov_thresh = kwargs.get("cov", 0) # minimum coverage value between intergenic regions for genes to belong to the same TU
	totalcov_thresh = kwargs.get("totalcov", 10) # total coverage between intergenic regions required for filtering
	nbconds_thresh = kwargs.get("nbcond", 3) # nb of conditions necessary

	conds = kwargs.get("conds", gen.TUs.keys())
	for cond in conds:
		newcond = "cov_"+cond
		gen.TUs[newcond] = {}
		for idx in gen.TUs[cond].keys():
			try:
				sTU = gen.TUs[cond][idx]
				newTU = []
				i = 0 ; j = 0
				newTU.append([sTU.genes[i]])
				if len(sTU.genes) > 1:				
					while i < len(sTU.genes)-1:
						candidate = sTU.genes[i+1]
						covs = [cov[2] for cov in sTU.intergenic_cov if sTU.genes[i] in cov and candidate in cov][0]
						count = 0
						for cov in covs:
							if np.sum(cov) >= totalcov_thresh:
								if np.shape(np.where(cov <= cov_thresh))[1] > 0:
									count += 1
						
						if count >= nbconds_thresh:
							newTU.append([candidate])
							j += 1
						else:
							newTU[j].append(candidate)
						
						i += 1
				for x in newTU:
					gen.TUs[newcond][gen.genes[x[0]].start] = TU(start=gen.genes[x[0]].start, stop = gen.genes[x[-1]].end, orientation = sTU.orientation, genes = x)
			except:
				pass

########################################################################
# PUTATIVE TSS AND TTS
########################################################################

def add_TSS_TTS_TU(gen,*arg, **kwargs):
	if not hasattr(gen, 'TUs'):
		gen.load_TU()
	if not hasattr(gen, 'TSSs'):
		gen.load_TSS()
	if not hasattr(gen, 'TTSs'):
		gen.load_TTS()

	condTSS = kwargs.get("condTSS", "raw_TSS")
	condTTSrhodpdt = kwargs.get("condTTS", "RhoTerm")
	condTTSrhoindpdt = kwargs.get("condTTS", "ARNold")
	condTU = kwargs.get("condTU", "operons")

	allTSS = gen.TSSs[condTSS].keys()
	allTTS = gen.TTSs[condTTSrhodpdt].keys() + gen.TTSs[condTTSrhoindpdt].keys()

	for idx in gen.TUs[condTU].keys():
		sTU = gen.TUs[condTU][idx]
		start = gen.TUs[condTU][idx].start
		stop = gen.TUs[condTU][idx].stop
		if sTU.orientation:
			TSSs = [TSS for TSS in allTSS if TSS >= start - 100 and TSS <= start - 10]
			TTSs = [TTS for TTS in allTTS if TTS >= stop + 10 and TTS <= stop + 100]
		else:
			TSSs = [TSS for TSS in allTSS if TSS <= start + 100 and TSS >= start + 10]
			TTSs = [TTS for TTS in allTTS if TTS >= stop - 100 and TTS <= stop - 10]

		gen.TUs[condTU][idx].add_TSS(TSSs)
		gen.TUs[condTU][idx].add_TTS(TTSs)

		# Mean 5'UTR: Genome-Wide Identification of Transcription Start Sites, Promoters and Transcription Factor Binding Sites in E. coli
		# Mean 3'UTR: AU-Rich Long 3' Untranslated Region Regulates Gene Expression in Bacteria


########################################################################
# QUANTITATIVE TU MAPPING
########################################################################

