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
	for s in [0,1]:
		lgenes[s] = []

	for gene in gen.genes.keys():
		try:
			g = gen.genes[gene]
			lgenes[g.strand].append((g.start,g.end,gene))
		except:
			pass

	for s in lgenes.keys():
		if s:
			lgenes[s] = sorted(lgenes[s])
		else:
			lgenes[s] = sorted(lgenes[s])[::-1]

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

########################################################################
# DIVERSE FUNCTIONS 
########################################################################

def export_TUs(gen):
	if not hasattr(gen, 'TUs'):
		gen.load_TU()

	conds = gen.TUs.keys()		
	conds = ['cov_ratio_corr_pred distance', 'cov_ratio_corr_operons','operons']

	for cond in conds:
		try:
			res = {}
			for idx in gen.TUs[cond].keys():
				gen.TUs[cond][idx].genes = [gen.genes[x].name for x in gen.TUs[cond][idx].genes]
				res[idx] = gen.TUs[cond][idx].__dict__
			df = pd.DataFrame.from_dict(res,orient='index')
			df.sort_index(inplace=True)
			df.to_csv(basedir+"data/"+gen.name+"/TU/res/"+cond+".csv",sep='\t',encoding='utf-8',columns=['start', 'stop', 'left', 'right','orientation','genes'], index=False)
		except:
			pass

def compare_TUs(gen):
	if not hasattr(gen, 'TUs'):
		gen.load_TU()

	conds = gen.TUs.keys()
	pairs = [(conds[i],conds[j]) for i in range(len(conds)) for j in range(i+1, len(conds))]
	res = {}
	# for c1,c2 in pairs:
	conds = ['cov_ratio_corr_pred distance', 'cov_ratio_corr_operons','operons']
	conds = ['pred distance','cov_ratio_corr_pred distance', 'cov_ratio_corr_operons','operons','corr_pred distance','cov_pred distance','ratio_pred distance']
	for cond in conds:
		res[cond] = {} ; res[cond]["all"] = 0
		for idx in gen.TUs[cond].keys():
			tu = gen.TUs[cond][idx]
			if len(tu.genes) not in res[cond].keys():
				res[cond][len(tu.genes)] = 0

			res[cond][len(tu.genes)] += 1
			res[cond]["all"] += 1

	print "Number of genes for each TUs",res
	df = pd.DataFrame.from_dict(res,orient='index')
	df.to_csv(basedir+"data/"+gen.name+"/TU/res/compare_lists.csv",sep='\t',encoding='utf-8')


		# res[(c1,c2)] = {}
		# starts1 = [gen.TUs[c1][tu].start for tu in gen.TUs[c1].keys()]
		# starts2 = [gen.TUs[c2][tu].start for tu in gen.TUs[c2].keys()]
		# stops1 = [gen.TUs[c1][tu].stop for tu in gen.TUs[c1].keys()]
		# stops2 = [gen.TUs[c2][tu].stop for tu in gen.TUs[c2].keys()]

		# res[(c1,c2)]["TUs"] = (len(starts1), len(starts2))
		# res[(c1,c2)]["start"] = len(set(starts1).intersection(set(starts2)))
		# res[(c1,c2)]["stop"] = len(set(stops1).intersection(set(stops2)))  

		# if c1 in ["cov_ratio_corr_operons","operons"] and c2 in ["cov_ratio_corr_operons","operons"]:
		# 	print c1,c2,res[(c1,c2)]
		# if c1 in ["cov_ratio_corr_pred distance","operons"] and c2 in ["cov_ratio_corr_pred distance","operons"]:
		# 	print c1,c2,res[(c1,c2)]


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
	
	if draw_distrib:
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
	min_cov = 20
	draw_distrib = kwargs.get("draw_distrib",True)

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
						pairs.append((g[i],g[i+1],[gen.cov_pos[c][g1.end:g2.start] for c in gen.cov_pos.keys() if np.sum(gen.cov_pos[c][g1.end:g2.start] > min_cov)]))
					else:
						pairs.append((g[i],g[i+1],[gen.cov_pos[c][g2.start:g1.end] for c in gen.cov_pos.keys() if np.sum(gen.cov_pos[c][g2.start:g1.end] > min_cov)]))

				elif not sTU.orientation:
					if g1.end - g2.start > 0:
						pairs.append((g[i],g[i+1],[gen.cov_neg[c][g2.start:g1.end] for c in gen.cov_neg.keys() if np.sum(gen.cov_neg[c][g2.start:g1.end] > min_cov)]))
					else:
						pairs.append((g[i],g[i+1],[gen.cov_neg[c][g1.end:g2.start] for c in gen.cov_neg.keys() if np.sum(gen.cov_neg[c][g1.end:g2.start] > min_cov)]))
					
			gen.TUs[cond][idx].add_intergenic_cov(pairs)

			if len(g) > 1:
				for p in pairs:
					for c in p[2]:
						all_cov += list(c)

		if draw_distrib:
			generate_hist(all_cov,cond)


	if draw_distrib:
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
						covs = [list(x) for x in [cov[c][g1.end:g2.start] for c in cov.keys() if np.sum(cov[c][g1.end:g2.start]) > min_cov]]

					else:					
						covs = [list(x) for x in [cov[c][g2.start:g1.end] for c in cov.keys() if np.sum(cov[c][g2.start:g1.end]) > min_cov]]
					
					for c in covs:
						all_cov += list(c)

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
	notlog = kwargs.get("log",False)

	def generate_hist(data,cond):
		weights = np.ones_like(data)/float(len(data))
		plt.hist(data, weights=weights,range = (0,1),bins=20,align = 'left')
		plt.xlabel("ratio expression") ; plt.ylabel("probability")
		plt.xticks(np.arange(0, 1.1, 0.1))
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

					
					pairs_expr.append((g1,g2,np.mean([min(v)/max(v) for v in values])))

				gen.TUs[cond][idx].add_expression_ratio(pairs_expr)
				if len(g) > 1:
					all_expr += [y[2] for y in pairs_expr]
			except:
				pass	

		if draw_distrib:
			generate_hist(all_expr,cond)

	if draw_distrib:
		all_expr = []
		g = gen.genes.keys()
		pairs = [(g[i],g[j]) for i in range(len(g)) for j in range(i+1, len(g))]
		pairs_expr = []
		# expression ratio among all genes
		for g1,g2 in pairs:
			try:
				if notlog:
					values = [[2**gen.genes[g1].expression[c],2**gen.genes[g2].expression[c]] for c in gen.genes_valid_expr[cond_RNAseq]["conditions"]]
				else:
					values = [[gen.genes[g1].expression[c],gen.genes[g2].expression[c]] for c in gen.genes_valid_expr[cond_RNAseq]["conditions"]]

				pairs_expr.append((g1,g2,np.mean([min(v)/max(v) for v in values])))		
			except:
				pass

		generate_hist([y[2] for y in pairs_expr],"all")






########################################################################
# FUNCTIONS FOR TU FILTERING
########################################################################

def filter_correlation_TU(gen,*arg, **kwargs):
	compute_TUs_corr(gen,draw_distrib = False)
	corr_thresh = kwargs.get("corr_thresh", 0.9) # minimum corr coeff for successive genes to belong to the same TU
	tolerance_thresh = kwargs.get("tolerance_thresh", 0.05) # minimum corr coeff for successive genes to belong to the same TU
	clustering = kwargs.get("clustering", False)
	for cond in gen.TUs.keys():
		newcond = "corr_"+cond
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
							if correlations_TU == []:
								if med > corr_thresh:
									newTU[j].append(candidate)
								else:
									newTU.append([candidate])			
									j += 1
								
							else:
								medTU = np.median(correlations_TU)
								if med >= (medTU - tolerance_thresh):
									newTU[j].append(candidate)
								else:
									newTU.append([candidate])			
									j += 1

						else:
							if med > corr_thresh:
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
	ratio_thresh = kwargs.get("ratio_thresh", 0.3) # minimum expression ratio value for successive genes to belong to the same TU
	tolerance_thresh = kwargs.get("tolerance_thresh", 0.1) # minimum corr coeff for successive genes to belong to the same TU
	clustering = kwargs.get("clustering", False)

	for cond in gen.TUs.keys():
		newcond = "ratio_"+cond
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
						ratios_TU = [ratio[2] for ratio in  sTU.expression_ratio for a,b in pairs if a in ratio and b in ratio]

						med = np.median(ratios)

						if clustering:
							if ratios_TU == []:
								if med > ratio_thresh:
									newTU[j].append(candidate)
								else:
									newTU.append([candidate])			
									j += 1
								
							else:
								medTU = np.median(ratios_TU)
								if med >= (medTU - tolerance_thresh):
									newTU[j].append(candidate)
								else:
									newTU.append([candidate])			
									j += 1

						else:

							if med > ratio_thresh:
								newTU[j].append(candidate)
							else:
								newTU.append([candidate])			
								j += 1
						
						i += 1
				for x in newTU:
					gen.TUs[newcond][gen.genes[x[0]].start] = TU(start=gen.genes[x[0]].start, stop = gen.genes[x[-1]].end, orientation = sTU.orientation, genes = x)
			except:
				pass

def filter_expr_ratio_corr_TU(gen,*arg, **kwargs):
	compute_TUs_expression_ratio(gen,draw_distrib = False)
	compute_TUs_corr(gen,draw_distrib = False)
	idx_thresh = kwargs.get("idx_thresh", 0.3) # minimum value for successive genes to belong to the same TU

	for cond in gen.TUs.keys():
		newcond = "ratio*corr_"+cond
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
						idxs = [ratio[2]*corr[2] for corr in sTU.correlation for ratio in sTU.expression_ratio for g in newTU[j] if g in corr and candidate in corr and g in ratio and candidate in ratio]
						if np.median(idxs) > idx_thresh:
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
	compute_TUs_cov(gen,draw_distrib = False)
	cov_thresh = kwargs.get("cov", 0) # minimum coverage value between intergenic regions for genes to belong to the same TU
	totalcov_thresh = kwargs.get("totalcov", 10) # total coverage between intergenic regions required for filtering
	nbconds_thresh = kwargs.get("nbcond", 3) # nb of conditions necessary

	for cond in gen.TUs.keys():
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
