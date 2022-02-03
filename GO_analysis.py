import sys
import os
import numpy as np
import pandas as pd
from globvar import *
import scipy.stats as stats
from statsmodels.stats.multitest import *
import Genome


def GO_enrichment(gen,*args,**kwargs):
	'''
	Computes GO enrichment analysis from lists of genes
	'''
	if not hasattr(gen,'GO'):
		load_GO(gen)
	cond = kwargs.get("cond","GO")

	files = os.listdir(basedir+"data/"+gen.name+"/GO_analysis/lists/")
	#files.remove("old_lists")
	# Lists of genes within GO enrichment analysis has to be computed
	for filename in files:
		with open(basedir+"data/"+gen.name+"/GO_analysis/lists/"+filename) as f:
			content = [x.strip() for x in f.readlines()] ; f.close()
		res = {}
		valid_genes = []
		for gene in content:
			try:
				g = gen.genes[gene]
				for term in g.GO[cond]:
					if term not in res.keys():
						res[term] = {}
						res[term]['list'] = []
						res[term]['input'] = 0
					
					if g.name not in res[term]['list']:
						res[term]['list'].append(g.name)
						res[term]['input'] += 1
						valid_genes.append(gene)
			except:
				pass

		g = len(content) # number of submitted genes
		k = len(list(set(valid_genes))) # selection = number of submitted genes annotated at least once
		
		l = [gen.GO[cond][a] for a in gen.GO[cond].keys()] 
		N = len(list(set([item for sublist in l for item in sublist]))) # total number of genes with some annotation in the category e.g. BP
		d = GO_dict(gen,cond)
		# selection of k balls in an urn containing m marked and n non-marked balls, and the observation that the selection contains x marked ball
		file = []
		for term in res.keys():
			try:
				m = len(gen.GO[cond][term]) # number of marked elements = genes annotated for selected term
				n = N - m # number of non marked elements = with some annotation but not to the selected GO term
				x = res[term]['input']# number of marked elements in selection

				pval = np.round(stats.hypergeom.sf(x-1,N,m,k),5) # pval = probability to observe at least x marked balls in the selection.
				if m > 4: # If there is at least 4 genes annotated to that term, otherwise it is not relevant
					file.append([term,d[term],", ".join(res[term]['list']),pval,x,m])
			except:
				pass
		
		df = pd.DataFrame(data=file,columns=['GO','Description','Genes','P-value','Count','Genome'])
		df['P-value adj (FDR)'] = fdrcorrection(df['P-value'])[1]
		df['%'] = (df['Count']/df['Genome'])*100
		df.sort_values(by=['P-value adj (FDR)'],ascending=True,inplace=True)
		#statsmodels.stats.multitest.multipletests
		df = df.round({'%':0,'P-value':4,'P-value adj (FDR)':4})
		df = df[['GO','Description','Genes','Count','Genome','%','P-value adj (FDR)','P-value']]
		df.to_csv(basedir+"data/"+gen.name+"/GO_analysis/res/"+filename[0:-4]+'_results.csv',sep='\t',index=False)

def GO_dict(gen,cond):
	'''
	Associates each GO term to its description
	'''
	d = {}
	name = "GOall.csv" if cond == "GO" else "domains_descr.csv"
	with open(basedir+"data/"+gen.name+"/GO_analysis/"+name) as f:
		skiphead = next(f)
		for line in f:
			line = line.strip().split('\t')
			d[line[0]] = line[2]
	f.close()
	return d