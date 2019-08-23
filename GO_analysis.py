import sys
import os
import numpy as np
import pandas as pd
from globvar import *
import scipy.stats as stats
from statsmodels.stats.multitest import *
import Genome

def load_GO(gen,*args,**kwargs):
	if not gen.genes: # if no genes loaded
		# try to load them
		gen.load_annotation()

	if os.path.exists(basedir+"data/"+gen.name+"/GO_analysis/GO.info"):
		with open(basedir+"data/"+gen.name+"/GO_analysis/GO.info","r") as f:
			skiphead = next(f) # skip head
			for header in f:
				header=header.strip()
				header=header.split('\t')
		 		file = header[0]
				tagcol = int(header[1])
				GOcol = int(header[2])
		f.close()
	else:
		print("No GO.info file, please create one")

	with open(basedir+"data/"+gen.name+"/GO_analysis/"+file, 'r') as f:
		header=next(f)
		for line in f:
			line = line.strip('\n').split('\t')
			try:
				gen.genes[line[tagcol]].GO = line[GOcol].split(',')
			except:
				pass    
	f.close()
	gen.GO = {}
	for gene in gen.genes.keys():
		try:
			g = gen.genes[gene]
			for term in g.GO:
				if term not in gen.GO.keys():
					gen.GO[term] = []
				gen.GO[term].append(gene)
		except:
			pass

def GO_enrichment(gen,*args,**kwargs):
    if not hasattr(gen,'GO'):
        load_GO(gen)
    files = os.listdir(basedir+"data/"+gen.name+"/GO_analysis/lists/")
    for filename in files:
		with open(basedir+"data/"+gen.name+"/GO_analysis/lists/"+filename) as f:
			content = [x.strip() for x in f.readlines()] ; f.close()
		res = {}
		valid_genes = []
		for gene in content:
			try:
				g = gen.genes[gene]
				for term in g.GO:
					if term not in res.keys():
						res[term] = {}
						res[term]['list'] = []
						res[term]['input'] = 0
					res[term]['list'].append(gene)
					res[term]['input'] += 1
					valid_genes.append(gene)
			except:
				pass

		g = len(content) # number of submitted genes
		k = len(list(set(valid_genes))) # selection = number of submitted genes annotated at least once
		
		l = [gen.GO[a] for a in gen.GO.keys()] 
		N = len(list(set([item for sublist in l for item in sublist]))) # total number of genes with some annotation in the category e.g. BP
		d = GO_dict(gen)
		# selection of k balls in an urn containing m marked and n non-marked balls, and the observation that the selection contains x marked ball
		file = []
		for term in res.keys():
			try:
				m = len(gen.GO[term]) # number of marked elements = genes annotated for selected term
				n = N - m # number of non marked elements = with some annotation but not to the selected GO term
				x = res[term]['input']# number of marked elements in selection

				pval = np.round(stats.hypergeom.sf(x-1,N,m,k),5) # pval = probability to observe at least x marked balls in the selection.
				if m > 4:
					file.append([term,d[term],pval,x,m])
			except:
				pass
		df = pd.DataFrame(data=file,columns=['GO','Description','P-value','Count','Genome'])
		df['P-value adj (FDR)'] = fdrcorrection(df['P-value'])[1]
		df['%'] = (df['Count']/df['Genome'])*100
		df.sort_values(by=['P-value adj (FDR)'],ascending=True,inplace=True)
		#statsmodels.stats.multitest.multipletests
		df = df.round({'%':0,'P-value':4,'P-value adj (FDR)':4})
		df = df[['GO','Description','Count','Genome','%','P-value adj (FDR)','P-value']]
		df.to_csv(basedir+"data/"+gen.name+"/GO_analysis/res/"+filename[0:-4]+'_results.csv',sep='\t',index=False)

def GO_dict(gen):
	d = {}
	with open(basedir+"data/"+gen.name+"/GO_analysis/GOall.csv") as f:
		skiphead = next(f)
		for line in f:
			line = line.strip().split(',')
			d[line[0]] = line[2]
	f.close()
	return d

# print'total number in population: ' + sys.argv[1]
# print 'total number with condition in population: ' + sys.argv[2]
# print 'number in subset: ' + sys.argv[3]
# print 'number with condition in subset: ' + sys.argv[4] 

# print 'p-value <= ' + sys.argv[4] + ': ' + str(stats.hypergeom.cdf(int(sys.argv[4]) ,int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])))
# print 'p-value >= ' + sys.argv[4] + ': ' + str(stats.hypergeom.sf(int(sys.argv[4]) - 1,int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])))
