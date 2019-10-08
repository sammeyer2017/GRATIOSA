import sys
import os
import numpy as np
import pandas as pd
from globvar import *
import scipy.stats as stats
from statsmodels.stats.multitest import *
import Genome

def load_GO(gen,*args,**kwargs):
	'''
	Loads file specified in GO.info to assign GO terms to genes
	One condition corresponds to one annotation system, e.g. GO, COG, but also domain assignment for domain enrichment
	'''
	if not hasattr(gen, "genes"): # if no genes loaded
		# try to load them
		gen.load_annotation()

	gen.GO = {}
	if os.path.exists(basedir+"data/"+gen.name+"/GO_analysis/GO.info"):
		with open(basedir+"data/"+gen.name+"/GO_analysis/GO.info","r") as f1:
			skiphead = next(f1) # skip head
			for header in f1:
				header=header.strip()
				header=header.split('\t')
				cond = header[0] # condition name
		 		file = header[1] # file containing assignment gene / functions
				tagcol = int(header[2])
				GOcol = int(header[3])
				gen.GO[cond] = {}

				with open(basedir+"data/"+gen.name+"/GO_analysis/"+file, 'r') as f2:
					header=next(f2)
					for line in f2: 
						line = line.strip('\n').split('\t')
						try:
							try:
								gen.genes[line[tagcol]].GO[cond] = line[GOcol].split(',')
							except:
								gen.genes[line[tagcol]].GO = {}
								gen.genes[line[tagcol]].GO[cond] = line[GOcol].split(',')
						except Exception as e:
							# print e
							pass   
				f2.close()
				for gene in gen.genes.keys():
					try:
						g = gen.genes[gene]
						for term in g.GO[cond]:
							if term not in gen.GO[cond].keys():
								gen.GO[cond][term] = []
							gen.GO[cond][term].append(gene)
					except Exception as e:
						# print e
						pass

		f1.close()
	else:
		print("No GO.info file, please create one")

 

def GO_enrichment(gen,*args,**kwargs):
	'''
	Computes GO enrichment analysis from lists of genes
	'''
	if not hasattr(gen,'GO'):
		load_GO(gen)
	cond = kwargs.get("cond","GO")

	files = os.listdir(basedir+"data/"+gen.name+"/GO_analysis/lists/")
	files.remove("old_lists")
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
					res[term]['list'].append(gene)
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

def GO_dict(gen,cond):
	d = {}
	name = "GOall.csv" if cond == "GO" else "domains_descr.csv"
	with open(basedir+"data/"+gen.name+"/GO_analysis/"+name) as f:
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
