import numpy as np
from matplotlib import patches
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import sqrt

def draw_melting_energy_circles(self,windows=500, increment=4):
	# generate melting energy circles 
	path = basedir+"data/"+self.name+"/fold_changes/melting_energy_circles-"+str(datetime.now())
	os.makedirs(path)
	if not self.list_cond_fc: # if no fc loaded
		try:
			self.load_fc()
		except:
			print 'Unable to load fc'
			sys.exit()
	if not self.seq:
		try:
			self.load_seq()
		except:
			print'Unable to load seq'
			sys.exit()

	for cond_fc in self.list_cond_fc:
		bins = [] # bins = windows of the genome : [start coordinate,end coordinate,melting energy]
		for i in range(1,self.length,increment): # create bins depending on windows size and increment value
			if (i+windows) <= self.length: # enough length to create a bin
				bins.append([i,i+windows,0])
			else: # i + windows > genome size, last bin of remaining length
				bins.append([i,i+(self.length-i),0])
		bins = np.array(bins) # convert to .npy
		for start,end,melting_energy in bins:
			seq = Seq(self.seq[start-1:end])
			melting_energy = mt.Tm_NN(seq)
		melting_energy = list(bins[:,2])
		cScale_fc = plt.get_cmap('seismic')
		# default values = smallest zscore, highest
		# normalisation of colours
		cNorm_fc  = colors.Normalize(vmin=min(melting_energy), vmax=max(melting_energy)) 
		# map which assigns a colour depending on value between vmin and vmax
		cMap_fc = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_fc) 
		# config, see globvar for more
		# init plots
		fig, ax = plt.subplots(1,1) ; fig.set_size_inches(fig_width, fig_height)
		plt.axis([0, fig_width, 0, fig_height]) ; ax.set_axis_off()

		angle = 360.0/len(melting_energy) # angle between two fragments
		i=0
		for value in melting_energy:
			# edgecolor = assign a colour depending on value using cMap
			# draw arc on circle
			arc = patches.Arc((center_x,center_y), radius, radius, angle=0,theta1=i, theta2=i+angle, edgecolor=cMap_fc.to_rgba(value), lw=10)
			ax.add_patch(arc)
			i+= angle

		cMap_fc._A = [] # fake array to print map
		plt.colorbar(cMap_fc).set_label("Melting energy")
		fig.savefig(path+"/melting_energy_circle-"+cond_fc+".png", format='png', dpi=400, transparent=True) # png (72,300,400 dpi) or svg		

def draw_expression_circles(self, *arg, **kwargs):
	# generate density circles based on FC 
	# opt arguments : colormap, vmin vmax (color normalisation)
	path = basedir+"data/"+self.name+"/fold_changes/circles-"+str(datetime.now())
	os.makedirs(path)
	if not self.list_cond_fc: # if no fc loaded
		try:
			self.load_fc()
		except:
			print 'Unable to load fc'
			sys.exit()

	for cond_fc in self.list_cond_fc:
		gen_states = compute_table_genes(self,cond_fc)
		bins = count_genes_in_windows(self, gen_states, windows=500, increment=4,cond_fc) 
		zscores = compute_zscores(bins)
		# Colormap for fold change
		colormap= kwargs.get('colormap','seismic') # default value seismic
		try:
			cScale_fc = plt.get_cmap(colormap)
		except:
			print 'Incorrect colormap, please check https://matplotlib.org/users/colormaps.html'
			print 'Loading default seismic'
			cScale_fc = plt.get_cmap('seismic')
		# default values = smallest zscore, highest
		vmin = kwargs.get('vmin', min(zscores))
		vmax = kwargs.get('vmax', max(zscores))
		# normalisation of colours
		cNorm_fc  = colors.Normalize(vmin=vmin, vmax=vmax) 
		# map which assigns a colour depending on value between vmin and vmax
		cMap_fc = cmx.ScalarMappable(norm=cNorm_fc, cmap=cScale_fc) 
		# config
		fig_width = 3;fig_height = 3 # figure dimensions
		center_x = fig_width/2;center_y = fig_height/2 # center circles
		radius = fig_width / 2 # radius circles
		# init plots
		fig, ax = plt.subplots(1,1) ; fig.set_size_inches(fig_width, fig_height)
		plt.axis([0, fig_width, 0, fig_height]) ; ax.set_axis_off()
		
		angle = 360.0/len(zscores) # angle between two fragments
		i=0
		for value in zscores:
			# edgecolor = assign a colour depending on value using cMap
			# draw arc on circle
			arc = patches.Arc((center_x,center_y), radius, radius, angle=0,theta1=i, theta2=i+angle, edgecolor=cMap_fc.to_rgba(value), lw=10)
			ax.add_patch(arc)
			i+= angle

		cMap_fc._A = [] # fake array to print map
		plt.colorbar(cMap_fc).set_label("Z-score")
		fig.savefig(path+"/circle-"+cond_fc+".png", format='png', dpi=400, transparent=True) # png (72,300,400 dpi) or svg

def compute_table_genes(self, cond): 
	# returns a npy where a row is a gene caracterised by a start pos and a gene state
	# gene is considered activated above a given fc, repressed below a given fc
	gen_states = []
	for gen in list(self.genes.keys()): # for each gene
		# if activated
		if self.genes[gen].all_fc[cond] >= fc_treshold_pos and self.genes[gen].all_pval[cond] <= pval_treshold:
			gen_states.append([self.genes[gen].start,1])
		# if repressed
		elif self.genes[gen].all_fc[cond] <= fc_treshold_neg and self.genes[gen].all_pval[cond] <= pval_treshold:
			gen_states.append([self.genes[gen].start,-1])
		# if not affected
		else:
			gen_states.append([self.genes[gen].start,0])
	
	gen_states = np.array(gen_states)
	return gen_states

def count_genes_in_windows(self, gen_states, windows, increment,cond):
    if not self.length: # if seq not loaded and thus genome length not available
	    try:
	        self.load_seq()
	    except:
	        print 'Error, unable to load sequence'
	        sys.exit()

	bins = [] # bins = windows of the genome : [start coordinate,end coordinate,nb of activated genes,nb of repressed genes,nb of genes not affected]
	for i in range(1,self.length,increment): # create bins depending on windows size and increment value
		if (i+windows) <= self.length: # enough length to create a bin
			bins.append([i,i+windows,0,0,0])
		else: # i + windows > genome size, last bin of remaining length
			bins.append([i,i+(self.length-i),0,0,0])
	bins = np.array(bins) # convert to .npy
	
	for start,state in gen_states: # reminder gene state : a row = beginning of gene, state (activated, repressed or not affected)
		if state == 1: # activated
		# test to which bins the gene belongs to, and add one to the nb of activated genes of these bins
			bins[np.where((bins[:,0] < start) & (bins[:,1] > start)),2] += 1
		elif state == -1: # repressed gene
			bins[np.where((bins[:,0] < start) & (bins[:,1] > start)),3] += 1
		elif state == 0: # not affected gene
			bins[np.where((bins[:,0] < start) & (bins[:,1] > start)),4] += 1
	return bins

def compute_zscores(bins):
	nb_exp = (np.sum(bins[:,2])) / (np.sum(bins[:,2]) + np.sum(bins[:,3]))
	zscores = []
	for start,end,nb_act,nb_repr,nb_null in bins:
		nb_obs = (nb_act)/(nb_act+nb_repr)
		zscore = 
		zscores.append((p_obs - ))



# proportion de g+ parmi g+ g- est la même que sur le génome = H0
# distribution de proba sous H0 = proportion attendue. pvalue = p(g+ >= 14)
# approx loi normale : B de param (p, n) = N (np, racine np(1-p)). Zscore alors = g+ obs - np / racine avec n (g+g-)