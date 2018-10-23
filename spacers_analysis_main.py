org = raw_input('Enter organism name : ')
from globvar import *
import Genome
try:
	gen = Genome.Genome(name=org)
	try:
		gen.load_fc_pval() ; print 'FC data successfully loaded : {}.'.format(gen.genes_valid.keys())
	except Exception as e:
		print 'Error while loading FC data :',e
	try:
		gen.load_TSS() ; print 'TSS data successfully loaded : {}.'.format(gen.TSSs.keys())
	except Exception as e:
		print 'Error while loading TSS data :',e
	try:
		gen.load_seq() ; print 'Seq data successfully loaded.'
	except Exception as e:
		print 'Error while loading seq data :',e

	TSScond = raw_input('Enter studied TSS condition : ')
	FCcond = raw_input('Enter studied FC condition : ')
	try:
		shift = input('compute_magic_prom : enter shift (default 0) : ')
		gen.compute_magic_prom(TSScond,shift=shift) ; print 'Magic prom successfully computed.'
	except Exception as e:
		print 'Error while computing magic prom :',e
	try:
		pval = input('compute_spacer_response : enter p-val threshold (default 0.05) : ')
		spacer_sigma = gen.compute_spacer_response(FCcond,TSScond,thresh_pval=pval) ; print 'Spacer response successfully computed.'
	except Exception as e:
		print 'Error while computing spacer response :',e
	try:
		res = raw_input('Still working on sigma factor {} and spacers {} ? : (yes/no) '.format(sigfactor,spacers))
		if res == 'no':
			sigfactor = raw_input('Please enter new sigma factor : ')
			spacers = input('Please enter new list of spacers (e.g. 15,16,17,18,19) : ')

		res = raw_input('Wanna draw results ? : (yes/no) ')
		if res == 'yes':
			params = input('Please enter parameters (FC threshold for actrep float default 0, aggregation True/False default True, print list of genes True/False default False) : ')
			import spacers
			spacers.draw_results(spacer_sigma, thresh_fc=params[0], aggregation=params[1], l_genes=params[2])
			print 'Results for spacer response successfully drawn.'
	except Exception as e:
		print 'Error while drawing spacer response results :',e
	
	try:
		res = raw_input('Wanna compare spacer response results to predicted results ? : (yes/no) ')
		if res == 'yes':
			gen.compute_magic_prom(TSScond,shift=2)
			validTSS = gen.valid_TSS(FCcond,TSScond,params[0])
			import spacers_modelling
			print 'Computing total theta and stiffness...'
			NPd1,ABCd1 = spacers_modelling.compute_angle_stiff(validTSS['act'])
			NPd2,ABCd2 = spacers_modelling.compute_angle_stiff(validTSS['rep'])
			params = input('Please enter parameters (var for model all/theta/kt/none, sig default -0.06, dsig default 0.03) : ')
			print 'Computing free energy and predicted log2FC...'
			res1 = spacers_modelling.compute_free_energy(NPd1,ABCd1,var=params[0],sig=params[1],dsig=params[2])
			res2 = spacers_modelling.compute_free_energy(NPd2,ABCd2,var=params[0],sig=params[1],dsig=params[2])
			spacers_modelling.MW_stats(res1,res2,'ABC')
	except Exception as e:
		print 'Error while comparing spacer response results to predicted results :',e

except Exception as e:
	print 'Error in organism',org,':',e