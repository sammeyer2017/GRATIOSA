import os,pysam, numpy as np, pandas as pd
from datetime import datetime

def process_bam_paired_end(bam_file): # /!\ paired-end only /!\ -> return fragments for each strand
	if not os.path.exists(bam_file+".bai"): # if index needed, created using samtools (.bai file)
		os.system("samtools index -b %s"%bam_file)
    
    bamfile = pysam.AlignmentFile(bam_file, "rb") # BAM opening, alignment file object
	print 'Header :',bamfile.header
	print 'Reads mapped :',bamfile.mapped
	print 'Reads without coordinate :',bamfile.nocoordinate
	print 'Reads not mapped :',bamfile.unmapped
	# /!\ 0-based coordinate system /!\
	# fastest way to build a numpy matrix -> store every coordinate of interest into lists, then merge into numpy array
	# lists storing coordinates of paired-end fragments for +/- strands
	Rpos_start = []
	Rpos_end = []
	Rneg_start = []
	Rneg_end = []
	# per default, only mapped reads to reference fetched
	for read in bamfile.fetch(): 
		# if correctly paired and R1/R2 in different strands and fragment size < 1000 kb
		if read.is_read1 and read.is_paired and not read.mate_is_unmapped and read.is_reverse != read.mate_is_reverse and abs(read.template_length) < 1000 :
			if read.is_reverse: 
				Rneg_start.append(read.reference_end)
				Rneg_end.append(read.reference_end + abs(read.template_length))
			elif not read.is_reverse:
				Rpos_start.append(read.reference_start)
				Rpos_end.append(read.reference_start + read.template_length)
	bamfile.close();
 
	# conversion step into numpy array
	Rpos_start = np.array(Rpos_start,dtype=int)
	Rpos_end = np.array(Rpos_end,dtype=int)
	Rneg_start = np.array(Rneg_start,dtype=int)
	Rneg_end = np.array(Rneg_end,dtype=int)

	Rpos = np.column_stack((Rpos_start,Rpos_end))
	Rneg = np.column_stack((Rneg_start,Rneg_end))

	# delete rows where one coordinate is missing
	Rpos = Rpos[np.isfinite(Rpos).all(axis=1)] 
	Rneg = Rneg[np.isfinite(Rneg).all(axis=1)]

	# if rnaseq_reads folder not exist
	if not os.path.exists(os.getcwd()+'/rnaseq_reads'):
		os.makedirs('rnaseq_reads')
	# if reads.info not exist
	if not os.path.exists(os.getcwd()+'/rnaseq_reads/reads.info'):
		file = open(os.getcwd()+'/rnaseq_reads/reads.info','w') 
		file.write('Condition\tReads file\tDate\tBAM file')
		file.close() 
	# save results ; .npz contains two .npy Rpos and Rneg
	file = open(os.getcwd()+'/rnaseq_reads/reads.info','a')
	file.write('\n'+bam_file[:-4]+'\t'+bam_file[:-4]+'_reads.npz\t'+str(datetime.now())+'\t'+bam_file)  
	file.close()	
	np.savez(os.getcwd()+'/rnaseq_reads/'+bam_file[:-4]+'_reads.npz', Rpos=Rpos, Rneg=Rneg)


# compute coverage from reads (.npz, one .npy per strand)
def cov_from_reads(npz_file, genome_length):
	# load npz file
	npzfile = np.load(os.getcwd()+'/rnaseq_reads/'+npz_file)
	Rpos = npzfile["Rpos"]
	Rneg = npzfile["Rneg"]
	# init cov
	cov_pos = np.zeros(genome_length, dtype=int)
	cov_neg = np.zeros(genome_length, dtype=int)
	# compute cov
	# on positive strand
	for start,end in Rpos:
		cov_pos[start:end+1] += 1
	# on negative strand
	for start,end in Rneg:
		cov_neg[start:end+1] += 1
	# if rnaseq_cov folder not exist		
	if not os.path.exists(os.getcwd()+'/rnaseq_cov'):
		os.makedirs('rnaseq_cov')
	# if cov.info not exist
	if not os.path.exists(os.getcwd()+'/rnaseq_cov/cov.info'):
		file = open(os.getcwd()+'/rnaseq_cov/cov.info','w') 
		file.write('Condition\tCov file\tDate\tReads file')
		file.close() 
	# save results
	file = open(os.getcwd()+'/rnaseq_cov/cov.info','a')
	file.write('\n'+npz_file[0:-10]+'\t'+npz_file[:-4]+'_cov.npz\t'+str(datetime.now())+'\t'+npz_file)  
	file.close()
	# save results ; .npz contains two .npy cov_pos and cov_neg
	np.savez(os.getcwd()+'/rnaseq_cov/'+npz_file[:-4]+'_cov.npz', cov_pos=cov_pos, cov_neg=cov_neg)

##### MAIN #####
# for a list of .bam
samples=["E%d"%x for x in range(1,15)]+["F%d"%x for x in range(1,9)+[13,14]]

for s in samples:
	process_bam_paired_end('%s.bam'%s)
    cov_from_reads('%s_reads.npz'%s, 4922802) # dickeya genome length
