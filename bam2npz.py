import pysam, numpy as np, pandas as pd, datetime

def process_bam_paired_end(bam_file): # /!\ paired-end only /!\ -> return fragments for each strand
	bamfile = pysam.AlignmentFile(bam_file, "rb") # BAM opening, alignment file object
	print 'Header :',bamfile.header
	print 'Reads mapped :',bamfile.mapped
	print 'Reads without coordinate :',bamfile.nocoordinate
	print 'Reads not mapped :',bamfile.unmapped
	# /!\ 0-based coordinate system /!\
	# fastest way to build a numpy matrix -> store every coordinate of interest into lists, then merge into numpy array
	# lists storing coordinates of R1 for +/- strands
	Rpos_start = []
	Rpos_end = []
	Rneg_start = []
	Rneg_end = []

	for read in bamfile.fetch(): # only mapped reads to reference fetched
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
	if not os.path.exists(os.getcwd()+'/rna_seq_reads'):
		os.makedirs('rna_seq_reads')
	if not os.path.exists(os.getcwd()+'/rna_seq_reads/reads.info'):
		file = open('reads.info','w') 
		file.write('BAM file\tDate\tReads file')
		file.close() 

	file = open('reads.info','a')
	file.write('\n'+bam_file+'\t'+datetime.now()+'\t'+bam_file[:-4]+'_reads.npz')  
	file.close()

	np.savez(os.getcwd()+'/rna_seq_reads/'+bam_file[:-4]+'_reads.npz', Rpos=Rpos, Rneg=Rneg)


# compute coverage from reads (.npz, one .npy per strand)
def cov_from_reads(npz_file, genome_length): # compute coverage from reads (.npz, one .npy per strand)
	if not os.path.exists(os.getcwd()+'/rna_seq_cov'):
		os.makedirs('rna_seq_cov')

	npzfile = np.load(os.getcwd()+'/rna_seq_reads/'+npz_file)
	Rpos = npzfile["Rpos"]
	Rneg = npzfile["Rneg"]
	cov_pos = np.zeros(genome_length, dtype=int)
	cov_neg = np.zeros(genome_length, dtype=int)
	for start,end in Rpos:
		cov_pos[start:end+1] += 1

	for start,end in Rneg:
		cov_neg[start:end+1] += 1

	if not os.path.exists(os.getcwd()+'/rna_seq_cov/cov.info'):
		file = open(os.getcwd()+'/rna_seq_cov/cov.info','w') 
		file.write('Reads file\tDate\tCov file')
		file.close() 

	file = open(os.getcwd()+'/rna_seq_cov/cov.info','a')
	file.write('\n'+npz_file+'\t'+datetime.now()+'\t'+npz_file[:-4]+'_cov.npz')  
	file.close()

	np.savez(os.getcwd()+'/rna_seq_cov/'+npz_file[:-4]+'_cov.npz', cov_pos=cov_pos, cov_neg=cov_neg)


##### MAIN #####
# process_bam_paired_end('E1.bam')
# cov_from_reads('E1_reads.npz', 4922802)


	
#grep -rl 'depth' ~/Documents/stage_GB5/test_grep | xargs sed -i 's/depth/cov/g'






# list of genes, start, end, FC
# P < 0.05, FC > 1.2

# elif value_type == 'stat':
#     # Multiply by 0.1 so that the
#     # values range from -0.5 to 0.5.
#     value = round(float(value)*0.1, 2))

#     # Get the absolute value because fold change is
#     # subjective. It doesn't matter which population
#     # expresses gene more as long as the difference is
#     # apparent.
#     value = abs(value)

#     r = 0 + value
#     g = 0.5 - value
#     b = 0.3
# nb of bp in each chrom
# BP = [('I' , 28185914), ('II' , 23295652), ('III' , 16798506),
#       ('IV' , 32632948), ('V' , 12251397), ('VI' , 17083675),
#       ('VII' , 27937443), ('VIII' , 19368704), ('IX' , 20249479),
#       ('X' , 15657440), ('XI' , 16706052), ('XII' , 18401067),
#       ('XIII' , 20083130), ('XIV' , 15246461), ('XV' , 16198764),
#       ('XVI' , 18115788), ('XVII' , 14603141), ('XVIII' , 16282716),
#       ('XIX' , 20240660), ('XX' , 19732071), ('XXI' ,11717487)]

# totalBP = 400788495

# def parse_stat_file(f):
#     """
#     Parses the stat file without headers and creates a dictionary:
#     {'chromosome number' : [(bp, log(stat)), (bp, log(stat)), 
#                             (bp, log(stat))]} 
#         * log(stat) ranges from -5 to 5.
#     """
#     stat_dict = {}

#     for line in f:
#         line = line.split('\t')
#         group = line[1]

#         if "group" in group:
#             group = group[5:]
#             BP = int(line[2])
#             stat = math.log(float(line[6]))
        
                
#             current_entry = (BP, stat)

                            
#             if group in stat_dict:
#                 stat_dict[group].append(current_entry)

#             else:
#                 stat_dict[group] = []
#                 stat_dict[group].append(current_entry)

#             stat_dict[group].sort()
            
# return stat_dict