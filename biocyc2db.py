import sys,os, pysam, numpy as np, pandas as pd, timeit

file = 'TU-with-box.txt'
names = 'Accession-1'
TSS = 'Absolute-Plus-1-Position'
strand = 'Sequence - coordinates of DNA region'
sigma = 'Binds-Sigma-Factor'
m10L = 'MINUS-10-LEFT'
m10R = 'MINUS-10-RIGHT'
m35L = 'MINUS-35-LEFT'
m35R = 'MINUS-35-RIGHT'

df = pd.read_table(file,sep="\t", usecols=[names,TSS,strand,sigma,m10L,m10R,m35L,m35R], \
	converters={names:str,TSS:str,strand:str,sigma:str,m10L:str,m10R:str,m35L:str,m35R:str})

df.replace('',np.nan,inplace=True) ; df.dropna(inplace=True)

def clean_df(df):
	'''
	Clean DF : for each row, replace //, delete spaces...
	'''
	def clean_row(row,row_accumulator,fake):
		new_row = row.to_dict()
		new_row[names] = row[names].replace("//",",").replace("\"","").replace(" ","")
		new_row[sigma] = row[sigma].replace("//",",").replace("\"","").replace(" ","")
		new_row[m10L] = row[m10L].replace("//",",").replace("\"","").replace(" ","")
		new_row[m10R] = row[m10R].replace("//",",").replace("\"","").replace(" ","")
		new_row[m35L] = row[m35L].replace("//",",").replace("\"","").replace(" ","")
		new_row[m35R] = row[m35R].replace("//",",").replace("\"","").replace(" ","")
		row_accumulator.append(new_row)
	
	new_rows = []
	df.apply(clean_row,axis=1,args=(new_rows,''))
	new_df = pd.DataFrame(new_rows)
	return new_df

df = clean_df(df)

def split_strand_df(df,target_column,separator):
	'''
	Starting from df and column containing genome coordinates
	and strand, extract only strand +/- 
	'''
	def split_strand_row(row,row_accumulator,target_column,separator):

		new_row = row.to_dict() # row to dict
		new_row[target_column] = row[target_column].split(separator)[2] # replace strand column
		row_accumulator.append(new_row)

	
	new_rows = []
	df.apply(split_strand_row,axis=1,args =(new_rows,target_column,separator)) # apply to all rows
	new_df = pd.DataFrame(new_rows)
	return new_df

df = split_strand_df(df,strand,' ')

df['sites'] = '' # sites of shape [m10l,m10r,m35l,m35r]

print df.head(10)

def split_prom_df(df,target_columns,separator):
	'''
	Starting from df and target columns : sigma factor, m10l,m10r,m35l,m35r,
	split rows in order to have one row : one sigma factor, one coordinate for each site
	'''
	def split_prom_row(row,row_accumulator,target_columns,separator):

		prom = target_columns[0] ; m10l = target_columns[1]
		m10r = target_columns[2] ; m35l = target_columns[3]
		m35r = target_columns[4]
		sr1 = row[prom].split(separator)
		sr2 = row[m10l].split(separator)
		sr3 = row[m10r].split(separator)
		sr4 = row[m35l].split(separator)
		sr5 = row[m35r].split(separator)
		try:

			while len(sr1) > len(sr2):
				sr2.append(sr2[-1])
			while len(sr1) > len(sr3):
				sr3.append(sr3[-1])
			while len(sr1) > len(sr4):
				sr4.append(sr4[-1])
			while len(sr1) > len(sr5):
				sr5.append(sr5[-1])


			for s1,s2,s3,s4,s5 in zip(sr1,sr2,sr3,sr4,sr5):
				new_row = row.to_dict()
				new_row[prom] = s1
				new_row[m10l] = s2
				new_row[m10r] = s3
				new_row[m35l] = s4
				new_row[m35r] = s5
				new_row['sites'] = s2+','+s3+','+s4+','+s5
				row_accumulator.append(new_row)

		except:
			print 'Weird line, skipped'

	new_rows = []
	df.apply(split_prom_row,axis=1,args =(new_rows,target_columns,separator))
	new_df = pd.DataFrame(new_rows)
	print 'Done'
	return new_df

df = split_prom_df(df,[sigma,m10L,m10R,m35L,m35R],',') ; print df.head(10)
df = df.loc[:,[TSS,strand,names,sigma,'sites']].rename(columns = {strand : 'strand', names : 'genes', sigma : 'sigma'})
df.to_csv('biocyc2db.csv',sep='\t',header=True, index=False, quotechar='"')