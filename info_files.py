import os

path_db = '/home/kwl06/Documents/stage_GB5/data/data/'
# path = os.environ['TOPO_DATABASE_PATH']+'data/'
for organism in os.listdir(path_db):
	path = path_db+organism
	write_fc_info_file(path)
	write_TSS_info_file(path)

def write_fc_info_file(path):
	try:	
		fcinfo = open(path+'/fold_changes/fc.info','w')
		fcinfo.write('Condition\tFilename\tTagColumn\tFCcolumn\tColumnSeparator\tStartline\tP-valueColumn\n')
		fcinfo.close()
	# for file in files:
	# 	fcinfo.write(file[:-4]+'\t'+file+'\t'+str(tagcolumn)+'\t'+str(fcolumn)+'\t'+columnseparator+'\t'+str(startline)+'\t'+str(pvaluecolumn)+'\n')
	except Exception as e:
		print 'Error for organism :',path
		print 'Error :',e
#write_fc_info_file(0,2,'\\t',2,6)

def write_TSS_info_file(path):
	try:
		info = open(path+'/TSS/TSS.info','w')
		info.write('Condition\tFilename\tTagColumn\tTSSColumn\tStartline\tColSeparator\tStrandColumn\tSigmaColumn\tSitesColumn\n')
		info.close()

	except Exception as e:
		print 'Error for organism :',path
		print 'Error :',e
#,tagcol,tsscol,colsep,startline,strandcol,sigmacol,sitescol
#info.write(tagcol+'\t',tsscol+'\t',colsep+'\t',startline+'\t',strandcol+'\t',sigmacol+'\t',sitescol+'\t')
