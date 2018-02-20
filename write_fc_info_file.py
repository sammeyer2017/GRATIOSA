import os 
def write_fc_info_file(tagcolumn,fcolumn, columnseparator, startline, pvaluecolumn):
	files = os.listdir(os.getcwd())
	files.remove('write_fc_info_file.py')
	fcinfo = open('fc.info','w')
	fcinfo.write('Condition\tFilename\tTagColumn\tFCcolumn\tColumnSeparator\tStartline\tP-valueColumn\n')
	for file in files:
		fcinfo.write(file[:-4]+'\t'+file+'\t'+str(tagcolumn)+'\t'+str(fcolumn)+'\t'+columnseparator+'\t'+str(startline)+'\t'+str(pvaluecolumn)+'\n')
	fcinfo.close()

write_fc_info_file(0,2,'\\t',2,6)