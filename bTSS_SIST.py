###################### EXTRAITS DE GENOME #######################

    """
    from TSS import TSS
    from TTS import TTS
    from TU import TU
    from btssfinder import *

    def predict_promoter_from_TSS(self,list_TSS,*args,**kwargs): #running bTSSfinder
        freedom = kwargs.get('free',0)
        nameOut = kwargs.get('out',list_TSS+'-btss-'+str(freedom))
        '''
        kwargs possible : out & free
        bTSSfinder need to be install on this computer
        for more informations you can read btssfinder.py
        obj is Genome object
        list_TSS is the TSS list eg. biocyc
        out is  the name of fasta file (whiout extension .fasta) or if don't exist will become file's name of all out
        free is number of additionnal base at the TSS region
        NEXT
        convert out-bTSSfinder.gff on basedir/data/[nom_liste_TSS]
        AND FINALY
        write the localization of new csv in file TSS.info to the next load.TSS()
        '''
        try:
            test = self.TSSs[list_TSS]
            test = self.seq[0:4]
        except:
            self.load_TSS()
            self.load_seq()

        run_btssfinder(self,list_TSS,nameOut,freedom)

        gff2csv(self,list_TSS,nameOut,freedom)
        TSSinfo = basedir+"data/"+self.name+"/TSS/TSS.info"
        if Path(TSSinfo).exists():
            exist = False
            f = open(TSSinfo,"r")
            for i in f.readlines():
                line = i.split('\t')
                if line[0] == nameOut:
                    exist = True
            f.close()
            if not exist:
                f = open(TSSinfo,"a")
                f.write(nameOut+'\t'+nameOut+'.csv'+'\t'+"2"+'\t'+"0"+'\t'+"2"+'\t'+"\\t"+'\t'+"1"+'\t'+"3"+'\t'+"4"+'\n')
                f.close()
        else:
            print("TSS info not found")
        print("Finished…"+'\n'+"Now, you can visualise file "+TSSinfo+" or you can just reload TSS list.")
            

    def load_SIST(self, start, end,*args, **kwargs):
        if not hasattr(self, 'SIST_profile'):
            self.SIST_profile={}
            self.load_seq()
        option = kwargs.get('option')
        if option:
            if option == 'A':
                self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq, option=option)
            elif option == 'Z':
                self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq, option=option)
            elif option == 'C':
                self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq, option=option)
            else:
                print("This option doesn't exist")
        else:
            self.SIST_profile=load_profile(basedir+"data/"+self.name+"/sequence.fasta", start, end, self.seq)

    """