# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
from useful_functions_Chipseq import *
from pathlib import Path
from globvar import *

class Chipseq:

    def __init__(self, name, *args, **kwargs):
        """ 
        Called when an Chipseq instance is created, 
        e.g. Chip = Chipseq.Chipseq("ecoli")
        Initialize the attribute name, e.g. self.name = "ecoli"
        """
        self.name = name

    def load_signal(self, cond="all"):
        '''
        load_signal imports a 1D distribution along the chromosome (typically a 
        CHIPSeq distribution) from a data file (typically a .bedGraph obtained 
        with bamCompare) containing :
        (1) bin starting positions, (2) bin ending positions and (3) signal in each bin.
        Creates 2 new attributes of the instance :
        + self.signal = {condition : array containing one value per genomic position}
        + self.signal_binsize = {condition : minimal binsize of the signal}.

        Requirements :
        The data importation requires an .info file that contains the columns 
        positions of each information in the data file, in the following order :
        [0] condition, [1] filename, [2] separator used in the data file, 
        [3] bin_start, [4] bin_end, [5] signal
        
        Option :
        A selection of one or several conditions (1 condition corresponds to 
        1 data file) can be made using 'cond' argument
        e.g. cond = ['cond1','cond2'] or cond = 'cond0', 
        by default : cond ='all' (all available signals are loaded).

        Output : 2 new attributes
        + self.signal = {cond:[signal per genomic position]}
        + signal_binsize[cond] = minimal binsize in the data file
        '''
        # gets the path to data and .info files
        path2dir = f"{basedir}data/{self.name}/chipseq/signals/"

        # converts from str to list (useful if only 1 condition is selected)
        if isinstance(cond, str):
            cond = [cond]

        # tries to open .info file
        list_cond_info = []  # list of the conditions described in the .info file
        if Path(f"{path2dir}signals.info").exists():
            with open(f"{path2dir}signals.info", "r") as f:
                header = next(f)

                # 1 line in the .info file corresponds to 1 condition 
                #(1 data file)
                for line in f:
                    line = line.strip().split('\t')
                    list_cond_info.append(line[0])

                    # For each selected condition :
                    if cond == ["all"] or line[0] in cond:
                        print('Loading condition', line[0])
                        # loads data
                        separator = line[2]
                        if separator == '\\t':
                            separator = '\t'
                        data = pd.read_csv(
                            path2dir + line[1], sep=separator, header=None)
                        signal = data.iloc[:, int(line[5])]

                        # computes binsize of each bin in the file
                        bs = data.iloc[:, int(line[4])] - \
                            data.iloc[:, int(line[3])]
                        bs.iloc[-1] = bs.iloc[-1]+1

                        # creates a new attribute to the Chipseq instance : 
                        # signal[cond]= signal per pos
                        if not hasattr(self, "signal"):
                            self.signal = {} 
                        self.signal[line[0]] = np.array(np.repeat(signal, bs))

                        # creates another attribute : signal_binsize[cond] = 
                        #minimal binsize in the data file
                        if not hasattr(self, "signal_binsize"):
                            self.signal_binsize = {}  
                        self.signal_binsize[line[0]] = min(bs)
            f.close()

        else:
            print(f"Unable to locate signals.info in {path2dir}")

        # warns that (one of) the selected condition(s) isn't described in the
        # .info file
        for c in cond:
            if c not in list_cond_info and c != "all":
                print(f"Please add {c} in signals.info")

    def load_binned_signal(self, binsize, cond="all"):
        '''
        If a data file containing the data for the chosen condition and binsize 
        exists, load_binned_signal loads these data using the load_signal method
        Else, load_binned signal performs the following steps: 
        1 - imports a 1D distribution along the chromosome (typically a CHIPSeq 
        distribution) using the load_signal method
        2 - performs the binning at the chosen binsize using the 
        useful_functions_Chipseq.binning
        3 - saves the binned data in a file 

        See ChipSeq.load_signal method for the data requirements.

        Arguments :
        + Binsize has to be given using 'binsize' argument
        + A selection of one or several conditions (1 condition corresponds to 
        1 data file) can be made using 'cond' argument
        e.g. cond = ['cond1','cond2'] or cond = 'cond0', 
        by default : cond ='all' (all available signals are loaded).

        Output : 
        a new attribute self.binned_signal = {cond_bin:[signal per genomic position]}   
        with cond_bin the condition name merged with the binsize (example: WT_bin200b)     
        '''
        # gets the path to data and .info files
        path2files = f"{basedir}data/{self.name}/chipseq/signals/"
        f_path = f"{path2files}/binned_data/"

        # gets the list of all existing conditions using .info file
        if Path(path2files + 'signals.info').exists() :
            if cond == "all":
                with open(path2files + 'signals.info', "r") as fi:
                    header = next(fi)
                    cond = []
                    for line in fi:  # for each condition
                        line = line.strip().split('\t')
                        cond.append(line[0])
                fi.close()

            # converts cond from str to list if only 1 condition is selected
            elif isinstance(cond, str):
                cond = [cond]

        else:
                print(f"Unable to locate signals.info in {path2files}")

        # For each selected condition :
        if not hasattr(self, "binned_signal"):
            self.binned_signal = {} # creates the new attribute 

        for c in cond:
            cond_name = f"{c}_bin{str(binsize)}b"
            print('Condition :' + cond_name)
            
            # if a file contains data for this condition and binsize 
            if Path(f"{f_path}{cond_name}.npy").exists():
                #adds the data to the binned_signal attribute
                print(f"loading existing file : {cond_name}.npy")
                binned_data = np.load(f"{f_path}{cond_name}.npy")
                self.binned_signal[cond_name] =np.array(np.repeat(binned_data, binsize))

            # if no file contains data for this condition and binsize
            else:
                print("performing the binning")
                # loads data
                self.load_signal(cond=c)
                # performs the binning using the 
                #useful_functions_Chipseq.binning function
                binned_data = binning(self.signal[c],binsize=binsize,
                    stat="mean").statistic
                # adds the data to the binned_signal attribute
                self.binned_signal[cond_name] = np.array(np.repeat(binned_data, binsize))

                # creates the "binned_data" directory if necessary
                Path(f_path).mkdir(parents=True, exist_ok=True)
                # saves the data as a .npy file
                np.save(f"{f_path}{cond_name}.npy", binned_data)

    def load_smoothed_signal(self, window, cond="all", *args, **kwargs):
        '''
        If a data file containing the data for the chosen condition and 
        smoothing exists, load_binned_signal loads these data using the 
        load_signal method
        Else, load_smoothed_signal performs the following steps: 
        1 - imports a 1D distribution along the chromosome (typically 
        a CHIPSeq distribution) using the load_signal method
        2 - performs the smoothing (moving average with the chosen 
        window size) using the useful_functions_Chipseq.smoothing
        3 - saves the smoothed data in a file 

        See Chipseq.load_signal method for the data requirements.
        
        Arguments :
        + Window size has to be given using 'window' argument
        + A selection of one or several conditions (1 condition corresponds to 
        1 data file) can be made using 'cond' argument
        e.g. cond = ['cond1','cond2'] or cond = 'cond0', 
        by default : cond ='all' (all available signals are loaded).

        Output : 
        a new attribute self.smoothed_signal = {cond_smoo:[signal per 
        genomic position]} with cond_smoo the condition name merged with 
        the smoothing window (example: WT_smooth200b)  
        '''
        # gets the path to data and .info files
        path2files = f"{basedir}data/{self.name}/chipseq/signals/"
        f_path = f"{path2files}/smoothed_data/"

        # gets the list of all existing conditions using .info file
        info_filename = kwargs.get('info_filename', 'signals.info')
        if Path(path2files + info_filename).exists():
            if cond == "all" : 
                with open(path2files + info_filename, "r") as fi:
                    header = next(fi)
                    cond = []
                    for line in fi:  # for each condition
                        line = line.strip().split('\t')
                        cond.append(line[0])
                fi.close()
            # converts cond from str to list if only 1 condition is selected
            elif isinstance(cond, str):
                cond = [cond]
        else:
            print("Unable to locate " + info_filename)
 
        if not hasattr(self, "smoothed_signal"):
            self.smoothed_signal = {} #creates the new attribute

        # For each selected condition :
        for c in cond:
            cond_name = f"{c}_smooth{str(window)}b"
            print('Condition :' + cond_name)

            # if a file contains data for this condition and window size
            if Path(f"{f_path}{cond_name}.npy").exists():
                # adds the data to the smoothed_signal attribute
                print(f"loading existing file : {cond_name}.npy")
                self.smoothed_signal[cond_name] = np.load(
                    f"{f_path}{cond_name}.npy")

            # if no file contains data for this condition and binsize
            else:
                print("performing the smoothing")
                # loads data
                self.load_signal(cond=c)
                # performs the smoothing using the 
                # useful_functions_Chipseq.smoothing function
                smooth_data = smoothing(
                    self.signal[c], window=window)
                # adds the data to the smoothed_signal attribute
                self.smoothed_signal[cond_name] = smooth_data
                # creates the "smoothed_data" directory if necessary 
                Path(f_path).mkdir(parents=True, exist_ok=True)
                # saves the data as a .npy file
                np.save(f"{f_path}{cond_name}.npy", smooth_data)

    def load_signals_average(self, list_cond, average_name, *args, **kwargs):
        '''
        Load_signals_average loads the average of signals replicates

        Requirements :
        + A selection of conditions that will be averaged has to be given
        with the "cond_name" argument, e.g. cond = ['cond1','cond2']
        + A name has to be specified for the output using the "average_name"
        argument
        + See Chipseq.load_signal method for the data requirements

        Option : 
        + Binning can be performed on the signals (before average) by specifying 
        data_treatment="binning" and using an additional argument "binsize"
        + Smoothing can be performed on the signals (before average) by 
        specifying data_treatment="smoothing" and using an additional argument 
        "window"
        + data_treatment = "binning", "smoothing" or None (None by default) 

    
        Output : 
        a new attribute self.signal_average[average_name] = [average signal 
        per genomic position]
        '''

        f_path = f"{basedir}data/{self.name}/chipseq/signals/average_data/"
        data_treatment = kwargs.get('data_treatment', None)
        if not hasattr(self, "signals_average"):
            self.signals_average = {}  
        # test if a file for this condition (cond_name) and this binsize
        # already exists

        if Path(f"{f_path}{average_name}.npy").exists():
            with open(f"{f_path}data.info", "r") as f:
                header = next(f).strip("\n").split("\t")
                # 1 line in the .info file corresponds to 1 condition 
                for line in f:
                    line = line.strip().split('\t')
                    if line[0] == average_name :
                        print("loading the file obtained with the following parameters:\n")
                        for i in np.arange(len(header)):
                            print(f"{header[i]}: {line[i]}\n")

                        print("Please change \"average_name\" to use other signals")
            self.signals_average[average_name] = np.load(
                f"{f_path}{average_name}.npy")

        else:
            print("performing the average")
            # get the list of replicates filenames (files containing the
            # experimental name given as argument)

            if data_treatment == "binning":
                print("data_treatment : binning")
                try : 
                    size = kwargs.get('binsize')  
                except :
                    sys.exit("please give a binsize") 
                self.load_binned_signal(size, cond=list_cond)
                data_to_av = [self.binned_signal[f"{i}_bin{str(size)}b"] for i in list_cond]

            elif data_treatment == "smoothing":
                print("data_treatment : smoothing")
                try : 
                    size = kwargs.get('window')  
                except :
                    sys.exit("please give a window size using \"window\" argument")
                self.load_smoothed_signal(size, cond=list_cond)
                data_to_av = [self.smoothed_signal[f"{i}_smooth{str(size)}b"] for i in list_cond]
            else:
                data_treatment = "No treatment"
                size = "NA"
                print("data_treatment : No treatment")
                self.load_signal(cond=list_cond)
                data_to_av = [self.signal[i] for i in list_cond]

            # performs the average between replicates
            average = np.mean(data_to_av, axis=0)       

            # adds the signals_average attribute to Chipseq instance
            self.signals_average[average_name] = average

            # saves the data as a .npy file
            Path(f_path).mkdir(parents=True, exist_ok=True)
            np.save(f"{f_path}{average_name}.npy", average)

            # saves information about the used replicates and treatment
            if not Path(f"{f_path}data.info").exists():
                file = open(f"{f_path}data.info",'w')
                file.write('Name\t Replicates\t Data treatment\t Size (window or bin, in b)')
                file.close()
            file = open(f"{f_path}data.info",'a')
            file.write(f"\n{average_name}\t{list_cond}\t{data_treatment}\t{size}")
            file.close()

#=============================================================================#        
    def load_peaks(self):  
        """ 
        load_peaks imports a list of peaks from a data file (typically 
        a .BED file of peaks obtained with MACS2) containing, at least,
        the start and end positions of peaks.

        Requirements:
        + The data importation requires a peaks.info file that contains the 
        column indices of each information in the data file and some additional 
        information, in the following order :
        [0] Condition, [1] Filename, [2] Startline, 
        [3] Separator, [4] StartCol, [5] StopCol
        + peaks.info and data file have to be in the /chipseq/peaks/ directory 
        
        Output: 
            self.peaks (dict.): new attribut of the Chipseq instance. 
                                self.peaks is a dictionary of shape
                                {cond:[(start,end)]}. It contains, for each 
                                condition, a list of tuples of shape (start,end).
                                One tuple represents one peak.
        """

        # gets the path to data and .info files
        path2dir = f"{basedir}data/{self.name}/chipseq/peaks/"

        # tries to open paks.info file
        if Path(f"{path2dir}peaks.info").exists():
            with open(f"{path2dir}peaks.info", "r") as f:
                skiphead = next(f)  
                # 1 line in peaks.info corresponds to 1 condition (1 data file)
                for line in f:
                    line = line.strip('\n').split('\t')

                    # loads data file information for this condition
                    cond = line[0]
                    path2file = f"{path2dir}{line[1]}"
                    startline, sep = int(line[2]), line[3]
                    start_col, end_col  = int(line[4]), int(line[5])

                    # creates the new attribute using the load_sites_cond function
                    # from useful_functions_Chipseq 
                    if not hasattr(self, "peaks"):
                        self.peaks = {}  # creates the new attribute 
                    self.peaks[cond] = load_sites_cond(path2file, startline, sep, 
                                                       start_col, end_col)
        else:
            print("No peaks.info, unable to load peaks")