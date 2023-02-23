# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
from useful_functions_Chipseq import *
from pathlib import Path
from globvar import *
from Genome import Genome


class Chipseq:

    def __init__(self, name, *args, **kwargs):
        """ 
        Called when an Chipseq instance is created,
        initialize the attribute name.
        Example:
        >>> ch = Chipseq.Chipseq("ecoli")
        >>> ch.name
        ecoli
        """
        self.name = name

    
    def load_signal(self, cond="all"):
        '''
        load_signal imports a 1D distribution along the chromosome (typically a 
        CHIPSeq distribution) from a data file (typically a .bedGraph obtained 
        with bamCompare) containing :
        (1) bin starting positions, (2) bin ending positions and (3) signal in each bin.

        Requirements :
        The data importation requires an signal.info file that contains the columns 
        positions of each information in the data file, in the following order :
        [0] condition, [1] filename, [2] separator used in the data file, 
        [3] bin_start, [4] bin_end, [5] signal
        
        Args:
            self (Chipseq instance)
            cond (list of str.): 
                    Selection of one or several conditions (1 condition corresponds 
                    to 1 data file). By default : cond ='all' ie all available 
                    signals are loaded. All selected conditions have to be 
                    listed in signal.info file.
        
        Outputs : creates or adds items to 2 Chipseq instance attributes
            self.signal (dict.) :
                   Dictionary of shape {condition : array containing one signal
                   value per genomic position}
            self.all_signals (dict.) :
                   Dictionary of shape {condition : array containing one signal
                   value per genomic position}
            The difference between these 2 attributes is that self.signal contains
            the signals loaded with this method only whereas self.all_signals
            contains the signals loaded with all Chipseq methods.

        Example:
        >>> ch = Chipseq.Chipseq("ecoli")
        >>> ch.load_signal()
        >>> ch.signal['Signal_Test']
        array([100,   0,   0, ..., 100, 100,   0])
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

                        # creates a new attribute to the Chipseq instance : 
                        # signal[cond]= signal per pos
                        if not hasattr(self, "signal"):
                            self.signal = {} 
                        self.signal[line[0]] = np.array(np.repeat(signal, bs))

                        if not hasattr(self, "all_signals"):
                            self.all_signals = {} 
                        self.all_signals[line[0]] = self.signal[line[0]] 
            f.close()

        else:
            print(f"Unable to locate signals.info in {path2dir}")

        # warns that some selected conditions aren't described in the
        # .info file
        for c in cond:
            if c not in list_cond_info and c != "all":
                print(f"Please add {c} in signals.info")


    def load_binned_signal(self, binsize, cond="all"):
        '''
        If a file containing the data for the chosen condition and binsize 
        exists, load_binned_signal loads these data using the load_signal method
        Else, load_binned signal performs the following steps: 
        1 - imports a 1D distribution along the chromosome (typically a CHIPSeq 
        distribution) using the load_signal method
        2 - performs the binning at the chosen binsize using the 
        useful_functions_Chipseq.binning
        3 - saves the binned data in a file 

        See ChipSeq.load_signal method for the data requirements.

        Args:
            binsize (int.): bin size 
            cond (list of str.): 
                    Selection of one or several conditions (1 condition corresponds 
                    to 1 data file). By default : cond ='all' ie all available 
                    signals are loaded. All selected conditions have to be 
                    listed in signal.info file.
        
        Outputs: creates or adds items to 3 Chipseq instance attributes
            self.signal (dict.) :
                   Dictionary of shape {condition : array containing one signal
                   value per genomic position (before binning)}
            self.binned_signal (dict.) :
                    Dictionary of shape {cond_bin : array containing one binned 
                    signal value per genomic position} with cond_bin the condition 
                    name merged with the bin size (example: WT_bin200b). 
            self.all_signals (dict.) :
                   Dictionary of shape {cond_bin : array containing one binned 
                    signal value per genomic position}
            The difference between binned_signal and all_signals attributes is that 
            binned_signal contains the signals loaded with this method only whereas 
            self.all_signals contains the signals loaded with all Chipseq methods.

        Example:
        >>> ch = Chipseq.Chipseq("ecoli")
        >>> ch.load_binned_signal(binsize=100,cond='Signal_Test')
        >>> ch.all_signals["Signal_Test_bin100b"]
        array([100,   100,   100, ..., 10, 10, 10])
        '''
        if not hasattr(self, "length"):
            gen = Genome(self.name)
            gen.load_seq() 
            self.length = gen.length

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
                self.binned_signal[cond_name] =np.array(np.repeat(binned_data, binsize))[:self.length]

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
                self.binned_signal[cond_name] = np.array(np.repeat(binned_data, binsize))[:self.length]

                # creates the "binned_data" directory if necessary
                Path(f_path).mkdir(parents=True, exist_ok=True)
                # saves the data as a .npy file
                np.save(f"{f_path}{cond_name}.npy", binned_data)

        if not hasattr(self, "all_signals"):
            self.all_signals = {} 
        self.all_signals[cond_name] = self.binned_signal[cond_name] 


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
        
        Args:
            window (int.): window size. The value of the smoothed signal 
                           of a position p is equal to the average of the 
                           signal between p - window/2 and p + window/2.
            cond (list of str.): 
                    Selection of one or several conditions (1 condition corresponds 
                    to 1 data file). By default : cond ='all' ie all available 
                    signals are loaded. All selected conditions have to be 
                    listed in signal.info file.


        Outputs: creates or adds items to 3 Chipseq instance attributes
            self.signal (dict.) :
                   Dictionary of shape {condition : array containing one signal
                   value per genomic position (before smoothing)}
            self.smoothed_signal (dict.) :
                    Dictionary of shape {cond_smoo : array containing one binned 
                    signal value per genomic position} with cond_smoo the condition 
                    name merged with the smoothing window size (example: 
                    WT_smooth200b).
            self.all_signals (dict.) :
                   Dictionary of shape {cond_smoo : array containing one binned 
                    signal value per genomic position}

            The difference between smoothed_signal and all_signals attributes is that 
            smoothed_signal contains the signals loaded with this method only whereas 
            self.all_signals contains the signals loaded with all Chipseq methods. 

        Example:
        >>> ch = Chipseq.Chipseq("ecoli")
        >>> ch.load_smoothed_signal(window=100,cond='Signal_Test')
        >>> ch.all_signals["Signal_Test_smooth100b"]
        array([100.1,   99.8,   98.8, ..., 10.1, 11.1, 10.8])
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

        if not hasattr(self, "all_signals"):
            self.all_signals = {} 
        self.all_signals[cond_name] = self.smoothed_signal[cond_name]

    def load_signals_average(self, list_cond, average_name, *args, **kwargs):
        '''
        Load_signals_average computes and loads the average of signals replicates.
        First, the function loads the signals (which must be listed in the signals.info 
        file in the /chipseq/signals/ folder) using the load_signal method. 
        It can then process the data using the load_binned_signal or 
        load_smoothed_signal methods. Finally, the average of these signals at each 
        genomic position of the genome is calculated. This average signal is assigned 
        to the Chipseq instance as signals_average attribute. 

        Required args. :
            list_cond (list of str.): 
                    Selection of conditions that will be averaged. All selected conditions 
                    have to be listed in signal.info file.
            average_name (str.): 
                    Name of the obtained signal

        Optionnal args: 
            data_treatment (str. "binning", "smoothing" or None):
                    Treatment to be applied to the different signals 
                    before averaging them (None by default).
            window (int.): 
                    window size used only if data_treatment = "smoothing"
                    The value of the smoothed signal of a position p is 
                    equal to the average of the signal between 
                    p - window/2 and p + window/2
            binsize (int.): bin size used only if data_treatment = "binning"

        See Chipseq.load_signal method for the data requirements

        Outputs: creates or adds items to the following Chipseq instance attributes
            self.signals_average (dict.) :
                    Dictionary of shape {average_name : array containing one averaged 
                    signal value per genomic position}.
            self.all_signals (dict.) :
                    Dictionary of shape {average_name : array containing one averaged 
                    signal value per genomic position}.
            The difference between signals_average and all_signals attributes is that 
            signals_average contains the signals loaded with this method only whereas 
            self.all_signals contains the signals loaded with all Chipseq methods. 

        Example:
        >>> ch = Chipseq.Chipseq("ecoli")
        >>> ch.load_signals_average(list_cond=["Signal1","Signal2"],
                                    average_name="Mean_signal",
                                    data_treatment = "smoothing",
                                    window=500)
        >>> ch.all_signals["Mean_signal"]
        array([100.1,   99.8,   98.8, ..., 10.1, 11.1, 10.8])
        >>> ch.signals_average["Mean_signal"]
        array([100.1,   99.8,   98.8, ..., 10.1, 11.1, 10.8])
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

        if not hasattr(self, "all_signals"):
            self.all_signals = {} 
        self.all_signals[average_name] = self.signals_average[average_name]

    def load_signal_per_genes(self, cond='all',window=0):
        """ 
        Load_signal_per_genes computes the mean signal for each Gene (i.e. 
        the mean signal between the gene start and gene end). Genomic signals
        have to be loaded (using load_signal,load_binned_signal, load_smoothed_signal
        or load_signals_average methods) before using load_signal_per_genes.

        Args: 
            self (Chipseq instance)
            cond (Optional [list of str.]): selection of one or several 
                                            conditions (1 condition corresponds
                                            to 1 data file).
                                            By default cond ='all' ie all 
                                            available loaded signals are used.
            window (Optional int.): to include the signal around the gene. The mean 
                                    signal will be computed between gene start - window
                                    and gene end + window.
                                    By default, window = 0
        
        Outputs: 
            self.genes[locus].signal (float.): new attribute of Gene instances 
                                               related to the Chipseq instance 
                                               given as argument.Contains the 
                                               Gene mean signal. 
            self.signals_gene (dict. of dict.): Dictionary containing one 
                                                subdictionary per condition 
                                                given in input. Each subdictionary
                                                contains the signal (value) for 
                                                each gene (locus_tag as key).
        N.B.: This method needs a genomic annotation. If no annotation is 
        loaded, the load_annotation method with the default "sequence.gff3" 
        file is computed. To use another annotation, please load an annotation to 
        your Transcriptome instance with the following commands before using 
        this method:
            >>> import Genome, Chipseq
            >>> ch = Chipseq.Chipseq("ecoli")
            >>> g = Genome.Genome(ch.name)
            >>> g.load_annotation(annot_file=chosen_file)
            >>> ch.genes = g.genes    

        Example:
        >>> ch = Chipseq.Chipseq("ecoli")
        >>> ch.load_signal()
        >>> ch.load_signal_per_genes()
        >>> ch.signals_gene["WT"]["b0984"]
        0.40348639903919425
        >>> ch.genes["b0984"].signal
        {'WT': 0.40348639903919425,
         'Signal_Test': 0.210505235030067}
        """

        if not hasattr(self, "signals_gene"):
            self.signals_gene = {}
        if not hasattr(self, "genes"): 
            gen = Genome(self.name)
            gen.load_annotation()
            self.genes = gen.genes
        
        if not hasattr(self, "signals_gene"):
            self.signals_gene = {cond_name: {}}

        if cond == "all" :
            cond =  list(self.all_signals.keys())
        elif isinstance(cond, str):
            cond = [cond]

        for c in cond :
            self.signals_gene[c] = {}
            for locus in self.genes.keys():
                g = self.genes[locus]
                s = np.mean(self.all_signals[c][g.left-window:g.right+window])
                self.genes[locus].add_signal(c, s)
                self.signals_gene[c][locus] = s
        print(f"Gene signals successfully computed for: {cond}")
      
    def load_peaks(self):  
        """ 
        load_peaks imports a list of peaks from a data file (typically 
        a .BED file of peaks obtained with MACS2) containing, at least,
        the start and end positions of peaks.
        
        Args: 
            self (Chipseq instance)
            
        Requirements:
            The data importation requires a peaks.info file that contains the 
                column indices of each information in the data file and some 
                additional information, in the following order :
                [0] Condition, [1] Filename, [2] Startline, 
                [3] Separator, [4] StartCol, [5] StopCol
            peaks.info and data file have to be in the /chipseq/peaks/ directory 
        
        Output: 
            self.peaks (dict.): new attribut of the Chipseq instance. 
                                self.peaks is a dictionary of shape
                                {cond:[(start,end)]}. It contains, for each 
                                condition, a list of tuples of shape (start,end).
                                One tuple represents one peak.
        """

        # gets the path to data and .info files
        path2dir = f"{basedir}data/{self.name}/chipseq/peaks/"

        # tries to open peaks.info file
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