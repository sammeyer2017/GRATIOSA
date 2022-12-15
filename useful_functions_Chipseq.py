#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from scipy import stats
#from useful_functions import *
#from globvar import *

#==============================================================================#

def binning(data, binsize,stat = "mean"):
    """ 
     #AJOUT COMMENTAIRE SUR POS BIN !!!!!!!!!!!!!!!!!!
   data binning at a chosen binsize
           N.B.: The position (pos) is the bin number. If data are binned at 2kb,
        the bin with number 10 corresponds to data between 20000 and 22000. ?????????????? A VERIFIER 
   input has to be an array 
   stat can be : mean, std, median, count, sum, min, max, user-defined function (see binned_statistic documentation)
    """
    #determines the bin coordinates
    l_data = len(data)
    bins = np.arange(0,l_data,binsize)
    bins = np.append(bins,l_data)

    #performs the binning 
    binned_data = stats.binned_statistic(np.arange(l_data),data,bins=bins,statistic=stat)
    
    return binned_data


def average_replicates(data,binsize = None) :
    """ 
    compute the average between replicates
    data is a list of np.arrays with same lengths
    """
    binned_data = []

    for d in data :
        d_new = binning(d,binsize,stat=stat).statistic
        binned_data.append(d_new)

    average = np.mean(binned_data,axis=0)

    return average


def smoothing(data,window) : 
    """ 
    Data smoothing by moving average with the window size given as argument
    """ 
    s_data = []
    i=0
    w = int(window/2)
    #circular DNA
    while i < w :
        window_average = np.mean(list(data[i-w:])+list(data[0:i+w+1]))
        s_data.append(window_average)
        i += 1
    while i < len(data) - w :
        window_average = np.mean(data[i-w:i+w+1])
        s_data.append(window_average)
        i += 1
    while i < len(data) :
        window_average = np.mean(list(data[i-w:])+list(data[0:i+w-len(data)+1]))
        s_data.append(window_average)
        i += 1

    return s_data

def load_sites_cond(path2file, startline, sep, start_col, end_col):
    ''' 
    '''
    peaks = []
    if sep == '\\t': sep = '\t' 

    with open(path2file, 'r') as f:
        i=1
        while i < startline:
            header=next(f)
            i+=1
        for line in f:
            line = line.strip('\n').split(sep)
            try: 
                peaks.append((int(line[start_col]),int(line[end_col])))
            except Exception as e:
                print(e)
    f.close()

    return peaks