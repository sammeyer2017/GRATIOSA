#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions called by Chipseq methods
"""

import numpy as np
from scipy import stats

def binning(data, binsize, stat="mean"):
    """
    Data binning at a chosen binsize.

    Args:
        data (array-like): A sequence of values to be binned.
        binsize (int.): The bin size
        stat (str.): The statistic to compute ('mean' by default). The 
                following statistics are available: mean, std, median, count, 
                sum, min, max. See scipy.stats.binned_statistic documentation 
                for more details)

    Returns:
        Array containing the values of the selected statistic in each bin.
    """
    # determines the bin coordinates
    l_data = len(data)
    bins = np.arange(0, l_data, binsize)
    bins = np.append(bins, l_data)

    # performs the binning
    binned_data = stats.binned_statistic(
        np.arange(l_data), data, bins=bins, statistic=stat)

    return binned_data


def smoothing(data, window):
    """
    Data smoothing by moving average of CIRCULAR data.

    Args:
        data (array-like): A sequence of values to be smoothed.
        window (int.): window size. The value of the smoothed signal of a 
                position p is equal to the average of the signal between 
                p - window/2 and p + window/2.

    Returns:
        Array containing the smoothed data.
    """
    s_data = []
    i = 0
    w = int(window / 2)
    # circular DNA
    while i < w:
        window_average = np.mean(list(data[i - w:]) + list(data[0:i + w + 1]))
        s_data.append(window_average)
        i += 1
    while i < len(data) - w:
        window_average = np.mean(data[i - w:i + w + 1])
        s_data.append(window_average)
        i += 1
    while i < len(data):
        window_average = np.mean(
            list(data[i - w:]) + list(data[0:i + w - len(data) + 1]))
        s_data.append(window_average)
        i += 1
    return s_data


def load_sites_cond(path2file, startline, sep, start_col, end_col, val_col=None):
    '''
    Called by load_peaks, allows the peaks data to be loaded by specifying
    files (typically a .BED file of peaks obtained with MACS2), and where 
    each information is (one column for each type of information).

    Args:
        path2file (str.): path to the file containing the data
        startline (int.): file start line
        sep (str.): file separator
        start_col (int.): index of the column containing the position of the
                beginning of the peak
        end_col (int.): index of the column containing the position of the 
                end of the peak
        val_col (int.): index of the column containing the position of the 
                value associated to the peak (default: None)

    Returns:
        List of tuples of shape (start,end). One tuple represents one peak.

    Note: 
        Column numbering starts at 0.
    '''
    peaks = {}
    if sep == '\\t':
        sep = '\t'

    with open(path2file, 'r') as f:
        i = 0
        while i < startline:
            header = next(f)
            i += 1
        for line in f:
            line = line.strip('\n').split(sep)
            try:
                if val_col != None :
                    peaks[(int(line[start_col]), int(line[end_col]))] = float(line[val_col])
                else : 
                    peaks[(int(line[start_col]), int(line[end_col]))] = None
            except Exception as e:
                print(e)
    f.close()

    return peaks