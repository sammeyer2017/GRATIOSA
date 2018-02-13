#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Vincent CABELI
"""

import sys
import os
import numpy as np
import scipy.spatial.distance
from scipy.optimize import curve_fit
from numpy.lib.stride_tricks import as_strided
from random import shuffle
from Domain import Domain


def periodic_dist(genes, max_size):
    """ Computes the periodic distance matrix of points on a circle of size
    max_size. 
    """
    pos_left = np.transpose([[x.left for x in genes]])
    pos_right = np.transpose([[x.right for x in genes]])
    pos_middle = np.transpose([[x.middle for x in genes]])
    distance_matrix = scipy.spatial.distance_matrix(pos_middle, pos_middle)
    #lambda u,v
    distance_matrix = np.where(distance_matrix > 0.5*max_size, 
        max_size-distance_matrix, distance_matrix)
    return distance_matrix


def change_dict_keys(genes_dict):
    """ Change a genes dictionnary's keys to the gene names.
    """
    for tag in list(genes_dict.keys()):
        gene = genes_dict[tag]
        genes_dict[gene.name] = gene
        del genes_dict[tag]
    return genes_dict


def gen_successive_domains(genes_dict, k):
    """ Generates overlapping domains of adjacent genes on the whole genome. It
    will return a list of n_genes-k domains.
    """
    sorted_genes = list(genes_dict.values())
    sorted_genes.sort(key = lambda x: x.start)
    domains = []
    for i in range(len(sorted_genes)-k):
        if all([hasattr(gene, 'mean_expression') for gene in sorted_genes[i:i+k]]):
            domains.append(Domain(sorted_genes[i:i+k]))
    return domains


def gen_random_domains(sorted_genes, n, domains=None, b=0.575):
    """ Generates a list of n random domains of adjacent genes. Uses an
    exponential law to determine the length of each domain, with the formula 
    for exponential decay exp(-b * x). Default = 0.575 as observed in the 
    synteny segments size distribution. Will determine a new b parameter by
    fitting an exponential curve if a domains list is specified.
    """
    n_genes = len(sorted_genes)
    # consider max 30 genes domains
    max_size = 30
    if domains:
        t = [len([domain for domain in domains if 
                  domain.number_of_genes==i]) for i in range(2, max_size)]
        popt, pcov = curve_fit(exp_to_fit, np.asarray(list(range(2, max_size))), t)
        b = popt[1]
    probs = np.exp(-b * np.asarray(list(range(max_size))))
    rand_domains = []
    rand_domains_ext = []
    for i in range(n):
        gene_pos = int(np.random.uniform(1,n_genes-max_size))
        size = np.max(np.where(np.random.rand() < probs)) + 2
        rand_domains.append(Domain(sorted_genes[gene_pos:gene_pos+size]))
        rand_domains_ext.append(Domain(sorted_genes[gene_pos-1:gene_pos+size+1]))
    return (rand_domains, rand_domains_ext)


def gen_random_segment_like_domains(sorted_genes, n, maps, domains,
                                    successions_types=None, 
                                    successions_types_counts=None):
    """ Generates a list of len(domains) random domains with the same number of
    convergent, divergent and same sense gene successions as those observed in
    the given list of domains. 
    """
    rand_domains = []
    rand_domains_ext = []
    if successions_types and successions_types_counts:
        for i, successions_type in enumerate(successions_types):
            n_genes = sum(successions_type)+1
            candidates = np.where([(list(row)==successions_type) for row in maps[n_genes-2]])[0]
            positions = np.random.choice(candidates, size=successions_types_counts[i])
            for pos in positions:
                rand_domains.append(Domain(sorted_genes[pos:pos+n_genes]))
                rand_domains_ext.append(Domain(sorted_genes[pos-1:pos+n_genes+1]))
        return (rand_domains, rand_domains_ext)
    else:
        for domain in domains:
            domain_successions = [count_successions([domain], np.asarray([True,True])) +
                                  count_successions([domain], np.asarray([False,False])),
                                  count_successions([domain], np.asarray([False,True])),
                                  count_successions([domain], np.asarray([True,False]))]
            candidates = np.where([(list(row)==domain_successions) for row in maps[domain.number_of_genes-2]])[0]
            pos = np.random.choice(candidates)
            rand_domains.append(Domain(sorted_genes[pos:pos+domain.number_of_genes]))
            rand_domains_ext.append(Domain(sorted_genes[pos-1:pos+domain.number_of_genes+1]))
        return (rand_domains, rand_domains_ext)


def count_successions_types(domains):
    seen_successions = []
    seen_successions_count = []
    for domain in domains:
        domain_successions = [count_successions([domain], np.asarray([True,True])) +
                              count_successions([domain], np.asarray([False,False])),
                              count_successions([domain], np.asarray([False,True])),
                              count_successions([domain], np.asarray([True,False]))]
        if domain_successions not in seen_successions:
            seen_successions.append(domain_successions)
            seen_successions_count.append(1)
        else:
            seen_successions_count[seen_successions.index(domain_successions)] += 1
    return (seen_successions, seen_successions_count)


def convergent_significance(genes_dict, domains, n):
    """ Generates n times random domains list of the same size as the list
    passed as argument. Returns a tuple with the number of observed convergences
    as well as the number of convergences observed only in direct domain 
    neighbours.
    """
    sorted_genes = list(genes_dict.values())
    sorted_genes.sort(key = lambda x: x.left)
    n_convergents = []
    diff = []
    # fit exp and get parameter b
    max_size = max([x.number_of_genes for x in domains])
    print('Fitting exponential...')
    domains_sizes = [len([domain for domain in domains if 
                     domain.number_of_genes==i]) for i in range(2, max_size)]
    popt, pcov = curve_fit(exp_to_fit, np.asarray(list(range(2, max_size))), domains_sizes)
    b = popt[1]
    print('Generating maps...')
    maps = gen_maps(sorted_genes, max_size)
    n_domains = len(domains)
    print('Generating domains...')
    (successions_types, successions_types_counts) = count_successions_types(domains)
    for i in range(n):
        print(i)
        #(d, d_ext) = gen_random_domains(sorted_genes, n_domains, b=b)
        (d, d_ext) = gen_random_segment_like_domains(sorted_genes, 1, maps, 
                                                     domains, 
                                                     successions_types,
                                                     successions_types_counts)
        n_convergents.append(count_successions(d, np.array([False, True])))
        diff.append(count_successions(d_ext, np.array([False,True])) - n_convergents[i])
    return (n_convergents, diff)


def count_successions(domains, succession):
    """ Counts the successions of 
    """
    n = 0
    b = np.size(succession)
    for domain in domains:
        for i in range(domain.number_of_genes-1):
            if np.all(succession == domain.genes_orientations[i:i+b]):
                n += 1
    return n
        

def exp_to_fit(x, a, b):     
    return a * np.exp(-b * x)


def gen_maps(sorted_genes, max_size):
    """ Generates a [max_size-1, len(sorted_genes)-max_size, 3]-shaped matrix 
    that contains a map of
    """
    maps = np.zeros((max_size, len(sorted_genes)-max_size, 3), dtype='int')
    for size in range(max_size):
        # allow for a 1-gene displacement for extended domains
        for pos in range(1, len(sorted_genes)-max_size-1):
            domain = Domain(sorted_genes[pos:pos+size+2])
            maps[size,pos,:] = [count_successions([domain], np.asarray([True,True])) + 
                                count_successions([domain], np.asarray([False,False])), 
                                count_successions([domain], np.asarray([False,True])),
                                count_successions([domain], np.asarray([True,False]))]
    return maps


def compare_domains_lists(domains1, domains2):
    """ Compares two domains lists and returns
    """
    max_size = np.max([x.number_of_genes for x in domains1] + 
                      [x.number_of_genes for x in domains2])
    equalities = np.zeros(max_size, dtype='int')
    for domain1 in domains1:
        for domain2 in domains2:
            if domain1 == domain2 or domain1.includes(domain2) or domain2.includes(domain1):
                equalities[domain1.number_of_genes-1] += 1
                break
    return equalities


def build_domains_list_from_matrix(genes_dict, domains_matrix):
    """ Builds a list of Domain objects from a square matrix that is returned
    by the Dickeya_parser.domains_matrix() function. 
    """
    if (len(np.shape(domains_matrix)) != 2) or (
        np.shape(domains_matrix)[0] != np.shape(domains_matrix)[1]):
        print('The matrix is not a square matrix, cannot build list of domains')
        return

    genes = list(genes_dict.values())
    domains = []
    i = 0
    j = 0
    current_domain = []
    seen_indices = []
    for i in range(np.shape(domains_matrix)[1]):
        if i in seen_indices:
            continue
        for j in range(i, np.shape(domains_matrix)[1]):
            if domains_matrix[i,j]:
                current_domain.append(genes[j])
                seen_indices.append(j)
        if current_domain:
            domains.append(Domain(current_domain))
            current_domain = []
    return domains


# ------------


from itertools import groupby

def cov_is_zero(cov,winlen=20,maxval=5):
    """ cov analysis: extract positions where cov is zero
    criterium: maximum value of winlen consecutive basepairs
    """
    wherezero = cov<=maxval
    a=[(k,len(list(g))) for k,g in groupby(wherezero)]
    q=np.zeros(len(wherezero),dtype=bool)
    index=0
    for (k,le) in a:
        if k==True and le>=winlen:
            q[index:index+le]=True
        index+=le
    return q
