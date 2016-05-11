
"""
    Prediction TODO list:
        1. Validate using the Danish 23andme genomes.
        2. Validate using 1K genomes.
        3. Report predicted percentile among Danish high-school students.
        4. Report predicted percentile among 1K genomes individuals.
        5. Support more genotype platform (and trait combinations).
"""

import logging
import itertools as it
import h5py
import scipy as sp
from scipy import stats
import pandas
import random
import gzip
import os
log = logging.getLogger(__name__)


                   
def predict(indiv_genot,trait_folder,
                   sex=None,pcs=None,**kwargs):
    """
    predict height of an individual.
    
    sex: 1 for male, and 2 for female
    """
    indiv_genot = os.path.abspath(indiv_genot)
    trait_folder = os.path.abspath(trait_folder)
    snp_weights_file = '%s/snp_weights.hdf5' % trait_folder
    prs_weigths_file = '%s/prs_weights.hdf5' % trait_folder

    log_extra = kwargs.get('log_extra',{'progress':0})
    partial_progress_inc = (100-log_extra['progress'])/22
    log.info('Starting prediction for %s' % trait_folder)

    h5f_ig = h5py.File(indiv_genot)
    swh5f = h5py.File(snp_weights_file,'r')
    prs = 0
    for chrom_i in range(1,23):
        chr_str = 'Chr%d'%chrom_i
        log_extra['progress']+=partial_progress_inc
        log.info('Computing weights for Chr %s'% chrom_i, extra=log_extra)
        snps = h5f_ig[chr_str]['snps'][...]
        snp_weights = swh5f[chr_str]['ldpred_betas'][...] #These are on a per-allele scale.
        prs += sp.dot(snp_weights,snps)
    h5f_ig.close()
    swh5f.close()
    pwh5f = h5py.File(prs_weigths_file,'r')
    if sex is not None:
        log.info('Calculating final prediction score for gender %s' % ('male' if sex == 1 else 'female') )
        weights = pwh5f['sex_adj']
        pred_phen = weights['Intercept'][...]+weights['ldpred_prs_effect'][...]*prs + weights['sex'][...]*sex
    else:
        log.info('Calculating final prediction score for unknown gender')
        weights = pwh5f['unadjusted']
        pred_phen = weights['Intercept'][...]+weights['ldpred_prs_effect'][...]*prs
    log.info('Finished prediction',extra=log_extra)
    return pred_phen



