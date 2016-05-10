
"""
    Prediction TODO list:
        1. Validate using the Danish 23andme genomes.
        2. Validate using 1K genomes.
        3. Report predicted percentile among Danish high-school students.
        4. Report predicted percentile among 1K genomes individuals.
        5. Support more genotype platform (and trait combinations).
"""

import logging
from plinkio import plinkfile
import itertools as it
import h5py
import scipy as sp
from scipy import stats
import pandas
import random
import gzip
log = logging.getLogger(__name__)


                   
def predict(indiv_genot,trait,genotype_version,
                   sex=None,pcs=None,**kwargs):
    """
    predict height of an individual.
    
    sex: 1 for male, and 2 for female
    """
    (snp_weights_file,prs_weigths_file) = get_weight_files(trait,genotype_version)

    log_extra = kwargs.get('log_extra',{'progress':0})
    partial_progress_inc = (100-log_extra['progress'])/22
    log.info('Starting prediction for %s' % trait,extra=log_extra)

    h5f_ig = h5py.File(indiv_genot)
    swh5f = h5py.File(snp_weights_file)
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
    
    pwh5f = h5py.File(prs_weigths_file)
    if sex is not None:
        weights = pwh5f['sex_adj']
        pred_phen = weights['Intercept'][...]+weights['ldpred_prs_effect'][...]*prs + weights['sex'][...]*sex
    else:
        weights = pwh5f['unadjusted']
        pred_phen = weights['Intercept'][...]+weights['ldpred_prs_effect'][...]*prs
    log.info('Finished prediction',extra=log_extra)
    return pred_phen
    

def get_weight_files(trait,genotype_version):
    """
    Here we simply need to use the right snp_weights_file for each (trait,genot_chip_version)
    """
    if trait=='height':
        if genotype_version=='23andmev4':
            snp_weights_file='/faststorage/project/TheHonestGene/prediction_data/weight_files/height_weights.hdf5',
            prs_weigths_file='/faststorage/project/TheHonestGene/prediction_data/weight_files/prs_weights.hdf5',
        else:
            raise NotImplementedError
    else:
        raise NotImplementedError
    
    return (snp_weights_file, prs_weigths_file)



