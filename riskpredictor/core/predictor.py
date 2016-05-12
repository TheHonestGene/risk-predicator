
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
import imputor
import sys


                   
def predict(indiv_genot,trait_folder,
                   sex=None,pcs=None,**kwargs):
    """
    predict height of an individual.
    
    sex: 1 for male, and 2 for female
    """
    indiv_genot = os.path.abspath(indiv_genot)
#     trait_folder = os.path.abspath(trait_folder)
    snp_weights_file = '%s/snp_weights.hdf5' % trait_folder
    prs_weigths_file = '%s/prs_weights.hdf5' % trait_folder
    print snp_weights_file

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
        log.info('Calculating final prediction score without gender information')
        weights = pwh5f['unadjusted']
        pred_phen = weights['Intercept'][...]+weights['ldpred_prs_effect'][...]*prs
    log.info('Finished prediction',extra=log_extra)
    return pred_phen




def validate_predictions(K=10):
    """
    Use the 23andme genomes to validate prediction, by pulling them through the pipeline 
    """
    #Pull out K individuals with phenotypes
    pred_res = pandas.read_csv('/home/bjarni/TheHonestGene/faststorage/prediction_data/weight_files/23andme_v4_height_prs.txt', 
                               skipinitialspace=True)
    pred_res = pred_res[:K]
    pred_phens = []
    for indiv_id in pred_res['IID']:
        input_file = '/home/bjarni/TheHonestGene/faststorage/prediction_data/23andme-genome/%s.genome'%indiv_id
        assert os.path.isfile(input_file), 'Unable to find file: %s'%input_file
        output_file = '/home/bjarni/TheHonestGene/faststorage/prediction_data/23andme-genomes_imputed/%s.genome.hdf5'%indiv_id
        if not os.path.isfile(output_file):
            args = {'input_file':input_file, 'output_file':output_file}
            imputor.parse_genotype(args)
        input_file = output_file
        output_file = '/home/bjarni/TheHonestGene/faststorage/prediction_data/23andme-genomes_imputed/%s.genome_converted.hdf5'%indiv_id
        if not os.path.isfile(output_file):
            args = {'input_file':input_file, 'output_file':output_file, 
                    'nt_map_file':'/home/bjarni/TheHonestGene/faststorage/data_for_pipeline/NT_DATA/23andme_v4_nt_map.pickled.new'}
            imputor.convert(args)
        input_file = output_file
        output_file = '/home/bjarni/TheHonestGene/faststorage/prediction_data/23andme-genomes_imputed/%s.genome_imputed.hdf5'%indiv_id
        if not os.path.isfile(output_file):
            args = {'genotype_file':input_file, 'output_file':output_file, 
                    'ld_folder':'/home/bjarni/TheHonestGene/faststorage/data_for_pipeline/LD_DATA/23andme_v4',
                    'validation_missing_rate':0.01,
                    'min_ld_r2_thres':0.05}
            imputor.impute(args)
        
        
        pred_phen = predict(output_file,'/home/bjarni/TheHonestGene/faststorage/data_for_pipeline/PRED_DATA/height')
        pred_phens.append(pred_phen)
    
    print sp.corrcoef(pred_phen,pred_res['true_phens'])
    
        
    
    
