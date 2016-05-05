"""

Risk prediction pre-calculation pipeline:
     0. Parse relevant summary statistics into HDF5 files, and prepare prediction genotype maps (for each supported genotype data format).
         a) Sum stats into HDF5
         b) LD ref file
         c) ...
     1. Identify the set of SNPs across: 
         a) the summary statistics 
         b) the LD-reference panel 
         c) the validation data set 
         d) the prediction genotype map.
     2. Run LDpred using the summary statistics and the LD-reference panel.  
     3. Validate using the danish High school students
     4. Set up prediction for a 23andme genotype, and store weights in a trait_file
         - weights for SNPs that we ignore will be set to 0.
         - 
    
When the resulting SNP weight files for each trait are available, we can use them to obtain the polygenic scores for each 
trait relatively easily using 

"""



import logging
#from plinkio import plinkfile
import itertools as it
import h5py
import scipy as sp
from scipy import stats
import random



def parse_BMI_HEIGHT():
    bmi_file = '/home/bjarni/TheHonestGene/faststorage/prediction_data/SNP_gwas_mc_merge_nogc.tbl.uniq'
    height_file = '/home/bjarni/TheHonestGene/faststorage/prediction_data/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt'
    KGpath = '/home/bjarni/TheHonestGene/faststorage/1Kgenomes/'
    comb_hdf5_file = '/home/bjarni/TheHonestGene/faststorage/prediction_data/HEIGHT_BMI.hdf5'
    bimfile = '/home/bjarni/TheHonestGene/faststorage/prediction_data/wayf.bim'
    parse_sum_stats(bmi_file,comb_hdf5_file,'BMI',KGpath,bimfile=bimfile)

#     parse_sum_stats(height_file,comb_hdf5_file,'height',KGpath,bimfile=bimfile)



                   
log = logging.getLogger(__name__)


def predict(genotype_file,trait,**kwargs):
    log_extra = kwargs.get('log_extra',{'progress':0})
    partial_progress_inc = (100-log_extra['progress'])/22
    log.info('Starting prediction for %s' % trait,extra=log_extra)
    for i in range(23):
        log_extra['progress']+=partial_progress_inc
        log.info('Computing weights for Chr%s'% i,extra=log_extra)
        import time
        time.sleep(1)
    log.info('Finished prediction',extra=log_extra)
    import random
    return random.random()



#-----------------------------Code for parsing summary statistics (of various formats) ----------------------------------------------
lc_2_cap_map = {'a':'A', 'c':'C', 'g':'G', 't':'T'}

headers = {'SSGAC1':['MarkerName', 'Effect_Allele', 'Other_Allele', 'EAF', 'Beta', 'SE', 'Pvalue'],
           'SSGAC2':['MarkerName', 'Effect_Allele', 'Other_Allele', 'EAF', 'OR', 'SE', 'Pvalue'],
           'CHIC':['SNP', 'CHR', 'BP', 'A1', 'A2', 'FREQ_A1', 'EFFECT_A1', 'SE', 'P'],
           'GCAN':['chromosome', 'position', 'SNP', 'reference_allele', 'other_allele', 'eaf', 'OR', 
                             'OR_se', 'OR_95L', 'OR_95U', 'z', 'p_sanger', '_-log10_p-value', 'q_statistic', 
                             'q_p-value', 'i2', 'n_studies', 'n_samples', 'effects'],
           'TESLOVICH':['MarkerName', 'Allele1', 'Allele2', 'Weight', 'GC.Zscore', 'GC.Pvalue', 'Overall', 'Direction'],
           'GIANT1':['MarkerName', 'Allele1', 'Allele2', 'FreqAllele1HapMapCEU', 'b', 'se', 'p', 'N'],
           'GIANT1b':['MarkerName', 'Allele1', 'Allele2', 'Freq.Allele1.HapMapCEU', 'b', 'SE', 'p', 'N'],
           'GIANT1c':['MarkerName', 'Chr', 'Pos', 'Allele1', 'Allele2', 'FreqAllele1HapMapCEU', 'b', 'se', 'p', 'N'],
           'GIANT2':['SNP', 'A1', 'A2', 'Freq1.Hapmap', 'b', 'se', 'p', 'N'],
           'MAGIC':['snp', 'effect_allele', 'other_allele', 'maf', 'effect', 'stderr', 'pvalue'],
           'CARDIoGRAM':['SNP', 'chr_pos_(b36)', 'reference_allele', 'other_allele', 'ref_allele_frequency', 'pvalue', 'het_pvalue', 'log_odds', 'log_odds_se', 'N_case', 'N_control', 'model'],
           'DIAGRAM':['SNP', 'CHROMOSOME', 'POSITION', 'RISK_ALLELE', 'OTHER_ALLELE', 'P_VALUE', 'OR', 'OR_95L', 'OR_95U', 'N_CASES', 'N_CONTROLS'],
           'TAG':['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A', 'FRQ_U', 'INFO', 'OR', 'SE', 'P'],
           'CD':['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A_5956', 'FRQ_U_14927', 'INFO', 'OR', 'SE', 'P', 'Direction', 'HetISqt', 'HetPVa'],
           'UC':['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A_6968', 'FRQ_U_20464', 'INFO', 'OR', 'SE', 'P', 'Direction', 'HetISqt', 'HetPVa'],
           'GEFOS':['chromosome', 'position', 'rs_number', 'reference_allele', 'other_allele', 'eaf', 'beta', 'se', 'beta_95L', 'beta_95U', 'z', 'p-value', '_-log10_p-value', 'q_statistic', 'q_p-value', 'i2', 'n_studies', 'n_samples', 'effects'],
           'RA':['SNPID','Chr','Position(hg19)','A1','A2','OR(A1)','OR_95%CIlow','OR_95%CIup','P-val'],
           'ASTHMA':['Chr', 'rs', 'position', 'Allele_1', 'Allele_2', 'freq_all_1_min', 'freq_all_1_max', 'OR_fix', 'ORl_fix', 'ORu_fix', 'P_fix'],
           'ICBP': ['ID', 'Analysis', 'ID', 'SNP', 'ID', 'P-value', 'Rank', 'Plot', 'data', 'Chr', 'ID', 'Chr', 'Position', 'Submitted', 'SNP', 'ID', 'ss2rs', 'rs2genome', 'Allele1', 'Allele2', 'Minor', 'allele', 'pHWE', 'Call', 'Rate', 'Effect', 'SE', 'R-Squared', 'Coded', 'Allele', 'Sample', 'size', 'Bin', 'ID']
           }

def get_sid_pos_map(sids, KGenomes_prefix):
    h5fn = '%ssnps.hdf5'%(KGenomes_prefix)
    h5f = h5py.File(h5fn,'r')
    sid_map = {}
    for chrom_i in range(1,23):
        cg = h5f['chrom_%d' % chrom_i]
        sids_1k = cg['sids'][...]
        sids_filter_1k = sp.in1d(sids_1k, sp.array(sids)) 
        common_sids = sids_1k[sids_filter_1k]
        common_positions = cg['positions'][sids_filter_1k]
        eur_mafs = cg['eur_mafs'][sids_filter_1k]
        for sid,pos,eur_maf in it.izip(common_sids,common_positions,eur_mafs):
            sid_map[sid]={'pos':pos, 'chrom':chrom_i, 'eur_maf':eur_maf}
    return sid_map


def parse_sum_stats(filename,
                    comb_hdf5_file,
                    ss_id,
                    KGpath,
                    bimfile=None,):
    """
    """
    h5f = h5py.File(comb_hdf5_file)
    if bimfile!=None:
        print 'Parsing SNP list'
        valid_sids = []
        print 'Parsing .bim file: %s'%bimfile
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                valid_sids.append(l[1])
        print 'Found %d SNPs in .bim file'%len(valid_sids)
        valid_sids = sp.array(valid_sids)
    chrom_dict = {}

    print 'Parsing SNP rsIDs from summary statistics.'
    sids = []
    with open(filename) as f:
        
        line = f.next()
        header = line.split()        
        if header==['hg19chrc', 'snpid', 'a1', 'a2', 'bp', 'info', 'or', 'se', 'p', 'ngt'] or header==headers['TAG'] or header==headers['CD'] or header==headers['UC'] or header==headers['ASTHMA']:
            for line in f:
                l = line.split()
                sids.append(l[1])
        elif header==['Chromosome', 'Position', 'MarkerName', 'Effect_allele', 'Non_Effect_allele', 'Beta', 'SE', 'Pvalue'] or header==headers['GCAN'] or header==headers['GEFOS'] or header==headers['ICBP']:
#             line_count =0
            for line in f:
#                 line_count +=1
#                 if line_count<100000:
                l = line.split()
                sids.append(l[2])
#                 elif random.random()<0.1:
#                     l = line.split()
#                     sids.append(l[2])
                        
        else:
            for line in f:
                l = line.split()
                sids.append(l[0])
                
    if bimfile!=None:
        valid_sids_filter = sp.in1d(valid_sids, sp.array(sids)) 
        sids = valid_sids[valid_sids_filter]

    print 'Retrieving 1K genomes positions.'
    sid_map = get_sid_pos_map(sids,KGpath)
    assert len(sid_map)>0, 'WTF?'

    print 'Parsing the file: %s' % filename
    with open(filename) as f:
        line = f.next()
        header = line.split()
        line_i = 0
        if header== ['snpid', 'hg18chr', 'bp', 'a1', 'a2', 'or', 'se', 'pval', 'info', 'ngt', 'CEUaf']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[7])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = -sp.log(float(l[5]))
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i
        
        elif header==['hg19chrc', 'snpid', 'a1', 'a2', 'bp', 'info', 'or', 'se', 'p', 'ngt']:    
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[1]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[8])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = -sp.log(float(l[6]))
                    if random.random()>0.5:
                        nt = [l[2], l[3]]
                    else:
                        nt = [l[3], l[2]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i
        
        elif header== ['snpid', 'hg18chr', 'bp', 'a1', 'a2', 'zscore', 'pval', 'CEUmaf']:     
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[6])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[5])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i
        
        elif header==['SNP', 'CHR', 'BP', 'A1', 'A2', 'OR', 'SE', 'P', 'INFO', 'EUR_FRQ']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[7])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = -sp.log(float(l[5]))
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i
        
        elif header==['Chromosome', 'Position', 'MarkerName', 'Effect_allele', 'Non_Effect_allele', 'Beta', 'SE', 'Pvalue']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[2]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[7])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[5])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i
        
        elif header ==headers['SSGAC1']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[6])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[4])
                    if random.random()>0.5:
                        nt = [l[1], l[2]]
                    else:
                        nt = [l[2], l[1]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i             
        elif header ==headers['SSGAC2']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[6])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = sp.log(float(l[4]))
                    if random.random()>0.5:
                        nt = [l[1], l[2]]
                    else:
                        nt = [l[2], l[1]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i             
        elif header ==headers['CHIC']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[8])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[8])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i             
        
        elif header==headers['GCAN']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[2]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[11])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = -sp.log(float(l[6]))
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i             
        
        elif header==headers['TESLOVICH']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[5])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[4])
                    if random.random()>0.5:
                        nt = [lc_2_cap_map[l[1]], lc_2_cap_map[l[2]]]
                    else:
                        nt = [lc_2_cap_map[l[2]], lc_2_cap_map[l[1]]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    chrom_dict[chrom]['weights'].append(int(float(l[3])))
                if line_i%100000==0:
                    print line_i                            
        
        elif header==headers['GIANT1'] or header==headers['GIANT1b'] or header==headers['GIANT2']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[6])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[4])
                    if random.random()>0.5:
                        nt = [l[1], l[2]]
                    else:
                        nt = [l[2], l[1]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    chrom_dict[chrom]['weights'].append(int(float(l[7])))
                if line_i%100000==0:
                    print line_i                                     
        elif header==headers['GIANT1c']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[8])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[6])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    chrom_dict[chrom]['weights'].append(int(float(l[9])))
                if line_i%100000==0:
                    print line_i   
        elif header==headers['MAGIC']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[6])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[4])
                    if random.random()>0.5:
                        nt = [lc_2_cap_map[l[1]], lc_2_cap_map[l[2]]]
                    else:
                        nt = [lc_2_cap_map[l[2]], lc_2_cap_map[l[1]]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['CARDIoGRAM']: 

            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[5])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[7])
                    if random.random()>0.5:
                        nt = [l[2], l[3]]
                    else:
                        nt = [l[3], l[2]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = float(l[9]) +float(l[10])
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['DIAGRAM']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[5])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = sp.log(float(l[6]))
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = float(l[9]) +float(l[10])
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['TAG']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[1]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[10])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[8])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['CD'] or header==headers['UC']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[1]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[10])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = sp.log(float(l[8]))
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['GEFOS']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[2]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[11])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[6])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)
#                     weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
#                     weight = (z/beta)**2
                    weight = float(l[16])  # Number of studies used.
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
#           'RA':['SNPID','Chr','Position(hg19)','A1','A2','OR(A1)','OR_95%CIlow','OR_95%CIup','P-val'],
#           'ASTHMA':['Chr', 'rs', 'position', 'Allele_1', 'Allele_2', 'freq_all_1_min', 'freq_all_1_max', 'OR_fix', 'ORl_fix', 'ORu_fix', 'P_fix'],
        elif header==headers['RA']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    raw_beta = sp.log(float(l[5]))
                    if raw_beta==0:
                        continue
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[8])
                    chrom_dict[chrom]['ps'].append(pval)
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['ASTHMA']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[1]
                d = sid_map.get(sid,None)
                if d is not None:
                    raw_beta = sp.log(float(l[7]))
                    pval = float(l[10])
                    if raw_beta==0 or pval == 0:
                        continue
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    chrom_dict[chrom]['ps'].append(pval)
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['ICBP']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[2]
                d = sid_map.get(sid,None)
                coded_allele = l[16]
                if d is not None and coded_allele in valid_nts:
#                     raw_beta = sp.log(float(l[7]))
                    pval = float(l[3])
                    if pval == 0:
                        continue
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    chrom_dict[chrom]['ps'].append(pval)
#                     if random.random()>0.5:
                    nt = [l[11], l[12]]
                    sign = 1
#                     else:
#                         sign = -1
                    if coded_allele==nt[1] or opp_strand_dict[coded_allele]==nt[1]:
                        nt = [l[12], l[11]]
                        sign = -1#*sign
#                     else:
#                         assert coded_allele==nt[0] or opp_strand_dict[coded_allele]==nt[0]
                    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sign * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = float(l[17])
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   

        else:
            raise Exception('Wrong or unknown file format')
    
        assert sp.all(sp.isreal(chrom_dict[1]['zs'])), 'WTF?'

    print 'SS file loaded, now sorting and storing in HDF5 file.'
    assert not ss_id in h5f.keys(), 'Summary stats with this name are already in the HDF5 file?'
    ssg = h5f.create_group(ss_id)
    num_snps = 0
    for chrom in chrom_dict.keys():
        print 'Parsed summary stats for %d SNPs on chromosome %d'%(len(chrom_dict[chrom]['positions']),chrom)
        sl = zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['ps'], chrom_dict[chrom]['zs'], chrom_dict[chrom]['eur_maf'], 
                 chrom_dict[chrom]['weights'])
        sl.sort()
        ps = []
        nts = []
        sids = []
        positions = []
        zs = []
        eur_mafs = []
        weights = []
        prev_pos = -1
        for pos, sid, nt, p, z, eur_maf, weight in sl:
            if pos == prev_pos:
                print 'duplicated position %d' % pos
                continue
            else:
                prev_pos = pos
            ps.append(p)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            zs.append(z)
            eur_mafs.append(eur_maf)
            weights.append(weight)
        g = ssg.create_group('chrom_%d' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('eur_mafs', data=eur_mafs)
        g.create_dataset('positions', data=positions)
        g.create_dataset('zs', data=zs)
        g.create_dataset('weights', data=weights)
        num_snps +=len(sids)
        h5f.flush()
    
    print 'In all, %d SNPs parsed from summary statistics file.'%num_snps



#----------------------------------------- Code for coordinating the datasets summary statistics (of various formats) ---------------------------------------------

def _get_chrom_dict_(loci, chromosomes):
    chr_dict = {}
    for chrom in chromosomes:
        chr_str = 'chrom_%d'%chrom
        chr_dict[chr_str] = {'sids':[],'snp_indices':[],'positions':[], 'nts':[]}
     
    for i, l in enumerate(loci):
        chrom = l.chromosome
        pos = l.bp_position
        chr_str = 'chrom_%d'%chrom
        chr_dict[chr_str]['sids'].append(l.name)
#         chr_dict[chr_str]['sids'].append('%d_%d'%(chrom,pos))
        chr_dict[chr_str]['snp_indices'].append(i)
        chr_dict[chr_str]['positions'].append(pos)
        chr_dict[chr_str]['nts'].append([l.allele1,l.allele2])
     
    print 'Genotype dictionary filled'
    return chr_dict


def coordinate_prediction_datasets(val_genotype_file = None,
                                    ref_genotype_file = None,
                                    sum_stat_file = None,
                                    pred_genotype_map=None,
                                    min_maf=0.01):
    """
    Coordinates the genotypes across 4 data sets: 
        a) the summary statistics 
        b) the LD-reference panel 
        c) the validation data set 
        d) the prediction genotype
    """
    
#   recode_dict = {1:'A', 2:'T', 3:'C', 4:'G'} #1K genomes recoding..
    print 'Coordinating datafiles w validation genotype file: %s \nref. genot. file: %s\nsum. stat. file: %s\npred. genotype map file:%s'%(val_genotype_file, ref_genotype_file, sum_stat_file, pred_genotype_map) 
    plinkf = plinkfile.PlinkFile(val_genotype_file)
    
    #Loads only the individuals... (I think?)
    samples = plinkf.get_samples()
    num_individs = len(samples)
    Y = [s.phenotype for s in samples]
    fids = [s.fid for s in samples]
    iids = [s.iid for s in samples]
    
    unique_phens = sp.unique(Y)
    if len(unique_phens)==1:
        print 'Unable to find phenotype values.'
        has_phenotype=False
    elif len(unique_phens)==2:
        cc_bins = sp.bincount(Y)
        assert len(cc_bins)==2, 'Problems with loading phenotype'
        print 'Loaded %d controls and %d cases'%(cc_bins[0], cc_bins[1])
        has_phenotype=True
    else:
        print 'Found quantitative phenotype values'
        has_phenotype=True

    #Figure out chromosomes and positions.  
    print 'Parsing validation genotype bim file'
    loci = plinkf.get_loci()
    plinkf.close()
    gf_chromosomes = [l.chromosome for l in loci] 

    chromosomes = sp.unique(gf_chromosomes)
    chromosomes.sort()
    
    chr_dict = _get_chrom_dict_(loci, chromosomes)

    print 'Parsing LD reference genotype bim file'
    plinkf_ref = plinkfile.PlinkFile(ref_genotype_file)
    loci_ref = plinkf_ref.get_loci()
    plinkf_ref.close()
    
    chr_dict_ref = _get_chrom_dict_(loci_ref, chromosomes)
#     chr_dict_ref = _get_chrom_dict_bim_(reference_genotype_file+'.bim', chromosomes)
    
    
    print 'Parsing the summary statistics'
    
    #Open HDF5 file and prepare out data
    assert not 'iids' in hdf5_file.keys(), 'Something is wrong with the HDF5 file?'
    if has_phenotype:
        hdf5_file.create_dataset('y', data=Y)
    
    hdf5_file.create_dataset('fids', data=fids)
    hdf5_file.create_dataset('iids', data=iids)
    ssf = hdf5_file['sum_stats']
    cord_data_g = hdf5_file.create_group('cord_data')

    maf_adj_risk_scores = sp.zeros(num_individs)
    num_common_snps = 0
    #corr_list = []
    
    tot_g_ss_nt_concord_count = 0
    tot_rg_ss_nt_concord_count = 0
    tot_g_rg_nt_concord_count = 0    
    tot_num_non_matching_nts = 0
   
    #Now iterate over chromosomes
    for chrom in chromosomes:
        ok_indices = {'g':[], 'rg':[], 'ss':[]}
        
        chr_str = 'chrom_%d'%chrom
        print 'Working on chromsome: %s'%chr_str
        
        chrom_d = chr_dict[chr_str]
        chrom_d_ref = chr_dict_ref[chr_str]
        try:
            ssg = ssf['chrom_%d' % chrom]
        except Exception, err_str:
            print err_str
            print 'Did not find chromsome in SS dataset.'
            print 'Continuing.'
            continue

        ssg = ssf['chrom_%d' % chrom]
        g_sids = chrom_d['sids']
        rg_sids = chrom_d_ref['sids']
        ss_sids = ssg['sids'][...]
        print 'Found %d SNPs in validation data, %d SNPs in LD reference data, and %d SNPs in summary statistics.'%(len(g_sids), len(rg_sids), len(ss_sids))
        common_sids = sp.intersect1d(ss_sids, g_sids)
        common_sids = sp.intersect1d(common_sids, rg_sids)
        print 'Found %d SNPs on chrom %d that were common across all datasets'%(len(common_sids), chrom)

        ss_snp_map = []
        g_snp_map = []
        rg_snp_map = []
        
        ss_sid_dict = {}
        for i, sid in enumerate(ss_sids):
            ss_sid_dict[sid]=i

        g_sid_dict = {}
        for i, sid in enumerate(g_sids):
            g_sid_dict[sid]=i

        rg_sid_dict = {}
        for i, sid in enumerate(rg_sids):
            rg_sid_dict[sid]=i
            
        for sid in common_sids:
            g_snp_map.append(g_sid_dict[sid])
        
        #order by positions
        g_positions = sp.array(chrom_d['positions'])[g_snp_map]
        order = sp.argsort(g_positions)
        #order = order.tolist()
        g_snp_map = sp.array(g_snp_map)[order]
        g_snp_map = g_snp_map.tolist()
        common_sids = sp.array(common_sids)[order]

        #Get the other two maps
        for sid in common_sids:
            rg_snp_map.append(rg_sid_dict[sid])
        
        for sid in common_sids:
            ss_snp_map.append(ss_sid_dict[sid])
            
        
        g_nts = sp.array(chrom_d['nts'])
        rg_nts = sp.array(chrom_d_ref['nts'])
        rg_nts_ok = sp.array(rg_nts)[rg_snp_map]
#         rg_nts_l = []
#         for nt in rg_nts_ok:
#             rg_nts_l.append([recode_dict[nt[0]],recode_dict[nt[1]]])
#         rg_nts_ok = sp.array(rg_nts_l)
        ss_nts = ssg['nts'][...]
        betas = ssg['betas'][...]
        log_odds = ssg['log_odds'][...]

        if 'freqs' in ssg.keys():
            ss_freqs = ssg['freqs'][...]

        g_ss_nt_concord_count = sp.sum(g_nts[g_snp_map] == ss_nts[ss_snp_map])/2.0
        rg_ss_nt_concord_count = sp.sum(rg_nts_ok == ss_nts[ss_snp_map])/2.0
        g_rg_nt_concord_count = sp.sum(g_nts[g_snp_map] == rg_nts_ok)/2.0
        print 'Nucleotide concordance counts out of %d genotypes: vg-g: %d, vg-ss: %d, g-ss: %d'%(len(g_snp_map),g_rg_nt_concord_count, g_ss_nt_concord_count, rg_ss_nt_concord_count)
        tot_g_ss_nt_concord_count += g_ss_nt_concord_count
        tot_rg_ss_nt_concord_count += rg_ss_nt_concord_count
        tot_g_rg_nt_concord_count += g_rg_nt_concord_count


        num_non_matching_nts = 0
        num_ambig_nts = 0


        #Identifying which SNPs have nucleotides that are ok..
        ok_nts = []
        for g_i, rg_i, ss_i in it.izip(g_snp_map, rg_snp_map, ss_snp_map):
            
            #To make sure, is the SNP id the same?
            assert g_sids[g_i]==rg_sids[rg_i]==ss_sids[ss_i], 'Some issues with coordinating the genotypes.'
            
            g_nt = g_nts[g_i]
            rg_nt = rg_nts[rg_i]
#             rg_nt = [recode_dict[rg_nts[rg_i][0]],recode_dict[rg_nts[rg_i][1]]]
            ss_nt = ss_nts[ss_i]

            #Is the nucleotide ambiguous.
            g_nt = [g_nts[g_i][0],g_nts[g_i][1]]
            if tuple(g_nt) in ambig_nts:
                num_ambig_nts +=1
                tot_num_non_matching_nts += 1                
                continue
            
            #First check if nucleotide is sane?
            if (not g_nt[0] in valid_nts) or (not g_nt[1] in valid_nts):
                num_non_matching_nts += 1
                tot_num_non_matching_nts += 1                
                continue
            
            os_g_nt = sp.array([opp_strand_dict[g_nt[0]], opp_strand_dict[g_nt[1]]])

            flip_nts = False
            if not ((sp.all(g_nt == ss_nt) or sp.all(os_g_nt == ss_nt)) and (sp.all(g_nt == rg_nt) or sp.all(os_g_nt == rg_nt))):
                if sp.all(g_nt == rg_nt) or sp.all(os_g_nt == rg_nt):
                    flip_nts = (g_nt[1] == ss_nt[0] and g_nt[0] == ss_nt[1]) or (os_g_nt[1] == ss_nt[0] and os_g_nt[0] == ss_nt[1])
                    #Try flipping the SS nt
                    if flip_nts:
                        betas[ss_i] = -betas[ss_i]                        
                        log_odds[ss_i] = -log_odds[ss_i]    
                        if 'freqs' in ssg.keys():
                            ss_freqs[ss_i] = 1-ss_freqs[ss_i]
                    else:
                        print "Nucleotides don't match after all?: g_sid=%s, ss_sid=%s, g_i=%d, ss_i=%d, g_nt=%s, ss_nt=%s" % \
                            (g_sids[g_i], ss_sids[ss_i], g_i, ss_i, str(g_nt), str(ss_nt))
                        num_non_matching_nts += 1
                        tot_num_non_matching_nts += 1
                        continue

                    
                else:
                    num_non_matching_nts += 1
                    tot_num_non_matching_nts += 1
                    continue
                    # Opposite strand nucleotides
            
           
            # everything seems ok.
            ok_indices['g'].append(g_i)
            ok_indices['rg'].append(rg_i)
            ok_indices['ss'].append(ss_i)

            ok_nts.append(g_nt)
#             if flip_nts:
#                 ok_nts.append([ss_nt[1],ss_nt[0]])
#             else:
#                 ok_nts.append(ss_nt)                

                        
        #print '%d SNPs in LD references to be flipped.'%((len(ref_snp_directions)-sp.sum(ref_snp_directions))/2.0)
        print '%d SNPs had ambiguous nucleotides.' % num_ambig_nts 
        print '%d SNPs were excluded due to nucleotide issues.' % num_non_matching_nts 
        print '%d SNPs were retained on chromosome %d.' % (len(ok_indices['g']), chrom)

        #Resorting by position
        positions = sp.array(chrom_d['positions'])[ok_indices['g']]
#         order = sp.argsort(positions)
#         sorted_positions = positions[order]
#         assert sp.all(sorted_positions==positions), 'Perhaps something is wrong here?'
#         ok_indices['g'] = list(sp.array(ok_indices['g'])[order])
#         ok_indices['ss'] = list(sp.array(ok_indices['ss'])[order])

        
        #Now parse SNPs ..
        snp_indices = sp.array(chrom_d['snp_indices'])
        snp_indices = snp_indices[ok_indices['g']] #Pinpoint where the SNPs are in the file.
        raw_snps,freqs = _parse_plink_snps_(genotype_file, snp_indices)
        
        snp_indices_ref = sp.array(chrom_d_ref['snp_indices'])
        snp_indices_ref = snp_indices_ref[ok_indices['rg']] #Pinpoint where the SNPs are in the file.
        raw_ref_snps, freqs_ref = _parse_plink_snps_(reference_genotype_file, snp_indices_ref)
        
        
        snp_stds_ref = sp.sqrt(2*freqs_ref*(1-freqs_ref)) 
        snp_means_ref = freqs_ref*2

        snp_stds = sp.sqrt(2*freqs*(1-freqs)) 
        snp_means = freqs*2
        
        betas = betas[ok_indices['ss']]  # * sp.sqrt(freqs * (1 - freqs))
        log_odds = log_odds[ok_indices['ss']]  # * sp.sqrt(freqs * (1 - freqs))

        ps = ssg['ps'][...][ok_indices['ss']]
        nts = sp.array(ok_nts)#[order]
        sids = ssg['sids'][...][ok_indices['ss']]

        #For debugging...
#         g_sids = sp.array(chrom_d['sids'])[ok_indices['g']]
#         rg_sids = sp.array(chrom_d_ref['sids'])[ok_indices['rg']]
#         ss_sids = ssg['sids'][...][ok_indices['ss']]
#         assert sp.all(g_sids==rg_sids) and sp.all(rg_sids==ss_sids), 'WTF?'
        
        #Check SNP frequencies..
        if check_mafs and 'freqs' in ssg.keys():
            ss_freqs = ss_freqs[ok_indices['ss']]
            freq_discrepancy_snp = sp.absolute(ss_freqs-(1-freqs))>0.15
            if sp.any(freq_discrepancy_snp):
                print 'Warning: %d SNPs were filtered due to high allele frequency discrepancy between summary statistics and validation sample'%sp.sum(freq_discrepancy_snp)
#                 print freqs[freq_discrepancy_snp]
#                 print ss_freqs[freq_discrepancy_snp]
                 
                #Filter freq_discrepancy_snps
                ok_freq_snps = sp.negative(freq_discrepancy_snp)
                raw_snps = raw_snps[ok_freq_snps]
                snp_stds = snp_stds[ok_freq_snps]
                snp_means = snp_means[ok_freq_snps]
                raw_ref_snps = raw_ref_snps[ok_freq_snps]
                snp_stds_ref = snp_stds_ref[ok_freq_snps]
                snp_means_ref = snp_means_ref[ok_freq_snps]
                freqs = freqs[ok_freq_snps]
                freqs_ref = freqs_ref[ok_freq_snps]
                ps = ps[ok_freq_snps]
                positions = positions[ok_freq_snps]
                nts = nts[ok_freq_snps]
                sids = sids[ok_freq_snps]
                betas = betas[ok_freq_snps]
                log_odds = log_odds[ok_freq_snps]
                #For debugging...
#         if sp.any(freq_discrepancy_snp):
#             g_sids = g_sids[ok_freq_snps]
#             rg_sids = rg_sids[ok_freq_snps]
#             ss_sids = ss_sids[ok_freq_snps]
#         assert sp.all(g_sids==rg_sids) and sp.all(rg_sids==ss_sids), 'WTF?'

        
        
        #Filter minor allele frequency SNPs.
        maf_filter = (freqs>min_maf)*(freqs<(1-min_maf))
        maf_filter_sum = sp.sum(maf_filter)
        n_snps = len(maf_filter)
        assert maf_filter_sum<=n_snps, "WTF?"
        if sp.sum(maf_filter)<n_snps:
            raw_snps = raw_snps[maf_filter]
            snp_stds = snp_stds[maf_filter]
            snp_means = snp_means[maf_filter]
            raw_ref_snps = raw_ref_snps[maf_filter]
            snp_stds_ref = snp_stds_ref[maf_filter]
            snp_means_ref = snp_means_ref[maf_filter]
            freqs = freqs[maf_filter]
            freqs_ref = freqs_ref[maf_filter]
            ps = ps[maf_filter]
            positions = positions[maf_filter]
            nts = nts[maf_filter]
            sids = sids[maf_filter]
            betas = betas[maf_filter]
            log_odds = log_odds[maf_filter]
#         if sp.sum(maf_filter)<n_snps:
#             g_sids = g_sids[maf_filter]
#             rg_sids = rg_sids[maf_filter]
#             ss_sids = ss_sids[maf_filter]
#         assert sp.all(g_sids==rg_sids) and sp.all(rg_sids==ss_sids), 'WTF?'
        
        
        
        maf_adj_prs = sp.dot(log_odds, raw_snps)
        if has_phenotype:
            maf_adj_corr = sp.corrcoef(Y, maf_adj_prs)[0, 1]
            print 'Log odds, per genotype PRS correlation w phenotypes for chromosome %d was %0.4f' % (chrom, maf_adj_corr)

        genetic_map = [] 
        if genetic_map_dir is not None:
            with gzip.open(genetic_map_dir+'chr%d.interpolated_genetic_map.gz'%chrom) as f:
                for line in f:
                    l = line.split()
                    if l[0] in sid_set:
                        genetic_map.append(l[0])
        
        
        print 'Now storing coordinated data to HDF5 file.'
        ofg = cord_data_g.create_group('chrom_%d' % chrom)
        ofg.create_dataset('raw_snps_val', data=raw_snps, compression='lzf')
        ofg.create_dataset('snp_stds_val', data=snp_stds)
        ofg.create_dataset('snp_means_val', data=snp_means)
        ofg.create_dataset('freqs_val', data=freqs)
        ofg.create_dataset('raw_snps_ref', data=raw_ref_snps, compression='lzf')
        ofg.create_dataset('snp_stds_ref', data=snp_stds_ref)
        ofg.create_dataset('snp_means_ref', data=snp_means_ref)
        ofg.create_dataset('freqs_ref', data=freqs_ref)
        ofg.create_dataset('nts', data=nts)
        ofg.create_dataset('ps', data=ps)
        ofg.create_dataset('positions', data=positions)
        ofg.create_dataset('sids', data=sids)
        if genetic_map_dir is not None:
            ofg.create_dataset('genetic_map', data=genetic_map)
        ofg.create_dataset('betas', data=betas)
        ofg.create_dataset('log_odds', data=log_odds)
        ofg.create_dataset('log_odds_prs', data=maf_adj_prs)
#         print 'Sum betas', sp.sum(betas ** 2)
        #ofg.create_dataset('prs', data=prs)
        
        
        #risk_scores += prs
        maf_adj_risk_scores += maf_adj_prs
        num_common_snps += len(betas)
        
    # Now calculate the prediction r^2
    if has_phenotype:
        maf_adj_corr = sp.corrcoef(Y, maf_adj_risk_scores)[0, 1]
        #print 'PRS correlation for the whole genome was %0.4f (r^2=%0.4f)' % (corr, corr ** 2)
        print 'Log odds, per PRS correlation for the whole genome was %0.4f (r^2=%0.4f)' % (maf_adj_corr, maf_adj_corr ** 2)
    print 'Overall nucleotide concordance counts: g_rg: %d, g_ss: %d, rg_ss: %d'%(tot_g_rg_nt_concord_count, tot_g_ss_nt_concord_count, tot_rg_ss_nt_concord_count)
    print 'There were %d SNPs in common' % num_common_snps
    print 'In all, %d SNPs were excluded due to nucleotide issues.' % tot_num_non_matching_nts 
    print 'Done!'




    