from __future__ import print_function, division
import argparse, subprocess, os
import pandas as pd
import scipy.stats as stats
import numpy as np
import sys

import metadata as sm
import phenotype as ph

def create_beta_and_profiles(s):
    for chrnum in s.chromosomes:
        betas = s.architecture.draw_2beta(chrnum)
        maf = s.dataset.frq_df(chrnum)['MAF'].values
        for beta_num in range(len(betas)): 
# TODO: instead of a list of dataframes, have a single DF that we select on to
# save out betas. More memory efficient that way.
            beta = betas[beta_num]
            beta = pd.concat([s.dataset.bim_df(chrnum)[['SNP','A1','A2']],
                beta[['BETA']]], axis=1)
            beta.to_csv(s.beta_filename(beta_num, chrnum), index = False, sep = '\t')

            if beta.shape[0] != len(maf):
                raise Exception("Length of betas and maf must match")

            v = np.var(beta['BETA'])
            beta['BETA'] /= (np.sqrt(2 * maf * (1-maf)) * \
                    np.sqrt(2749160 * v / s.h2g))
                # convert beta to plink's per-allele scale

            # write sparse beta as well (plink works faster if beta is explicitly sparse)
            sparsebeta = beta.loc[beta['BETA'] != 0]
            sparsebeta.to_csv(s.sparsebeta_filename(beta_num, chrnum), index = False, sep = '\t')
            cmd = ['plink',
                    '--bfile', s.dataset.bfile(chrnum),
                    '--allow-no-sex',
                    '--score', s.sparsebeta_filename(beta_num, chrnum), '1', '2', '4',
                    'sum',
                    'center',
                    '--out', s.chr_filestem(beta_num, chrnum)]
            #print 'executing ' + ' '.join(cmd)
            os.system(' '.join(cmd))

def make_noiseless_pheno(s, beta_num):
    def get_profile(fname):
        df = pd.read_csv(fname, header=0, delim_whitespace=True)
        df.drop(['PHENO', 'CNT', 'CNT2'], axis=1, inplace=True)
        return df

    print('merging phenotypes. chr', s.chromosomes[0])
    phenotype = get_profile(s.noiselessYchr_filename(beta_num, s.chromosomes[0]))
    phenotype.rename(columns={'SCORESUM' : 'PHENO'}, inplace=True)
    for chrnum in s.chromosomes[1:]:
        print('merging phenotypes. chr', chrnum)
        profile = get_profile(s.noiselessYchr_filename(beta_num, chrnum))
        phenotype = pd.merge(phenotype, profile, how='inner', on=['FID', 'IID'])
        phenotype['PHENO'] += phenotype['SCORESUM']
        phenotype.drop(['SCORESUM'], axis=1, inplace=True)
    print('variance of noiseless phenotype is:',
            np.var(phenotype.PHENO), 'It should be', s.h2g)

    phenotype.to_csv(s.noiselessY_filename(beta_num),
            index=False,
            sep='\t')
    return phenotype


# add noise to resulting phenotype
def add_noise_and_save(s, beta_num, phenotype):
    sigma2e = 1-s.h2g
    print('adding noise. sigma2e =', sigma2e)
    phenotype['PHENO'] += np.sqrt(sigma2e) * np.random.randn(len(phenotype))
    phenotype.to_csv(s.noisyY_filename(beta_num),
            index=False,
    sep='\t')

# call plink to compute sumstats
def make_qassoc(s, beta_num):
    # compute one set of sumstats per chromosome
    for chrnum in s.chromosomes:
        print('computing sumstats for chr', chrnum)
        cmd = ['plink',
                '--bfile', s.dataset.bfile(chrnum),
                '--pheno', s.noisyY_filename(beta_num),
                '--allow-no-sex',
                '--assoc',
                '--out', s.chr_filestem(beta_num, chrnum)]
        print('executing', ' '.join(cmd))
        os.system(' '.join(cmd))

def make_sumstats(s, beta_num):
    # create one large df from the qassoc files.
    print('merging chromosomes of sumstats')
    sumstats = []
    for chrnum in s.chromosomes:
        print(chrnum)
        qassoc = pd.read_csv(s.chr_filestem(beta_num, chrnum)+'.qassoc',
                             header=0, delim_whitespace=True)
        qassoc['A1'] = s.dataset.bim_df(chrnum)['A1']
        qassoc['A2'] = s.dataset.bim_df(chrnum)['A2']
        qassoc['Z'] = qassoc['BETA'] / qassoc['SE']
        sumstats.append(qassoc[['SNP', 'A1', 'A2', 'Z', 'NMISS']])

    sumstats = pd.concat(sumstats, axis=0)
    sumstats.rename(columns = {'NMISS': 'N'}, inplace = True)
    sumstats['A1'] = 'A'
    sumstats['A2'] = 'G'

    if np.sum(np.isnan(sumstats['Z'])) > 0 or np.sum(np.isinf(sumstats['Z'])) > 0:
        print('ERROR: some summary statistics were either inf or nan. aborting')
        print(np.where(np.isnan(sumstats['Z']))[0])
        print(np.where(np.isinf(sumstats['Z']))[0])
        return

#TODO: decide if we need print snps; ask yakir about the file
#    # filter down to the snps that we want
#    if s.print_snps != 'none':
#        print('filtering to', s.print_snps)
#        rsids = pd.read_csv(s.print_snps, header=None, names=['SNP'])
#        sumstats = pd.merge(sumstats, rsids, on='SNP', how='inner')
#    print(len(sumstats), 'snps left')


    # save the big sumstats file
    print('saving result to', s.sumstats_filename(beta_num))
    sumstats.to_csv(s.sumstats_file(beta_num, mode='w'), sep='\t', index=False)

def pdmatrix(tau1 = None, tau2 = None, ups = None):
    if tau1 is None:
        tau1 = np.square(np.random.normal(0, 1))
    if tau2 is None:
        tau2 = np.square(np.random.normal(0, 1))
    if ups is None:
        ups = np.random.uniform(-np.sqrt(tau1*tau2), np.sqrt(tau1*tau2))
    if np.abs(ups) > np.sqrt(tau1 * tau2):
        raise Exception("Matrix is not positive definite")
    else:
        return(np.matrix([[tau1, ups], [ups, tau2]]))

if __name__ == '__main__':
    arrayind = sys.argv[2] if len(sys.argv) > 2 else 0
    s = sm.Simulation.from_json(sys.argv[1], arrayind)
    create_beta_and_profiles(s)
    for i in range(2):
        pheno = make_noiseless_pheno(s, i)
        add_noise_and_save(s, i, pheno)
        make_qassoc(s, i)
        make_sumstats(s, i)

        # If we're running low on space or need to do a lot of simulations,
        # we'll have to delete our intermediate files
        if s.low_space:
            s.remove_intermediate_files(beta_num)


