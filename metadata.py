#class Dataset(prd.Dataset):
from __future__ import print_function, division
import gzip
import pandas as pd
from pyutils import fs, memo

from phenotype import Architecture
from gprim.dataset import Dataset as gpds

class Dataset(gpds):
    def __init__(self, path, name):
        gpds.__init__(self, path)
        self.name = name

class Simulation(object):
    def __init__(self, dataset, architecture, h2g, chromosomes, name, low_space = True):
        self.dataset = dataset
        self.architecture = architecture
        self.h2g = h2g
        self.chromosomes = chromosomes
        self.name = name
        self.low_space = low_space

#    @property
#    @memo.memoized
#    def architecture(self):
#        self.__architecture = Architecture(self.annot_files, self.vcov_effects)
#        #self.__architecture.set_pheno_var(self.h2g, self.chromosomes)
#        return self.__architecture
    #TODO: import architecture, check against this

    def root_folder(self, create=True):
        #NB: used to use a paths string here, see yakir's code for reference
        folder = '{}/{}/'.format(self.dataset.name, self.name)
        if create:
            fs.makedir(folder)
        return folder
    def beta_folder(self, beta_num, create=True):
        folder = '{}{}/'.format(self.root_folder(), beta_num)
        if create:
            fs.makedir(folder)
        return folder
    def beta_filename(self, beta_num, chrnum):
        return '{}{}.beta.gz'.format(
            self.beta_folder(beta_num), chrnum)
    def beta_file(self, beta_num, chrnum, mode='r'):
        return gzip.open(self.beta_filename(beta_num, chrnum), mode=mode)
    def sparsebeta_filename(self, beta_num, chrnum, mode='r'):
        return '{}{}.betanz'.format(
            self.beta_folder(beta_num), chrnum)
    def sparsebeta_file(self, beta_num, chrnum, mode='r'):
        return open(self.sparsebeta_filename(beta_num, chrnum), mode=mode)

    def chr_filestem(self, beta_num, chrnum):
        return '{}{}'.format(self.beta_folder(beta_num), chrnum)
    def noiselessYchr_filename(self, beta_num, chrnum):
        return self.chr_filestem(beta_num, chrnum) + '.profile'
    def noiselessY_filename(self, beta_num):
        return self.beta_folder(beta_num) + 'noiseless.pheno'
    def noisyY_filename(self, beta_num):
        return self.beta_folder(beta_num) + 'noisy.pheno'
    def sumstats_filename(self, beta_num):
        return self.beta_folder(beta_num) + 'all.sumstats.gz'
    def sumstats_file(self, beta_num, mode='r'):
        return gzip.open(self.sumstats_filename(beta_num), mode=mode)

    
    
