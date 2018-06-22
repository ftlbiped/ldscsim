#class Dataset(prd.Dataset):
from __future__ import print_function, division
import gzip
import pandas as pd
import numpy as np
from pyutils import fs, memo
import json
import os

from phenotype import Architecture
from gprim.dataset import Dataset as gpds
import paths

class Dataset(gpds):
    def __init__(self, path, name):
        gpds.__init__(self, path)
        self.name = name

class Simulation(object):
    def __init__(self, dataset, architecture, h2g, chromosomes, name, low_space = True,
            arrayind = 0):
        self.dataset = dataset
        self.architecture = architecture
        self.h2g = h2g
        self.chromosomes = chromosomes
        self.name = name
        self.low_space = low_space
        self.arrayind = arrayind

    @classmethod
    def from_json(cls, name, arrayind = 0):
        config = '{}{}.json'.format(paths.simconfig, name)
        with open(config, 'r') as f:
            cfg = json.load(f)

        dataset = Dataset(path = cfg['dataset_path'], 
                name = cfg['dataset_name'])

        chrnums = range(1, 23) if cfg['chromosomes'] == 'all' else cfg['chromosomes']

        # Vcov effects are stored in config as the per-category contribution to
        # VCOV, but architecture understands on a per-SNP scale.
        vcov_effects = [np.matrix(x) for x in cfg['vcov'].values()]

        if type(cfg['annot']) is str:
            raise TypeError("Did you mean to feed a list of annotations?")
        
        new = cls(
                dataset = Dataset(path = cfg['dataset_path'], 
                    name = cfg['dataset_name']),
                architecture = Architecture(annot_files = cfg['annot'],
                    vcov_effects = vcov_effects),
                h2g = cfg['h2g'],
                chromosomes = chrnums,
                name = name,
                arrayind = arrayind)
        return new
                     
#    @property
#    @memo.memoized
#    def architecture(self):
#        self.__architecture = Architecture(self.annot_files, self.vcov_effects)
#        #self.__architecture.set_pheno_var(self.h2g, self.chromosomes)
#        return self.__architecture
    #TODO: import architecture, check against this

    def root_folder(self, create=True):
        # NB: used to use a paths string here, see yakir's code for reference
        # folder = '{}/{}/'.format(self.dataset.name, self.name)
        folder = '{}{}/{}/'.format(paths.simsumstats, self.name, self.arrayind)
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

    def remove_intermediate_files(self, beta_num):
        #import gzip, shutil
        print('compressing/deleting intermediate files')
        fnames = []
        for chrnum in self.chromosomes:
            # # make a gzipped copy of the betanz files before we delete them
            # fname = s.chr_filestem(beta_num, chrnum)+'.betanz'
            # with open(fname, 'rb') as f_in, gzip.open(fname+'.gz', 'wb') as f_out:
            #     shutil.copyfileobj(f_in, f_out)

            # Yakir's original filetypes
            filetypes = ['betanz','log','nosex','profile','qassoc','assoc.linear']
            
            # We're extra low on space, because we want to do a lot of
            # simulations, so we're going to delete everything except the
            # sumstats files.
            filetypes += ['beta.gz', 'nopred']

            for suff in filetypes:
                fname = self.chr_filestem(beta_num, chrnum)+'.'+suff
                fnames.append(fname)

        fnames.append(self.noiselessY_filename(beta_num))
        fnames.append(self.noisyY_filename(beta_num))

        for fname in fnames:
            print('removing', fname)
            if os.path.exists(fname): os.remove(fname)

    
    
