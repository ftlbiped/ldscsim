from __future__ import print_function, division
import numpy as np
import pandas as pd
from gprim.annotation import Annotation

class Architecture(object):
    def __init__(self, annot_files, vcov_effects):
        #self.name = name
        self.vcov_effects = vcov_effects
        #TODO: type checking, ensure vcov is positive definite matrix
        
        #TODO: manually add in all the arguments here

        self.annot_files = annot_files
        self.annotations = {}
        for annot_file in self.annot_files:
            self.annotations[annot_file] = Annotation(annot_file)

        #TODO: ensure vcov dimensions match annotations

#TODO: this doesn't work because gprim.annotation assumes sannot instead of annot.
#      So talk with Yakir, sort this out, and make draw betas more robust
#    def __find_annot(self, n, chrnum, signed=True):
#        matches = [a for a in self.annotations.values() if n in a.names(chrnum)]
#        if signed:
#            return matches[0].sannot_df(chrnum)[n].values # there should only be one match
#        else:
#            return matches[0].annot_df(chrnum)[n].values # there should only be one match

#    def draw_beta():
        #for drawing a single beta

#    def draw_2beta(self, chrnum):
#        #draw a pair of betas
#        result = self.annotations.values()[0].annot_df(chrnum).copy()
#        result = result[['SNP', 'A1', 'A2']]
#        result['BETA'] = 0
#
#        for vcov in 
    
    #TODO: this is still sloppy and cutting corners. Fix this.
    def draw_2beta(self, chrnum):
        annot_df = self.annotations.values()[0].annot_df(chrnum)
        cat = annot_df.iloc[:, np.arange(4, annot_df.shape[1])].as_matrix()
        cat = [[i for i, x in enumerate(x) if x == 1] for x in cat]
        vcov = [np.sum([self.vcov_effects[i] for i in x], 0) for x in cat]
        betas = [np.random.multivariate_normal([0,0], vc) for vc in vcov] 

        #result = self.annotations.values()[0].annot_df(chrnum).copy()

        #TODO: return single dataframe instead of list of DFs
        def beta_df(i):
            beta = np.array([b[i] for b in betas])
            beta = pd.DataFrame(beta)
            beta.columns = ['BETA']
            #beta = pd.concat([result[['SNP']], beta], axis=1)
            return beta

        return [beta_df(i) for i in range(2)]








            


