#!/usr/bin/env python

from boruta import BorutaPy
import timeit
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
import pandas as pd
import os
import pickle
from pathlib import Path
import argparse
from joblib import Parallel, delayed, parallel_backend

class borutaObj:
    def __init__(self, form='csv', percent=100, alpha=0.05, v=0, max_iter=2500, repeat=1):
        self.form = form
        self.percent = percent
        self.alpha = alpha
        self.v = v
        self.max_iter = max_iter
        self.repeat = repeat
        
    def XY_index(self, X, Y, form):
        if form=='csv':
            X=pd.read_csv(X, sep='\t', index_col=0)
            Y=pd.read_csv(Y, sep='\t', index_col=0)
            Y.columns = ['RESPONSE']
            X = Y.join(X)
            Y = X['RESPONSE']
            X = X.drop(columns=['RESPONSE'])
        elif form =='pickle':
            #with open(X, 'rb') as filehandle:
            #    X=pickle.load(filehandle)
            X=pd.read_pickle(X) #, sep='\t', index_col='STRAIN'
            Y=pd.read_csv(Y, sep=' ', header=None, index_col=0).drop(columns=[1,3])        
            #Y= pd.read_pickle(Y) #, sep='\t', index_col='STRAIN'
            Y.columns=['RESPONSE']
            X = Y.join(X)
            Y = X['RESPONSE']
            X = X.drop(columns=['RESPONSE'])
        if X.isnull().all().all():
             print("Your indices do not match between X and Y!")
        # if X is all nan, something is wrong
        return X, Y
    def BorutaSelect(self, mod, percent, alpha, X, Y, Borut_Est, v, max_iter, rand): 
        y=Y.values # NEED .ravel() IF Y IS DATAFRAME INSTEAD OF SERIES
        x=X.values
        feat_selector = BorutaPy(mod, n_estimators=Borut_Est, alpha=alpha, verbose=v, random_state=rand, perc=percent, max_iter=max_iter)
        feat_selector.fit(x, y)
        tmp=pd.DataFrame(feat_selector.support_)
        best = tmp[tmp[0]==True].reset_index()
        FEATS=best['index'].to_list()
        X_filtered = X.iloc[:,FEATS]
        Important =X_filtered.columns.to_list()
        
        tmp2=pd.DataFrame(feat_selector.support_weak_)
        notbest = tmp2[tmp2[0]==True].reset_index()
        maybeFEATS=notbest['index'].to_list()
        X_maybe = X.iloc[:,maybeFEATS]
        MAYBE =X_maybe.columns.to_list()
        FEATURES = pd.DataFrame({'SELECTED': pd.Series(Important, dtype=object), 'TENTATIVE': pd.Series(MAYBE, dtype=object)})
           
        return FEATURES
    
    def BorutaOrganize(self, FEATURES, repeat):
        GENE_SELECT = pd.concat(FEATURES)
        GENE_SELECT_DF=pd.DataFrame(GENE_SELECT.SELECTED.value_counts())
        GENE_SELECT_DF['PERCENT_SELECT'] = GENE_SELECT_DF.SELECTED/repeat
        
        GENE_TENTATIVE_DF=pd.DataFrame(GENE_SELECT.TENTATIVE.value_counts())
        GENE_TENTATIVE_DF['PERCENT_SELECT'] = GENE_TENTATIVE_DF.TENTATIVE/repeat
        return GENE_SELECT_DF, GENE_TENTATIVE_DF
