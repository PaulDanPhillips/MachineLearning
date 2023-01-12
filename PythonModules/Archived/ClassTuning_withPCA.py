# import timeit
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
# from xgboost import XGBClassifier
import pandas as pd
import numpy as np
import os
# import pickle
from pathlib import Path
from sklearn.model_selection import GridSearchCV
import argparse

class Tuning:
    def __init__(self, n_splits, repeat, n_jobs, Eparams={'C': np.logspace(-3,2,5), 'l1_ratio': np.linspace(0,1,4)}, LRparams={'C': np.logspace(0,1,4), 'penalty':['l1', 'l2']}, GBparams = {'learning_rate' : np.linspace(0.0001,0.01,5), 'n_estimators' : [int(x) for x in np.linspace(50,200,5)], 'max_depth' : [int(x) for x in np.linspace(2,20,3)], 'max_features' : [ 'sqrt', 'log2', None]}, RFparams = {'max_features' : [ 'sqrt', 'log2', None], 'criterion' : ['gini', 'entropy'], 'min_samples_split' : [int(x) for x in np.linspace(2,10,3)], 'max_depth': [50, 500, 1000, 2000, None], 'n_estimators':[int(x) for x in np.linspace(50, 100,3)]}, enet=LogisticRegression(penalty= 'elasticnet', solver='saga', max_iter=1e6, class_weight='balanced'), lr = LogisticRegression(solver='saga', max_iter=1e7, class_weight='balanced'), gb = GradientBoostingClassifier(), rf = RandomForestClassifier(class_weight='balanced')):
        # self.Xpath = Xpath
        # self.Ypath = Ypath Xpath, Ypath,
        self.n_splits = n_splits
        self.repeat = repeat
        self.n_jobs = n_jobs
        self.Eparams = Eparams
        self.LRparams = LRparams
        self.GBparams = GBparams
        self.RFparams = RFparams
        self.enet = enet
        self.lr = lr
        self.gb = gb
        self.rf = rf


    def XY_index(self, X, Y, form='csv'): # Xpath, Ypath,
        if form=='csv':
            X=pd.read_csv(X, sep='\t', index_col=0)
            Y=pd.read_csv(Y, sep='\t', index_col=0)
            Y.columns = ['RESPONSE']
            X = Y.join(X)
            Y = X['RESPONSE']
            X = X.drop(columns=['RESPONSE'])
        elif form =='pickle':
            X=pd.read_pickle(X) #, sep='\t', index_col='STRAIN'
            # Y= pd.read_pickle(Y) #, sep='\t', index_col='STRAIN'
            Y=pd.read_csv(Y, sep=' ', header=None, index_col=0).drop(columns=[1,3])
            Y.columns=['RESPONSE']
            X = Y.join(X)
            Y = X['RESPONSE']
            X = X.drop(columns=['RESPONSE'])
        if X.isnull().all().all():
             print("Your indices do not match between X and Y!")
        # if X is all nan, something is wrong
        return X, Y
    
    #### How to orgainze the dataframe from score_list
    def Tuning_Kfold_classify(self, X, Y, repeat, n_splits, score, mod, hyperparameters, n_jobs):
        dfL = []
        for i in range(0,repeat):
            skf = StratifiedKFold(n_splits=n_splits, random_state=i, shuffle=True) 
            boosted_grid = GridSearchCV(mod, hyperparameters, scoring=score, cv=skf, verbose=0, refit='AUC', error_score=np.nan, n_jobs=n_jobs, return_train_score=True)
            grid_fit = boosted_grid.fit(X, Y)
            DF = pd.DataFrame(grid_fit.cv_results_)
            DF['Iteration'] = i
            dfL.append(DF)
        DFall = pd.concat(dfL)
        return DFall

#### Need to nest a for loop above to extract training and testing data.
# https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.StratifiedKFold.html