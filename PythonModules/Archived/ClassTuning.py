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
from sklearn.model_selection import RandomizedSearchCV
import argparse
from time import time
import inspect

class Tuning:
    def __init__(self, n_splits=None, repeat=None, n_jobs=None, class_weight=None, Eparams={'C': np.logspace(-3,2,5), 'l1_ratio': np.linspace(0,1,4)}, LRparams={'C': np.logspace(0,1,4), 'penalty':['l1', 'l2']}, GBparams = {'learning_rate' : np.linspace(0.0001,0.01,5), 'n_estimators' : [int(x) for x in np.linspace(50,200,5)], 'max_depth' : [int(x) for x in np.linspace(2,20,3)], 'max_features' : [ 'sqrt', 'log2', None]}, RFparams = {'max_features' : [ 'sqrt', 'log2', None], 'criterion' : ['gini', 'entropy'], 'min_samples_split' : [int(x) for x in np.linspace(2,10,3)], 'max_depth': [50, 500, 1000, 2000, None], 'n_estimators':[int(x) for x in np.linspace(50, 100,3)]}, enet=LogisticRegression(penalty= 'elasticnet', solver='saga', max_iter=1e6), lr = LogisticRegression(solver='saga', max_iter=1e7), gb = GradientBoostingClassifier(), rf = RandomForestClassifier()):
        # self.Xpath = Xpath
        # self.Ypath = Ypath Xpath, Ypath,
        self.n_splits = n_splits
        self.repeat = repeat
        self.n_jobs = n_jobs
        self.class_weight=class_weight
        self.Eparams = Eparams
        self.LRparams = LRparams
        self.GBparams = GBparams
        self.RFparams = RFparams
        self.enet = enet
        self.lr = lr
        self.gb = gb
        self.rf = rf
        
    def _timer_func(func):
    # This function shows the execution time of 
    # the function object passed
        def wrap_func(*args, **kwargs):
            t1 = time()
            result = func(*args, **kwargs)
            t2 = time()
            argnames = func.__code__.co_varnames[:func.__code__.co_argcount]
            func_args = inspect.signature(func).bind(*args, **kwargs).arguments
            # func_args_str = ", ".join(map("{0[0]} = {0[1]!r}".format, func_args.items()))
            # print(f"{func.__module__}.{func.__qualname__} ( {func_args_str} from {func.__name__!r} executed in {(t2-t1):.4f}s)")

            print(f'Function args: {func_args} from {func.__name__!r} executed in {(t2-t1):.4f}s')
            return result
        return wrap_func

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
        if self.n_splits==None:
            for i in Y.unique():
                if Y.value_counts()[i]>=12:
                    self.n_splits=10
                elif Y.value_counts()[i]>=5:
                    self.n_splits=5
                elif Y.value_counts()[i]>=3:
                    self.n_splits=3
                else:
                    print("You may want to consider using subset splitting or LOOCV")
                    self.n_splits=2
        else:
            self.n_splits=self.n_splits
        if self.class_weight==None:
            for i in Y.unique():
                if Y.value_counts()[i]/len(Y) <= 0.2:
                    self.class_weight='balanced'
                else:
                    self.class_weight=None
        else:
            self.class_weight=self.class_weight
        return X, Y, self.n_splits, self.class_weight
    
    #### How to orgainze the dataframe from score_list
    @_timer_func
    def Tuning_GridSearch_classify(self, X, Y, repeat, n_splits, score, mod, hyperparameters, n_jobs):
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
    
    @_timer_func
    def Tuning_RandomSearch_classify(self, X, Y, repeat, n_splits, score, mod, hyperparameters, n_iter, n_jobs):
        dfL = []
        for i in range(0,repeat):
            skf = StratifiedKFold(n_splits=n_splits, random_state=i, shuffle=True) 
            boosted_grid = RandomizedSearchCV(mod, hyperparameters, n_iter=n_iter, scoring=score, cv=skf, verbose=0, refit='AUC', error_score=np.nan, n_jobs=n_jobs, return_train_score=True)
            grid_fit = boosted_grid.fit(X, Y)
            DF = pd.DataFrame(grid_fit.cv_results_)
            DF['Iteration'] = i
            dfL.append(DF)
        DFall = pd.concat(dfL)
        return DFall
    