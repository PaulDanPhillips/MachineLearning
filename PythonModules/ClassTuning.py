# import timeit
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import KFold
import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from ClassIndexing import Indexing

class Tuning(Indexing):
    """
    Tuning of machine learning models

    Attributes
    ----------
    n_splits: int
        The K in K-fold cross-validation
    repeat: int
        The number of repetions to perform on different k splits
    n_jobs: int
        The number of cores to use to parallelize
    class_weight: str or None
        assigns the class wheight to either None or 'balanced'


    Methods
    ----------
    _timer_func
        A timer decorator function
    XY_index
        Aligns predictor and repsonse indices and calculates optimal n_splits and class_weight
    Tuning_RandomSearch_classify
        Using RandomSearch to preliminarily investigate large hyperparameter space to find subset of that space to evaluate further with GridSearch.
    Tuning_GridSearch_classify
        Using GridSearch to thoroughly investigate hyperparameter space to find optimal model for the given dataset.

    """
    def __init__(self, n_splits=None, repeat=None, n_jobs=None, class_weight=None):
        self.n_splits = n_splits
        self.repeat = repeat
        self.n_jobs = n_jobs
        self.class_weight = class_weight

    @Indexing._timer_func
    def Tuning_RandomSearch_classify(self, X, Y, repeat, n_splits, scorer, mod, hyperparameters, n_iter, n_jobs, stratify=True):
        """
        Using RandomSearch to preliminarily investigate large hyperparameter space to find subset of that space to evaluate further with GridSearch.

        Parameters
        ----------
        X : dataframe
            Predictor or Feature dataframe
        Y : series
            Repsonse variables
        repeat : int
            Number of stratified k-fold cross validations to perform
        n_splits: int
            number of splits in k-fold cross validation
        scorer : str
            Metric of which to compare model performance to calculate subsequent feature importances
        mod : model object
            scikit learn style model (e.g. Random Forest, XGBoost, logistic, or linear).
        hyperparameters : dict
            A dictionary of all hyperparameter keys and associated values desired to be tuned.
        n_iter: int
            The number of hyperparameters tested per each k-fold cross-validation splitting
        n_jobs : int
            The number of cores to utilize in parallelization
        """
        print("Number of k-fold cross-validations run is: " + str(repeat))
        print("Number of hyperparamter iterations tested is: " + str(n_iter))
        dfL = []
        for i in range(0,repeat):
            if stratify==True:
                cv = StratifiedKFold(n_splits=n_splits, random_state=i, shuffle=True) 
            else:
                cv = KFold(n_splits=n_splits, random_state=i, shuffle=True)
            boosted_grid = RandomizedSearchCV(mod, hyperparameters, n_iter=n_iter, scoring=scorer, cv=cv, verbose=0, refit=True, error_score=np.nan, n_jobs=n_jobs, return_train_score=True)
            grid_fit = boosted_grid.fit(X, Y)
            DF = pd.DataFrame(grid_fit.cv_results_)
            DF['Iteration'] = i
            dfL.append(DF)
        DFall = pd.concat(dfL)
        return DFall

    @Indexing._timer_func
    def Tuning_GridSearch_classify(self, X, Y, repeat, n_splits, scorer, mod, hyperparameters, n_jobs, stratify=True):
        """
        Using GridSearch to thoroughly investigate hyperparameter space to find optimal model for the given dataset..

        Parameters
        ----------
        X : dataframe
            Predictor or Feature dataframe
        Y : series
            Repsonse variables
        repeat : int
            Number of stratified k-fold cross validations to perform
        n_splits: int
            number of splits in k-fold cross validation
        scorer : str
            Metric of which to compare model performance to calculate subsequent feature importances
        mod : model object
            scikit learn style model (e.g. Random Forest, XGBoost, logistic, or linear).
        hyperparameters : dict
            A dictionary of all hyperparameter keys and associated values desired to be tuned.
        n_jobs : int
            The number of cores to utilize in parallelization
        """
        print("Number of repeats run is: " + str(repeat))
        dfL = []
        for i in range(0,repeat):
            if stratify==True:
                cv = StratifiedKFold(n_splits=n_splits, random_state=i, shuffle=True) 
            else:
                cv = KFold(n_splits=n_splits, random_state=i, shuffle=True)
            boosted_grid = GridSearchCV(mod, hyperparameters, scoring=scorer, cv=cv, verbose=0, refit=True, error_score=np.nan, n_jobs=n_jobs, return_train_score=True)
            grid_fit = boosted_grid.fit(X, Y)
            DF = pd.DataFrame(grid_fit.cv_results_)
            DF['Iteration'] = i
            dfL.append(DF)
        DFall = pd.concat(dfL)
        return DFall