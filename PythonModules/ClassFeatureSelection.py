from boruta import BorutaPy
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.inspection import permutation_importance
from ClassIndexing import Indexing

class FeatSelection(Indexing):
    """
    Selecting and Ranking of Features.

    Attributes
    ----------
    percent: int
        percent stregth of shadow features in Boruta.
    alpha: float
        threshold to reject corrected p-values.
    v: int
        level of verbosity in output
    max_iter: int
        The number of maximum iterations
    repeat:
        The number of repetitions to perform with a different random-seed.

    Methods
    ----------
    _timer_func
        A timer decorator function
    XY_index
        Aligns predictor and repsonse indices and calculates optimal n_splits and class_weight
    BorutaSelect
        Using Boruta to select most important features to Tree based model.
    BorutaOrganize
        Oranize the Boruta output into Selected and Tentative Dataframes with counts and percent representations.
    RegRegress
        Calculates the coefficient strength from a regularized linear or logistic regression model.
    PermImp
        Calculates the permutation feature importance of each feature in the dataframe.
    VIF
        Calculates the variance inflation factor (Should NOT be used with large datasets)
    IterVIF
        Uses above VIF and a selected threshold to drop any features with a VIF above.
    """
    def __init__(self, percent=100, alpha=0.05, v=0, max_iter=2500, repeat=1):
        self.percent = percent
        self.alpha = alpha
        self.v = v
        self.max_iter = max_iter
        self.repeat = repeat

    # @Indexing._timer_func # Not sure why this doesn't work. Maybe simply changing varaible names in _timer_func?
    def BorutaSelect(self, mod, percent, alpha, X, Y, Borut_Est, v, max_iter, rand): 
        """Selecting most influential Features using the Boruta algorihtm.
        
        Parameters
        ----------
        mod : model object
            scikit learn style tree based model (e.g. random forest or XGBoost).
        percent : int [0-100]
            percent stregth of shadow features in Boruta
        alpha : float
            threshold to reject corrected p-values.
        X : dataframe
            Predictor or Feature dataframe
        Y : series
            Repsonse variables
        Borut_Est : int
            Number of trees in the tree based model
        v : int
            level of verbosity in output
        max_iter : int
            The number of maximum iterations
        rand : int
            random state of model

        Returns
        ----------
        FEATURES: dataframe
            All Boruta selected features over all repetitions.
        """

        y=Y.values
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
        """
        Oranize the Boruta output into Selected and Tentative Dataframes with counts and percent representations.

        Parameters
        ----------
        FEATURES : dataframe
            dataframe of all Boruta selected features.
        repeat : int
            The number of repetitions Boruta was run

        Returns
        ----------
        GENE_SELECT_DF : dataframe
            Each Boruta selected feature with number of repetitions it was selected and percent selected.
        GENE_TENTATIVE_DF : dataframe
            Each feature neither boruta selected or removed with number of repetitions it was selected and percent selected.
        """
        GENE_SELECT = pd.concat(FEATURES)
        GENE_SELECT_DF=pd.DataFrame(GENE_SELECT.SELECTED.value_counts())
        GENE_SELECT_DF['PERCENT_SELECT'] = GENE_SELECT_DF.SELECTED/repeat
        GENE_TENTATIVE_DF=pd.DataFrame(GENE_SELECT.TENTATIVE.value_counts())
        GENE_TENTATIVE_DF['PERCENT_SELECT'] = GENE_TENTATIVE_DF.TENTATIVE/repeat
        return GENE_SELECT_DF, GENE_TENTATIVE_DF
    
    def RegRegress(self, mod, X, Y, rand):
        """
        Calculates the coefficient strength from a regularized linear or logistic regression model.

        Parameters
        ----------
        mod : model object
            scikit learn style penalized regresion based model (e.g. logistic or linear).
        X : dataframe
            Predictor or Feature dataframe
        Y : series
            Repsonse variables
        rand : int
            random state of model

        Returns
        ----------
        Coefs : dataframe
            sorted dataframe of model coeficients
        """
        mod.random_state=rand
        mod.fit(X,Y)
        Coefs = pd.DataFrame(mod.coef_)
        Coefs.columns = X.columns
        Coefs=Coefs.loc[:, (Coefs != 0).any(axis=0)].T.rename(columns={0:"VarStrength"})
        Coefs = Coefs['VarStrength'].sort_values(ascending=False, key=abs)
        return Coefs
    
    # @Indexing._timer_func
    def PermImp(self, mod, X, Y, internalRep, scorer, rand):
        """
        Calculates the permutation feature importance of each feature in the dataframe.

        Parameters
        ----------
        mod : model object
            scikit learn style model (e.g. Random Forest, XGBoost, logistic, or linear).
        X : dataframe
            Predictor or Feature dataframe
        Y : series
            Repsonse variables
        internalRep : int
            Number of of times to permute a single feature
        scorer : str
            Metric of which to compare model performance to calculate subsequent feature importances
        rand : int
            random state of model

        Returns
        ----------
        Importances : dataframe
            All non-zero feature importances.
        """

        mod.random_state=rand
        mod.fit(X, Y)
        PI = permutation_importance(mod, X, Y, n_repeats=internalRep, scoring=scorer, random_state=rand)
        Importances = pd.DataFrame(PI.importances_mean)
        Importances.index = X.columns.to_list()
        Importances=Importances[Importances[0]!=0].copy()
        return Importances

    def FeatImp(self, mod, X, Y, rand):
        """
        Calculates the permutation feature importance of each feature in the dataframe.

        Parameters
        ----------
        mod : model object
            scikit learn style model (e.g. Random Forest, XGBoost, logistic, or linear).
        X : dataframe
            Predictor or Feature dataframe
        Y : series
            Repsonse variables
        rand : int
            random state of model

        Returns
        ----------
        Importances : dataframe
            All non-zero feature importances.
        """

        mod.random_state=rand
        mod.fit(X, Y)
        FeatImp = mod.feature_importances_
        Importances = pd.DataFrame(FeatImp)
        Importances.index = X.columns.to_list()
        Importances=Importances[Importances[0]!=0].copy()
        return Importances

    def VIF(self, X):
        """
        Calculates the variance inflation factor (Should NOT be used with large datasets)

        Parameters
        ----------
        X : dataframe
            Predictor or Feature dataframe

        Returns
        ----------
        vif : series
            VIF values for each feature in the X dataframe.
        """
        corr_mat = np.array(X.corr())
        inv_corr_mat = np.linalg.inv(corr_mat)
        vif=pd.Series(np.diag(inv_corr_mat), index=X.columns)
        return vif
    
    def IterVIF(self, X, value):
        """
        Uses above VIF and a selected threshold to drop any features with a VIF above.

        Parameters
        ----------
        X : dataframe
            Predictor or Feature dataframe
        value: int
            Threshold of VIF to remove all features in X greater than.
        
        Returns
        ----------
        X : dataframe
            A new X with less features that have been downsized based on their VIF values
        """
        while True:
            vif = self.VIF(X)
            # print(vif)
            if vif.max() >=value:
                X.drop(columns=[vif.idxmax()], inplace=True)
                # print(X.columns)
            else:
                break
        return X # return vif or X?