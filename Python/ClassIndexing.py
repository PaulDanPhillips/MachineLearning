from time import time
import pandas as pd

class Indexing:
    """
    Ensuring Predictor and Response Dataframes contain the same index.

    Attributes
    ----------
    n_splits: int
        number of splits in k-fold cross validation
    class_weight: str or None
        assigns the class wheight to either None or 'balanced'

    Methods
    ----------
    _timer_func
        A timer decorator function
    XY_index
        Aligns predictor and repsonse indices and calculates optimal n_splits and class_wheight
    """
    def __init__(self, n_splits=None, class_weight=None):
        self.n_splits = n_splits
        self.class_weight = class_weight

    def _timer_func(func):
        """
        A timer decorator function

        Parameters
        ----------
        func : function
            A function whose runtime will be evaluated
        
        Returns
        ----------
        Func_result : output
            Whatever the output is for the function being timed
        wrap_func : func
            decorators do not alter the calling signature or return value of function being wrapped.
            DOES THIS RESET THE WRAP FUNCTION?
        """
    # This function shows the execution time of 
        def wrap_func(*args, **kwargs):
            t1 = time()
            Func_result = func(*args, **kwargs)
            t2 = time()
            # argnames = func.__code__.co_varnames[:func.__code__.co_argcount]
            print(f'Function args: {func.__name__!r} executed in {(t2-t1):.4f}s')
            return Func_result
        return wrap_func

    def XY_index(self, X, Y, BacGWASim=False): # Xpath, Ypath,
        """
        Aligns predictor and repsonse indices and calculates optimal n_splits and class_wheight

        Parameters
        ----------
        X : file
            Predictor or Feature dataframe
        Y : file
            Repsonse variables
        BacGWASim : bool
            If the data is the raw data from BacGWASim or not

        Returns
        ----------
        X : dataframe
            X predictor dataframe that has been joined to the repsonse variable based on index values and split back to just the predictor features
        Y : series
            Y response series that has been joined to the predictor variable based on index values and split back to just the response features
        n_splits: int
            If n_splits not deliberately defined, number of splits in k-fold cross validation
        class_weight: str or None
            Ic class_weight not deliberatley defined, assigns the class wheight to either None or 'balanced'
        """
        if BacGWASim==False:
            X=pd.read_csv(X, sep='\t', index_col=0)
            Y=pd.read_csv(Y, sep='\t', index_col=0)
            Y.columns = ['RESPONSE']
            X = Y.join(X)
            Y = X['RESPONSE']
            X = X.drop(columns=['RESPONSE'])
        elif BacGWASim ==True:
            X=pd.read_pickle(X) #, sep='\t', index_col='STRAIN'
            # Y= pd.read_pickle(Y) #, sep='\t', index_col='STRAIN'
            Y=pd.read_csv(Y, sep=' ', header=None, index_col=0).drop(columns=[1,3])
            Y.columns=['RESPONSE']
            X = Y.join(X)
            Y = X['RESPONSE']
            X = X.drop(columns=['RESPONSE'])
        if X.isnull().all().all():
             print("Your indices do not match between X and Y! Please use either raw BacGWASim files or tab delimeted txt files")
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