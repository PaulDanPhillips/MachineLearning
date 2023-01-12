import pandas as pd
import os
import re
# import glob
import matplotlib.pyplot as plt
import numpy as np
from plotnine import *
from plotnine import ggplot
from scipy import stats
import seaborn as sns
import pickle


class Param_analysis:
    """
    Analyzing the tuning of machine learning models

    Attributes
    ----------
    path : str
        working directory

    Methods
    ----------
    RandomSearchAnalysis
        Investigate which set of features was found most often in the top scoring models
    RandSearchParamGrid
        A dictionary of hyperparamters from randomsearch that will be used in gridsearch.
    Best_mean
        Generating the a smaller dataframe of only the best scoring models.
    Param_figs
        Faceted and colored boxpots
    KfoldPairedTtest
        calculating pvalues from a kfold corrected ttest.
    """
    # def __init__(self, path):
    #     self.path = path
    def __init__(self):
            pass

    def RandomSearchAnalysis(self, DF, percBest, RoundTestScore, TestScoreColumnName, figTitle='Random Forest Parameter Count \n of Top scorers from RandomSearchCV', figname="RandSearchFIG.png"):
        """
        Investigate which set of features was found most often in the top scoring models

        Parameters
        ----------
        DF : dataframe
            The output dataframe from the Tuning_RandomSearch_classify function
        percBest : float
            The percent of best scorers to visualize
        RoundTestScore : int
            The number of significant figures to round the TestScoreColumnName to. This is important for the figure as without rounding the axis labels will be unruly.
        TestScoreColumnName : str
            The scoring column in DF that should be analyzed (should either be mean_test_score, or mean_train_score)
        figTitle : str
            The title displayed on the figure
        figname : str
            The filename to save the figure as

        Returns
        ----------
        countdata list
            list of dataframes that contain the counts of how often each hyperparameter was found in the top scorers
        """
        #### Read in Data, collect highest scoring 25%, and collect only hyperparameter columns and mean_test_score(will only work if one scoring metric is used)
        tmp = int(len(DF)*percBest)
        DF[TestScoreColumnName] = round(DF[TestScoreColumnName],RoundTestScore)
        filter_col = [col for col in DF if col.startswith('param_')]
        filter_col.append(TestScoreColumnName)
        DF = DF[filter_col]
        colnames = []
        for item in filter_col:
            new = re.sub('param_', '', item)
            colnames.append(new)
        DF.columns = colnames
        DFbest = DF.nlargest(tmp, TestScoreColumnName)

        #### Count the number of times each value of each hyperparamter is seen in the highest scoring 25%.
        DFcount = DFbest.copy()
        DFcount.max_depth.replace(np.nan, 'None', regex=True, inplace=True)
        DFcount.fillna('None', inplace=True)
        DFcount.max_depth = DFcount.max_depth.astype('category')
        DFcount.mean_test_score = round(DFcount.mean_test_score,2)
        countdata = []
        for col in DFcount:
            DFcount[col] = DFcount[col].sort_values(ascending=False)
            countdata.append(DFcount[col].value_counts(dropna=False).reset_index())

        #### Plotting
        fig, axes = plt.subplots(nrows=int((len(countdata)+1)/2), ncols=2)
        Flat_ax = axes.flatten()
        tmpest = zip(countdata, Flat_ax)
        for i in tmpest:
            sns.barplot(data=i[0], x=i[0].iloc[0:,0], y=i[0].iloc[0:,1], order=i[0]['index'], ax=i[1])
            i[1].set_title(str(i[0].columns[1]))
            i[1].set_ylabel('') 
            i[1].set_xlabel('')
            i[1].tick_params(rotation=45)
            i[1].set_title(str(i[0].columns[1]))
        fig.suptitle(figTitle)
        fig.tight_layout() 
        fig.savefig(figname, dpi=300)
        return countdata


    def RandSearchParamGrid(self, countdata, num_opt=3, outdir=os.getcwd(), outfile='RFGridParamSpace.pickle'):
        """
        Generate a dictionary of hyperparamters from Tuning_RandomSearch_classify that will be used in gridsearch.

        Parameters
        ----------
        countdata : list
            list of dataframes that contain the counts of how often each hyperparameter was found in the top scorers
        num_opt : int
            The number of the most frequently found hyperparameters to utilize downstream in gridsearch tuning.
        outdir : str
            Where to save the dictionary as a pickle file
        outfile : str
            The filename to save the dictionary of hyperparameters
        """

        RFGridParamSpace = {}
        # for i in range(len(countdata)-1):
        #     if len(countdata[i])<=num_opt:
        #         tmp = countdata[i]['index'].to_list()
        #         RFGridParamSpace[countdata[i].columns[1]] = tmp
        #     elif len(countdata[i])>=num_opt+1:
        #         tmp = countdata[i]['index'].to_list()[0:num_opt]
        #         RFGridParamSpace[countdata[i].columns[1]] = tmp
        # for i in range(len(RFGridParamSpace['max_features'])):
        #     if RFGridParamSpace['max_features'][i]=='None':
        #         RFGridParamSpace['max_features'][i] = None
        #     elif RFGridParamSpace['max_features'][i] != 'sqrt' and RFGridParamSpace['max_features'][i] != 'log2':
        #         RFGridParamSpace['max_features'][i] = float(RFGridParamSpace['max_features'][i])
        # for i in range(len(RFGridParamSpace['max_depth'])):
        #     if RFGridParamSpace['max_depth'][i]=='None':
        #         RFGridParamSpace['max_depth'][i] = None
        # print(RFGridParamSpace)
        # os.chdir(outdir)   
        # with open(outfile, 'wb') as handle:
        #     pickle.dump(RFGridParamSpace, handle, protocol=pickle.HIGHEST_PROTOCOL)
        # return RFGridParamSpace
        for i in range(len(countdata)-1):
            if len(countdata[i])<=num_opt:
                tmp = countdata[i]['index'].to_list()
                RFGridParamSpace[countdata[i].columns[1]] = tmp 
            elif len(countdata[i])>=num_opt+1:
                tmp = countdata[i]['index'].to_list()[0:num_opt]
                RFGridParamSpace[countdata[i].columns[1]] = tmp 
        if 'max_features' in RFGridParamSpace:### Need a control here if max_features and max_depth not hyperparameters in model.
            for i in range(len(RFGridParamSpace['max_features'])):
                if RFGridParamSpace['max_features'][i]=='None':
                    RFGridParamSpace['max_features'][i] = None
                elif RFGridParamSpace['max_features'][i] != 'sqrt' and RFGridParamSpace['max_features'][i] != 'log2':
                    RFGridParamSpace['max_features'][i] = float(RFGridParamSpace['max_features'][i])
        if 'max_depth' in RFGridParamSpace:
            for i in range(len(RFGridParamSpace['max_depth'])):
                if RFGridParamSpace['max_depth'][i]=='None':
                    RFGridParamSpace['max_depth'][i] = None
        print(RFGridParamSpace)
        os.chdir(outdir)   
        with open(outfile, 'wb') as handle:
            pickle.dump(RFGridParamSpace, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return RFGridParamSpace


    def Best_mean(self, df, TestMetric, n_best):
        """
        Generating the a smaller dataframe of only the best scoring models.

        Parameters
        ----------
        df : dataframe
            The dataframe output from Tuning_GridSearch_classify function
        TestMetric : str
            The column name you would like analyzed
        n_best : int
            The number of best scoring models you would like to collect

        Returns
        ----------
        DF_BEST : dataframe
            The best scoring parameters found
        """
        DF_BEST = df.groupby('params')[TestMetric].mean().nlargest(n_best)
        df.set_index('params', inplace=True)
        DF_BEST = df[df.index.isin(DF_BEST.index)].copy()
        return DF_BEST

    def Param_figs(self, BEST_params, TestMetric, Expname, fill_param=None, facet_param=None): 
        """
        Generating faceted and colored boxpots

        Parameters
        ----------
        BEST_params : dataframe
            The selected amount of best scoring models from Best_mean function
        TestMetric : str
            The column name you would like analyzed
        Expname : str
            The experiment name that you would like to see in the figure title
        fill_param : str
            The column to color all of the box plots by
        facet_param : str
            The column to facet around
        
        Returns
        ----------
        Best_Box : plotnine object
            The faceted and colored boxplot displaying the best scorers from the Tuning_GridSearch_classify function
        """
        if fill_param==None and facet_param != None:
            BEST_params.reset_index(inplace=True)
            Best_Box=(ggplot(BEST_params) 
                + aes(x='params', y=TestMetric) #, fill='test' , fill= 'Model'
                # + scale_fill_manual(values=color_dict)
                + geom_boxplot() #show_legend=False
                # + facet_grid(str(facet_param1)+' ~ '+str(facet_param2), scales='free_x')
                + facet_wrap(facet_param, scales='free_x')
                # + theme(axis_text_x=element_blank())
                + theme(axis_text_x = element_text(angle = 45, hjust =1))
                + labs(title= str(Expname)+' Best Scoring Parmeters for '+str(TestMetric), fill=fill_param)
                )
            return Best_Box #, Best_worst_Box

        elif facet_param==None and fill_param != None:
            BEST_params.reset_index(inplace=True)
            Best_Box=(ggplot(BEST_params) 
                + aes(x='params', y=TestMetric, fill='factor('+str(fill_param)+')') #, fill='test' , fill= 'Model'
                # + scale_fill_manual(values=color_dict)
                + geom_boxplot() #show_legend=False
                # + facet_grid(str(facet_param1)+' ~ '+str(facet_param2), scales='free_x')
                # + facet_wrap(facet_param, scales='free_x')
                # + theme(axis_text_x=element_blank())
                + theme(axis_text_x = element_text(angle = 45, hjust =1))
                + labs(title= str(Expname)+' Best Scoring Parmeters for '+str(TestMetric), fill=fill_param)
                )
            return Best_Box #, Best_worst_Box

        elif facet_param==None and fill_param == None:
            BEST_params.reset_index(inplace=True)
            Best_Box=(ggplot(BEST_params) 
                + aes(x='params', y=TestMetric) #, fill='test' , fill= 'Model'
                # + scale_fill_manual(values=color_dict)
                + geom_boxplot() #show_legend=False
                # + facet_grid(str(facet_param1)+' ~ '+str(facet_param2), scales='free_x')
                # + facet_wrap(facet_param, scales='free_x')
                # + theme(axis_text_x=element_blank())
                + theme(axis_text_x = element_text(angle = 45, hjust =1))
                + labs(title= str(Expname)+' Best Scoring Parmeters for '+str(TestMetric), fill=fill_param)
                )
            return Best_Box #, Best_worst_Box
        else:
            BEST_params.reset_index(inplace=True)
            Best_Box=(ggplot(BEST_params) 
                + aes(x='params', y=TestMetric, fill='factor('+str(fill_param)+')') #, fill='test' , fill= 'Model'
                # + scale_fill_manual(values=color_dict)
                + geom_boxplot() #show_legend=False
                # + facet_grid(str(facet_param1)+' ~ '+str(facet_param2), scales='free_x')
                + facet_wrap(facet_param, scales='free_x')
                # + theme(axis_text_x=element_blank())
                + theme(axis_text_x = element_text(angle = 45, hjust =1))
                + labs(title= str(Expname)+' Best Scoring Parmeters for '+str(TestMetric), fill=fill_param)
                )
            return Best_Box #, Best_worst_Box    

    # def Param_figs(self, BEST_params, TestMetric, Expname, fill_param, facet_param): 
    #     """
    #     Generating faceted and colored boxpots

    #     Parameters
    #     ----------
    #     BEST_params : dataframe
    #         The selected amount of best scoring models from Best_mean function
    #     TestMetric : str
    #         The column name you would like analyzed
    #     Expname : str
    #         The experiment name that you would like to see in the figure title
    #     fill_param : str
    #         The column to color all of the box plots by
    #     facet_param : str
    #         The column to facet around
        
    #     Returns
    #     ----------
    #     Best_Box : plotnine object
    #         The faceted and colored boxplot displaying the best scorers from the Tuning_GridSearch_classify function
    #     """
    #     BEST_params.reset_index(inplace=True)
    #     Best_Box=(ggplot(BEST_params) 
    #         + aes(x='params', y=TestMetric, fill='factor('+str(fill_param)+')') #, fill='test' , fill= 'Model'
    #         # + scale_fill_manual(values=color_dict)
    #         + geom_boxplot() #show_legend=False
    #         # + facet_grid(str(facet_param1)+' ~ '+str(facet_param2), scales='free_x')
    #         + facet_wrap(facet_param, scales='free_x')
    #         # + theme(axis_text_x=element_blank())
    #         + theme(axis_text_x = element_text(angle = 45, hjust =1))
    #         + labs(title= str(Expname)+' Best Scoring Parmeters for '+str(TestMetric), fill=fill_param)
    #         )
    #     return Best_Box #, Best_worst_Box

    def KfoldPairedTtest(self, r, k, n2, n1, a, b, score, ABX):
        '''  
        Perform repeated K-fold corrected paired T-test

        Parameters
        ----------
        r : int
            Number of repetitions.
        k : int
            Number of folds in cross-validation.
        n2 : int
            Number of total observations in testing folds int((1/k)*nsample*r)
        n1 : int
            Number of total observations in training folds int(((k-1)/k)*nsample*r)
        a : DataFrame
            1st DataFrame containing metric desired to be compared.
        b : DataFrame
            2nd DataFrame containing metric desired to be compared.
        score : str
            Column name of metric desired to be compared between a and b. The name of the columns needs to match in both a and b.
        ABX : str
            Used for Labeling purposes by titling the QQ-plots from InMatrixCompare function.

        Returns
        ----------
        t : float
            Tcrit of differences of metrics from a and b.
        pval: float
            Associated p-value calculated from Tcrit.
        fig : png
            QQ-plot on the differnce values between a and b.
        shp : Shapiro-Wilkes
            Tuple of W-statistic and associated p-value from Shapiro-Wilkes normality test.

        '''
        a.reset_index(inplace=True, drop=True)
        b.reset_index(inplace=True, drop=True)
        x = a[str(score)]-b[str(score)]
        x = x.drop_duplicates(keep='first')
        m = (1/((r*k)))*sum(x)
        if m == 0:
            t = 0
            pval = stats.t.sf(np.abs(t), (k*r)-1)*2 
            fig = np.nan #plt.figure()
            # ax = fig.add_subplot(111)
            # qq = stats.probplot(x, dist="norm", plot=ax)
            # ax.set_title(str(ABX)+' QQ-Plot')
            
            # shp = np.nan
        else:
            tmp = x-m
            tmp = tmp.pow(2)
            sigma2 = (1/(((k*r))-1))*sum(tmp)
            denom = (((1/((k*r))))+((n2/n1)))*sigma2
            t = m/(denom**(1/2))
            pval = stats.t.sf(np.abs(t), (k*r)-1)*2 
            # fig = plt.figure()
            # ax = fig.add_subplot(111)
            # qq = stats.probplot(x, dist="norm", plot=ax)
            # ax.set_title(str(ABX)+' QQ-Plot')
            # shp = stats.shapiro(x)

        return t, pval         