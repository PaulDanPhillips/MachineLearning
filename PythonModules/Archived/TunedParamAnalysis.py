import pandas as pd
import os
import re
import numpy as np
# import glob
# import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import *
from plotnine import ggplot
from scipy import stats
import pickle

class Param_analysis:
    def __init__(self, path):
        self.path = path
    
    
    def RandomSearchAnalysis(self, DF, percBest, RoundTestScore, TestScoreColumnName, outdir=os.getcwd(), figTitle='Random Forest Parameter Count \n of Top scorers from RandomSearchCV', figname="RandSearchFIG.png"):

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
        fig, axes = plt.subplots(nrows=int(len(countdata)/2), ncols=2)
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
        os.chdir(outdir)
        fig.savefig(figname, dpi=300)
        return countdata
    def RandSearchParamGrid(self, countdata, outdir=os.getcwd(), outfile='RFGridParamSpace.pickle'):
        RFGridParamSpace = {}
        for i in range(len(countdata)-1):
            if len(countdata[i])<=3:
                tmp = countdata[i]['index'].to_list()
                RFGridParamSpace[countdata[i].columns[1]] = tmp
            elif len(countdata[i])>=4:
                tmp = countdata[i]['index'].to_list()[0:3]
                RFGridParamSpace[countdata[i].columns[1]] = tmp
        for i in range(len(RFGridParamSpace['max_features'])):
            if RFGridParamSpace['max_features'][i]=='None':
                RFGridParamSpace['max_features'][i] = None
        for i in range(len(RFGridParamSpace['max_depth'])):
            if RFGridParamSpace['max_depth'][i]=='None':
                RFGridParamSpace['max_depth'][i] = None
        print(RFGridParamSpace)
        
        with open('RFGridParamSpace.pickle', 'wb') as handle:
            pickle.dump(RFGridParamSpace, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return RFGridParamSpace
    
    def File_input(self, path, LR_name=None, Enet_name=None, RF_name=None, GB_name=None):
        os.chdir(path)
        if LR_name!=None:
            LR = pd.read_csv(LR_name, sep='\t').rename(columns={'Unnamed: 0':'Iteration'})
            LR['Model'] = 'Logistic Regression'
        else:
            LR = 'Empty'
        
        if Enet_name!=None:
            Enet = pd.read_csv(Enet_name, sep='\t').rename(columns={'Unnamed: 0':'Iteration'})
            Enet['Model'] = 'Elastic Net'
        else:
            Enet = 'Empty'
        
        if RF_name!=None:
            RF = pd.read_csv(RF_name, sep='\t').rename(columns={'Unnamed: 0':'Iteration'})
            RF['Model'] = 'Random Forest'
        else:
            RF = 'empty'
        
        if GB_name!=None:
            GB = pd.read_csv(GB_name, sep='\t').rename(columns={'Unnamed: 0':'Iteration'})
            GB['Model'] = 'GBoost'
        else:
            GB = 'empty'
        return LR, Enet, RF, GB
    
    
    def BestNworst_mean(self, df, metric, n_bestNworst):
        DF_Best = df.groupby('params')[metric].mean().nlargest(n_bestNworst)
        DF_worst = df.groupby('params')[metric].mean().nsmallest(n_bestNworst)
        df2 = df.copy()
        df2.set_index('params', inplace=True)
        # DF_Best=df.Mean_Test_AUC.nlargest(n_bestNworst)
        DF_BEST = df2[df2.index.isin(DF_Best.index)]
        # DF_worst=df.Mean_Test_AUC.nsmallest(n_bestNworst)
        DF_WORST = df2[df2.index.isin(DF_worst.index)]
        # DF_ALL = pd.concat([DF_BEST, DF_WORST])
        return DF_BEST, DF_WORST
        
    def Param_figs(self, BEST_params, metric, Expname): #SHRUNK
        BEST_params.reset_index(inplace=True)
        Best_Box=(ggplot(BEST_params) #SHRUNK
            + aes(x='params', y=metric) #, fill='test' , fill= 'Model'
            # + scale_fill_manual(values=color_dict)
            + geom_boxplot() #show_legend=False
            # + theme(axis_text_x=element_blank())
            + theme(axis_text_x = element_text(angle = 45, hjust =1))
            + labs(title= str(Expname)+' Best Scoring Parmeters for '+str(metric))
            )
        # Best_worst_Box=(ggplot(SHRUNK) #SHRUNK
        #     + aes(x='Params', y='Mean_Test_AUC', fill= 'Model') #, fill='test'
        #     # + scale_fill_manual(values=color_dict)
        #     + geom_boxplot() #show_legend=False
        #     # + theme(axis_text_x=element_blank())
        #     + theme(axis_text_x = element_text(angle = 45, hjust =1))
        #     + labs(title=str(sim_num)+' Best & Worst Scoring Parmeters for AUC')
        #     )
        return Best_Box #, Best_worst_Box

    def KfoldPairedTtest(self, r, k, n2, n1, a, b, score, ABX):
        '''
        Description
        --------------------------------------------------------------------------    
        Perform repeated K-fold corrected paired T-test

        Parameters
        --------------------------------------------------------------------------
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
        --------------------------------------------------------------------------
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
            
            shp = np.nan
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
            shp = stats.shapiro(x)

        return t, pval, shp                 
        
# Sim_Exp1 = Param_analysis('H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/SIMULATED_DATA_ANALYSIS/250_Heff_Hcorr/TUNING/OUT')

# Enet, RF, GB = Sim_Exp2.File_input(Sim_Exp2.path, 'Sim10_Enet_results.txt', 'Sim10_RF_results.txt', 'Sim10_GB_results.txt')

# ENET_ALL, ENET_BEST= Sim_Exp2.BestNworst_mean(Enet, 5)
# RF_ALL, RF_BEST= Sim_Exp2.BestNworst_mean(RF, 5)
# GB_ALL, GB_BEST= Sim_Exp2.BestNworst_mean(GB, 5)

# AllModelAllParams = pd.concat([Enet, RF, GB])
# BestAndWorst = pd.concat([ENET_ALL, GB_ALL, RF_ALL])

# BEST = pd.concat([ENET_BEST, RF_BEST, GB_BEST])
# BestNWorst=AllModelAllParams[AllModelAllParams.index.isin(BestAndWorst.index)]
# BEST.reset_index(inplace=True)
# BestNWorst.reset_index(inplace=True)
# BestNWorst = BestNWorst.loc[BestNWorst['Mean_Test_AUC']>0.5]
# Bbox, BWbox = Sim_Exp2.Param_figs(BEST, BestNWorst, 'Sim1')

# BEST.Params = BEST.Params.replace(':', '=', regex=True)
# BEST.Params = BEST.Params.replace('{', '', regex=True)
# BEST.Params = BEST.Params.replace('}', '', regex=True)
# BEST.Params = BEST.Params.replace('\'', '', regex=True)


# Sim_Exp1 = Param_analysis('H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/SIMULATED_DATA_ANALYSIS/250_Heff_Mcorr/TUNING/OUT')

# Enet, RF, GB = Sim_Exp1.File_input(Sim_Exp1.path, 'Sim10_Enet_results.txt', 'Sim10_RF_results.txt', 'Sim10_GB_results.txt')

# ENET_ALL, ENET_BEST= Sim_Exp1.BestNworst_mean(Enet, 5)
# RF_ALL, RF_BEST= Sim_Exp1.BestNworst_mean(RF, 5)
# GB_ALL, GB_BEST= Sim_Exp1.BestNworst_mean(GB, 5)

# AllModelAllParams = pd.concat([Enet, RF, GB])
# BestAndWorst = pd.concat([ENET_ALL, GB_ALL, RF_ALL])

# BEST = pd.concat([ENET_BEST, RF_BEST, GB_BEST])
# BestNWorst=AllModelAllParams[AllModelAllParams.index.isin(BestAndWorst.index)]
# BEST.reset_index(inplace=True)
# BestNWorst.reset_index(inplace=True)
# BestNWorst = BestNWorst.loc[BestNWorst['Mean_Test_AUC']>0.5]
# Bbox, BWbox = Sim_Exp1.Param_figs(BEST, BestNWorst, 'Sim1')

# BEST.Params = BEST.Params.replace(':', '=', regex=True)
# BEST.Params = BEST.Params.replace('{', '', regex=True)
# BEST.Params = BEST.Params.replace('}', '', regex=True)
# BEST.Params = BEST.Params.replace('\'', '', regex=True)