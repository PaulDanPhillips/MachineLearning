import pandas as pd
import os
import glob
import numpy as np
from plotnine import *
from plotnine import ggplot
from pathlib import Path
import csv



class Simulation_Analysis:
    def __init__(self, cntrl_path, pyseer_path, boruta_path, regress_path):
        self.cntrl_path = cntrl_path
        self.pyseer_path = pyseer_path
        self.boruta_path = boruta_path
        self.regress_path = regress_path
        # self.total = sum(1 for line in open(total))

    def control_gene_lst(self, cntrl_path):
    
        os.chdir(cntrl_path)
        all_files=sorted(glob.glob("*.par"))# add sorted around all glob.glob
        cntrl_lst_DF = []

        for filename in all_files:
            df1 = pd.read_csv(filename, sep='\t')
            df1['Simulation'] = str(filename)
            df1['Simulation'] = df1['Simulation'].replace("\.par", '', regex=True)
            df1.set_index('QTL', inplace=True)
            cntrl_lst_DF.append(df1)
        # order = [0, 2, 3, 4, 5, 6, 7, 8, 9, 1] # Needed to be re-ordered as the gene list by simulation number
        # cntrl_lst_DF = [cntrl_lst_DF[i] for i in order]    
        cntrlDF =pd.concat(cntrl_lst_DF).sort_values('Simulation')
        grouped_cntrl = cntrlDF.groupby('Simulation')
        cntrl_lst_DF = []
        for i in grouped_cntrl:
            cntrl_lst_DF.append(i[1])
        
        return cntrl_lst_DF


    def totalFeats_lst(self, vcfpath, header=13):
        os.chdir(vcfpath)
        all_files=sorted(glob.glob('*.vcf')) # add sorted around all glob.glob
        TotalFeatsLST=[]
        simLST = []
        TotalFeatsDFlst = []
        for filename in all_files:
            print(str(filename))
            with open(filename) as filehandle:
                TotalFeats = sum(1 for line in filehandle)-header
            TotalFeatsLST.append(TotalFeats)
            simLST.append(str(filename).replace("_sims.vcf", ''))
        df2 = pd.DataFrame({'SimID':simLST, 'TotalFeats':TotalFeatsLST})
        TotalFeatsDFlst.append(df2)
        TotalFeatsDF = pd.concat(TotalFeatsDFlst).sort_values('SimID')
        grouped = TotalFeatsDF.groupby('SimID')
        TotalFeatsOrder = []
        for i in grouped:
            TotalFeatsOrder.append(i[1])
        return TotalFeatsOrder


    def confusion_org(self, Confusion):
        Means_DF =pd.DataFrame(Confusion.groupby('Threshold')['PosPredVal'].mean()).rename(columns={'PosPredVal':'PosPredVal_mean'})
        Vars_DF = pd.DataFrame(Confusion.groupby('Threshold')['PosPredVal'].var()).rename(columns={'PosPredVal':'PosPredVal_variance'})
        
        STATS_DF = Means_DF.join(Vars_DF)
        THRESH_DF = STATS_DF.join(Confusion.set_index('Threshold')).reset_index()
        THRESH_DF['PosPredVal_meanRound'] = THRESH_DF['PosPredVal_mean'].round(4)
        THRESH_DF['PosPredVal_varianceRound'] = "+/-" + THRESH_DF['PosPredVal_variance'].round(5).astype(str)
        return THRESH_DF


    def ConfusionTableExport(self, BorConfusion, PysConfusion, RegConfusion):
        CONFUSION = pd.concat([BorConfusion, PysConfusion, RegConfusion])
        CONFUSION=CONFUSION.round({'TruePos_mean':1, 'FalsePos_mean':1, 'TrueNeg_mean':1, 'FalseNeg_mean':1,'PosPredVal_mean':3, 'TPR_mean':3, 'TNR_mean':3})
        print(CONFUSION.to_markdown())
        return CONFUSION
# =============================================================================
# PYSEER ANALYSIS
# =============================================================================

# class pyseerAnalysis(Simulation_Analysis):
    def pyseer_lst(self, pyseer_path):
        os.chdir(pyseer_path)
        #### SORTED CHANGED ORDER. THE FILES NEED TO BE RENAMED WITH SIM# AT END
        all_files=sorted(glob.glob("*_pyseer.txt"))# add sorted around all glob.glob
        pyseer_lst_DF = []
    
        for filename in all_files:
            df = pd.read_csv(filename, sep='\t')
            df['Simulation'] = str(filename)
            df['Simulation'] = df['Simulation'].replace("_pyseer\.txt", '', regex=True)
            df['PlotValue'] = -(np.log(df['lrt-pvalue'])/np.log(10))
            df['variant'] = df['variant'].replace('_', ':', regex=True)
            df[['chrom', 'loci', 'SNP:SNP']] = df['variant'].str.split(':', 2, expand=True)
            df.set_index('variant', inplace=True)
            cols = [ 'af', 'filter-pvalue', 'lrt-pvalue',
            'beta', 'beta-std-err', 'intercept', 'notes', 'chrom', 'loci', 'SNP:SNP', 'Simulation', 'PlotValue']
            df = df[cols]
            df = df.sort_values('PlotValue', ascending=False)
            pyseer_lst_DF.append(df)
    
        # order = [0, 2, 3, 4, 5, 6, 7, 8, 9, 1] # Needed to be re-ordered as the gene list by simulation number
        # pyseer_lst_DF = [pyseer_lst_DF[i] for i in order]
        
        pyseerDF= pd.concat(pyseer_lst_DF).sort_values('Simulation')
        grouped_pyseer = pyseerDF.groupby('Simulation')
        pyseer_lst_DF = []
        for i in grouped_pyseer:
            pyseer_lst_DF.append(i[1])
        return pyseer_lst_DF


    def pyseer_select(self, pyseer_lst_df, cntrl_lst_DF):
        Pyseer_Select_lst = []
        for i in range(len(pyseer_lst_df)):
            Selected = cntrl_lst_DF[i].join(pyseer_lst_df[i], rsuffix='_pyseer').sort_values('PlotValue', ascending=False) # lsuffix='_cntrl'
            Pyseer_Select_lst.append(Selected)
        return Pyseer_Select_lst
    
    
    def pyseer_confusion(self,Pyseer_lst, Control_lst, Pyseer_thresh, Pyseer_Select_lst, Sim1TotalFeats, simExp=True):
        sim_lst = []
        FPlst = []
        TPlst = []
        FNlst = []
        PPVlst = []
        thresh = []
        TNlst = []
        # TN = len(Pyseer_lst[0])-len(Control_lst[0])
        for j in Pyseer_thresh:
            for i in range(0,len(Pyseer_lst)):
                # TN = len(Pyseer_lst[i])-len(Control_lst[i]) pyseer removes some as first step.
                if simExp==True:
                    TN = Sim1TotalFeats[i].iloc[0,1]
                elif simExp==False:
                    TN=Sim1TotalFeats
                else:
                    TN="BROKEN"
                sim = Pyseer_lst[i].Simulation[0]
                total_select = len(Pyseer_lst[i][Pyseer_lst[i]['lrt-pvalue']<=j])
                TP = len(Pyseer_Select_lst[i][Pyseer_Select_lst[i]['lrt-pvalue']<=j])
                FN = len(Control_lst[0])-TP
                FP = total_select-TP
                PPV = TP/(TP+FP)
                FPlst.append(FP)
                TPlst.append(TP)
                FNlst.append(FN)
                PPVlst.append(PPV)
                thresh.append(j)
                sim_lst.append(sim)
                TNlst.append(TN)
        Confusion = pd.DataFrame({'SimID':sim_lst, 'Threshold':thresh, 'TruePos':TPlst, 'FalsePos':FPlst, 'TrueNeg':TNlst, 'FalseNeg':FNlst, 'PosPredVal':PPVlst})
        
        Confusion['TPR'] = Confusion.TruePos/(Confusion.TruePos+Confusion.FalseNeg)
        Confusion['TNR'] = Confusion.TrueNeg/(Confusion.TrueNeg+Confusion.FalsePos)
        
        MeanConfusion=pd.DataFrame(Confusion.groupby('Threshold').mean()).rename(columns={0:'Mean'}).join(pd.DataFrame(Confusion.groupby('Threshold').std()).rename(columns={0:'STD'}), lsuffix='_mean', rsuffix='_STD')
        MeanConfusion['Method'] = 'pyseer'
        
        return Confusion, MeanConfusion
    
    
    def ThresholdStats_pyseer(self, PyseerConf_forBox):
        COUNT=PyseerConf_forBox[['TruePos', 'FalsePos']].stack().reset_index().set_index('level_0')
        COUNT = COUNT.join(PyseerConf_forBox).rename(columns={'level_1':'Group',0:'Count'})
        COUNT['SimID']=COUNT.SimID.replace('Sim','', regex=True)
        COUNT['SimID']=COUNT.SimID.astype(int)
        Count2 = COUNT.sort_values(['Group', 'Count'], ascending=[True, False])
        SimOrder = Count2.sort_values('SimID')['SimID'].tolist()
        Count2['Simulation'] = pd.Categorical(Count2['SimID'], categories=pd.unique(SimOrder), ordered=True)
        Count2['Group'] = Count2['Group'].replace('NUM_', '', regex=True)
        Count3 = pd.DataFrame(Count2.groupby(['Group', 'Threshold'])['Count'].mean()).rename(columns={'Count':'Mean'})
        Count4 = pd.DataFrame(Count2.groupby(['Group', 'Threshold'])['Count'].std()).rename(columns={'Count':'StDeviation'})
        
        Avgs_pyseer = Count3.join(Count4).reset_index()
        Avgs_pyseer['Method'] = 'PYSEER'
        return Avgs_pyseer
    

    def pyseer_scatter(self, Pyseer_lst, Pyseer_Select_lst, facet=True): #, facet =True

        All_Pyseer = pd.concat(Pyseer_lst)
        Correct_Pyseer = pd.concat(Pyseer_Select_lst)
        Correct_Pyseer['Control'] = 'Control'
        cols= Correct_Pyseer.columns
        All_Pyseer_plot = All_Pyseer.join(Correct_Pyseer, how='outer', rsuffix='_Drop')
        All_Pyseer_plot=All_Pyseer_plot[cols]
        All_Pyseer_plot['Control']=All_Pyseer_plot['Control'].fillna('Non-Control')
        All_Pyseer_plot=All_Pyseer_plot.sort_values('loci')
        # Facet_scatter=Sim1.pyseer_scatter(All_Pyseer_plot)
        
        Facet_scatter = (ggplot(All_Pyseer_plot)
            + aes(x='loci', y='PlotValue', color='Control')
            + geom_point()
            # + facet_wrap('Simulation')
            + labs(title='Facet Manhattan Plots per Simulation', x='Loci', y='-log(lrt-pval)/log(10)')
            + theme(axis_text_x=element_blank())
        )
        if facet == False:
            Facet_scatter=Facet_scatter
        elif facet == True:
            Facet_scatter += facet_wrap('Simulation')
            # Bar1 += scale_fill_manual(values = color_dict)
        else:
            Facet_scatter='Something is broken'
        # scatter= scatter+facet_wrap('Simulation')
        return Facet_scatter
        
# =============================================================================
# BORUTA ANALYSIS
# =============================================================================
# class bourtaAnalysis(Simulation_Analysis):
    def Boruta_lst(self, boruta_path):
        os.chdir(boruta_path)
        all_files=sorted(glob.glob('*_Selected.txt')) # *_SELECTION*.txt
        Boruta_lst_DF=[]
    
        for filename in all_files:
            df1 = pd.read_csv(filename, sep='\t').rename(columns={'Unnamed: 0':'GENE', '0':'ITER_SELECTED'})
            df1['Simulation'] = str(filename)
            df1['Simulation'] = df1['Simulation'].replace("_Selected\.txt", '', regex=True) #_GENE_SELECTION_DF\.txt
            df1['Simulation'] = df1['Simulation'].replace("[A-Za-z]+_", '', regex=True)
            df1['Simulation'] = df1['Simulation'].replace("\.txt", '', regex=True)
            df1.set_index('GENE', inplace=True)
            Boruta_lst_DF.append(df1)
            
        #### Manually reodering in python to ensure SimIDs match
        BorutaDF = pd.concat(Boruta_lst_DF).sort_values('Simulation')
        grouped_Boruta = BorutaDF.groupby('Simulation')
        Boruta_lst_DF = []
        for i in grouped_Boruta:
            Boruta_lst_DF.append(i[1])
        return Boruta_lst_DF
    
    def BorutaCorrectSelect(self, Boruta_lst_DF, cntrl_lst_DF):
        Boruta_Select_lst = []
        for i in range(len(Boruta_lst_DF)):
            Selected = cntrl_lst_DF[i].join(Boruta_lst_DF[i], how='inner', rsuffix='_Boruta').sort_values('PERCENT_SELECT', ascending=False) # lsuffix='_cntrl'
            Boruta_Select_lst.append(Selected)
        return Boruta_Select_lst
    

    
    def Boruta_Confusion(self, Sim1TotalFeats, Control_lst, Bor_thresh, Boruta_lst, Boruta_Correct_lst, simExp=True):
        thresh = []
        sim_lst = []
        FPlst = []
        TPlst = []
        FNlst = []
        PPVlst = []
        # TN = total-len(Control_lst[0])
        TNlst = []
        for j in Bor_thresh:
            for i in range(0,len(Boruta_lst)):
                if simExp==True:
                    TN = Sim1TotalFeats[i].iloc[0,1] # Better to Join on SimID
                elif simExp==False:
                    TN = Sim1TotalFeats
                else:
                    TN = "BROKEN"
                sim = Boruta_lst[i].Simulation[0]
                total_select=len(Boruta_lst[i][Boruta_lst[i].PERCENT_SELECT>=j])
                TP=len(Boruta_Correct_lst[i][Boruta_Correct_lst[i].PERCENT_SELECT>=j])
                FN = len(Control_lst[0])-TP
                FP = total_select-TP
                PPV = TP/(TP+FP)
                FPlst.append(FP)
                TPlst.append(TP)
                FNlst.append(FN)
                PPVlst.append(PPV)
                sim_lst.append(sim)
                thresh.append(j)
                TNlst.append(TN)
        Confusion = pd.DataFrame({'SimID':sim_lst, 'Threshold':thresh, 'TruePos':TPlst, 'FalsePos':FPlst, 'TrueNeg':TNlst, 'FalseNeg':FNlst, 'PosPredVal':PPVlst})
            
        Confusion['TPR'] = Confusion.TruePos/(Confusion.TruePos+Confusion.FalseNeg)
        Confusion['TNR'] = Confusion.TrueNeg/(Confusion.TrueNeg+Confusion.FalsePos)
        MeanConfusion=pd.DataFrame(Confusion.groupby('Threshold').mean()).rename(columns={0:'Mean'}).join(pd.DataFrame(Confusion.groupby('Threshold').std()).rename(columns={0:'STD'}), lsuffix='_mean', rsuffix='_STD')
        MeanConfusion['Method'] = 'boruta'
        return Confusion, MeanConfusion

    def ThresholdStats_boruta(self, BorutaConf_forBox, percent=100, alpha=0.05):
        COUNT=BorutaConf_forBox[['TruePos', 'FalsePos']].stack().reset_index().set_index('level_0')
        COUNT = COUNT.join(BorutaConf_forBox).rename(columns={'level_1':'Group',0:'Count'})
        COUNT['SimID']=COUNT.SimID.replace('Sim','', regex=True)
        COUNT['SimID']=COUNT.SimID.astype(int)
        Count2 = COUNT.sort_values(['Group', 'Count'], ascending=[True, False])
        SimOrder = Count2.sort_values('SimID')['SimID'].tolist()
        Count2['Simulation'] = pd.Categorical(Count2['SimID'], categories=pd.unique(SimOrder), ordered=True)
        Count2['Group'] = Count2['Group'].replace('NUM_', '', regex=True)
        Count3 = pd.DataFrame(Count2.groupby(['Group', 'Threshold'])['Count'].mean()).rename(columns={'Count':'Mean'})
        Count4 = pd.DataFrame(Count2.groupby(['Group', 'Threshold'])['Count'].std()).rename(columns={'Count':'StDeviation'})
        
        Avgs_boruta = Count3.join(Count4).reset_index()
        Avgs_boruta['Method'] = 'BORUTA'
        Avgs_boruta['Percent'] = int(percent)
        Avgs_boruta['alpha'] = float(alpha)
        return Avgs_boruta


# =============================================================================
# LR SELECTION
# =============================================================================
# class PenalizedRegressAnalysis(Simulation_Analysis):
    def Regress_lst(self, regress_path):
        os.chdir(regress_path)
        all_files=sorted(glob.glob('*_LRSelected.txt')) # add sorted around all glob.glob
        Regress_lst=[]
        
        for filename in all_files:
            df1 = pd.read_csv(filename, sep='\t').rename(columns={'Unnamed: 0':'GENE', '0':'ITER_SELECTED'})
            df1['Simulation'] = str(filename)
            df1['Simulation'] = df1['Simulation'].replace("_LRSelected\.txt", '', regex=True)
            df1['Simulation'] = df1['Simulation'].replace("[A-Za-z]+_", '', regex=True)
            df1['Simulation'] = df1['Simulation'].replace("\.txt", '', regex=True)
            df1.set_index('GENE', inplace=True)
            Regress_lst.append(df1)
            df1['PercentSelect'] = df1.Count/1000
            
            df1['RelativePercentStength']= np.nan
            for i in range(len(df1)):
                df1.iloc[i,4] = df1.iloc[i,0]/df1.iloc[0,0]
            
        #### Manually reodering in python to ensure SimIDs match
        RegressDF = pd.concat(Regress_lst).sort_values('Simulation')
        grouped_Boruta = RegressDF.groupby('Simulation')
        Regress_lst = []
        for i in grouped_Boruta:
            Regress_lst.append(i[1])
        return Regress_lst


    def Regress_select(self, Regress_lst, Control_lst):
        Regress_Select_lst = []
        for i in range(len(Regress_lst)):
            Selected = Control_lst[i].join(Regress_lst[i], how='inner', rsuffix='_Regress').sort_values('RelativePercentStength', ascending=False) # lsuffix='_cntrl'
            Regress_Select_lst.append(Selected)
        return Regress_Select_lst


    def Regress_Confusion(self, Sim1TotalFeats, Control_lst, Regress_thresh, Regress_lst, Regress_Select_lst, simExp=True):
        thresh = []
        sim_lst = []
        FPlst = []
        TPlst = []
        FNlst = []
        PPVlst = []
        # TN = total-len(Control_lst[0])
        TNlst = []
        for j in Regress_thresh:
            for i in range(0,len(Regress_lst)):
                if simExp==True:
                    TN = Sim1TotalFeats[i].iloc[0,1] # Better to Join on SimID
                elif simExp==False:
                    TN = Sim1TotalFeats
                else:
                    TN = "BROKEN"
                sim = Regress_lst[i].Simulation[0]
                total_select=len(Regress_lst[i][Regress_lst[i].RelativePercentStength>=j])
                TP=len(Regress_Select_lst[i][Regress_Select_lst[i].RelativePercentStength>=j])
                FN = len(Control_lst[0])-TP
                FP = total_select-TP
                PPV = TP/(TP+FP)
                FPlst.append(FP)
                TPlst.append(TP)
                FNlst.append(FN)
                PPVlst.append(PPV)
                sim_lst.append(sim)
                thresh.append(j)
                TNlst.append(TN)
        Confusion = pd.DataFrame({'SimID':sim_lst, 'Threshold':thresh, 'TruePos':TPlst, 'FalsePos':FPlst, 'TrueNeg':TNlst, 'FalseNeg':FNlst, 'PosPredVal':PPVlst})
            
        Confusion['TPR'] = Confusion.TruePos/(Confusion.TruePos+Confusion.FalseNeg)
        Confusion['TNR'] = Confusion.TrueNeg/(Confusion.TrueNeg+Confusion.FalsePos)
        MeanConfusion=pd.DataFrame(Confusion.groupby('Threshold').mean()).rename(columns={0:'Mean'}).join(pd.DataFrame(Confusion.groupby('Threshold').std()).rename(columns={0:'STD'}), lsuffix='_mean', rsuffix='_STD')
        MeanConfusion['Method'] = 'Regress'
        return Confusion, MeanConfusion

    def ThresholdStats_regress(self, Regress_Confusion):
        COUNT=Regress_Confusion[['TruePos', 'FalsePos']].stack().reset_index().set_index('level_0')
        COUNT = COUNT.join(Regress_Confusion).rename(columns={'level_1':'Group',0:'Count'})
        COUNT['SimID']=COUNT.SimID.replace('Sim','', regex=True)
        COUNT['SimID']=COUNT.SimID.astype(int)
        Count2 = COUNT.sort_values(['Group', 'Count'], ascending=[True, False])
        SimOrder = Count2.sort_values('SimID')['SimID'].tolist()
        Count2['Simulation'] = pd.Categorical(Count2['SimID'], categories=pd.unique(SimOrder), ordered=True)
        Count2['Group'] = Count2['Group'].replace('NUM_', '', regex=True)
        Count3 = pd.DataFrame(Count2.groupby(['Group', 'Threshold'])['Count'].mean()).rename(columns={'Count':'Mean'})
        Count4 = pd.DataFrame(Count2.groupby(['Group', 'Threshold'])['Count'].std()).rename(columns={'Count':'StDeviation'})
        
        Avgs_regress = Count3.join(Count4).reset_index()
        Avgs_regress['Method'] = 'LOGIT'
        return Avgs_regress


# =============================================================================
# FIGURE GENERATION
# =============================================================================
# class FigGen(Simulation_Analysis):

    def Comparison_box2(self, BorutaConf_forBox, PyseerConf_forBox, RegressConf_forBox, expname, BorIter):
        PyseerConf_forBox['Method'] = 'PYSEER'
        BorutaConf_forBox['Method'] = 'BORUTA'
        RegressConf_forBox['Method'] = 'LOGIT'
        Thresh_concat = pd.concat([PyseerConf_forBox, BorutaConf_forBox, RegressConf_forBox ])
        # Thresh_concat = pd.concat([PyseerConf_forBox, BorutaConf_forBox])
        Threshold_Box_comparison=(ggplot(Thresh_concat) #SHRUNK
            + aes(x='factor(Threshold)', y='PosPredVal') #, fill='factor(Threshold)'
            + geom_boxplot() #show_legend=False
            + labs(title= str(expname)+' Boruta vs. Pyseer Significant Threshold Analysis', x='Boruta Threshold of %Selected per '+str(BorIter)+' Iterations \nLogit Threshold of %Varibale Strength where the stongest variable is set to 100% \nPyseer Thresholds of lrt-pvalue \n(LOGIT = l1 penalized logistic regression)', y='Positive Predictive Value (n=10)') #, fill='Threshold'
            + theme(axis_text_x = element_text(angle = 45, hjust =1))
            # + geom_text(aes(y = 1.15, label="Means +/- Variances"), size=10)
            + geom_text(aes(x='factor(Threshold)', y = 1.1, label='factor(PosPredVal_meanRound)'), size=7, angle=45) #  min(Thresh_concat.PosPredVal) - 0.001, color='factor(Threshold)' inside aes
            + geom_text(aes(x='factor(Threshold)', y = 1.05, label='factor(PosPredVal_varianceRound)'), size=5, angle=45) # min(Thresh_concat.PosPredVal)- 0.005, color='factor(Threshold)')
            # + ylim(0.0, 0.2)
            + facet_wrap('Method', scales='free_x')
            + theme(subplots_adjust={'wspace': 0.25})
            # + annotate('text', x='factor(Threshold)', y = 'MEANS', label='MEANS')
        )
        
        return Threshold_Box_comparison


    def Comparison_BarCounts(self, Avgs_pyseer, Avgs_boruta, Avgs_regress, expname, BorIter):

        CountConcat = pd.concat([Avgs_pyseer, Avgs_boruta, Avgs_regress])
        
        AverageCounts=(ggplot(CountConcat)
            + aes(x='factor(Threshold)', y='Mean', fill='Group', position_stack=False)
            + geom_col(stat='identity', position='dodge') #position_dodge()
            # + geom_text(aes(label='Count'), position='stack', size=9, va='center')
            # + geom_text(aes(label='Count'), position=position_stack(vjust=0.5), size=9, va='center')
            + geom_text(aes(label='Mean'), position=position_dodge(width=0.9), size=8, va='bottom')
            + geom_errorbar(aes(x="factor(Threshold)", ymin="Mean-StDeviation",ymax="Mean+StDeviation"), width=0.1, position=position_dodge(0.9))
            + theme(axis_text_x = element_text(angle = 45, hjust =1))
            + labs(title=str(expname)+' Average Correct and Incorrect Gene Selections', y='Average Number of Genes Selected (n=10)', x='Boruta Threshold of %Selected per '+str(BorIter)+' Iterations \nLogit Threshold of %Varibale Strength where the stongest variable is set to 100% \nPyseer Thresholds of lrt-pvalue \n(LOGIT = l1 penalized logistic regression)', fill='Genes Selected')
            # + facet_wrap('Threshold')
            + facet_wrap('Method', scales='free_x')
            # + ylim(0, 250)
        )
        return AverageCounts


    def Boruta_Param_Boxplots(self, Bourta_param_DF, expname):
        Threshold_Box_comparison=(ggplot(Bourta_param_DF) #SHRUNK
            + aes(x='factor(THRESH)', y='RATIO') #, fill='factor(THRESH)'
            + geom_boxplot() #show_legend=False
            + labs(title=str(expname)+' Boruta Percent and alpha analysis', x='Boruta Threshold of %Selected per 250 Iterations', y='Positive Predictive Value (n=10)') #, fill='Threshold'
            + theme(axis_text_x = element_text(angle = 45, hjust =1))
            + geom_text(aes(x='factor(THRESH)', y = -0.05, label='factor(MeanRound)'), size=7) # , color='factor(THRESH)' inside aes
            + geom_text(aes(x='factor(THRESH)', y = -0.07, label='factor(VARRound)'), size=6) # color='factor(THRESH)')
            # + ylim(0.0, 0.2)
            + facet_grid('Percent ~ alpha', scales='free_x')
            # + annotate('text', x='factor(THRESH)', y = 'MEANS', label='MEANS')
        )
        return Threshold_Box_comparison
    
    def Boruta_Param_Barplots(self, CountConcat, expname):
        AverageCounts=(ggplot(CountConcat)
            + aes(x='factor(THRESH)', y='Mean', fill='Group', position_stack=False)
            + geom_col(stat='identity', position='dodge') #position_dodge()
            # + geom_text(aes(label='Count'), position='stack', size=9, va='center')
            # + geom_text(aes(label='Count'), position=position_stack(vjust=0.5), size=9, va='center')
            + geom_text(aes(label='Mean'), position=position_dodge(width=0.9), size=8, va='bottom')
            + geom_errorbar(aes(x="factor(THRESH)", ymin="Mean-StDeviation",ymax="Mean+StDeviation"), width=0.1, position=position_dodge(0.9))
            + theme(axis_text_x = element_text(angle = 45, hjust =1))
            + labs(title=str(expname)+' Average Correct and Incorrect Gene Selections', y='Average Number of Genes Selected (n=10)', x='Boruta Threshold of %Selected per 250 Iterations', fill='Genes Selected')
            # + facet_wrap('THRESH')
            # + facet_wrap('Method', scales='free_x')
            + facet_grid('Percent ~ alpha', scales='free_x')
            # + ylim(0, 250)
        )
        return AverageCounts


# =============================================================================
# GOOD CONFUSION MATRIX
# =============================================================================
    # def plot_confusion_matrix(cm,
    #                           target_names,
    #                           title='Confusion matrix',
    #                           cmap=None,
    #                           normalize=True):
    #     """
    #     given a sklearn confusion matrix (cm), make a nice plot
    
    #     Arguments
    #     ---------
    #     cm:           confusion matrix from sklearn.metrics.confusion_matrix
    
    #     target_names: given classification classes such as [0, 1, 2]
    #                   the class names, for example: ['high', 'medium', 'low']
    
    #     title:        the text to display at the top of the matrix
    
    #     cmap:         the gradient of the values displayed from matplotlib.pyplot.cm
    #                   see http://matplotlib.org/examples/color/colormaps_reference.html
    #                   plt.get_cmap('jet') or plt.cm.Blues
    
    #     normalize:    If False, plot the raw numbers
    #                   If True, plot the proportions
    
    #     Usage
    #     -----
    #     plot_confusion_matrix(cm           = cm,                  # confusion matrix created by
    #                                                               # sklearn.metrics.confusion_matrix
    #                           normalize    = True,                # show proportions
    #                           target_names = y_labels_vals,       # list of names of the classes
    #                           title        = best_estimator_name) # title of graph
    
    #     Citiation
    #     ---------
    #     http://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html
    
    #     """
    #     import matplotlib.pyplot as plt
    #     import numpy as np
    #     import itertools
    
    #     accuracy = np.trace(cm) / np.sum(cm).astype('float')
    #     misclass = 1 - accuracy
    
    #     if cmap is None:
    #         cmap = plt.get_cmap('Blues')
    
    #     plt.figure(figsize=(8, 6))
    #     plt.imshow(cm, interpolation='nearest', cmap=cmap)
    #     plt.title(title)
    #     plt.colorbar()
    
    #     if target_names is not None:
    #         tick_marks = np.arange(len(target_names))
    #         plt.xticks(tick_marks, target_names, rotation=45)
    #         plt.yticks(tick_marks, target_names)
    
    #     if normalize:
    #         cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    
    
    #     thresh = cm.max() / 1.5 if normalize else cm.max() / 2
    #     for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
    #         if normalize:
    #             plt.text(j, i, "{:0.4f}".format(cm[i, j]),
    #                      horizontalalignment="center",
    #                      color="white" if cm[i, j] > thresh else "black")
    #         else:
    #             plt.text(j, i, "{:,}".format(cm[i, j]),
    #                      horizontalalignment="center",
    #                      color="white" if cm[i, j] > thresh else "black")
    
    
    #     plt.tight_layout()
    #     plt.ylabel('True label')
    #     plt.xlabel('Predicted label\naccuracy={:0.4f}; misclass={:0.4f}'.format(accuracy, misclass))
    #     plt.show()