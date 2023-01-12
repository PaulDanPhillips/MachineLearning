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
    
    # def totalFeats_lst(self, pickle_path, simExp=True, vcfFile=None, header=7):
    #     if simExp==True:
    #         os.chdir(pickle_path)
    #         all_files=sorted(glob.glob('*.pickle')) # add sorted around all glob.glob
    #         TotalFeats=[]
    #         print('Through naming files:')
    #         print(all_files)
    #         for filename in all_files:
    #             df1 = pd.read_pickle(filename)
    #             print(filename)
    #             df1['Simulation'] = str(filename)
    #             df1['Simulation'] = df1['Simulation'].replace("_sims\.pickle", '', regex=True)
    #             df2 = pd.DataFrame({'SimID':df1['Simulation'],'TotalFeats':len(df1.columns)-1})
    #             df2.set_index('SimID', inplace=True)
    #             df2 = pd.DataFrame(df2.reset_index().iloc[0,:]).T
    #             TotalFeats.append(df2)
    #             del df1
    #         print('Through for loop')
    #         print(TotalFeats)
    #         TotalFeatsDF = pd.concat(TotalFeats).sort_values('SimID')
    #         grouped = TotalFeatsDF.groupby('SimID')
    #         TotalFeatsOrder = []
    #         for i in grouped:
    #             TotalFeatsOrder.append(i[1])
            
    #     elif simExp == False:
    #         TotalFeatsOrder = sum(1 for line in open(vcfFile))-header
        
    #     else:
    #         TotalFeatsOrder = 'Broken'
    #     return TotalFeatsOrder
    
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



class pyseerAnalysis(Simulation_Analysis):
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
    
    def pyseer_calcs(self, Pyseer_Select_lst, pyseer_lst_DF):
        nuisance=[]
        right = []
        thresholds = []
        pyseer_analysis_lst = []
        sim_lst = []
        for i in range(0, len(Pyseer_Select_lst)):
            for j in range(0, len(Pyseer_Select_lst[i])):
                thresh = Pyseer_Select_lst[i].iloc[j,-1]
                sim = Pyseer_Select_lst[i].iloc[j,-2]
                sim_lst.append(sim)
                thresholds.append(thresh)
                nuisance.append(len(pyseer_lst_DF[i][pyseer_lst_DF[i]['PlotValue']>=thresh]))
                right.append(len(Pyseer_Select_lst[i][Pyseer_Select_lst[i]['PlotValue']>=thresh]))
                pyseerDF = pd.DataFrame({'NumCorrect':right, 'NumSelected':nuisance, 'THRESHOLD':thresholds, 'Simulation':sim_lst})
                pyseerDF['RATIO'] = pyseerDF['NumCorrect']/pyseerDF['NumSelected'] # Make to FPR = FP/FP+TN
                pyseerDF['FPR'] = pyseerDF['NumSelected']/(pyseerDF['NumCorrect']+pyseerDF['NumSelected']) # type1 error
                pyseerDF = pyseerDF.sort_values('RATIO', ascending=False)
                pyseer_analysis_lst.append(pyseerDF)
        PYSEER = pd.concat(pyseer_analysis_lst)
        PYSEER['ALL_SIMS']='ALL_SIMS'
        return PYSEER
    
    def pyseer_threshold_calcs(self, Pyseer_thresh, pyseer_lst, Pyseer_Select_lst):
        ratio_lst = []
        sim_lst = []
        thresh = []
        total_select = []
        total_correct = []
        for j in Pyseer_thresh:
            for i in range(0,len(pyseer_lst)):
                total_select_thresh=len(pyseer_lst[i][pyseer_lst[i]['lrt-pvalue']<=j])
                corect_select_thresh=len(Pyseer_Select_lst[i][Pyseer_Select_lst[i]['lrt-pvalue']<=j])
                ratio = corect_select_thresh/total_select_thresh
                sim = pyseer_lst[i].Simulation[0]
                total_select.append(total_select_thresh)
                total_correct.append(corect_select_thresh)
                sim_lst.append(sim)
                ratio_lst.append(ratio)
                thresh.append(j)
        
            DATA = pd.DataFrame({'SIM':sim_lst, 'THRESH':thresh, 'NUM_Correct':total_correct, 'NUM_Incorrect':total_select, 'RATIO':ratio_lst})
            DATA['ALL_SIMS'] = 'Over 10 Simulations'
        Means_DF =pd.DataFrame(DATA.groupby('THRESH')['RATIO'].mean()).rename(columns={'RATIO':'MEANS'})
        Vars_DF = pd.DataFrame(DATA.groupby('THRESH')['RATIO'].var()).rename(columns={'RATIO':'VARIANCES'})
        
        STATS_DF = Means_DF.join(Vars_DF)
        THRESH_DF = STATS_DF.join(DATA.set_index('THRESH')).reset_index()
        THRESH_DF['MeanRound'] = THRESH_DF['MEANS'].round(4)
        THRESH_DF['VARRound'] = "+/-" + THRESH_DF['VARIANCES'].round(5).astype(str)
        return THRESH_DF
    
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
        
        # Separated_Boxes=(ggplot(PYSEER) #SHRUNK
        #     + aes(x='Simulation', y='RATIO') #, fill='test'
        #     # + scale_fill_manual(values=color_dict)
        #     + geom_boxplot() #show_legend=False
        #     + labs(title='Pyseer Analysis Plot of Every Arbitrary p-value Threshold by All Simulations', x='Separated by Simulation Experiment', y='#Correct/#Selected') #
        #     # + theme(axis_text_x=element_blank())
        #     + theme(axis_text_x = element_text(angle = 45, hjust =1))
        # )
        # Overall_Box=(ggplot(PYSEER) #SHRUNK
        #     + aes(x='ALL_SIMS', y='RATIO') #, fill='test'
        #     + geom_boxplot() #show_legend=False
        #     + labs(title='Pyseer Analysis Plot of Every Arbitrary p-value Threshold', x='Over 10 Simulations', y='Proportion Correct Selected') 
        #     # + theme(axis_text_x = element_text(angle = 45, hjust =1))
        #     + theme(axis_text_x=element_blank())
        # )
        
        # Thresh_pyseer=(ggplot(THRESH_pyseer_DF) #SHRUNK
        #     + aes(x='factor(THRESH)', y='RATIO') #, fill='factor(THRESH)'
        #     + geom_boxplot() #show_legend=False
        #     + labs(title='Pyseer Significant Threshold Analysis', x='Thresholds of lrt-pvalue', y='Proportion Correct Selected') #, fill='Threshold'
        #     + theme(axis_text_x = element_text(angle = 45, hjust =1))
        #     + geom_text(aes(x='factor(THRESH)', y = 0.02, label='factor(MeanRound)'), size=9) # , color='factor(THRESH)' inside aes
        #     + geom_text(aes(x='factor(THRESH)', y = 0.015, label='factor(VARRound)'), size=8) # color='factor(THRESH)')
        #     + ylim(0.0, 0.2)
        #     # + annotate('text', x='factor(THRESH)', y = 'MEANS', label='MEANS')
        #     )
        # COUNT=THRESH_pyseer_DF[['NUM_Correct', 'NUM_Incorrect']].stack().reset_index().set_index('level_0')
        # COUNT = COUNT.join(THRESH_pyseer_DF).rename(columns={'level_1':'Group',0:'Count'})
        # COUNT['SIM']=COUNT.SIM.replace('Sim','', regex=True)
        # COUNT['SIM']=COUNT.SIM.astype(int)
        # Count2 = COUNT.sort_values(['Group', 'Count'], ascending=[True, False])
        # SimOrder = Count2.sort_values('SIM')['SIM'].tolist()
        # Count2['Simulation'] = pd.Categorical(Count2['SIM'], categories=pd.unique(SimOrder), ordered=True)
        # Count2['Group'] = Count2['Group'].replace('NUM_', '', regex=True)
        
        # Stacked_Bar=(ggplot(Count2)
        #     + aes(x='Simulation', y='Count', fill='Group', position_stack=False)
        #     + geom_col(stat='identity', position='stack') #position_dodge()
        #     + geom_text(aes(label='Count'), position=position_stack(vjust=0.5), size=9, va='center')
        #     + theme(axis_text_x = element_text(angle = 45, hjust =1))
        #     + labs(title='Number of Correct Selected Genes from Pyseer', y='Number of Genes Selected', x='Simulated Repetition', fill='Genes Selected')
        #     + facet_wrap('THRESH')
        #     + ylim(0, 1250)
        #     )
        
        # return Separated_Boxes, Overall_Box, Facet_scatter , Thresh_pyseer, Stacked_Bar
# =============================================================================
# MAYBE BREAK INTO TWO CLASSES
# =============================================================================
class bourtaAnalysis(Simulation_Analysis):
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

    # def boruta_select(self, Boruta_lst_DF, cntrl_lst_DF):
    #     Boruta_Select_lst = []
    #     for i in range(len(Boruta_lst_DF)):
    #         Selected = cntrl_lst_DF[i].join(Boruta_lst_DF[i], rsuffix='_Boruta').sort_values('PERCENT_SELECT', ascending=False) # lsuffix='_cntrl'
    #         Boruta_Select_lst.append(Selected)
    #     return Boruta_Select_lst
    
    def boruta_select(self, Boruta_lst_DF, cntrl_lst_DF):
        Boruta_Select_lst = []
        for i in range(len(Boruta_lst_DF)):
            Selected = cntrl_lst_DF[i].join(Boruta_lst_DF[i], how='inner', rsuffix='_Boruta').sort_values('PERCENT_SELECT', ascending=False) # lsuffix='_cntrl'
            Boruta_Select_lst.append(Selected)
        return Boruta_Select_lst
        
    def Boruta_calcs(self, Boruta_SELECT_lst, Boruta_lst_DF):
        nuisance=[]
        right = []
        thresholds = []
        Boruta_analysis_lst = []
        sim_lst = []
        for i in range(0, len(Boruta_SELECT_lst)):
            for j in range(0, len(Boruta_SELECT_lst[i])):
                thresh = Boruta_SELECT_lst[i].iloc[j,-2]
                sim = Boruta_SELECT_lst[i].iloc[j,-1]
                sim_lst.append(sim)
                thresholds.append(thresh)
                nuisance.append(len(Boruta_lst_DF[i][Boruta_lst_DF[i]['PERCENT_SELECT']>=thresh]))
                right.append(len(Boruta_SELECT_lst[i][Boruta_SELECT_lst[i]['PERCENT_SELECT']>=thresh]))
                BorutaDF = pd.DataFrame({'NumCorrect':right, 'NumSelected':nuisance, 'THRESHOLD':thresholds, 'Simulation':sim_lst})
                BorutaDF['RATIO'] = BorutaDF['NumCorrect']/BorutaDF['NumSelected'] # Make to FPR = FP/FP+TN
                BorutaDF['FPR'] = BorutaDF['NumSelected']/(BorutaDF['NumCorrect']+BorutaDF['NumSelected']) # type1 error
                BorutaDF = BorutaDF.sort_values('RATIO', ascending=False)
                Boruta_analysis_lst.append(BorutaDF)
        BORUTA = pd.concat(Boruta_analysis_lst)
        BORUTA['ALL_SIMS']='ALL_SIMS'
        return BORUTA
    
    
    def Boruta_Threshold_calcs(self, Bor_thresh, Boruta_lst_100, Boruta_SELECT_lst_100):
        ratio_lst = []
        sim_lst = []
        thresh = []
        total_select = []
        total_correct = []
        for j in Bor_thresh:
            for i in range(0,len(Boruta_lst_100)):
                total_select_thresh=len(Boruta_lst_100[i][Boruta_lst_100[i].PERCENT_SELECT>=j])
                corect_select_thresh=len(Boruta_SELECT_lst_100[i][Boruta_SELECT_lst_100[i].PERCENT_SELECT>=j])
                ratio = corect_select_thresh/total_select_thresh
                sim = Boruta_lst_100[i].Simulation[0]
                total_select.append(total_select_thresh)
                total_correct.append(corect_select_thresh)
                sim_lst.append(sim)
                ratio_lst.append(ratio)
                thresh.append(j)
        
        DATA = pd.DataFrame({'SIM':sim_lst, 'THRESH':thresh, 'NUM_Correct':total_correct, 'NUM_Incorrect':total_select,'RATIO':ratio_lst})
        DATA['ALL_SIMS'] = 'Over 10 Simulations'
        
        Means_DF =pd.DataFrame(DATA.groupby('THRESH')['RATIO'].mean()).rename(columns={'RATIO':'MEANS'})
        Vars_DF = pd.DataFrame(DATA.groupby('THRESH')['RATIO'].var()).rename(columns={'RATIO':'VARIANCES'})
        
        STATS_DF = Means_DF.join(Vars_DF)
        THRESH_DF = STATS_DF.join(DATA.set_index('THRESH')).reset_index()
        THRESH_DF['MeanRound'] = THRESH_DF['MEANS'].round(4)
        THRESH_DF['VARRound'] = "+/-" + THRESH_DF['VARIANCES'].round(4).astype(str)
        return THRESH_DF  
    
    
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



# =============================================================================
# LR SELECTION
# =============================================================================
class PenalizedRegressAnalysis(Simulation_Analysis):
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

    # def boruta_select(self, Boruta_lst_DF, cntrl_lst_DF):
    #     Boruta_Select_lst = []
    #     for i in range(len(Boruta_lst_DF)):
    #         Selected = cntrl_lst_DF[i].join(Boruta_lst_DF[i], rsuffix='_Boruta').sort_values('PERCENT_SELECT', ascending=False) # lsuffix='_cntrl'
    #         Boruta_Select_lst.append(Selected)
    #     return Boruta_Select_lst
    
    def Regress_select(self, Regress_lst, Control_lst):
        Regress_Select_lst = []
        for i in range(len(Regress_lst)):
            Selected = Control_lst[i].join(Regress_lst[i], how='inner', rsuffix='_Regress').sort_values('RelativePercentStength', ascending=False) # lsuffix='_cntrl'
            Regress_Select_lst.append(Selected)
        return Regress_Select_lst
        
    # def Regress_calcs(self, Regress_Select_lst, Regress_lst):
    #     nuisance=[]
    #     right = []
    #     thresholds = []
    #     Boruta_analysis_lst = []
    #     sim_lst = []
    #     for i in range(0, len(Regress_Select_lst)):
    #         for j in range(0, len(Regress_Select_lst[i])):
    #             thresh = Regress_Select_lst[i].iloc[j,-2]
    #             sim = Regress_Select_lst[i].iloc[j,-1]
    #             sim_lst.append(sim)
    #             thresholds.append(thresh)
    #             nuisance.append(len(Regress_lst[i][Regress_lst[i]['PERCENT_SELECT']>=thresh]))
    #             right.append(len(Regress_Select_lst[i][Regress_Select_lst[i]['PERCENT_SELECT']>=thresh]))
    #             BorutaDF = pd.DataFrame({'NumCorrect':right, 'NumSelected':nuisance, 'THRESHOLD':thresholds, 'Simulation':sim_lst})
    #             BorutaDF['RATIO'] = BorutaDF['NumCorrect']/BorutaDF['NumSelected'] # Make to FPR = FP/FP+TN
    #             BorutaDF['FPR'] = BorutaDF['NumSelected']/(BorutaDF['NumCorrect']+BorutaDF['NumSelected']) # type1 error
    #             BorutaDF = BorutaDF.sort_values('RATIO', ascending=False)
    #             Boruta_analysis_lst.append(BorutaDF)
    #     BORUTA = pd.concat(Boruta_analysis_lst)
    #     BORUTA['ALL_SIMS']='ALL_SIMS'
    #     return BORUTA
    
    
    # def Regress_Threshold_calcs(self, Regress_thresh, Regress_lst, Regress_Select_lst):
    #     ratio_lst = []
    #     sim_lst = []
    #     thresh = []
    #     total_select = []
    #     total_correct = []
    #     for j in Regress_thresh:
    #         for i in range(0,len(Regress_lst)):
    #             total_select_thresh=len(Regress_lst[i][Regress_lst[i].PERCENT_SELECT>=j])
    #             corect_select_thresh=len(Regress_Select_lst[i][Regress_Select_lst[i].PERCENT_SELECT>=j])
    #             ratio = corect_select_thresh/total_select_thresh
    #             sim = Regress_lst[i].Simulation[0]
    #             total_select.append(total_select_thresh)
    #             total_correct.append(corect_select_thresh)
    #             sim_lst.append(sim)
    #             ratio_lst.append(ratio)
    #             thresh.append(j)
        
    #     DATA = pd.DataFrame({'SIM':sim_lst, 'THRESH':thresh, 'NUM_Correct':total_correct, 'NUM_Incorrect':total_select,'RATIO':ratio_lst})
    #     DATA['ALL_SIMS'] = 'Over 10 Simulations'
        
    #     Means_DF =pd.DataFrame(DATA.groupby('THRESH')['RATIO'].mean()).rename(columns={'RATIO':'MEANS'})
    #     Vars_DF = pd.DataFrame(DATA.groupby('THRESH')['RATIO'].var()).rename(columns={'RATIO':'VARIANCES'})
        
    #     STATS_DF = Means_DF.join(Vars_DF)
    #     THRESH_DF = STATS_DF.join(DATA.set_index('THRESH')).reset_index()
    #     THRESH_DF['MeanRound'] = THRESH_DF['MEANS'].round(4)
    #     THRESH_DF['VARRound'] = "+/-" + THRESH_DF['VARIANCES'].round(4).astype(str)
    #     return THRESH_DF  
    
    
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


        
    # def boruta_figs(self, BORUTA100, THRESH_Bor_DF):
    #     Separated_Boxes=(ggplot(BORUTA100) #SHRUNK
    #         + aes(x='Simulation', y='RATIO') #, fill='test'
    #         # + scale_fill_manual(values=color_dict)
    #         + geom_boxplot() #show_legend=False
    #         + labs(title='', x='') #
    #         # + theme(axis_text_x=element_blank())
    #         + theme(axis_text_x = element_text(angle = 45, hjust =1))
    #         + labs(title='Boruta Analysis Plot of Every Arbitrary Selection% Threshold by All Simulations')
    #     )
    #     Overall_Box=(ggplot(BORUTA100) #SHRUNK
    #         + aes(x='ALL_SIMS', y='RATIO') #, fill='test'
    #         + geom_boxplot() #show_legend=False
    #         + labs(title='Boruta Analysis Plot of Every Arbitrary Selection% Threshold', x='Over 10 Simulations', y='Proportion Correct Selected') 
    #         # + theme(axis_text_x = element_text(angle = 45, hjust =1))
    #         + theme(axis_text_x=element_blank())
    #     )
    #     Thresh_Bor=(ggplot(THRESH_Bor_DF) #SHRUNK
    #         + aes(x='factor(THRESH)', y='RATIO') #, fill='factor(THRESH)'
    #         + geom_boxplot() #show_legend=False
    #         + labs(title='Boruta Significant Threshold Analysis', x='Threshold of %Selected per 250 Iterations', y='Proportion Correct Selected') #, fill='Threshold'
    #         + theme(axis_text_x = element_text(angle = 45, hjust =1))
    #         + geom_text(aes(x='factor(THRESH)', y = 0.0, label='factor(MeanRound)'), size=10) # , color='factor(THRESH)' inside aes
    #         + geom_text(aes(x='factor(THRESH)', y = -0.025, label='factor(VARRound)'), size=9) # color='factor(THRESH)')
    #         + ylim(-0.05, 0.25)
    #         # + annotate('text', x='factor(THRESH)', y = 'MEANS', label='MEANS')
    #     )
        
        
    #     COUNT=THRESH_Bor_DF[['NUM_Correct', 'NUM_Incorrect']].stack().reset_index().set_index('level_0')
    #     COUNT = COUNT.join(THRESH_Bor_DF).rename(columns={'level_1':'Group',0:'Count'})
    #     COUNT['SIM']=COUNT.SIM.replace('Sim','', regex=True)
    #     COUNT['SIM']=COUNT.SIM.astype(int)
    #     Count2 = COUNT.sort_values(['Group', 'Count'], ascending=[True, False])
    #     SimOrder = Count2.sort_values('SIM')['SIM'].tolist()
    #     Count2['Simulation'] = pd.Categorical(Count2['SIM'], categories=pd.unique(SimOrder), ordered=True)
    #     Count2['Group'] = Count2['Group'].replace('NUM_', '', regex=True)
    #     Stacked_Bar=(ggplot(Count2)
    #         + aes(x='Simulation', y='Count', fill='Group', position_stack=False)
    #         + geom_col(stat='identity', position='stack') #position_dodge()
    #         + geom_text(aes(label='Count'), position=position_stack(vjust=0.5), size=9, va='center')
    #         + theme(axis_text_x = element_text(angle = 45, hjust =1))
    #         + labs(title='Number of Correct Selected Genes from Boruta', y='Number of Genes Selected', x='Simulated Repetition', fill='Genes Selected')
    #         + facet_wrap('THRESH')
    #         + ylim(0, 1500)
    #     )
        

    #     Count3 = pd.DataFrame(Count2.groupby(['Group', 'THRESH'])['Count'].mean()).rename(columns={'Count':'Mean'})
    #     Count4 = pd.DataFrame(Count2.groupby(['Group', 'THRESH'])['Count'].std()).rename(columns={'Count':'StDeviation'})
        
    #     Avgs = Count3.join(Count4).reset_index()
    #     AverageCounts=(ggplot(Avgs)
    #         + aes(x='factor(THRESH)', y='Mean', fill='Group', position_stack=False)
    #         + geom_col(stat='identity', position='dodge') #position_dodge()
    #         # + geom_text(aes(label='Count'), position='stack', size=9, va='center')
    #         # + geom_text(aes(label='Count'), position=position_stack(vjust=0.5), size=9, va='center')
    #         + geom_text(aes(label='Mean'), position=position_dodge(width=0.9), size=8, va='bottom')
    #         + geom_errorbar(aes(x="factor(THRESH)", ymin="Mean-StDeviation",ymax="Mean+StDeviation"), width=0.1, position=position_dodge(0.9))
    #         + theme(axis_text_x = element_text(angle = 45, hjust =1))
    #         + labs(title='Average Correct and Incorrect Gene Selections', y='Average Number of Genes Selected', x='Threshold', fill='Genes Selected')
    #         # + facet_wrap('THRESH')
    # #         # + ylim(0, 250)
    # #     )
        
            
    #     return Separated_Boxes, Overall_Box, Thresh_Bor, Stacked_Bar, AverageCounts
class FigGen(Simulation_Analysis):
    def confusion_forPlot(self, Confusion):
        Means_DF =pd.DataFrame(Confusion.groupby('Threshold')['PosPredVal'].mean()).rename(columns={'PosPredVal':'PosPredVal_mean'})
        Vars_DF = pd.DataFrame(Confusion.groupby('Threshold')['PosPredVal'].var()).rename(columns={'PosPredVal':'PosPredVal_variance'})
        
        STATS_DF = Means_DF.join(Vars_DF)
        THRESH_DF = STATS_DF.join(Confusion.set_index('Threshold')).reset_index()
        THRESH_DF['PosPredVal_meanRound'] = THRESH_DF['PosPredVal_mean'].round(4)
        THRESH_DF['PosPredVal_varianceRound'] = "+/-" + THRESH_DF['PosPredVal_variance'].round(5).astype(str)
        return THRESH_DF
    
    def Comparison_box2(self, BorutaConf_forBox, PyseerConf_forBox, RegressConf_forBox, expname, BorIter):
        PyseerConf_forBox['Method'] = 'PYSEER'
        BorutaConf_forBox['Method'] = 'BORUTA'
        RegressConf_forBox['Method'] = 'LOGIT'
        Thresh_concat = pd.concat([PyseerConf_forBox, BorutaConf_forBox, RegressConf_forBox ])
        # Thresh_concat = pd.concat([PyseerConf_forBox, BorutaConf_forBox])
        Threshold_Box_comparison=(ggplot(Thresh_concat) #SHRUNK
            + aes(x='factor(Threshold)', y='PosPredVal') #, fill='factor(Threshold)'
            + geom_boxplot() #show_legend=False
            + labs(title=' Boruta vs. Pyseer Significant Threshold Analysis', x='Boruta Threshold of %Selected per  Iterations \nLogit Threshold of %Varibale Strength where the stongest variable is set to 100% \nPyseer Thresholds of lrt-pvalue \n(LOGIT = l1 penalized logistic regression)', y='Positive Predictive Value (n=10)') #, fill='Threshold'
            + theme(axis_text_x = element_text(angle = 45, hjust =1))
            # + geom_text(aes(y = 1.15, label="Means +/- Variances"), size=10)
            + geom_text(aes(x='factor(Threshold)', y = 1.1, label='factor(PosPredVal_meanRound)'), size=7, angle=45) #  min(Thresh_concat.PosPredVal) - 0.001, color='factor(Threshold)' inside aes
            + geom_text(aes(x='factor(Threshold)', y = 1.05, label='factor(PosPredVal_varianceRound)'), size=5, angle=45) # min(Thresh_concat.PosPredVal)- 0.005, color='factor(Threshold)')
            # + ylim(0.0, 0.2)
            + facet_wrap('Method', scales='free_x')
            + theme(subplots_adjust={'wspace': 0.25})
            # + annotate('text', x='factor(Threshold)', y = 'MEANS', label='MEANS')
        )
        # Threshold_Box_comparison=(ggplot(Thresh_concat) #SHRUNK
        #     + aes(x='factor(Threshold)', y='PosPredVal') #, fill='factor(Threshold)'
        #     + geom_boxplot() #show_legend=False
        #     + labs(title=str(expname)+' Boruta vs. Pyseer Significant Threshold Analysis', x='Boruta Threshold of %Selected per '+str(BorIter)+' Iterations & Pyseer Thresholds of lrt-pvalue', y='Positive Predictive Value (n=10)') #, fill='Threshold'
        #     + theme(axis_text_x = element_text(angle = 45, hjust =1))
        #     + geom_text(aes(x='factor(Threshold)', y = min(Thresh_concat.PosPredVal) - 0.001, label='factor(PosPredVal_meanRound)'), size=7) # , color='factor(Threshold)' inside aes
        #     + geom_text(aes(x='factor(Threshold)', y = min(Thresh_concat.PosPredVal)- 0.005, label='factor(PosPredVal_varianceRound)'), size=5) # color='factor(Threshold)')
        #     # + ylim(0.0, 0.2)
        #     + facet_wrap('Method', scales='free_x')
        #     # + annotate('text', x='factor(Threshold)', y = 'MEANS', label='MEANS')
        # )
        return Threshold_Box_comparison
    
    # def Comparison_box(self, THRESH_pyseer_DF, THRESH_Bor_DF, expname, BorIter):
    #     THRESH_pyseer_DF['Method'] = 'PYSEER'
    #     THRESH_Bor_DF['Method'] = 'BORUTA'
    #     Thresh_concat = pd.concat([THRESH_pyseer_DF, THRESH_Bor_DF])
                   
    #     Threshold_Box_comparison=(ggplot(Thresh_concat) #SHRUNK
    #         + aes(x='factor(THRESH)', y='RATIO') #, fill='factor(THRESH)'
    #         + geom_boxplot() #show_legend=False
    #         + labs(title=str(expname)+' Boruta vs. Pyseer Significant Threshold Analysis', x='Boruta Threshold of %Selected per '+str(BorIter)+' Iterations & Pyseer Thresholds of lrt-pvalue', y='Positive Predictive Value (n=10)') #, fill='Threshold'
    #         + theme(axis_text_x = element_text(angle = 45, hjust =1))
    #         + geom_text(aes(x='factor(THRESH)', y = 0.02, label='factor(MeanRound)'), size=7) # , color='factor(THRESH)' inside aes
    #         + geom_text(aes(x='factor(THRESH)', y = 0.015, label='factor(VARRound)'), size=5) # color='factor(THRESH)')
    #         # + ylim(0.0, 0.2)
    #         + facet_wrap('Method', scales='free_x')
    #         # + annotate('text', x='factor(THRESH)', y = 'MEANS', label='MEANS')
    #     )
    #     return Threshold_Box_comparison
    
    def ThresholdStats_pyseer(self, THRESH_pyseer_DF):
        COUNT=THRESH_pyseer_DF[['NUM_Correct', 'NUM_Incorrect']].stack().reset_index().set_index('level_0')
        COUNT = COUNT.join(THRESH_pyseer_DF).rename(columns={'level_1':'Group',0:'Count'})
        COUNT['SIM']=COUNT.SIM.replace('Sim','', regex=True)
        COUNT['SIM']=COUNT.SIM.astype(int)
        Count2 = COUNT.sort_values(['Group', 'Count'], ascending=[True, False])
        SimOrder = Count2.sort_values('SIM')['SIM'].tolist()
        Count2['Simulation'] = pd.Categorical(Count2['SIM'], categories=pd.unique(SimOrder), ordered=True)
        Count2['Group'] = Count2['Group'].replace('NUM_', '', regex=True)
        Count3 = pd.DataFrame(Count2.groupby(['Group', 'THRESH'])['Count'].mean()).rename(columns={'Count':'Mean'})
        Count4 = pd.DataFrame(Count2.groupby(['Group', 'THRESH'])['Count'].std()).rename(columns={'Count':'StDeviation'})
        
        Avgs_pyseer = Count3.join(Count4).reset_index()
        Avgs_pyseer['Method'] = 'PYSEER'
        return Avgs_pyseer
    
    def ThresholdStats_boruta(self, THRESH_Bor_DF, percent=100, alpha=0.05):
        COUNT=THRESH_Bor_DF[['NUM_Correct', 'NUM_Incorrect']].stack().reset_index().set_index('level_0')
        COUNT = COUNT.join(THRESH_Bor_DF).rename(columns={'level_1':'Group',0:'Count'})
        COUNT['SIM']=COUNT.SIM.replace('Sim','', regex=True)
        COUNT['SIM']=COUNT.SIM.astype(int)
        Count2 = COUNT.sort_values(['Group', 'Count'], ascending=[True, False])
        SimOrder = Count2.sort_values('SIM')['SIM'].tolist()
        Count2['Simulation'] = pd.Categorical(Count2['SIM'], categories=pd.unique(SimOrder), ordered=True)
        Count2['Group'] = Count2['Group'].replace('NUM_', '', regex=True)
        Count3 = pd.DataFrame(Count2.groupby(['Group', 'THRESH'])['Count'].mean()).rename(columns={'Count':'Mean'})
        Count4 = pd.DataFrame(Count2.groupby(['Group', 'THRESH'])['Count'].std()).rename(columns={'Count':'StDeviation'})
        
        Avgs_boruta = Count3.join(Count4).reset_index()
        Avgs_boruta['Method'] = 'BORUTA'
        Avgs_boruta['Percent'] = int(percent)
        Avgs_boruta['alpha'] = float(alpha)
        return Avgs_boruta
    
    def ThresholdStats_regress(self, Regress_Confusion):
        COUNT=Regress_Confusion[['TruePos', 'FalsePos']].stack().reset_index().set_index('level_0')
        COUNT = COUNT.join(THRESH_Bor_DF).rename(columns={'level_1':'Group',0:'Count'})
        COUNT['SIM']=COUNT.SIM.replace('Sim','', regex=True)
        COUNT['SIM']=COUNT.SIM.astype(int)
        Count2 = COUNT.sort_values(['Group', 'Count'], ascending=[True, False])
        SimOrder = Count2.sort_values('SIM')['SIM'].tolist()
        Count2['Simulation'] = pd.Categorical(Count2['SIM'], categories=pd.unique(SimOrder), ordered=True)
        Count2['Group'] = Count2['Group'].replace('NUM_', '', regex=True)
        Count3 = pd.DataFrame(Count2.groupby(['Group', 'THRESH'])['Count'].mean()).rename(columns={'Count':'Mean'})
        Count4 = pd.DataFrame(Count2.groupby(['Group', 'THRESH'])['Count'].std()).rename(columns={'Count':'StDeviation'})
        
        Avgs_regress = Count3.join(Count4).reset_index()
        Avgs_regress['Method'] = 'LOGIT'
        return Avgs_regress
    
    def Comparison_BarCounts(self, Avgs_pyseer, Avgs_boruta, expname):

        CountConcat = pd.concat([Avgs_pyseer, Avgs_boruta])
        
        AverageCounts=(ggplot(CountConcat)
            + aes(x='factor(THRESH)', y='Mean', fill='Group', position_stack=False)
            + geom_col(stat='identity', position='dodge') #position_dodge()
            # + geom_text(aes(label='Count'), position='stack', size=9, va='center')
            # + geom_text(aes(label='Count'), position=position_stack(vjust=0.5), size=9, va='center')
            + geom_text(aes(label='Mean'), position=position_dodge(width=0.9), size=8, va='bottom')
            + geom_errorbar(aes(x="factor(THRESH)", ymin="Mean-StDeviation",ymax="Mean+StDeviation"), width=0.1, position=position_dodge(0.9))
            + theme(axis_text_x = element_text(angle = 45, hjust =1))
            + labs(title=str(expname)+' Average Correct and Incorrect Gene Selections', y='Average Number of Genes Selected (n=10)', x='Boruta Threshold of %Selected per 250 Iterations & Pyseer Thresholds of lrt-pvalue', fill='Genes Selected')
            # + facet_wrap('THRESH')
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
    
    # def ConfusionTableExport(self, MeanBorConfusion, MeanPysConfusion, tablepath, filename):
    #     CONFUSION = pd.concat([MeanBorConfusion,MeanPysConfusion])

    #     cols = CONFUSION.columns.tolist()
    #     cols = cols[-1:] + cols[:-1]
    #     cols = cols[:8]
    #     CONFUSION=CONFUSION[cols]
    #     CONFUSION=CONFUSION.round({'TruePos_mean':1, 'FalsePos_mean':1, 'TrueNeg_mean':0, 'FalseNeg_mean':1,'PosPredVal_mean':3, 'TPR_mean':3, 'TNR_mean':3})
    #     # , 'TruePos_STD':2, 'FalsePos_STD':2, 'TrueNeg_STD':2, 'FalseNeg_STD':2, 'PosPredVal_STD':2,'TPR_STD':2, 'TNR_STD':2
    #     CONFUSION.reset_index(inplace=True)
    #     CONFUSION.loc[-1] = '-----'
    #     CONFUSION.index = CONFUSION.index+1
    #     CONFUSION = CONFUSION.sort_index()
    #     CONFUSION[CONFUSION.columns[0]] = '|' + CONFUSION[CONFUSION.columns[0]].astype(str)
    #     CONFUSION[CONFUSION.columns[-1]] = CONFUSION[CONFUSION.columns[-1]].astype(str) + '|'
    #     cols = CONFUSION.columns.tolist()
    #     cols[0] = '|' + cols[0]
    #     cols[-1] = cols[-1] + '|'
    #     CONFUSION.columns=cols
    #     os.chdir(tablepath)
    #     CONFUSION.to_csv(filename,sep='|', escapechar=' ', index=False, quoting=csv.QUOTE_NONE)
    
    def ConfusionTableExport(self, BorConfusion, PysConfusion):
        CONFUSION = pd.concat([BorConfusion, PysConfusion])
        CONFUSION=CONFUSION.round({'TruePos_mean':1, 'FalsePos_mean':1, 'TrueNeg_mean':1, 'FalseNeg_mean':1,'PosPredVal_mean':3, 'TPR_mean':3, 'TNR_mean':3})
        print(CONFUSION.to_markdown())
        return CONFUSION
    
    
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