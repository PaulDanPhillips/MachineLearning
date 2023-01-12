#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
from pathlib import Path
from math import isnan
import pickle
import argparse

class SNPnoSNP:
    # needs to be annotated nasp for SNPnoSNP generation, but not for snpFreqs
    # needs to be phenotype file without a header for below constructor
    def __init__(self, annotated_nasp, phenotype): # needs to be annotated nasp for SNPnoSNP generation, but not for snpFreqs
        self.annotated_nasp = pd.read_csv(annotated_nasp, sep='\t')
        self.phenotype = pd.read_csv(phenotype, sep='\t', header=None).rename(columns={0:'STRAIN',1:'Res_Sus'}).set_index('STRAIN') #.set_index('STRAIN', inplace=True)
    
    def IndexOrg(self,annotated_nasp):
        naspT=annotated_nasp.T # Transpose
        naspT.columns=naspT.iloc[0] #Make column names the first row of the DF
        # naspT=naspT.iloc[2:] #Remove first two rows (column names and reference)
        naspT = naspT.drop(['LocusID', 'Reference']) #probably best to drop by index name
        naspT.index.rename('STRAIN', inplace=True)
        return naspT
    def GatherPheno(self, phenotypeFile, phenoSelect): # add phenotype value variable (self, phenotypeFile, phenotypeSelect)
        # collects only phenotype of interest.
        # phenotype['logical'] = np.logical_or.reduce([phenotype['Res_Sus']!=1])
        # drug_sus=phenotype.loc[phenotype.logical,:]
        # drug_sus=drug_sus.reset_index()
        # drug_sus=drug_sus.set_index('STRAIN')
        # phenotypeFile.set_index('STRAIN', inplace=True)
        drug_sus = phenotypeFile.loc[phenotypeFile['Res_Sus']==phenoSelect]#.set_index('STRAIN')
        return drug_sus
    def SNPfreqs(self, drug_sus, naspT):
        drug_susNASP=drug_sus.join(naspT)
        drug_susNASP = drug_susNASP.drop(columns=['Res_Sus']) #['index', 'Res_Sus', 'logical']
        drug_susNASP= drug_susNASP.apply(pd.value_counts)
        drug_susNASP=(100. * drug_susNASP/drug_susNASP.sum()).round(2)
        return drug_susNASP
    def susSNPcalls(self, drug_susNASP):
        # df2 = df.copy() # Make a copy because otherwise this rewrites the original df
        # for col in df.columns: # sort each column by max. Do I need for loop?
        #     df[col] = df[col].idxmax()
        # df = df.drop_duplicates(keep='first').T
        # df.columns=['susSNPcall']
        # df.index.rename('LocusID', inplace=True)
        df = pd.DataFrame(drug_susNASP.idxmax())
        df.columns=['susSNPcall']
        df.index.rename('LocusID', inplace=True)
        return df
    def SNPcalls(self, annotated_nasp, PanSusVect):
        nasp = annotated_nasp.loc[:, :'#SNPcall'].iloc[:,:-1]
        nasp['ncbi_id']= annotated_nasp['ncbi_id'].copy()
        nasp=nasp.set_index('LocusID')
        naspPanSus=PanSusVect.join(nasp).reset_index()
        # nasp['Loci'] = nasp['ncbi_id'] +'_' + nasp['LocusID']
        NASP2=pd.DataFrame()
        for i in range(0,len(naspPanSus.columns)):
            SNPnoSNP = np.where(naspPanSus.iloc[:,i]==naspPanSus['susSNPcall'],0,1)
            NASP2[i]=SNPnoSNP.tolist()
        NASP2.columns = naspPanSus.columns
        return NASP2
    def SNPcount(self, NASP2, annotated_nasp):
        NASP2['ncbi_id']= annotated_nasp['ncbi_id']
        cols = list(NASP2)
        cols.insert(0, cols.pop(cols.index('ncbi_id')))
        NASP2 = NASP2.loc[:, cols]
        NASP2_sum=NASP2.groupby('ncbi_id').sum()
        NASP2_sum=NASP2_sum.drop(columns=['LocusID', 'susSNPcall', 'Reference'])
        SNPcount=NASP2_sum.T
        return SNPcount