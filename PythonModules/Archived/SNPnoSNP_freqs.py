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
    def __init__(self, annotated_nasp, phenotype=None): # needs to be annotated nasp for SNPnoSNP generation, but not for snpFreqs
        self.annotated_nasp = pd.read_csv(annotated_nasp, sep='\t')
        if phenotype != None:
            self.phenotype = pd.read_csv(phenotype, sep='\t')#.rename(columns={0:'STRAIN',1:'DRUGconc'})#.set_index('STRAIN') #.set_index('STRAIN', inplace=True) , header=None
            self.phenotype.columns=['STRAIN', 'DRUGconc']
            self.phenotype.set_index('STRAIN',inplace=True)
    def IndexOrg(self,annotated_nasp):
        naspT=annotated_nasp.T # Transpose
        naspT.columns=naspT.iloc[0] #Make column names the first row of the DF
        # naspT=naspT.iloc[2:] #Remove first two rows (column names and reference)
        naspT = naspT.drop(['LocusID', 'Reference']) #probably best to drop by index name
        naspT.index.rename('STRAIN', inplace=True)
        return naspT
    def GatherPheno(self, phenotypeFile, phenoSelect): # add phenotype value variable (self, phenotypeFile, phenotypeSelect)
        drug_sus = phenotypeFile.loc[phenotypeFile['DRUGconc']==phenoSelect]#.set_index('STRAIN')
        return drug_sus
    def SNPfreqs(self, drug_sus, naspT):
        drug_susNASP=drug_sus.join(naspT)
        drug_susNASP = drug_susNASP.drop(columns=['DRUGconc']) #['index', 'DRUGconc', 'logical']
        drug_susNASP= drug_susNASP.apply(pd.value_counts)
        drug_susNASP=(100. * drug_susNASP/drug_susNASP.sum()).round(2)
        return drug_susNASP
    def susSNPcalls(self, drug_susNASP):
        df = pd.DataFrame(drug_susNASP.idxmax())
        df.columns=['susSNPcall']
        df.index.rename('LocusID', inplace=True)
        return df
    def SNPcalls(self, annotated_nasp, PanSusVect):
        nasp = annotated_nasp.loc[:, :'#SNPcall'].iloc[:,:-1]
        nasp['ncbi_id']= annotated_nasp['ncbi_id'].copy()
        nasp=nasp.set_index('LocusID')
        naspPanSus=PanSusVect.join(nasp).reset_index()
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