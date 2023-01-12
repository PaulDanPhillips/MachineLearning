import pandas as pd
import random
import numpy as np
import os

# os.chdir('H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAS_FileFormatting')
# # LSBSR = pd.read_csv('LSBSR_indexed_RAW.txt', sep='\t')
# LSBSR = pd.read_csv('LSBSR_Balanced_SIM1.txt', sep='\t', index_col=0)
# # Plink_resp = pd.read_csv('Plink_response.txt', sep='\t')
# # plink = pd.read_csv('Plink_predictor.tsv', sep='\t', nrows=400)
# # scoary = pd.read_csv('Scoary_predictor.csv')
# # scoary_resp = pd.read_csv('Scoary_response.csv')
# BacGWASim = pd.read_pickle('BacGWASim_predictor.pickle')
# BacGWASim_resp = pd.read_csv('BacGWASim_response.phen', header=None, sep='\s+')
# # AMC = pd.read_csv('AMCnew.txt', sep='\t')
# # AMC['STRAIN2'] = AMC.STRAIN.str.split('_',expand=True)[0]
# # AMC.STRAIN2.replace('L1','',regex=True, inplace=True)
# # AMC.set_index('STRAIN2', inplace=True)

# treewas = pd.read_csv('TreeWAS_predictor.txt', sep='\t', nrows=100)
# treewas_response = pd.read_csv('TreeWAS_response.txt', sep='\t')
# # treewas_response['NEWid'] = treewas_response['0'].str.split('_', expand=True)[0]
# # treewas_response.NEWid.replace('L1', '', regex=True, inplace=True)
# # treewas_response.set_index('NEWid', inplace=True)

# SimResp = pd.read_csv('BalancedRandomResponse.txt', sep='\t', index_col=0)
# # for col in BacGWASim.columns:
#     # print(col[-3])
#     BacGWASim.loc[BacGWASim[col]==0, col]=col[-3]
#     BacGWASim.loc[BacGWASim[col]==1, col]=col[-1]
    
# BacGWASim.loc[BacGWASim['1:2:A:G']==0]
# BacGWASim2 = pd.read_csv('BacGWASim_predictor.vcf', skiprows=13, sep='\t')
# scoaryPred = pd.read_csv('H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAStutorials/Scoary/scoary/exampledata/Gene_presence_absence.csv')
# tmp=scoaryPred.columns
# tmp=tmp[1:14]
# tmp = tmp.to_list()
def BacGWASim_toScoary(features, resp, SimNum, path, filename, parentPath):
    os.chdir(parentPath)
    scoary_empty = ['Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Genome fragment', 'Order within fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']
    scoary=features.copy().T
    # scoary.set_index('Gene', inplace=True)
    for item in scoary_empty:
        scoary.insert(0,item,np.nan)
    scoary_resp = resp.copy()
    scoary_resp.columns = ['Name', 'Name_copy', 'BacGWASim_'+str(SimNum)]
    scoary_resp.set_index('Name', inplace=True)
    scoary_resp = scoary_resp[['BacGWASim_'+str(SimNum)]]
    scoary_resp.replace(1,0,regex=True, inplace=True)
    scoary_resp.replace(2,1,regex=True, inplace=True)

    if not os.path.isdir(str(path)):
        os.makedirs(str(path))
    os.chdir(str(path))
    
    scoary.reset_index().to_csv(filename+'_'+SimNum+'_FEATURES.csv')
    scoary_resp.reset_index().to_csv(filename+'_'+SimNum+'_RESPONSE.csv', index=False)
    # return #scoary.reset_index(), scoary_resp.reset_index()

def BacGWASim_toPlink(features, resp, SimNum, path, filename, parentPath):
    os.chdir(parentPath)
    # plink --no-fid --no-parents --no-sex --map3
    
    # plink_PED_required1 = ['Family ID', 'Individual ID', 'Paternal ID', 'Maternal ID', 'Sex (1=male; 2=female; other=unknown)', 'Phenotype']
    # plink_PED_required2 = ['Family ID', 'Paternal ID', 'Maternal ID', 'Sex (1=male; 2=female; other=unknown)']
    plink_PED = resp.copy()
    # plink_resp.insert(0,'Family ID', np.nan)
    plink_PED.columns = ['Individual ID', 'Individual ID_copy', 'Phenotype']
    plink_PED.drop('Individual ID_copy', axis=1, inplace=True)
    # for item in plink_PED_required2:
    #     plink_PED.insert(0, item, np.nan)
        
    Plink_PED_feats = features.copy()
    Plink_PED_feats.replace(1, '2 2', regex=True, inplace=True)
    Plink_PED_feats.replace(0, '1 1', regex=True, inplace=True)
    Plink_PED_feats.columns = ['']*len(Plink_PED_feats.columns)
    Plink_PED_feats.reset_index(inplace=True)
    Plink_PED_final = plink_PED.merge(Plink_PED_feats, left_on='Individual ID', right_on='index')
    Plink_PED_final.drop('index',axis =1, inplace=True)
    # Plink_feats.replace(1,2,regex=True, inplace=True)
    # Plink_feats.replace(0,1,regex=True, inplace=True)
    # tmp2= Plink_PED.merge(tmp, left_on='Individual ID', right_on='index')
    
    
    # plink_MAP_required = ['chromosome (1-22, X, Y or 0 if unplaced)', 'rs# or snp identifier', 'Base-pair position (Bp units']
    Plink_features = features.copy()
    Plink_MAP = pd.DataFrame({'rs# or snp identifier':Plink_features.columns.to_list()})
    Plink_MAP[['chromosome', 'Base-pair position (bp units)', 'RefAllele', 'SNPAllele']]= Plink_MAP['rs# or snp identifier'].str.split(':', expand=True)
    Plink_MAP.drop(['RefAllele', 'SNPAllele'], axis=1, inplace=True)
    cols = ['chromosome', 'rs# or snp identifier', 'Base-pair position (bp units)']
    Plink_MAP=Plink_MAP[cols]
    
    if not os.path.isdir(str(path)):
        os.makedirs(str(path))
    os.chdir(str(path))
    Plink_PED_final.to_csv(filename+'_'+SimNum+".ped", sep='\t', index=False, header=False)
    Plink_MAP.to_csv(filename+'_'+SimNum+'.map', sep='\t', index=False, header=False)
    
    # plink_resp.set_index(0, inplace=True)
    # plink_resp = plink_resp[[2]]
    # plink = features.copy()
    # for col in plink.columns:
    #     plink.loc[plink[col]==0, col]=col[-3]
    #     plink.loc[plink[col]==1, col]=col[-1]
    # plink_resp = resp.copy()
    # plink_resp.set_index(0, inplace=True)
    # plink_resp = plink_resp[[2]]
    # return #Plink_PED_final, Plink_MAP

def BacGWASim_toTreeWAS(features, resp, SimNum, path, filename, parentPath):
    os.chdir(parentPath)
    TreeWAS= features.copy()
    TreeWAS.reset_index()
    TreeWAS_resp = resp.copy()
    TreeWAS_resp.columns = ['Name', 'Name_copy', 'BacGWASim_'+str(SimNum)]
    TreeWAS_resp.set_index('Name', inplace=True)
    TreeWAS_resp = TreeWAS_resp[['BacGWASim_'+str(SimNum)]]
    TreeWAS_resp.replace(1,0,regex=True, inplace=True)
    TreeWAS_resp.replace(2,1,regex=True, inplace=True)
    
    if not os.path.isdir(str(path)):
        os.makedirs(str(path))
    os.chdir(str(path))
    TreeWAS.to_csv(filename+'_'+SimNum+"_FEATURES.txt", sep='\t')
    TreeWAS_resp.to_csv(filename+'_'+SimNum+'_RESPONSE.txt', sep='\t')
    # return #TreeWAS, TreeWAS_resp

def BacGWASim_toBugWAS(DF):
    pass

def BacGWASim_toFastLMM(features, resp, SimNum, path, filename, parentPath):
    os.chdir(parentPath)
    plink_PED = resp.copy()
    plink_PED.columns = ['Individual ID', 'Individual ID_copy', 'Phenotype']
    plink_PED.drop('Individual ID_copy', axis=1, inplace=True)

        
    Plink_PED_feats = features.copy()
    Plink_PED_feats.replace(1, '2 2', regex=True, inplace=True)
    Plink_PED_feats.replace(0, '1 1', regex=True, inplace=True)
    Plink_PED_feats.columns = ['']*len(Plink_PED_feats.columns)
    Plink_PED_feats.reset_index(inplace=True)
    Plink_PED_final = plink_PED.merge(Plink_PED_feats, left_on='Individual ID', right_on='index')
    Plink_PED_final.drop('index',axis =1, inplace=True)

    Plink_features = features.copy()
    Plink_MAP = pd.DataFrame({'rs# or snp identifier':Plink_features.columns.to_list()})
    Plink_MAP[['chromosome', 'Base-pair position (bp units)', 'RefAllele', 'SNPAllele']]= Plink_MAP['rs# or snp identifier'].str.split(':', expand=True)
    Plink_MAP.drop(['RefAllele', 'SNPAllele'], axis=1, inplace=True)
    cols = ['chromosome', 'rs# or snp identifier', 'Base-pair position (bp units)']
    Plink_MAP=Plink_MAP[cols]

    fastLMM_resp = resp.copy()
    fastLMM_resp.columns = ['Name', 'Name_copy', 'Response']
    fastLMM_resp.set_index('Name', inplace=True)
    ## Be sure to write the file without header
    
    if not os.path.isdir(str(path)):
        os.makedirs(str(path))
    os.chdir(str(path))
    Plink_PED_final.to_csv(filename+'_'+SimNum+".ped", sep='\t', index=False, header=False)
    Plink_MAP.to_csv(filename+'_'+SimNum+'.map', sep='\t', index=False, header=False)
    fastLMM_resp.to_csv(filename+'_'+SimNum+'_PHENOTYPE.txt', sep='\t', header=False)
    # return #Plink_PED_final, Plink_MAP, fastLMM_resp


###### LSBSR
def LSBSRsim_toScoary(features, resp, SimNum, path, filename):
    scoary_resp = resp.copy()
    
    scoary_empty = ['Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 'Avg sequences per isolate', 'Genome fragment', 'Order within fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']
    scoary=features.copy().T
    scoary.index.names = ['Gene']
    # scoary.set_index('Gene', inplace=True)
    for item in scoary_empty:
        scoary.insert(0,item,np.nan)
        
    scoary_resp.to_csv(filename+'_RESPONSE.csv')
    scoary.to_csv(filename+'_'+SimNum+'.csv')

def LSBSRsim_toPlink(feats, resp, SimNum, path, filename, parentPath):
    os.chdir(parentPath)
    Plink_PED_feats = feats.copy()
    plink_resp = resp.copy()
    plink_resp.replace(1, 2, regex=True, inplace=True)
    plink_resp.replace(0, 1, regex=True, inplace=True)
    Plink_PED_feats.replace(1, '2 2', regex=True, inplace=True)
    Plink_PED_feats.replace(0, '1 1', regex=True, inplace=True)
    Plink_PED_feats = plink_resp.join(Plink_PED_feats)
    
    Plink_features = feats.copy()
    Plink_MAP = pd.DataFrame({'rs# or snp identifier':Plink_features.columns.to_list()})
    Plink_MAP.reset_index(inplace=True)
    Plink_MAP['index'] = Plink_MAP['index']+1
    Plink_MAP[['chromosome']] = 1
    cols = ['chromosome', 'rs# or snp identifier', 'index']
    Plink_MAP=Plink_MAP[cols]
    
    if not os.path.isdir(str(path)):
        os.makedirs(str(path))
    os.chdir(str(path))
    Plink_PED_feats.reset_index().to_csv(filename+'_'+SimNum+".ped", sep='\t', index=False, header=False)
    Plink_MAP.to_csv(filename+'_'+SimNum+'.map', sep='\t', index=False, header=False)
    # return #Plink_PED_feats.reset_index(), Plink_MAP

def LSBSRsim_toTreeWAS(DF, resp, SimNum, indexDF, path, filename, parentPath):
    os.chdir(parentPath)
    indexDF['NEWid'] = indexDF['0'].str.split('_', expand=True)[0]
    indexDF.NEWid.replace('L1', '', regex=True, inplace=True)
    indexDF.set_index('NEWid', inplace=True)
    
    TreeWAS = DF.copy()
    TreeWAS = indexDF.join(TreeWAS)
    TreeWAS.set_index('0', inplace=True)
    TreeWAS.drop('RESPONSE', axis=1, inplace=True)
    TreeWAS.reset_index(inplace=True)
    TreeWAS_resp = resp.copy()
    TreeWAS_resp=TreeWAS_resp.join(indexDF)
    TreeWAS_resp.set_index('0', inplace=True)
    TreeWAS_resp = TreeWAS_resp[['response']]
    TreeWAS_resp.reset_index(inplace=True)
    TreeWAS_resp.columns=['Strains', str(SimNum)]
    
    if not os.path.isdir(str(path)):
        os.makedirs(str(path))
    os.chdir(str(path))
    TreeWAS.to_csv(filename+'_'+SimNum+"_FEATURES.txt", sep='\t', index=False)
    TreeWAS_resp.to_csv(filename+'_RESPOSNE.txt', sep='\t', index=False)
    # return #TreeWAS, TreeWAS_resp

def LSBSRsim_toBugWAS(DF):
    pass

def LSBSRsim_toFastLMM(feats, resp, SimNum, path, filename, parentPath):
    os.chdir(parentPath)
    Plink_PED_feats = feats.copy()
    plink_resp = resp.copy()
    plink_resp.replace(1, 2, regex=True, inplace=True)
    plink_resp.replace(0, 1, regex=True, inplace=True)
    Plink_PED_feats.replace(1, '2 2', regex=True, inplace=True)
    Plink_PED_feats.replace(0, '1 1', regex=True, inplace=True)
    Plink_PED_feats = plink_resp.join(Plink_PED_feats)
    
    Plink_features = feats.copy()
    Plink_MAP = pd.DataFrame({'rs# or snp identifier':Plink_features.columns.to_list()})
    Plink_MAP.reset_index(inplace=True)
    Plink_MAP['index'] = Plink_MAP['index']+1
    Plink_MAP[['chromosome']] = 1
    cols = ['chromosome', 'rs# or snp identifier', 'index']
    Plink_MAP=Plink_MAP[cols]
    
    fastLMM_resp = resp.copy()
    # fastLMM_resp.columns = ['Name', 'Name_copy', 'Response']
    fastLMM_resp['index_copy'] = fastLMM_resp.index.to_list()
    fastLMM_resp= fastLMM_resp.reset_index().set_index('index_copy')
    # fastLMM_resp.set_index('Name', inplace=True)
    
    if not os.path.isdir(str(path)):
        os.makedirs(str(path))
    os.chdir(str(path))
    Plink_PED_feats.reset_index().to_csv(filename+'_'+SimNum+".ped", sep='\t', index=False, header=False)
    Plink_MAP.to_csv(filename+'_'+SimNum+'.map', sep='\t', index=False, header=False)
    fastLMM_resp.to_csv('PHENOTYPE.txt', sep='\t', header=False)
    # return Plink_PED_feats.reset_index(), Plink_MAP

###### NASP
# def NASP_toScoary(DF):
#     pass

# def NASP_toPlink(DF):
#     pass

# def NASP_toTreeWAS(DF):
#     pass

def NASP_toBugWAS(DF):
    pass

def NASP_toFastLMM(DF):
    pass


# =============================================================================
# BacGWASim 
# =============================================================================
#### PLINK
# plink_feats_BacGWASim, plink_resp_BacGWASim = BacGWASim_toPlink(BacGWASim, BacGWASim_resp)
# plink_feats_BacGWASim.to_csv("H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAS_FileFormatting/TESTcases/plink_feats_BacGWASim.ped", sep='\t', index=False, header=False)
# plink_resp_BacGWASim.to_csv('H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAS_FileFormatting/TESTcases/plink_resp_BacGWASim.map', sep='\t', index=False, header=False)

# #### SCOARY
# scoary_feats_BacGWASim, scoary_resp_BacGWASim = BacGWASim_toScoary(BacGWASim,  BacGWASim_resp, 1)
# scoary_feats_BacGWASim.to_csv('H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAS_FileFormatting/TESTcases/scoary_feats_BacGWASim.csv')
# scoary_resp_BacGWASim.to_csv('H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAS_FileFormatting/TESTcases/scoary_resp_BacGWASim.csv', index=False)

# #### TreeWAS
# TreeWAS_feats_BacGWASim, TreeWAS_resp_BacGWASim = BacGWASim_toTreeWAS(BacGWASim, BacGWASim_resp, 'Heff_Hcorr_1')
# TreeWAS_feats_BacGWASim.to_csv("H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAS_FileFormatting/TESTcases/TreeWAS_feats_BacGWASim.txt", sep='\t')
# TreeWAS_resp_BacGWASim.to_csv('H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAS_FileFormatting/TESTcases/TreeWAS_resp_BacGWASim.txt', sep='\t')


# # =============================================================================
# # LSBSR
# # =============================================================================
# #### PLINK
# plink_feats_LSBSR, plink_resp_LSBSR = LSBSRsim_toPlink(LSBSR, SimResp)
# plink_feats_LSBSR.to_csv("H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAS_FileFormatting/TESTcases/plink_feats_LSBSR.ped", sep='\t', index=False, header=False)
# plink_resp_LSBSR.to_csv('H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAS_FileFormatting/TESTcases/plink_resp_LSBSR.map', sep='\t', index=False, header=False)

# #### SCOARY

# #### TreeWAS
# TreeWAS_feats_LSBSR, TreeWAS_resp_LSBSR = LSBSRsim_toTreeWAS(LSBSR, SimResp, 'LSBSR_Balanced_Sim1', treewas_response)
# TreeWAS_feats_LSBSR.to_csv("H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAS_FileFormatting/TESTcases/TreeWAS_feats_LSBSR.txt", sep='\t', index=False)
# TreeWAS_resp_LSBSR.to_csv('H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/GWAS_FileFormatting/TESTcases/TreeWAS_resp_LSBSR.txt', sep='\t', index=False)
