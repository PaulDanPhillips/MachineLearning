import pandas as pd
import random
import numpy as np
import os
import pickle
from scipy.stats import chi2_contingency

# =============================================================================
# Simulate Genes
# =============================================================================

def chi2_Sim_Test(filename, SimNum, CNTRLs, ABX):
    CNTRL_cols = [col for col in CNTRLs.columns if 'CNTRL' in col]
    plst=[]
    crssTb = []
    for col in CNTRL_cols:
        CT2 = pd.crosstab(CNTRLs[ABX], CNTRLs[col])
        stat, p, dof, expected = chi2_contingency(CT2)
        plst.append(p)
        crssTb.append(CT2)
    for i in range(len(plst)):
        print('CNTRL_'+str(i+1)+ ' chi2 pval is: ' + str(plst[i]))
    with open(filename+'_'+SimNum+'CHI2Vals.pickle', 'wb') as file:
        pickle.dump(plst, file)
    with open(filename+'_'+SimNum+'CROSStab.pickle', 'wb') as file:
        pickle.dump(crssTb, file)
    # return plst, crssTb

def geneSim(filename, SimNum, responseDF, inputDF, NumCol, PhenoToONE, TotalAssoc=False, numNonAssoc=None, fracAsocc=None):
    #May want to split this up into CNTRL_(1-5) and CNTRL_(6-10)
    predictorDF = inputDF.copy()
    print(len(predictorDF.columns))
    ResponseCOPY = responseDF.copy()
    col = ResponseCOPY.columns

    for i in range(1, int(NumCol)+1):
        ResponseCOPY['CNTRL_'+str(i)]=np.where(ResponseCOPY[col] ==PhenoToONE, 1, 0)

    if TotalAssoc==True:
        CNTRL_cols = [col for col in ResponseCOPY.columns if 'CNTRL' in col]
        responseDFcntrl=ResponseCOPY[CNTRL_cols].copy()
        for item in responseDFcntrl:
            Extract = responseDFcntrl[item]
            place = random.randint(0,(len(predictorDF.columns)))
            predictorDF.insert(place,item, Extract)
        print(len(predictorDF.columns))
        chi2_Sim_Test(filename, SimNum, responseDF.join(responseDFcntrl), 'response')
        predictorDF.to_csv(filename+'_'+SimNum+'.txt', sep='\t')
        return responseDF.join(responseDFcntrl)#, predictorDF
    else:
        CNTRL_cols = [col for col in ResponseCOPY.columns if 'CNTRL' in col]
        responseDFcntrl=ResponseCOPY[CNTRL_cols].copy()
        
        for col in responseDFcntrl:
            idx = responseDFcntrl.loc[responseDFcntrl[col] == 1].sample( n=numNonAssoc, frac=fracAsocc).index
            responseDFcntrl.loc[idx, col]=0
        for item in responseDFcntrl:
            Extract = responseDFcntrl[item]
            place = random.randint(0,(len(predictorDF.columns)))
            predictorDF.insert(place,item, Extract)
        print(len(predictorDF.columns))
        chi2_Sim_Test(filename, SimNum,responseDF.join(responseDFcntrl), 'response')
        predictorDF.to_csv(filename+'_'+SimNum+'.txt', sep='\t')
        return responseDF.join(responseDFcntrl)#, predictorDF
    
    

# =============================================================================
# Convert LSBSR for various GWAS programs
# =============================================================================
def LSBSRsim_toPyseer(features, resp, SimNum, path, filename):
    pyseer_feats = features.T.copy()
    pyseer_feats.index.names = ['Gene']
    pyseer_feats.to_csv(filename+'_'+SimNum+'.Rtab', sep='\t')
    pyseer_resp = resp.copy()
    pyseer_resp.to_csv(filename+'RESPONSE.txt', sep='\t')
    
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
