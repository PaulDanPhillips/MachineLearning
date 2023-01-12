import pandas as pd
import numpy as np
import os

class IToL_Generation:
    def __init__(self, mlst_file):
        self.mlst_file = pd.read_csv(mlst_file, sep='\t')
        
    def TEXT(self, outpath, filename):
        os.chdir(outpath)
        with open(str(filename), "w") as f:
            f.write("DATASET_TEXT")
            f.write("\nSEPARATOR COMMA")
            f.write("\nDATASET_LABEL,MLST Designation")
            f.write("\nCOLOR,#ff00ff")
            f.write("\nMARGIN,0")
            f.write("\nSHOW_INTERNAL,0")
            f.write("\nLABEL_ROTATION,0")
            f.write("\nALIGN_TO_TREE,0")
            f.write("\nSIZE_FACTOR,1")
            f.write("\nDATA\n")
            
    def COLOR_STRIP(self, outpath, filename):
        os.chdir(outpath)
        with open(str(filename), "w") as f:
            f.write("DATASET_COLORSTRIP")
            f.write("\nSEPARATOR COMMA")
            f.write("\nDATASET_LABEL,MLST Designation")
            f.write("\nCOLOR,#ff0000")
            f.write("\nSTRIP_WIDTH,25")
            f.write("\nMARGIN,0")
            f.write("\nBORDER_WIDTH,1")
            f.write("\nBORDER_COLOR,#000")
            f.write("\nSHOW_INTERNAL,0")
            f.write("\nLEGEND_TITLE,Strain Type Designation")
            f.write("\nLEGEND_SHAPES,1,1,1,1,1,1")
            f.write("\nLEGEND_COLORS,#FF0000,#0000FF,#33FF00,#660099,#CC6666,#FF9933")
            f.write("\nLEGEND_LABELS,109,56,266,36,132,novel")
            f.write("\nLEGEND_SHAPE_SCALES,0.7,0.7,0.7,0.7,0.7,0.7")
            f.write("\nDATA\n")
            
    def ReLABEL(self, outpath, filename):
        os.chdir(outpath)
        with open(str(filename), "w") as f:
            f.write("LABELS")
            f.write("\nSEPARATOR COMMA")
            f.write("\nDATA\n")
            
    def ColorRange(self, outpath, filename):
        os.chdir(outpath)
        with open(str(filename), "w") as f:
            f.write("TREE_COLORS")
            f.write("\nSEPARATOR COMMA")
            f.write("\nTREE_COLORS_LABEL,Origin Country")
            f.write("\nLEGEND_TITLE,Origin Country")
            f.write("\nLEGEND_SHAPES,0,0")
            f.write("\nLEGEND_COLORS,#eeffee,#ddddff")
            f.write("\nLEGEND_LABELS,Australia,Thailand")
            f.write("\nLEGEND_SHAPE_SCALES,1,1")
            f.write("\nDATA\n")
    def MLST_ORDER(self,mlst_file):
        mlst_file.all(axis ='columns')
        lst = ['N','X', 'T']
        DF=mlst_file[~mlst_file[['ace','gltB', 'gmhD', 'lepA', 'lipA', 'narK', 'ndh']].isin(lst)]
        DF_DROP =DF[DF.isnull().all(axis=1)].drop(columns=['genome', 'ST'])
        
        DF = DF.dropna(how='all', axis=0) #rows
        DF = DF.dropna(how='all', axis=1) #columns
        GenomeDF = mlst_file[['genome','ST']]
        DFall =DF.join(GenomeDF)

        DF_possibles =DFall[DFall.isnull().any(axis=1)]
        DF_sure = DFall.dropna()
        DF_DROP = DF_DROP.join(GenomeDF)
        
        cols = ['genome', 'ST', 'ace', 'gltB', 'gmhD', 'lepA', 'lipA', 'narK', 'ndh']
        # DF_possibles = DF_possibles[cols]
        DF_sure = DF_sure[cols]
        # DF_DROP = DF_DROP[cols]
        return DF_sure
    def MLST_TEXT(self, DF_sure, outpath, filename):
        DF_mlst = DF_sure[['genome', 'ST']].copy() #.set_index('genome').
        DF_mlst[['Strain', 'readID', 'SeqID', 'final', 'assembly']]=DF_mlst['genome'].str.split('_', expand=True)
        
        DF_mlst['TreeNodes']= DF_mlst['Strain']+'_' +DF_mlst['readID']+'_'+DF_mlst['SeqID']+'__pre-aligned_pre-called'
        
        dflst=[]
        for i in DF_mlst['ST'].unique():
            df = DF_mlst2.where(DF_mlst2['ST']==str(i)).dropna()
            dflst.append(df)
            
        ### for text labes
        color_lst = ['#FF0000', '#0000FF', '#33FF00', '#660099', '#663366', '#CC6666', '#FF9933']
        for i in range(0, len(dflst)):
            dflst[i]['position'] = '-1'
            dflst[i]['color'] = color_lst[i]
            dflst[i]['style'] = 'bold'
            dflst[i]['size_factor'] = '1.1'
            dflst[i]['rotation'] = '0'
        
        DF_mlst3 = pd.concat(dflst)
        os.chdir(otpath)
        DF_mlst3.to_csv(str(filename), mode='a', sep=',', header=False)


import seaborn as sns
palette = sns.color_palette(None, 3)
