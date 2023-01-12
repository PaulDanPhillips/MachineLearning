import pandas as pd
import os
import numpy as np
import scipy.stats as stats
# from sklearn.datasets import load_iris
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold, train_test_split
from ClassTuning import Tuning
import matplotlib.pyplot as plt
import seaborn as sns
# from statsmodels.stats.outliers_influence import variance_inflation_factor
# from statsmodels.tools.tools import add_constant
# from sklearn.preprocessing import StandardScaler
import TunedParamAnalysis as PMA
from pathlib import Path

# =============================================================================
# Preliminary Analysis
# =============================================================================

titanic = pd.read_csv('Titanic_train.csv', index_col=0)

titanic.drop(columns=['Cabin'], inplace=True)
titanic.dropna(inplace=True)

random_state = np.random.RandomState(42)

Y = titanic['Survived'].copy()
X = titanic.drop(columns=['Survived', 'Name', 'Ticket'])
X.dtypes

# X=pd.get_dummies(data=X, columns=['Sex', 'Embarked'])
X['Sex'] = X['Sex'].map({'male':0, 'female':1})
X['Embarked'] = X['Embarked'].map({'S':0, 'C':1, 'Q':2})


######## Preliminary Figs
df = titanic.copy()
df['Age'] =pd.cut(df['Age'], range(0,90,10))
df['Fare'] =pd.cut(df['Fare'], range(0,520,20))
df.drop(columns=['Ticket'], inplace=True)
val_counts = []
for col in df.columns:
    tmp = pd.DataFrame(df[col].value_counts())
    if len(tmp) < len(titanic):
        val_counts.append(tmp)
        tmp.plot.bar()


sns.pairplot(titanic.loc[:,titanic.dtypes != 'object'])

corr = titanic.loc[:,titanic.dtypes != 'object'].corr()
sns.heatmap(corr, xticklabels=corr.columns, yticklabels=corr.columns, cmap=sns.diverging_palette(220, 10, as_cmap=True))

# =============================================================================
# Tuning with RandomSeachCV
# =============================================================================

GBparam_space = {'learning_rate': stats.uniform(0.001, 0.2),
                'n_estimators': [50, 100, 200],
                'subsample': stats.uniform(0.8, 0.2),
                'min_samples_split': [int(x) for x in np.linspace(2,20,5)],
                'max_depth':[int(x) for x in np.linspace(3,50,5)],
                'max_features':['sqrt', 'log2', 0.1*len(X.columns)]}


max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
max_depth.append(None)
RFparam_space= {'bootstrap': [True, False],
  'max_depth': max_depth,
  'max_features': ['sqrt', 'log2'],
  'min_samples_leaf': [1, 2, 4],
  'min_samples_split': [2, 5, 10, 20],
  'n_estimators': [int(x) for x in np.linspace(start = 200, stop = 2000, num = 7)]}


LRparam_space= {'penalty': ['l1', 'l2'],
                  'C': np.logspace(4,-4,25)}


ENETparam_space= {'C': np.logspace(4,-4,25),
              'l1_ratio':np.linspace(0.1, 0.9, 9)}


os.chdir('H:/My Drive/Machine_Learning/TutorialDataSets/Titanic/RandomSearchTuned')
enet=LogisticRegression(penalty= 'elasticnet', solver='saga', max_iter=1e6)
lr = LogisticRegression(solver='saga', max_iter=1e7)
gb = GradientBoostingClassifier()
rf = RandomForestClassifier()


TuneObj = Tuning(n_splits=10, repeat=5,  n_jobs=4)
RF = TuneObj.Tuning_RandomSearch_classify(X=X, Y=Y, repeat=5, n_splits=10, score='accuracy', mod=rf, hyperparameters=RFparam_space, n_iter=100, n_jobs=4)
RF.to_csv('Titanic_RFrandSearch.txt', sep='\t')

LR = TuneObj.Tuning_RandomSearch_classify(X=X, Y=Y, repeat=5, n_splits=10, score='accuracy', mod=lr, hyperparameters=LRparam_space, n_iter=10, n_jobs=4)
LR.to_csv('Titanic_ENETrandSearch.txt', sep='\t')

ENET = TuneObj.Tuning_RandomSearch_classify(X=X, Y=Y, repeat=5, n_splits=10, score='accuracy', mod=enet, hyperparameters=ENETparam_space, n_iter=25, n_jobs=4)
ENET.to_csv('Titanic_ENETrandSearch.txt', sep='\t')

GB = TuneObj.Tuning_RandomSearch_classify(X=X, Y=Y, repeat=5, n_splits=10, score='accuracy', mod=gb, hyperparameters=GBparam_space, n_iter=50, n_jobs=4)
GB.to_csv('Titanic_GBrandSearch.txt', sep='\t')


# =============================================================================
# RandomSearch Analysis
# =============================================================================
Titan = PMA.Param_analysis('H:/My Drive/Machine_Learning/TutorialDataSets/Titanic/RandomSearchTuned')

LR, Enet, RF, GB = Titan.File_input(Titan.path, LR_name='Titanic_LRrandSearch.txt', Enet_name='Titanic_ENETrandSearch.txt', RF_name='Titanic_RFrandSearch.txt', GB_name='Titanic_GBrandSearch.txt')

LRbest, LRworst = Titan.BestNworst_mean(LR, 'mean_test_score', 10)
LRbest = LRbest.loc[LRbest['mean_test_score']>0.5]
LRbox =  Titan.Param_figs(LRbest, 'mean_test_score', 'Titanic Random LR Tuning')
LRparams = pd.Series(LRbest.params.unique())
# LRgridParams={'penalty':['l1', 'l2'],
#               'C':np.linspace(3,0.1, 5)
#     }

Enetbest, Enetworst = Titan.BestNworst_mean(Enet, 'mean_test_score', 10)
Enetbest = Enetbest.loc[Enetbest['mean_test_score']>0.5]
Enetbox =  Titan.Param_figs(Enetbest, 'mean_test_score', 'Titanic Random Enet Tuning')
Enetparams = pd.Series(Enetbest.params.unique())
# EnetgridParams={'C':np.linspace(5,0.05, 5),
#                 'l1_ratio': np.linspace(0.7, 0.9, 5)
#     }


GBbest, GBworst = Titan.BestNworst_mean(GB, 'mean_test_score', 10)
GBbest = GBbest.loc[GBbest['mean_test_score']>0.5]
GBbox =  Titan.Param_figs(GBbest, 'mean_test_score', 'Titanic Random GB Tuning')
GBparams = pd.Series(GBbest.params.unique())
# GBgridParams={'learning_rate':np.linspace(0.001, 0.05, 5),
#               'max_depth':np.linspace(3, 75, 7),
#               'max_features':['sqrt', 'log2', 0.5, 0.7],
#               'min_samples_split': [int(x) for x in np.linspace(10,25, 5)],
#               'n_estimators':[50, 200, 500],
#                'subsample':np.linspace(0.8, 0.97, 5)
#     }


RFbest, RFworst = Titan.BestNworst_mean(RF, 'mean_test_score', 10)
RFbest = RFbest.loc[RFbest['mean_test_score']>0.5]
RFbox =  Titan.Param_figs(RFbest, 'mean_test_score', 'Titanic Random RF Tuning')
RFparams = pd.Series(RFbest.params.unique())
n_estimators=RFbest.param_n_estimators.mean()
min_samples_split = RFbest.param_min_samples_split.mean()
min_samples_leaf = RFbest.param_min_samples_leaf.mean()
max_depth = RFbest.param_max_depth.mean()
# RFgridParams={'n_estimators':[750, 1200, 1400, 1500],
#               'min_samples_split':[10, 12, 20, 22],
#               'min_samples_leaf':[1,2,4],
#               'max_depth':[80, 85, 90, None]
#     }


# imgpath= Path('H:/My Drive/Machine_Learning/Projects_by_Organism/Bpseudo/DTRA_AMR/CORRECTED/TUNING/NewNASP/OUT/ParamAnalysis')
# RFbox.save(imgpath/'NASP-CAZm_RFTuningAUC.png', height=5, width=5, dpi=300)
# # RFbox2.save(imgpath/'NASP-CAZm_RFTuningAUC_BA.png', height=5, width=5, dpi=300)

# =============================================================================
# Tuning with GridSearchCV
# =============================================================================
TuneObj = Tuning(n_splits=10, repeat=5,  n_jobs=4)

lr = LogisticRegression(solver='saga', max_iter=1e7)
LRgridParams={'penalty':['l1', 'l2'],
              'C':np.linspace(3,0.1, 5)
    }
LR = TuneObj.Tuning_GridSearch_classify(X=X, Y=Y, repeat=5, n_splits=10, score='accuracy', mod=lr, hyperparameters=LRgridParams, n_jobs=5)


enet=LogisticRegression(penalty= 'elasticnet', solver='saga', max_iter=1e6)
EnetgridParams={'C':np.linspace(5,0.05, 5),
                'l1_ratio': np.linspace(0.7, 0.9, 5)
    }
ENET = TuneObj.Tuning_GridSearch_classify(X=X, Y=Y, repeat=5, n_splits=10, score='accuracy', mod=enet, hyperparameters=EnetgridParams, n_jobs=5) 

rf = RandomForestClassifier(max_features='log2', bootstrap=True)
RFgridParams={'n_estimators':[750, 1200, 1400, 1500],
              'min_samples_split':[10, 12, 20, 22],
              'min_samples_leaf':[1,2,4],
              'max_depth':[80, 85, 90, None]
    }
RF = TuneObj.Tuning_GridSearch_classify(X=X, Y=Y, repeat=5, n_splits=10, score='accuracy', mod=rf, hyperparameters=RFgridParams, n_jobs=5)

gb = GradientBoostingClassifier()
GBgridParams={'learning_rate':np.linspace(0.01, 0.1, 4),
              'max_depth':np.linspace(10, 50, 3),
              'max_features':['sqrt', 'log2', 0.7],
              'min_samples_split': [int(x) for x in np.linspace(10,25, 3)],
              'n_estimators':[50, 200, 250],
               'subsample':np.linspace(0.8, 0.97, 3)
    }
GB = TuneObj.Tuning_GridSearch_classify(X=X, Y=Y, repeat=5, n_splits=10, score='accuracy', mod=gb, hyperparameters=GBgridParams, n_jobs=3)

os.chdir('H:/My Drive/Machine_Learning/TutorialDataSets/Titanic/GridSearchTuned')
# LR.to_csv('LRGridTune.txt', sep='\t')
# ENET.to_csv('ENETGridTune.txt', sep='\t')
GB.to_csv('GBGridTune.txt', sep='\t')
# RF.to_csv('RFGridTune.txt', sep='\t')
