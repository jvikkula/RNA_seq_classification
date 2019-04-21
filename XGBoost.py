from xgboost import XGBClassifier
import xgboost as xgb
import pandas as pd
import numpy as np
import sklearn
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
import time

def load_datasets(*args):
    """args: paths to csv files"""
    return {key:pd.read_csv(path, sep='\t') for key, path in enumerate(args)}

def concat_datasets(*kwargs):
    """kwargs: a dictionary of the dataframes, keys are integers from 0...n"""
    return pd.concat(kwargs[0].values(), axis=1).transpose()

def create_labels(*kwargs):
    """kwargs: a dictionary of the dataframes, keys are integers from 0...n"""
    df_dict = kwargs[0]
    labels = np.empty(0, dtype=int)
    for i in range(len(df_dict)):
        frame = df_dict[i]
        labels_set = np.empty(len(frame.columns), dtype=int)
        labels_set.fill(i+1)
        labels = np.concatenate((labels, labels_set), axis=0)

    return labels

def hyper_param_gsearch(param_dict):
    """args: dictionary of parameter options"""
    xgb_classifier = XGBClassifier(learning_rate=0.1,
                                   n_estimators=140,
                                   max_depth=5,
                                   min_child_weight=1,
                                   gamma=0,
                                   subsample=0.8,
                                   colsample_bytree=0.8,
                                   objective='count:poisson',
                                   nthread=4,
                                   scale_pos_weight=1,
                                   seed=27)

    gsearch = GridSearchCV(estimator=xgb_classifier,
                           param_grid=param_dict,
                           n_jobs=4,
                           iid=False,
                           cv=5)
    return gsearch

if __name__=="__main__":
    path_ductal = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/brca/brca_rsem_unprocessed/brca_ductal_rsem_unprocessed.csv"
    path_lobular = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/brca/brca_rsem_unprocessed/brca_lobular_rsem_unprocessed.csv"
    datasets = load_datasets(path_ductal, path_lobular)
    full_data = concat_datasets(datasets)
    labels = create_labels(datasets)
    X_train, X_test, y_train, y_test = train_test_split(full_data, labels, test_size=0.33, random_state=1234)

    param_dict = {
        'learning_rate':list(np.linspace(0.04, 0.1, 5)),
        'max_depth':[4,6,8,10],
        'min_child_weight':range(1,10,1),
        'sub_sample':list(np.linspace(0.5, 1, 5)),
        'colsample_bytree':list(np.linspace(0.4, 1, 5))
    }
    gsearch = hyper_param_gsearch(param_dict)
    xgb_train_start = time.perf_counter()
    gsearch.fit(X_train, y_train)
    xgb_train_end = time.perf_counter()
    xgb_training_time = xgb_train_end-xgb_train_start

    print(gsearch.best_params_)
    print(gsearch.best_score_)
    print("Time consumed for training model: %6.5f mins" % (xgb_training_time/60))
