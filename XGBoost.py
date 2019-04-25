from xgboost import XGBClassifier
import xgboost as xgb
import pandas as pd
import numpy as np
import seaborn as sn
import sklearn
from xgboost import plot_importance
from matplotlib import pyplot
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import confusion_matrix
from collections import Counter
import time
import os
import sys

###############################################################################
# Paths to datasets

# # UNPROCESSED DATASETS
path_ductal = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/brca/brca_rsem_unprocessed/brca_ductal_rsem_unprocessed.csv"
path_lobular = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/brca/brca_rsem_unprocessed/brca_lobular_rsem_unprocessed.csv"
path_colon_adeno = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/coadread/coadread_rsem_unprocessed/colon_adenocarcinoma_rsem_unprocessed.csv"
path_colon_mucinous = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/coadread/coadread_rsem_unprocessed/colon_mucinous_adenocarcinoma_rsem_unprocessed.csv"
path_rectal_adeno = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/coadread/coadread_rsem_unprocessed/rectal_adenocarcinoma_rsem_unprocessed.csv"
path_kidney_chromo = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/kipan/kipan_rsem_unprocessed/kidney_chromophobe_rsem_processed.csv"
path_kidney_clear_cell = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/kipan/kipan_rsem_unprocessed/kidney_clear_cell_renal_carcinoma_rsem_processed.csv"
path_kidney_renal_cell = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/kipan/kipan_rsem_unprocessed/kidney_papillary_renal_cell_carcinoma_rsem_processed.csv"
path_lgg_astro = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/lgg/lgg_rsem_unprocessed/astrocytoma_rsem_unprocessed.csv"
path_lgg_oligoastro = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/lgg/lgg_rsem_unprocessed/oligoastrocytoma_rsem_unprocessed.csv"
path_lgg_oligodendro = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/lgg/lgg_rsem_unprocessed/oligodendroglioma_rsem_unprocessed.csv"
path_meso_biphasic = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/meso/meso_rsem_unprocessed/biphasic_mesothelioma_rsem_unprocessed.csv"
path_meso_epithelioid = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/meso/meso_rsem_unprocessed/epithelioid_mesothelioma_rsem_unprocessed.csv"
path_sarc_leiomyosarcoma = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/sarc/sarc_rsem_unprocessed/leiomyosarcoma_rsem_unprocessed.csv"
path_sarc_liposarcoma = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/sarc/sarc_rsem_unprocessed/dedifferentiated_liposarcoma_rsem_unprocessed.csv"
path_sarc_pleomorphic = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/sarc/sarc_rsem_unprocessed/pleomorphic_mfh_rsem_unprocessed.csv"
path_sarc_myxofibrosarcoma = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/sarc/sarc_rsem_unprocessed/myxofibrosarcoma_rsem_unprocessed.csv"
path_sarc_un_pleomorphic = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/sarc/sarc_rsem_unprocessed/undifferentiated_pleomorphic_sarcoma_rsem_unprocessed.csv"
path_sarc_nerve_sheath = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/sarc/sarc_rsem_unprocessed/malignant_peripheral_nerve_sheath_tumors_rsem_unprocessed.csv"
path_thca_classical = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/thca/thca_rsem_unprocessed/thyroid_papillary_carcinoma_classical_rsem_unprocessed.csv"
path_thca_follicular = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/thca/thca_rsem_unprocessed/thyroid_papillary_carcinoma_follicular_rsem_unprocessed.csv"
path_thca_tall_cell = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/thca/thca_rsem_unprocessed/thyroid_papillary_carcinoma_tall_cell_rsem_unprocessed.csv"
path_thym_a = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/thym/thym_rsem_unprocessed/thymoma_type_a_rsem_unprocessed.csv"
path_thym_ab = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/thym/thym_rsem_unprocessed/thymoma_type_ab_rsem_unprocessed.csv"
path_thym_b1 = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/thym/thym_rsem_unprocessed/thymoma_type_b1_rsem_unprocessed.csv"
path_thym_b2 = "/~/LST/rnaseqV2data/unprocessed_rnaseqV2/thym/thym_rsem_unprocessed/thymoma_type_b2_rsem_unprocessed.csv"
path_thym_b3 = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/thym/thym_rsem_unprocessed/thymoma_type_b3_rsem_unprocessed.csv"
path_thym_c = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/thym/thym_rsem_unprocessed/thymoma_type_c_rsem_unprocessed.csv"
path_ucec_endometrial = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/ucec/ucec_rsem_unprocessed/endometrioid_endometrial_adenocarcinoma_rsem_unprocessed.csv"
path_ucec_mixed_serous = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/ucec/ucec_rsem_unprocessed/mixed_serous_and_endometrioid_rsem_unprocessed.csv"
path_ucec_serous_endo = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/ucec/ucec_rsem_unprocessed/serous_endometrial_adenocarcinoma_rsem_unprocessed.csv"
path_ucs_mixed = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/usc/ucs_rsem_unprocessed/malignant_mixed_mullerian_tumor_rsem_unprocessed.csv"
path_ucs_heterologous = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/usc/ucs_rsem_unprocessed/mmmt_heterologous_type_rsem_unprocessed.csv"
path_ucs_homologous = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/usc/ucs_rsem_unprocessed/mmmt_homologous_type_rsem_unprocessed.csv"

# PROCESSED DATASETS
# path_ductal = "~/LST/rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/brca_ductal_rsem_processed.csv"
# path_lobular = "~/LST/rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/brca_lobular_rsem_processed.csv"
# path_colon_adeno = "~/LST/rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed/colon_adenocarcinoma_rsem_processed.csv"
# path_colon_mucinous = "~/LST/rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed/colon_mucinous_adenocarcinoma_rsem_processed.csv"
# path_rectal_adeno = "~/LST/rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed/rectal_adenocarcinoma_rsem_processed.csv"
# path_kidney_chromo = "~/LST/rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_chromophobe_rsem_processed.csv"
# path_kidney_clear_cell = "~/LST/rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_clear_cell_renal_carcinoma_rsem_processed.csv"
# path_kidney_renal_cell = "~/LST/rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_papillary_renal_cell_carcinoma_rsem_processed.csv"
# path_lgg_astro = "~/LST/rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/astrocytoma_rsem_processed.csv"
# path_lgg_oligoastro = "~/LST/rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/oligoastrocytoma_rsem_processed.csv"
# path_lgg_oligodendro = "~/LST/rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/oligodendroglioma_rsem_processed.csv"
# path_meso_biphasic = "~/LST/rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/biphasic_mesothelioma_rsem_processed.csv"
# path_meso_epithelioid = "~/LST/rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/epithelioid_mesothelioma_rsem_processed.csv"
# path_sarc_leiomyosarcoma = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/leiomyosarcoma_rsem_processed.csv"
# path_sarc_liposarcoma = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/dedifferentiated_liposarcoma_rsem_processed.csv"
# path_sarc_pleomorphic = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/pleomorphic_mfh_rsem_processed.csv"
# path_sarc_myxofibrosarcoma = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/myxofibrosarcoma_rsem_processed.csv"
# path_sarc_un_pleomorphic = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/undifferentiated_pleomorphic_sarcoma_rsem_processed.csv"
# path_sarc_nerve_sheath = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/malignant_peripheral_nerve_sheath_tumors_rsem_processed.csv"
# path_thca_classical = "~/LST/rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_classical_rsem_processed.csv"
# path_thca_follicular = "~/LST/rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_follicular_rsem_processed.csv"
# path_thca_tall_cell = "~/LST/rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_tall_cell_rsem_processed.csv"
# path_thym_a = "~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_a_rsem_processed.csv"
# path_thym_ab = "~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_ab_rsem_processed.csv"
# path_thym_b1 = "~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_b1_rsem_processed.csv"
# path_thym_b2 = "/~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_b2_rsem_processed.csv"
# path_thym_b3 = "~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_b3_rsem_processed.csv"
# path_thym_c = "~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_c_rsem_processed.csv"
# path_ucec_endometrial = "~/LST/rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/endometrioid_endometrial_adenocarcinoma_rsem_processed.csv"
# path_ucec_mixed_serous = "~/LST/rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/mixed_serous_and_endometrioid_rsem_processed.csv"
# path_ucec_serous_endo = "~/LST/rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/serous_endometrial_adenocarcinoma_rsem_processed.csv"
# path_ucs_mixed = "~/LST/UCS-datasets/processed_rnaseqV2/malignant_mixed_mullerian_tumor_rsem_processed.csv"
# path_ucs_heterologous = "~/LST/UCS-datasets/processed_rnaseqV2/mmmt_heterologous_type_rsem_processed.csv"
# path_ucs_homologous = "~/LST/UCS-datasets/processed_rnaseqV2/mmmt_homologous_type_rsem_processed.csv"

# PROCESSED LOG2 DATASETS
# path_ductal = "~/LST/rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed_log2/brca_ductal_rsem_processed_log2.csv"
# path_lobular = "~/LST/rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed_log2/brca_lobular_rsem_processed_log2.csv"
# path_colon_adeno = "~/LST/rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed_log2/colon_adenocarcinoma_rsem_processed_log2.csv"
# path_colon_mucinous = "~/LST/rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed_log2/colon_mucinous_adenocarcinoma_rsem_processed_log2.csv"
# path_rectal_adeno = "~/LST/rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed_log2/rectal_adenocarcinoma_rsem_processed_log2.csv"
# path_kidney_chromo = "~/LST/rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed_log2/kidney_chromophobe_rsem_processed_log2.csv"
# path_kidney_clear_cell = "~/LST/rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed_log2/kidney_clear_cell_renal_carcinoma_rsem_processed_log2.csv"
# path_kidney_renal_cell = "~/LST/rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed_log2/kidney_papillary_renal_cell_carcinoma_rsem_processed_log2.csv"
# path_lgg_astro = "~/LST/rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed_log2/astrocytoma_rsem_processed_log2.csv"
# path_lgg_oligoastro = "~/LST/rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed_log2/oligoastrocytoma_rsem_processed_log2.csv"
# path_lgg_oligodendro = "~/LST/rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed_log2/oligodendroglioma_rsem_processed_log2.csv"
# path_meso_biphasic = "~/LST/rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed_log2/biphasic_mesothelioma_rsem_processed_log2.csv"
# path_meso_epithelioid = "~/LST/rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed_log2/epithelioid_mesothelioma_rsem_processed_log2.csv"
# path_sarc_leiomyosarcoma = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed_log2/leiomyosarcoma_rsem_processed_log2.csv"
# path_sarc_liposarcoma = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed_log2/dedifferentiated_liposarcoma_rsem_processed_log2.csv"
# path_sarc_pleomorphic = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed_log2/pleomorphic_mfh_rsem_processed_log2.csv"
# path_sarc_myxofibrosarcoma = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed_log2/myxofibrosarcoma_rsem_processed_log2.csv"
# path_sarc_un_pleomorphic = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed_log2/undifferentiated_pleomorphic_sarcoma_rsem_processed_log2.csv"
# path_sarc_nerve_sheath = "~/LST/rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed_log2/malignant_peripheral_nerve_sheath_tumors_rsem_processed_log2.csv"
# path_thca_classical = "~/LST/rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed_log2/thyroid_papillary_carcinoma_classical_rsem_processed_log2.csv"
# path_thca_follicular = "~/LST/rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed_log2/thyroid_papillary_carcinoma_follicular_rsem_processed_log2.csv"
# path_thca_tall_cell = "~/LST/rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed_log2/thyroid_papillary_carcinoma_tall_cell_rsem_processed_log2.csv"
# path_thym_a = "~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed_log2/thymoma_type_a_rsem_processed_log2.csv"
# path_thym_ab = "~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed_log2/thymoma_type_ab_rsem_processed_log2.csv"
# path_thym_b1 = "~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed_log2/thymoma_type_b1_rsem_processed_log2.csv"
# path_thym_b2 = "/~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed_log2/thymoma_type_b2_rsem_processed_log2.csv"
# path_thym_b3 = "~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed_log2/thymoma_type_b3_rsem_processed_log2.csv"
# path_thym_c = "~/LST/rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed_log2/thymoma_type_c_rsem_processed_log2.csv"
# path_ucec_endometrial = "~/LST/rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed_log2/endometrioid_endometrial_adenocarcinoma_rsem_processed_log2.csv"
# path_ucec_mixed_serous = "~/LST/rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed_log2/mixed_serous_and_endometrioid_rsem_processed_log2.csv"
# path_ucec_serous_endo = "~/LST/rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed_log2/serous_endometrial_adenocarcinoma_rsem_processed_log2.csv"
# path_ucs_mixed = "~/LST/rnaseqV2data/processed_rnaseqV2/usc/ucs_rsem_processed_log2/malignant_mixed_mullerian_tumor_rsem_processed_log2.csv"
# path_ucs_heterologous = "~/LST/rnaseqV2data/processed_rnaseqV2/usc/ucs_rsem_processed_log2/mmmt_heterologous_type_rsem_processed_log2.csv"
# path_ucs_homologous = "~/LST/rnaseqV2data/processed_rnaseqV2/usc/ucs_rsem_processed_log2/mmmt_homologous_type_rsem_processed_log2.csv"

###############################################################################
###############################################################################
# Data reading/writing
def get_file_paths():
    # All the paths to a list
    rootdir = '/home/leos/LST/rnaseqV2data/'
    paths = []
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            paths.append(os.path.join(subdir, file))
    print(paths)
    return paths

# Load the data from the paths given as args
def load_datasets(*args):
    """args: paths to csv files"""
    return {key:pd.read_csv(path, sep='\t') for key, path in enumerate(args)}

# Concatenate the distinct datasets in the order they are given as args
def concat_datasets(*kwargs):
    """kwargs: a dictionary of the dataframes, keys are integers from 0...n"""
    return pd.concat(kwargs[0].values(), axis=1).transpose()

# Create labels for the datasets
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

###############################################################################
###############################################################################
# Hyperparameter grid search
def hyper_param_gsearch(X_test, y_test):
    """args: test data and labels"""
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

    grid_params = {
        'learning_rate':list(np.linspace(0.04, 0.1, 5)),
        'max_depth':[4,6,8,10],
        'min_child_weight':range(1,10,1),
        'sub_sample':list(np.linspace(0.5, 1, 5)),
        'colsample_bytree':list(np.linspace(0.4, 1, 5))
    }

    gsearch = GridSearchCV(estimator=xgb_classifier,
                           param_grid=grid_params,
                           fit_params=fit_params,
                           n_jobs=4,
                           iid=False,
                           cv=3)
    return gsearch

# Run Hyperparameter grid search
def run_grid_search(X_test, y_test):
    """args: test data and labels""""
    gsearch = hyper_param_gsearch(X_test, y_test)
    xgb_train_start = time.perf_counter()
    gsearch.fit(X_test, y_test)
    xgb_train_end = time.perf_counter()
    xgb_training_time = xgb_train_end-xgb_train_start
    print(gsearch.best_params_)
    print(gsearch.best_score_)
    print("Time consumed for training model: %6.5f mins" % (xgb_training_time/60))

###############################################################################
###############################################################################
# Model fitting and stratified k-fold cross validation
def fit_single_model(model, X_train, y_train, eval_set, multi=False):
    """args: model = xgb-classifier object,
             X_train, y_train = train data and labels,
             eval_set = a list of test and train data/label tuples
             multi = flag to set the classification to multiclass or binary
    """
    # eval_metric based on binary/multiclass classification paradigm
    eval_metric=["merror", "mlogloss"] if multi else ["error", "logloss"]

    model.fit(X_train,
              y_train,
              early_stopping_rounds=20,
              eval_metric=eval_metric,
              eval_set=eval_set,
              verbose=False)

    # make predictions for test data
    y_pred = model.predict(X_test)
    predictions = [round(value) for value in y_pred]

    # evaluate predictions
    accuracy = accuracy_score(y_test, predictions)
    print("Accuracy: %.2f%%" % (accuracy * 100.0))
    print(Counter(predictions))
    return predictions, accuracy


def cross_validate(model, full_data, labels, fit_params, multi=False):
    """args: model=the classifier object,
             full_data=full count data df,
             labels=labels of the full data,
             fit_params=dictionary of model fitting parameters,
             multi=flag to set the classification to multiclass or binary
    """
    # eval_metric based on binary/multiclass classification paradigm
    eval_metric=["merror", "mlogloss"] if multi else ["error", "logloss"]

    kfold = StratifiedKFold(n_splits=5, random_state=1234)
    results = cross_val_score(model, full_data, labels, cv=kfold, fit_params=fit_params)
    y_pred = cross_val_predict(model, full_data, labels, cv=kfold, fit_params=fit_params)

    print("Accuracy: %.2f%% std: (%.2f%%)" % (results.mean()*100, results.std()*100))
    return results, y_pred
###############################################################################
###############################################################################
# PLOT FEATURE IMPORTANCES AND CLASSIFICATION ERRORS
def plot_errors(model, save=False, filename=None, multi=False):
    """args: model=classifier object,
             save=flag to save the image,
             filename=image filename
    """

    if multi:
        error='merror'
        logloss='mlogloss'
    else:
        error='error'
        logloss='mlogloss'

    # retrieve performance metrics
    results = model.evals_result()
    epochs = len(results['validation_0'][error])
    x_axis = range(0, epochs)

    # plot log loss
    fig, ax = pyplot.subplots()
    ax.plot(x_axis, results['validation_0'][logloss], label='Train')
    ax.plot(x_axis, results['validation_1'][logloss], label='Test')
    ax.legend()
    pyplot.ylabel('Log Loss')
    pyplot.title('XGBoost Log Loss')
    pyplot.show()

    # plot classification error
    fig, ax = pyplot.subplots()
    ax.plot(x_axis, results['validation_0'][error], label='Train')
    ax.plot(x_axis, results['validation_1'][error], label='Test')
    ax.legend()
    pyplot.ylabel('Classification Error')
    pyplot.title('XGBoost Classification Error')
    pyplot.show()

    if save:
        fig.savefig(filename, bbox_inches='tight')

# PLOT THE 100 MOST IMPORTANT FEATURES OF ONE SINGLE BOOSTER
def plot_importances(model save=False, filename=None):
    """args: model=classifier object,
             save=flag to save the image,
             filename=image filename
    """
    pyplot.rcParams['figure.figsize'] = [18, 18]
    fig = plot_importance(model, importance_type='gain', max_num_features=100, height=1)

    if save:
        fig.savefig(filename, bbox_inches='tight')

    pyplot.show()

# PLOT CONFUSION MATRIX
def plot_confusion(labels, prediction, save=False, filename=None):
    """args: labels=labels of the data, prediction=redicted labels, save=flag to save the image, filename=image filename"""
    conf = confusion_matrix(labels, y_pred)
    columns = ['class {}'.format(i) for i in np.unique(labels)]
    df_cm = pd.DataFrame(conf, index=columns, columns=columns)
    plot = sn.heatmap(df_cm, annot=True, fmt="d")

    if save:
        plot.savefig(filename)

    pyplot.show()

###############################################################################
###############################################################################
if __name__=="__main__":

    paths = get_file_paths()

    # Uncomment and comment accordingly
    # Binary datasets
    datasets = load_datasets(path_ductal, path_lobular)
    #datasets = load_datasets(path_kidney_clear_cell, path_kidney_renal_cell)
    #datasets = load_datasets(path_meso_biphasic, path_meso_epithelioid)
    multi = False

    # Multiclass datasets
    #datasets = load_datasets(path_colon_adeno, path_colon_mucinous, path_rectal_adeno)
    #datasets = load_datasets(path_kidney_chromo, path_kidney_clear_cell, path_kidney_renal_cell)
    #datasets = load_datasets(path_lgg_astro, path_lgg_oligoastro, path_lgg_oligodendro)
    #datasets = load_datasets(path_sarc_leiomyosarcoma, path_sarc_liposarcoma, path_sarc_pleomorphic, path_sarc_myxofibrosarcoma, path_sarc_un_pleomorphic, path_sarc_nerve_sheath)
    #datasets = load_datasets(path_thca_classical, path_thca_follicular, path_thca_tall_cell)
    #datasets = load_datasets(path_thym_a, path_thym_ab, path_thym_b1, path_thym_b2, path_thym_b3, path_thym_c)
    #datasets = load_datasets(path_ucec_endometrial, path_ucec_mixed_serous, path_ucec_serous_endo)
    #datasets = load_datasets(path_ucs_mixed, path_ucec_heterologous, path_ucec_homologous)
    #multi = True

    full_data = concat_datasets(datasets)
    labels = create_labels(datasets)
    X_train, X_test, y_train, y_test = train_test_split(full_data, labels, test_size=0.33, random_state=1234)

    # evaluation sets for model fitting
    eval_set = [(X_train, y_train), (X_test, y_test)]

    # Fitting parameters for gsearch and cv
    fit_params={"early_stopping_rounds":10,
                "eval_metric":"error",
                "eval_set":eval_set}

    ###########################################################
    # USE THIS ONLY FOR GRID SEARCHING THE BEST HYPERPARAMETERS
    # run_grid_search(X_test, y_test) # takes a long time
    ###########################################################

    # Create a XGB model based on the best parameters
    model = XGBClassifier(learning_rate=0.1,
                          n_estimators=200,
                          max_depth=4,
                          min_child_weight=1,
                          gamma=0,
                          subsample=0.75,
                          colsample_bytree=0.75,
                          objective='count:poisson',
                          nthread=4,
                          scale_pos_weight=1,
                          seed=27)

    predictions, accuracy = fit_single_model(model, X_train, y_train, eval_set, multi=multi)
    k_fold_accuracies, k_fold_predictions = cross_validate(model, full_data, labels, fit_params, multi=multi)
    plot_importances(model)
    plot_confusion(model)
    plot_errors(model, multi=multi)
