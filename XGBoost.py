from xgboost import XGBClassifier
import xgboost as xgb
import pandas as pd
import numpy as np
import sklearn
from xgboost import plot_importance
from matplotlib import pyplot
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from collections import Counter
import time

###############################################################################
# Paths to datasets

# UNPROCESSED DATASETS
path_ductal = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/brca/brca_rsem_unprocessed/brca_ductal_rsem_unprocessed.csv"
path_lobular = "~/LST/rnaseqV2data/unprocessed_rnaseqV2/brca/brca_rsem_unprocessed/brca_lobular_rsem_unprocessed.csv"
path_colon_adeno = "/home/leos/LST/rnaseqV2data/unprocessed_rnaseqV2/coadread/coadread_rsem_unprocessed/colon_adenocarcinoma_rsem_unprocessed.csv"
path_colon_mucinous = "/home/leos/LST/rnaseqV2data/unprocessed_rnaseqV2/coadread/coadread_rsem_unprocessed/colon_mucinous_adenocarcinoma_rsem_unprocessed.csv"
path_rectal_adeno = "/home/leos/LST/rnaseqV2data/unprocessed_rnaseqV2/coadread/coadread_rsem_unprocessed/rectal_adenocarcinoma_rsem_unprocessed.csv"

###############################################################################
###############################################################################
# Data reading/writing
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

###############################################################################
###############################################################################
# Hypeeparameter grid search
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
def fit_single_model(model, X_train, y_train, eval_set):
    model.fit(X_train,
              y_train,
              early_stopping_rounds=20,
              eval_metric=["error", "logloss"],
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


def cross_validate(model, full_data, labels, fit_params):
    kfold = StratifiedKFold(n_splits=5, random_state=1234)
    results = cross_val_score(model, full_data, labels, cv=kfold, fit_params=fit_params)
    print("Accuracy: %.2f%% std: (%.2f%%)" % (results.mean()*100, results.std()*100))
    return results
###############################################################################
###############################################################################
# PLOT FEATURE IMPORTANCES AND CLASSIFICATION ERRORS
def plot_errors(model, save=False, filename=None):
    """args: model = classifier object, save=flag to save the image, filename=image filename"""
    # retrieve performance metrics
    results = model.evals_result()
    epochs = len(results['validation_0']['error'])
    x_axis = range(0, epochs)

    # plot log loss
    fig, ax = pyplot.subplots()
    ax.plot(x_axis, results['validation_0']['logloss'], label='Train')
    ax.plot(x_axis, results['validation_1']['logloss'], label='Test')
    ax.legend()
    pyplot.ylabel('Log Loss')
    pyplot.title('XGBoost Log Loss')
    pyplot.show()

    # plot classification error
    fig, ax = pyplot.subplots()
    ax.plot(x_axis, results['validation_0']['error'], label='Train')
    ax.plot(x_axis, results['validation_1']['error'], label='Test')
    ax.legend()
    pyplot.ylabel('Classification Error')
    pyplot.title('XGBoost Classification Error')
    pyplot.show()

    if save:
        fig.savefig(filename, bbox_inches='tight')

def plot_importances(model save=False, filename=None):
    pyplot.rcParams['figure.figsize'] = [18, 18]
    fig = plot_importance(model, importance_type='gain', max_num_features=100, height=1)

    if save:
        fig.savefig(filename, bbox_inches='tight')
###############################################################################
###############################################################################
if __name__=="__main__":
    datasets = load_datasets(path_ductal, path_lobular)
    full_data = concat_datasets(datasets)
    labels = create_labels(datasets)
    X_train, X_test, y_train, y_test = train_test_split(full_data, labels, test_size=0.33, random_state=1234)

    # Fitting parameters for gsearch and cv
    fit_params={"early_stopping_rounds":20,
                "eval_metric" : "error",
                "eval_set" : [[X_test, y_test]]}

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

    # evaluation sets for model fitting
    eval_set = [(X_train, y_train), (X_test, y_test)]

    predictions, accuracy = fit_single_model(model, X_train, y_train, eval_set)
    k_fold_accuracies = cross_validate(model, full_data, labels, fit_params)
