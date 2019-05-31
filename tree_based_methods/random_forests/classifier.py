#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 01 15:16:10 2019

@author: Shazam Kasher

RNA Classification

"""
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier

def read_split_csv(train_path, test_path, labels_path):

    """ Read and split data into train and val using iloc

     Args:
         args.train_path: Takes path to train_data
         args.test_path: Takes path to test_data
         args.labels_path: Takes path to train_labels
    
    Returns: 
     
        train_data, val_data, test_data, train_labels and val_labels
    """
    
    # Read the CSV file and split using iloc   
    df = pd.read_csv(train_path, header=None)
    test_df = pd.read_csv(test_path, header=None)
    labels = pd.read_csv(labels_path, header=None)

    # Split data and labels into training and validation   
    train_df = df.iloc[0:2750,:]
    val_df = df.iloc[2750:3750,:]
    train_labels = labels.iloc[0:2750,:]
    val_labels = labels.iloc[2750:3750,:]
    
    return train_df, val_df, test_df, train_labels, val_labels

def class_sample_count(predictions):
    
    """ Count the number of -1 and 1 labels in data

     Args:
         predictions: Takes validaion or test predictions
    
    Returns: 
     
        Total count of -1 or 1 labels
    """
    
    positive = 0
    negative = 0

    for i in range(len(predictions)):
        if predictions[i] == -1:
            negative += 1
        else:
            positive += 1
    return negative, positive

def plot(data, labels, plot_name=""):
    
    """ Plot training, validation and test data

     Args:
         data: PCA dataframes
         labels: labels of training, validation and test data
         plot_name: Name of the plot
    
    Returns: 
     
        Scatter plots
    """
    
    label_colors = ["red" if x == -1 else "blue" for x in labels.values]

    plt.figure()
    plt.scatter(
            x=data[:,0], 
            y=data[:,1], 
            alpha=0.4, 
            c=label_colors)
    plt.title(plot_name)
    plt.show()
    
def pca(train_df, val_df, test_df):
    
    """ Apply PCA to training, Validation and test data

     Args:
         train_df: Training dataframe
         val_df: Validation dataframe
         test_df: Test dataframe
    
    Returns: 
     
        PCA of training, validation and test dataframes
    """
    
    n_components = 2
    pca = PCA(n_components=n_components)
    train_pca = pca.fit_transform(train_df.copy())
    val_pca = pca.transform(val_df.copy())
    test_pca = pca.transform(test_df.copy())
    variance = pca.explained_variance_ratio_
    
    print("-"*40)
    print("Percentage of PCA Variances: {0}".format(sum(variance)*100))
    print("-"*40)
    
    return train_pca, val_pca, test_pca
    
def classifier(train_pca, val_pca, test_pca):
    
    """ Run the classifier on the PCA data and predict test data labels

     Args:
         train_pca: PCA train data
         val_pca: PCA validation data
         test_pca: PCA test data
    
    Returns: 
     
        Predicted labels for test data as numpy array
    """
    
    clf = RandomForestClassifier(n_estimators=100, 
                                 max_depth=None, 
                                 min_samples_split=2, 
                                 random_state=0)
    clf.fit(train_pca, train_labels.values.ravel())
    
    accuracy = clf.score(val_pca, val_labels)
    val_predictions = clf.predict(val_pca)
    negative, positive = class_sample_count(val_predictions)
    
    print("-"*40)
    print("Classifier Accuracy: {0}".format(accuracy) + "\n")
    print("Validation Set Predictions")
    print("Number of Data Points for Class -1: {0}".format(negative))
    print("Number of Data Points for Class 1: {0}".format(positive))
    print("-"*40)
    
    test_labels = clf.predict(test_pca)
    # Convert label numpy array to dataframe
    test_labels_df = pd.DataFrame(test_labels)
    
    negative, positive = class_sample_count(test_labels)        
    print("-"*40)
    print("Test Set Predictions")
    print("Number of Data Points for Class -1: {0}".format(negative))
    print("Number of Data Points for Class 1: {0}".format(positive))
    print("-"*40)

    return test_labels_df

if __name__=="__main__":    
    
    # Provide data paths
    train_data_path = "/home/shazam/Desktop/Data/train_data.csv"
    test_data_path = "/home/shazam/Desktop/Data/test_data.csv"
    train_labels_path = "/home/shazam/Desktop/Data/train_labels.csv"
    
    # Read CSV data and split into training and validation
    train_df, val_df, test_df, train_labels, val_labels = read_split_csv(
            train_data_path, 
            test_data_path, 
            train_labels_path)
        
    # Run PCA on data
    train_pca, val_pca, test_pca = pca(train_df, val_df, test_df)
        
    # Run classifier on PCA data and predict labels
    test_labels_df = classifier(train_pca, val_pca, test_pca)
        
    # Save as csv file in working directory
    test_labels_df.to_csv("test_labels.csv", sep=",", index=False)
    
    # Show scatter plots
    plot(train_pca, train_labels, plot_name="PCA Training Data")
    plot(val_pca, val_labels, plot_name="PCA Validation Data")
    plot(test_pca, test_labels_df, plot_name="Prediction Labels on Test Data")


    
