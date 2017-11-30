# Author: Yan Kai
# This script is for training a model from the training data.

import re,os,sys
from optparse import OptionParser
import sklearn
from sklearn import svm
from sklearn.metrics import *
from sklearn.cross_validation import *
from sklearn import preprocessing
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier, RandomTreesEmbedding
from sklearn.metrics import classification_report
from sklearn.externals import joblib
import numpy as np
import pandas as pd
from scipy import interp
import matplotlib.pyplot as plt

def get_true_and_proba(estimator, X, y, n_folds, cv):
    ys = [[],[]]
    for train_idx, valid_idx in cv:
        clf = estimator
        if isinstance(X, np.ndarray):
            clf.fit(X[train_idx], y[train_idx])
            prob = clf.predict_proba(X[valid_idx])
            cur_prob = prob[:,1]
        elif isinstance(X, pd.DataFrame):
            clf.fit(X.iloc[train_idx, :], y[train_idx]) 
            prob = clf.predict_proba(X.iloc[valid_idx, :])
            cur_prob = prob[:,1]
        else:
            raise Exception('Only numpy array and pandas DataFrame as types of X are supported')

        ys[0] = ys[0]+list(y[valid_idx])
        ys[1] = ys[1]+list(cur_prob)
    return ys

def fit_and_score_CV(estimator, X, y, n_folds=10, stratify=True):
    
    if not stratify:
        cv_arg = sklearn.cross_validation.KFold(y.size, n_folds)
    else:
        cv_arg = sklearn.cross_validation.StratifiedKFold(y, n_folds)
    

    ys = get_true_and_proba(estimator, X, y, n_folds, cv_arg)    
    return ys

def get_pr(clf,data):
    X = data.iloc[:,4:]
    Y = np.array(data['response'])
    X = X.as_matrix()
    y_true, y_pred = Y, clf.predict(X)
    probas = clf.predict_proba(X)
    precisions, recalls, thresholds = precision_recall_curve(y_true, probas[:,1])
    au_pr = average_precision_score(y_true, probas[:,1], average='micro')
    return ((recalls,precisions,au_pr))



def main(argv):
    parser = OptionParser()
    parser.add_option("-t", "--train", action="store", type="string", dest="train", metavar="<file>", help="The path of the training data")
    parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="The complete path for the resulting model and relevant results")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 4:
        parser.print_help()
        sys.exit(1)

    
    print 'Reading in features...'
    data = pd.read_table(opt.train)
    features = data.columns.values[4:]
    
    
    # Evaluating whether the data are balanced..
    num_pos = data[data['response'] == 1].shape[0]
    num_neg = data[data['response'] == 0].shape[0]
    n2p = float(num_neg)/float(num_pos)
    print 'Negative samples were '+str(n2p)+' times more than the positive loops...'
    
    X = data.iloc[:,4:].as_matrix()
    Y = np.array(data['response'])

    # Use Random Forest classifier to predict interactions
    random_state = np.random.RandomState(0)
    n_estimator = 100
    rf = RandomForestClassifier(max_depth = 36,max_features=18, n_estimators=n_estimator, random_state = random_state,n_jobs=-1, class_weight={0:1,1:n2p}) 

    print 'Training model...'
    ys = fit_and_score_CV(rf, X, Y, n_folds=10, stratify=True)         
    precisions, recalls, thresholds = precision_recall_curve(ys[0], ys[1])
    au_pr = average_precision_score(ys[0], ys[1], average='micro')
    plt.figure()
    plt.plot(recalls, precisions, lw=2, label="AU-PR=%0.2f"%(au_pr))
    plt.legend(loc='lower left',fontsize=14)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('PR curve',fontsize=16)
    plt.savefig(opt.output+'/PRC.png', dpi=300)
    print 'Precision-Recall Curve generated and saved...'


    fpr, tpr, thresholds = roc_curve(ys[0], ys[1])
    au_roc = auc(fpr, tpr)
    plt.figure()
    plt.plot(fpr, tpr, lw=2, label="AU-ROC=%0.2f"%(au_roc))
    plt.legend(loc='lower right',fontsize=14)
    plt.xlabel('1-Specificity')
    plt.ylabel('Sensitivity')
    plt.title('ROC curve',fontsize=16)
    plt.savefig(opt.output+'/ROC.png', dpi=300)
    print 'Receiver Operating Characteristic Curve generated and saved...'

    print 'Using the entire training data to generate a model ...'
    loop_model = rf.fit(X, Y)
    model = str(opt.output)+"/model.joblib.pkl"
    _ = joblib.dump(loop_model, model, compress=9)
    print 'Trained model saved.'


if __name__ == "__main__":
	main(sys.argv)
