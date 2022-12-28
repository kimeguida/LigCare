

# Copyright (c) 2022 LIT
# For academic use only.
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.




import numpy as np
import random
import operator
import argparse
import pickle
import json
import copy



from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, cross_validate
from sklearn import metrics
from sklearn.naive_bayes import GaussianNB
import xgboost as xgb

#import matplotlib as mpl
#import matplotlib.pyplot as plt

import iodesc, iomol2




parser = argparse.ArgumentParser()
parser.add_argument('-t', '--trainingdata', type=str, required=True,
            help='features_split.json')
parser.add_argument('-d', '--descdir', type=str, required=True,
            help='directory of descriptors')
parser.add_argument('-clf', '--classifier', type=str, required=True,
            choices=['rf', 'xgboost'],
            help='ML classifier, rf: random forest or xgboost: XGBoost')
args = parser.parse_args()





def list_pdb(training_data_, ph4_list_=None):
    pdb_entries = []
    if ph4_list_ is None:
        ph4_list_ = list(training_data_.keys())
    for ph4 in ph4_list_:
        for pdb, idx in training_data_[ph4]['0']:
            pdb_entries.append(pdb)
        for pdb, idx in training_data_[ph4]['1']:
            pdb_entries.append(pdb)

    return set(pdb_entries)





def load_desc(pdb_entries_, desc_dir_):
    descriptors = {}
    for pdb in pdb_entries_:
        ########################### content of descriptors ###########################
        #iodesc.write_desc(indexes, [query_ph4_fingerprint, query_bsa_fingerprint, query_dist_centroids,
        #                    neigh_count_fingerprint, neigh_bsa_fingerprint, bsa_desc,
        #                    interactions_desc, interactions_type_desc, interactions_type_all_desc],
        #                    args.output)


        # query_ph4_fingerprint : 8 bins --> [0-7]
        # positions: CA': 0, 'CZ': 1, 'N': 2, 'NZ': 3, 'O': 4, 'OD1': 5, 'OG': 6, 'DU': 7
        
        # query_bsa_fingerprint: 12 bins --> [8-19]
        # query_dist_centroids: 9 bins --> [20-28]
        # positions: all, CA': 0, 'CZ': 1, 'N': 2, 'NZ': 3, 'O': 4, 'OD1': 5, 'OG': 6, 'DU': 7

        # neigh_count_fingerprint: 8*3 = 24 bins --> [29-52]
        # neigh_bsa_fingerprint: (12*8 = 96) * 3 = 288 bins --> [53-340] for all neighbors,
        #                                                   --> [53-148] for first neighbor circle,
        #                                                   --> [149-244] for second neighbor circle,
        #                                                   --> [245-340] for third neighbor circle

        # buriedness: 1 bin --> [341-341]
        
        # interactions_desc: 3 bins --> [342-344]
        # interactions_type_desc: 8 bins --> [345-352]
        # interactions_type_all_desc: 8 bins --> [353-360]


        #iodesc.write_desc(indexes, [query_ph4_fingerprint, query_bsa_fingerprint, query_dist_centroids,
        #                    neigh_count_fingerprint, neigh_bsa_fingerprint,
        #                    [[buriedness[p]] for p in range(len(buriedness))],
        #                    interactions_desc, interactions_type_desc, interactions_type_all_desc],
        #                    args.output)

        
        #indexes, cav_desc = iodesc.read_desc(f'{desc_dir_}/{pdb}_desc.npy', [[8, 332],])
        #indexes, cav_desc = iodesc.read_desc(f'{desc_dir_}/{pdb}_desc.npy', [[8, 19], [21, 28], [45, 140]])
        #indexes, cav_desc = iodesc.read_desc(f'{desc_dir_}/{pdb}_desc.npy', [[53,148], [341,341], [342,344], [345,352]])
        ##############################################################################

        # selection of descriptors
        indexes, cav_desc = iodesc.read_desc(f'{desc_dir_}/{pdb}_desc.npy',
                                            [[53,148], [341,341], [345,352]])
        
        for idx in indexes:
            descriptors[(pdb, idx)] = cav_desc[idx]
    
    return descriptors



def write_report(cv_results_, model_results, naive_bayes_, ofile_):
    with open(ofile_, 'w') as of:
        of.write('### Cross validation accuracy\n')
        of.write('{:<10}{:<10}\n'.format('train', 'test'))
        for i in range(len(cv_results_['test_accuracy'])):
            of.write('{:<10}{:<10}\n'.format(round(cv_results_['train_accuracy'][i], 4),
                                            round(cv_results_['test_accuracy'][i], 4),))
        of.write('\n### Model accuracy\n')
        of.write('{:<10}{:<10}\n'.format('train', 'test'))
        of.write('{:<10}{:<10}\n'.format(*model_results))

        of.write('\n### Naive Bayes\n')
        of.write('Successful predictions: {} %\n'.format(naive_bayes_))





def train_rf(ph4_, train_features_, train_labels_, test_features_, test_labels_, n_estimators_=50):

    # Gaussian classifier
    #n_estimators = int(np.sqrt(train_features.shape[0]))
    #n_estimators = 100
    print('n_estimators', n_estimators_)
    
    # cross validation
    print('######################## cross validation')
    clf=RandomForestClassifier(n_estimators=n_estimators_)
    # default n split: sqrt(n), n = train_features.shape[1]
    cv_results = cross_validate(clf, train_features_, train_labels_, cv=5, scoring=['accuracy'], return_train_score=True)
    print(cv_results)
    del clf
    

    # train the model using the training sets
    print('######################## training')
    clf=RandomForestClassifier(n_estimators=n_estimators_)
    clf.fit(train_features_, train_labels_)

    train_pred=clf.predict(train_features_)
    train_acc = round(metrics.accuracy_score(train_labels_, train_pred), 4)
    print("Accuracy train:", train_acc)


    test_pred=clf.predict(test_features_)
    test_acc = round(metrics.accuracy_score(test_labels_, test_pred), 4)
    print("Accuracy external test:", test_acc)
    
    success = (test_labels_ == test_pred).sum()
    print("Number of correct predictions out of a total {} points : {} ({} %)".format(
                    test_features_.shape[0], success, round(100*success/test_features_.shape[0], 2))
         )

    #print(clf.feature_importances_)
    #feature_importances = [[f, imp] for f, imp in zip(range(333), clf.feature_importances_)]
    #feature_importances = list(sorted(feature_importances, key=operator.itemgetter(1), reverse=True))
    #print(feature_importances)

    # saving model
    print("...................saving model")
    with open(f'{ph4_}.model', 'wb') as mf:
        pickle.dump(clf, mf)

    del clf


    print('######################## naive_bayes')
    gnb = GaussianNB()
    test_pred = gnb.fit(train_features_, train_labels_).predict(test_features_)
    success = (test_labels_ == test_pred).sum()
    proportion = round(100*success/test_features_.shape[0], 2)
    print("Number of correct predictions out of a total {} points : {} ({} %)".format(
                    test_features_.shape[0], success, proportion)
         )


    write_report(cv_results, (train_acc, test_acc), proportion, ph4_+'.report')

    # testing
    #np.random.shuffle(train_labels)
    #print(train_labels)
   
    #print(train_features[:, 12])
    #np.random.shuffle(train_features[:, 12])
    #print(train_features[:, 12])
    #random.seed(42)


    #indexes_to_shuffle = random.sample(range(325), 300)
    #print(indexes_to_shuffle)
    #indexes_to_shuffle = [12, 13,]
    #for i in indexes_to_shuffle:
    #    np.random.shuffle(train_features[:, i])


    #s_features = copy.deepcopy(train_features[:, 12])
    #print(s_features)
    #print(s_features.flatten())
    #shuffled_features = np.random.shuffle(train_features[:, 12].flatten())
    #print(shuffled_features)
    #random.seed(42)
    #for i in range(len(shuffled_features)):
    #    train_features[i, 12] = shuffled_features[i]





def train_xgboost(ph4_, train_features_, train_labels_, test_features_, test_labels_,):

    
    # cross validation
    print('######################## cross validation')
    clf=xgb.XGBClassifier()
    # default n split: sqrt(n), n = train_features.shape[1]
    cv_results = cross_validate(clf, train_features_, train_labels_, cv=5, scoring=['accuracy'], return_train_score=True)
    print(cv_results)
    del clf
    

    # train the model using the training sets
    print('######################## training')
    clf=xgb.XGBClassifier()
    clf.fit(train_features_, train_labels_)

    train_pred=clf.predict(train_features_)
    train_acc = round(metrics.accuracy_score(train_labels_, train_pred), 4)
    print("Accuracy train:", train_acc)


    test_pred=clf.predict(test_features_)
    test_acc = round(metrics.accuracy_score(test_labels_, test_pred), 4)
    print("Accuracy external test:", test_acc)
    
    success = (test_labels_ == test_pred).sum()
    print("Number of correct predictions out of a total {} points : {} ({} %)".format(
                    test_features_.shape[0], success, round(100*success/test_features_.shape[0], 2))
         )

    #print(clf.feature_importances_)
    #feature_importances = [[f, imp] for f, imp in zip(range(333), clf.feature_importances_)]
    #feature_importances = list(sorted(feature_importances, key=operator.itemgetter(1), reverse=True))
    #print(feature_importances)

    # saving model
    print("...................saving model")
    with open(f'{ph4_}.model', 'wb') as mf:
        pickle.dump(clf, mf)

    del clf


    print('######################## naive_bayes')
    gnb = GaussianNB()
    test_pred = gnb.fit(train_features_, train_labels_).predict(test_features_)
    success = (test_labels_ == test_pred).sum()
    proportion = round(100*success/test_features_.shape[0], 2)
    print("Number of correct predictions out of a total {} points : {} ({} %)".format(
                    test_features_.shape[0], success, proportion)
         )


    write_report(cv_results, (train_acc, test_acc), proportion, ph4_+'.report')




    # testing
    #np.random.shuffle(train_labels)
    #print(train_labels)
   
    #print(train_features[:, 12])
    #np.random.shuffle(train_features[:, 12])
    #print(train_features[:, 12])
    #random.seed(42)


    #indexes_to_shuffle = random.sample(range(325), 300)
    #print(indexes_to_shuffle)
    #indexes_to_shuffle = [12, 13,]
    #for i in indexes_to_shuffle:
    #    np.random.shuffle(train_features[:, i])


    #s_features = copy.deepcopy(train_features[:, 12])
    #print(s_features)
    #print(s_features.flatten())
    #shuffled_features = np.random.shuffle(train_features[:, 12].flatten())
    #print(shuffled_features)
    #random.seed(42)
    #for i in range(len(shuffled_features)):
    #    train_features[i, 12] = shuffled_features[i]




def train(ph4_, data_, descriptors_, classifier_):


    if classifier_ not in ['xgboost', 'rf']:
        print('select appropriate classifier name: xgboost, rf')
        return


    # reading and spliting data into train and external test set
    positives = data_[ph4_]['1']
    negatives = data_[ph4_]['0']

    positives_train, positives_test = train_test_split(positives,
                                        test_size=0.25, random_state=42)
    negatives_train, negatives_test = train_test_split(negatives,
                                        test_size=0.25, random_state=42)

    train_features = []
    train_labels = []
    test_features = []
    test_labels = []


    for pdb_idx in positives_train:
        train_labels.append(1)
        train_features.append(descriptors_[tuple(pdb_idx)])

    for pdb_idx in negatives_train:
        train_labels.append(0)
        train_features.append(descriptors_[tuple(pdb_idx)])


    for pdb_idx in positives_test:
        test_labels.append(1)
        test_features.append(descriptors_[tuple(pdb_idx)])

    for pdb_idx in negatives_test:
        test_labels.append(0)
        test_features.append(descriptors_[tuple(pdb_idx)])


    train_features = np.array(train_features)
    train_labels = np.array(train_labels)
    test_features = np.array(test_features)
    test_labels = np.array(test_labels)
    #print(train_labels)

    if classifier_ == 'xgboost':
        print('########################################## training XGB')
        train_xgboost(ph4_, train_features, train_labels, test_features, test_labels)

    if classifier_ == 'rf':
        print('########################################## training Random Forest')
        train_rf(ph4_, train_features, train_labels, test_features, test_labels, n_estimators_=100)







print('########################################## loading data')

with open(args.trainingdata, 'r') as f:
    training_data = json.load(f)

for ph4 in ['NZ', 'OG', 'OD1', 'N', 'O', 'CA', 'CZ']:
#for ph4 in ['NZ',]:

    pdb_entries = list_pdb(training_data, [ph4])
    print(len(pdb_entries), 'pdb entries')
    
    print('########################################## loading descriptors for', ph4)
    descriptors = load_desc(pdb_entries, args.descdir)

    print('########################################## training ', ph4)
    train(ph4, training_data, descriptors, classifier_=args.classifier)


    del descriptors, pdb_entries

