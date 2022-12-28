


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
from collections import defaultdict

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, cross_validate
from sklearn import metrics
import json


import matplotlib as mpl
import matplotlib.pyplot as plt


import iodesc, iomol2



parser = argparse.ArgumentParser()
parser.add_argument('-v', '--validfile', type=str, required=True,
            help='external_test_pdb.json')
parser.add_argument('-c', '--cavdir', type=str, required=True,
            help='directory of cavities')
parser.add_argument('-d', '--descdir', type=str, required=True,
            help='directory of descriptors ')
parser.add_argument('--models', type=str, required=True, nargs='+',
            help='model')
parser.add_argument('--modelnames', type=str, required=True, nargs='+',
            help='model names in the same order')
parser.add_argument('-o', '--odir', type=str, required=True,
            help='directory of pruned cavities ')
args = parser.parse_args()





models = {}
for name, model in zip(args.modelnames, args.models):
    models[name] = pickle.load(open(model, 'rb'))


with open(args.validfile, 'r') as f:
    validation_data = json.load(f)


#for pdb in validation_data[:50]:
for pdb in validation_data:

    try:
        #print('................. predicting', pdb)
        coords, ph4_features = iomol2.read_cav_mol2(f'{args.cavdir}/{pdb}_cavityALL.mol2')
        _, desc = iodesc.read_desc(f'{args.descdir}/{pdb}_desc.npy', [[53,148], [341,341], [345,352]])
                                                                     # same as used in training
    except:
        print(f'problem encountered while reading: {args.cavdir}/{pdb}_cavityALL.mol2, {args.descdir}/{pdb}_desc.npy')
        continue

    groups_indexes = defaultdict(list)
    groups = defaultdict(list)
    for i, ph4 in enumerate(ph4_features):
        if ph4 == 'DU':
            continue
        groups_indexes[ph4].append(i)
        groups[ph4].append(desc[i])

    pred_coordinates_ph4 = []
    kept_indexes = []
    for ph4 in groups:
        valid_pred=models[ph4].predict(groups[ph4])
        #print(ph4, valid_pred)
        for label, idx in zip(valid_pred, groups_indexes[ph4]):
            if label == 1:
                pred_coordinates_ph4.append(list(coords[idx])+[ph4])
                kept_indexes.append(idx)

    ofile = f'{args.odir}/{pdb}_cavityALL_pharm.mol2'
    iomol2.write_mol2(ofile, pred_coordinates_ph4, macromol_="SMALL")

    
