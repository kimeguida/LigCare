




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



import argparse
import numpy as np
from sklearn.neighbors import NearestNeighbors
import json

import iodesc, iomol2

from collections import defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns




plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
mpl.rcParams["figure.dpi"] = 150



parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cavdir', type=str, required=True,
            help='directory of cavities data')
parser.add_argument('-pr', '--preddir', type=str, required=True,
            help='directory of pruned cavities')
parser.add_argument('-ap', '--applicationset', type=str, required=True,
            help='external_test_pdb.json')
parser.add_argument('-d', '--descdir', type=str, required=True,
            help='directory of descriptors')
parser.add_argument('-l', '--labeldir', type=str, required=True,
            help='directory of labels')
parser.add_argument('-an', '--annotation', type=str, required=True,
            help='scpdb_annotation.tsv')
args = parser.parse_args()





def xtract_annotation(file_):
    annotations = {}
    with open(file_, 'r') as f:
        data = f.read().split('\n')
        del data[-1]
        del data[0]
    for l in data:
        cols = l.split('\t')
        pdb = cols[0]
        uniprot = cols[2]
        name = cols[3]
        annotations[pdb] = (uniprot, name)
    return annotations




annotations = xtract_annotation(args.annotation)

with open(args.applicationset, 'r') as f:
    entries = json.load(f)


stats = {}
uniprot_ids = []
names = []

for pdb in entries:
    uniprot_ids.append(annotations[pdb][0])
    names.append(annotations[pdb][1])

    try:
        cav_labels = iodesc.read_labels(f'{args.labeldir}/{pdb}_labels.npy')
    except:
        print(f'problem encountered while reading: {args.labeldir}/{pdb}_labels.npy')
        continue
    final_labels = {}
    tn = 0
    tp = 0
    fn = 0
    fp = 0

    for idx, l1, l2, l3, l4 in cav_labels:
        idx = int(idx)
        if l2 == 1 and l3 == 1:
            label = 1
        else:
            label = 0
        final_labels[idx] = label

    try:
        cavall_coords, _ = iomol2.read_cav_mol2(f'{args.cavdir}/{pdb}_cavityALL.mol2')
    except:
        print(f'problem encountered while reading: {args.cavdir}/{pdb}_cavityALL.mol2')
        continue
    try:
        cavpred_coords, _ = iomol2.read_cav_mol2(f'{args.preddir}/{pdb}_cavityALL_pharm.mol2')
    except:
        print(f'problem encountered while reading: {args.preddir}/{pdb}_cavityALL_pharm.mol2')
        continue

    #if len(cavpred_coords) == 0:
    #    print('EMPTY')
    distances = []
    indices = []
    if len(cavpred_coords) > 0:
        neigh = NearestNeighbors(n_neighbors=len(cavpred_coords), radius=0.1,
                               algorithm='ball_tree').fit(cavpred_coords)
        distances, indices = neigh.radius_neighbors(cavall_coords)

    for i, n in enumerate(indices):
        #print(n.shape)
        if n.size == 0:
            if final_labels[i] == 0:
                tn += 1
            elif final_labels[i] == 1:
                fn += 1
        elif len(list(n.flat)) == 1:
            if final_labels[i] == 0:
                fp += 1
            elif final_labels[i] == 1:
                tp += 1

    size_cavall = len(list(cavall_coords))
    size_cavpred = len(list(cavpred_coords))

    proportion_kept = round(size_cavpred/size_cavall*100, 2)

    if tn+fp == 0:
        specificity = None
    else:
        specificity = round(tn/(tn+fp), 2)

    if tp+fn == 0:
        sensitivity = None
    else:
        sensitivity = round(tp/(tp+fn), 2)


    if sensitivity is None or specificity is None:
        balanced_accuracy = None
    else:
        balanced_accuracy = round((specificity + sensitivity)/2, 2)

    #print(pdb, size_cavall, size_cavpred, proportion_kept, tn, tp, fn, fp, specificity, sensitivity, balanced_accuracy)


    stats[pdb] = [size_cavall, size_cavpred, proportion_kept, tn, tp, fn, fp, specificity, sensitivity, balanced_accuracy]
    
    #if sensitivity >= 0.9 and specificity >= 0.9:
    #    print(pdb)

    #if sensitivity == 1:
    #    print(pdb)

    if specificity == 1:
        print(pdb)



with open('ext_application_1000.tsv', 'w') as of:
    of.write('pdb\tsize_cavity_all\tsize_prediction\tproportion_kept\t'
        'true_negative\ttrue_positive\tfalse_negative\tfalse_positive\t'
        'specificity\tsensitivity\tbalanced_accuracy\tuniprot_ac\tname\n')
    for pdb in stats:
        of.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(pdb, *stats[pdb], *annotations[pdb]))

        




all_specificity = []
all_sensitivity = []
all_baccuracy = []

for pdb in stats:
    if stats[pdb][9] is None:
        continue
    all_sensitivity.append(stats[pdb][8])
    all_specificity.append(stats[pdb][7])
    all_baccuracy.append(stats[pdb][9])


print(len(all_specificity))




all_proportion_kept = []
all_size_cavall = []
for pdb in stats:
    
    all_proportion_kept.append(stats[pdb][2])
    all_size_cavall.append(stats[pdb][0])



print('Uniprot', len(set(uniprot_ids)))
print('name', len(set(names)))
for n in set(names):
    #print(names.count(n))
    if names.count(n) >= 10:
        print(names.count(n), n)








g1 = sns.JointGrid(marginal_ticks=True, ratio=2)
sns.scatterplot(x=all_sensitivity, y=all_specificity, ax=g1.ax_joint, alpha=0.5, color=plt.cm.get_cmap('RdPu')(0.7))
sns.distplot(all_sensitivity, ax=g1.ax_marg_x, kde=False, bins=np.arange(0, 1+0.05, 0.05),
        color=plt.cm.get_cmap('RdPu')(0.7), hist_kws={'edgecolor':'k'})
sns.distplot(all_specificity, ax=g1.ax_marg_y, vertical=True, kde=False, bins=np.arange(0, 1+0.05, 0.05),
        color=plt.cm.get_cmap('RdPu')(0.7), hist_kws={'edgecolor':'k'})
g1.ax_marg_x.set_xlim(-0.05, 1.05)
#g1.ax_marg_x.set_yticks([0, 100, 200, 300])

g1.ax_marg_y.set_ylim(-0.05, 1.05)
#g1.ax_marg_y.set_xticks([0, 100, 200, 300])

g1.ax_joint.set_ylabel('Specificity')
g1.ax_joint.set_xlabel('Sensitivity')
g1.ax_marg_x.set_ylabel('Count')
g1.ax_marg_y.set_xlabel('Count')

plt.tight_layout()



g2 = sns.JointGrid(marginal_ticks=True, ratio=2)
sns.scatterplot(x=all_size_cavall, y=all_proportion_kept, ax=g2.ax_joint, alpha=0.5, color=plt.cm.get_cmap('RdPu')(0.7))
sns.distplot(all_size_cavall, ax=g2.ax_marg_x, kde=False, bins=np.arange(50, 900, 50),
        color=plt.cm.get_cmap('RdPu')(0.7), hist_kws={'edgecolor':'k'})
sns.distplot(all_proportion_kept, ax=g2.ax_marg_y, vertical=True, kde=False, bins=np.arange(0, 100, 5),
        color=plt.cm.get_cmap('RdPu')(0.7), hist_kws={'edgecolor':'k'})
#g2.ax_marg_x.set_xlim(-0.05, 1.05)
#g2.ax_marg_x.set_yticks([0, 100, 200, 300])

#g2.ax_marg_y.set_ylim(-0.05, 1.05)
#g2.ax_marg_y.set_xticks([0, 50, 100,])

g2.ax_joint.set_ylabel('Proportion kept (%)')
g2.ax_joint.set_xlabel("# points in 'cavity ALL'")
g2.ax_marg_x.set_ylabel('Count')
g2.ax_marg_y.set_xlabel('Count')
plt.tight_layout()



plt.figure(3, (4, 4))
g3 = sns.distplot(all_baccuracy, bins=np.arange(0, 1+0.05, 0.05), kde=False,
                color=plt.cm.get_cmap('RdPu')(0.7), hist_kws={'edgecolor':'k'})
g3.set_xlabel('Balanced accuracy')
g3.set_ylabel('Count')
print(g3)
plt.tight_layout()

plt.show()


