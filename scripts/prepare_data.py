

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




# splitting traning and test sets
# splitting modes: random, target classes, time,


import numpy as np
import random
import operator
import argparse
import os
from collections import defaultdict
import json

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, cross_validate
from sklearn import metrics

import iodesc, iomol2








def extract_entries(file_):
    with open(file_, 'r') as f:
        entries = f.read().split('\n')
    del entries[-1]
    entries = [e.split('/')[-1] for e in entries]
    return entries



def xtract_point_ph4(point_ph4_desc_):

    # same positions from iodesc
    
    positions = {'CA': 0,
                 'CZ': 1,
                 'N': 2,
                 'NZ': 3,
                 'O': 4,
                 'OD1': 5,
                 'OG': 6,
                 'DU': 7,
                 }

    # unique ph4 will be assigned even if multiple in desc

    if point_ph4_desc_[positions['CZ']] == 1:
        ph4 = 'CZ'
    elif point_ph4_desc_[positions['CA']] == 1:
        ph4 = 'CA'
    elif point_ph4_desc_[positions['OD1']] == 1:
        ph4 = 'OD1'
    elif point_ph4_desc_[positions['NZ']] == 1:
        ph4 = 'NZ'
    elif point_ph4_desc_[positions['OG']] == 1:
        ph4 = 'OG'
    elif point_ph4_desc_[positions['O']] == 1:
        ph4 = 'O'
    elif point_ph4_desc_[positions['N']] == 1:
        ph4 = 'N'
    elif point_ph4_desc_[positions['DU']] == 1:
        ph4 = 'DU'

    return ph4      




def count_unique_pdb(list_):
    pdbs = [t[0] for t in list_]
    counter  = {}
    for pdb in pdbs:
        counter[pdb] = counter.get(pdb, 0) + 1
    #print(counter)
    counts = {i:list(counter.values()).count(i) for i in range(1, 1+max(list(counter.values())))}
    
    return counts



def xtract_annotation(file_):
    missing_annotation = []
    pdb_kw = {}
    uniprot_ac = {}
    complexes = defaultdict(list)
    rscb_infos = {}
    lig_het = {}
    with open(file_, 'r') as f:
        for l in f:
            if l.startswith('pdb_nsite'):
                continue
            cols = l.split('\n')[0].split('\t')
            if 'None' in cols:
                missing_annotation.append(cols[0])
                continue
            pdb_nsite = cols[0]
            het = cols[1]
            ac = cols[2]
            name = cols[4]
            keywords = cols[5].split(';')
            keywords = [kw.strip() for kw in keywords]
            resolution = float(cols[6])
            depo_date = cols[7]
            release_date = cols[8]

            pdb_kw[pdb_nsite] = keywords
            uniprot_ac[pdb_nsite] = ac
            complexes[(het, ac)].append(pdb_nsite)
            rscb_infos[pdb_nsite] = (resolution, depo_date, release_date)
            lig_het[pdb_nsite] = het

    return pdb_kw, uniprot_ac, complexes, rscb_infos, lig_het, missing_annotation




def remove_redundant(complexes_, rscb_infos_):
    kept = []
    for het, ac in complexes_:
        if ac != 'None':
            if len(complexes_[(het, ac)]) > 1:
                #get resol
                pdb_resol = [(e, rscb_infos_[e][0]) for e in complexes_[(het, ac)] if rscb_infos_[e][0] != 'None']
                #print(pdb_resol)
                #print(sorted(pdb_resol, key=operator.itemgetter(1)))
                selected = sorted(pdb_resol, key=operator.itemgetter(1))[0][0]
                kept.append(selected)
            else:
                kept.append(complexes_[(het, ac)][0])
    return kept




def count_overlap_entries(set_1_, set_2_, uniprot_ac_, het_):
    uniprot_ac_entries_1 = set([uniprot_ac_[e] for e in set_1_])
    uniprot_ac_entries_2 = set([uniprot_ac_[e] for e in set_2_])

    overlap_ac = [e for e in uniprot_ac_entries_1 if e in uniprot_ac_entries_2]
    
    het_entries_1 = set([het_[e] for e in set_1_])
    het_entries_2 = set([het_[e] for e in set_2_])
    
    overlap_het = [e for e in het_entries_1 if e in het_entries_2]

    counts = {'Uniprot_AC_overlap': (len(uniprot_ac_entries_1), len(uniprot_ac_entries_2), len(overlap_ac)),
                'HET_overlap': ((len(het_entries_1), len(het_entries_2), len(overlap_het)))}
    
    return counts



if __name__ == '__main__':
    




    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=str, required=True,
                help='file list of cavities IDs, eg. 2rh1_1')
    parser.add_argument('-d', '--descdir', type=str, required=True,
                help='path to descriptors directory, data/desc')
    parser.add_argument('-l', '--labeldir', type=str, required=True,
                help='path to descriptors directory, data/labels')
    parser.add_argument('-a', '--annotation', type=str, required=True,
                help='annonation file: e.g. scpdb_annotation.tsv')
    parser.add_argument('--split', type=str, required=True,
                choices=['random', 'time', 'gpcr', 'kinase', 'protease', 'nuclear_receptor', 'druglike'],
                help='data spliting scheme')
    args = parser.parse_args()




    pdb_kw, uniprot_ac, complexes, rscb_infos, lig_het, missing_annotation = xtract_annotation(args.annotation)
    entries = extract_entries(args.file)
    non_redundant_entries = remove_redundant(complexes, rscb_infos)

    print('total entries:', len(entries))
    print('non redundant:', len(non_redundant_entries))
    #print(non_redundant_entries)
    print("missing_annotation:", len(missing_annotation))




    random.seed(42)


    # split data into model and application set

    # random split
    if args.split == 'random':
        sample_size = 1000
        external_test_pdb = random.sample(non_redundant_entries, sample_size)

    # time split
    elif args.split == 'time':
        year_thr = 2020
        external_test_pdb = [e for e in non_redundant_entries if int(rscb_infos[e][2].split('-')[0]) >= year_thr]

    # target classes split
    elif args.split == 'gpcr':
        external_test_pdb = [e for e in non_redundant_entries if 'G-protein coupled receptor' in pdb_kw[e]]

    elif args.split == 'kinase':
        external_test_pdb = [e for e in non_redundant_entries if 'Kinase' in pdb_kw[e]]

    elif args.split == 'nuclear_receptor':
        external_test_pdb = [e for e in non_redundant_entries if 'Nucleus' in pdb_kw[e] and 'Receptor' in pdb_kw[e]]

    elif args.split == 'protease':
        external_test_pdb = []
        for e in non_redundant_entries:
            for kw in pdb_kw[e]:
                if 'protease' in kw or 'Protease' in kw:
                    external_test_pdb.append(e)



    train_test_pdb = [e for e in non_redundant_entries if e not in set(external_test_pdb)]

    data_overlap_counts = count_overlap_entries(train_test_pdb, external_test_pdb, uniprot_ac, lig_het)



    print("######## external test set")
    print('Total:', len(external_test_pdb))
    #print(external_test_pdb)


    print("######## train test set")
    print('Total:', len(train_test_pdb))



    print('########### OVERLAP ############')
    print(data_overlap_counts)



    features = {'CA': {0:[], 1:[]},
                'CZ': {0:[], 1:[]},
                'O': {0:[], 1:[]},
                'N': {0:[], 1:[]},
                'OD1': {0:[], 1:[]},
                'OG': {0:[], 1:[]},
                'NZ': {0:[], 1:[]},
                'DU': {0:[], 1:[]},
                }


    for i, pdb in enumerate(train_test_pdb):
        if not os.path.isfile(f'{args.labeldir}/{pdb}_labels.npy') or\
             not os.path.isfile(f'{args.descdir}/{pdb}_desc.npy'):
             continue

        cav_labels = iodesc.read_labels(f'{args.labeldir}/{pdb}_labels.npy')
        _, cav_desc = iodesc.read_desc(f'{args.descdir}/{pdb}_desc.npy', [[0, 7]]) # 0-7 ph4 encoding
        # positives : near ligand of same ph4 and interacting according to IChem
        for idx, l1, l2, l3, l4 in cav_labels:
            idx = int(idx)
            if l2 == 1 and l3 == 1:
                label = 1
            else:
                label = 0
            ph4 = xtract_point_ph4(cav_desc[idx])
            features[ph4][label].append((pdb, idx))
        if not (i % 3000):
            print('progress', i)



    print("######## before balancing labels")
    for ph4 in features:
        print(ph4, '| 0:', len(features[ph4][0]), '| 1:', len(features[ph4][1]))

    # balancing positive and negative labels
    print("######## counter #points/PDB after balancing. contribution : number of PDB concerned")
    selected_features = {}
    for ph4 in features:
        if ph4 != 'DU':
            positives_sample = features[ph4][1]
            negatives_sample = random.sample(features[ph4][0], len(features[ph4][1]))
            counts = count_unique_pdb(negatives_sample)
            print('0 |', ph4, counts)
            counts = count_unique_pdb(positives_sample)
            print('1 |', ph4, counts)
            selected_features[ph4] = {1:positives_sample, 0:negatives_sample}
            #print(ph4, negatives_sample)


    print("######## after balancing labels")
    for ph4 in selected_features:
        print(ph4, '| 0:', len(selected_features[ph4][0]), '| 1:', len(selected_features[ph4][1]))




    with open('features_split.json', 'w') as of:
        json.dump(selected_features, of)


    with open('external_test_pdb.json', 'w') as of:
        json.dump(external_test_pdb, of)
        


    #print(external_test_pdb)