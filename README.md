# LigCare

## Prediction of interacting points
### 1. Calculation of cavity descriptors
```python scripts/cavity_descriptors.py -c data/cavities/2rh1_1_cavityALL.mol2 -p data/proteins/2rh1_1_protein.mol2 -o data/desc/2rh1_1_desc.npy -obsa bsa/2rh1_1_bsa2.tsv
```
If BSA data already available
```python scripts/cavity_descriptors.py -c data/cavities/2rh1_1_cavityALL.mol2 -p data/proteins/2rh1_1_protein.mol2 -o data/desc/2rh1_1_desc.npy -ibsa bsa/2rh1_1_bsa2.tsv
```
### 2. Calculation of cavity labels
Ligand ph4 features and interactions with the target:

```python scripts/ligand_to_ph4.py --database --duplicate -l data/ligands.list -idir data/ligands -odir data/ligph4```

```IChem ints data/proteins/2rh1_1_protein.mol2 data/ligands/2rh1_1_ligand.mol2```
python scripts/cavity_labels.py -c data/cavities/2rh1_1_cavityALL.mol2 -lph data/ligph4/2rh1_1_ligph4.mol2 -int data/ints/2rh1_1_ints.mol2 -o data/labels/2rh1_1_labels.npy

Compute labels
python scripts/cavity_labels.py -c data/cavities/2rh1_1_cavityALL.mol2 -lph data/ligph4/2rh1_1_ligph4.mol2 -int data/ints/2rh1_1_ints.mol2 -o data/labels/2rh1_1_labels.npy


### 3. Prepare training test data
After generation for all the database, split into balanced training and tests sets
```python scripts/prepare_data.py -f data/scpdb_2022.list -d data/desc -l data/labels -a data/scpdb_2022_annotation.tsv --split random > prepare_data.log```
outputs:
- `features_split.json`: training and test sets of points
- `external_test_pdb.json`: leave-out entire cavities, for application
- `prepare_data.log`: report, statistics of data splitting into positive, negative and application set


### 4. train
```python ../sharing/scripts/train.py -t data/split_train_test/random/features_split.json -d data/descriptors -clf rf```
outputs:
- `<model>.report.json`: training and test sets of points
- `external_test_pdb.json`: leave-out entire cavities, for application
- `prepare_data.log`: report


### 5. Apply prediction on cavities (leave-out) and generate pruned cavities:
```for f in random kinase gpcr nuclear_receptor protease time; do
  echo predincting for $f xx;
  python scripts/application.py -v data/split_train_test/$f/external_test_pdb.json -c data/scpdb_2022/cavities_ALL --models models/rf/$f/CA.model models/rf/$f/CZ.model models/rf/$f/O.model models/rf/$f/OG.model models/rf/$f/OD1.model models/rf/$f/N.model models/rf/$f/NZ.model --modelnames CA CZ O OG OD1 N NZ -d data/descriptors -o data/pred_pharm/rf/$f/
done
```
Output:
- pruned cavities
- Figures: balanced accuracy, specificity, sensitivity, degree of prunning,


### 6. Analyze prediction on pruned cavities:
Example for the prediction of GPCR cavities 
python ../sharing/scripts/statistics_application.py -c data/cavities_ALL -pr pred_pharm/rf/gpcr/ -ap ../sharing/data/split_train_test/gpcr/external_test_pdb.json -d ../data/scpdb_2022/descriptors -l ../data/scpdb_2022/labels_using_ligph4_v5_duplicated -an scpdb_2022_annotation_plus.tsv



## Graph-based positioning of ligands to protein cavities:
python ../sharing/scripts/ligcare_graph.py -c 1zzl_1_cavityALL.mol2 -p 1zzl_1_protein.mol2 -l 1zzl_1_ligph4.mol2 -lig 1zzl_1_ligand.mol2
