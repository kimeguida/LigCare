# LigCare: Processing, comparison and alignment of protein pharmacophoric points to ligand pharmacophores
[![Generic badge](https://img.shields.io/badge/version-0.1.0-orange.svg)](https://shields.io/)

The goal of this project is to investigate the comparability of protein pockets represented as negative image to ligands. <br>

Machine learning models were developed to automatically learn important features and simplify the cloud of points representing a protein binding site, without a priori knowledge of ligand binding area. Alignment of ligands to protein sites were performed by point cloud registration or by searching graph isomorphism. <br>

Please note that this repository is in beta version.


## Requirements
* [IChem](http://bioinfo-pharma.u-strasbg.fr/labwebsite/download.html)
* [Conda](https://docs.conda.io/en/latest/miniconda.html)
* Conda environment: [ligcare.yml](https://github.com/kimeguida/LigCare/blob/main/ligcare.yml)
* Dataset: [sc-PDB](http://bioinfo-pharma.u-strasbg.fr/scPDB/) v.2022 (to be published data)


## Content
* scripts
* data/: input and output data organization with examples
* models: storage of trained models (example), all trained models are in models.tgz archive
	

## Help and issues
`<script>.py --help` for options. <br>
To signal issues: https://github.com/kimeguida/LigCare/issues




## Prediction of interacting points in protein pockets

### 1. Computation of cavity descriptors
```python scripts/cavity_descriptors.py -c data/cavities/2rh1_1_cavityALL.mol2 -p data/proteins/2rh1_1_protein.mol2 -o data/desc/2rh1_1_desc.npy -obsa bsa/2rh1_1_bsa.tsv``` <br>

Computing the buriedness BSA is computationally clostly. If BSA data (.bsa files) are already available:
```python scripts/cavity_descriptors.py -c data/cavities/2rh1_1_cavityALL.mol2 -p data/proteins/2rh1_1_protein.mol2 -o data/desc/2rh1_1_desc.npy -ibsa bsa/2rh1_1_bsa.tsv```

### 2. Computation of cavity labels

Ligand pharmacophoric features (ph4) and interactions with the target:

```python scripts/ligand_to_ph4.py --database --duplicate -l data/ligands.list -idir data/ligands -odir data/ligph4```

```IChem ints data/proteins/2rh1_1_protein.mol2 data/ligands/2rh1_1_ligand.mol2```

```python scripts/cavity_labels.py -c data/cavities/2rh1_1_cavityALL.mol2 -lph data/ligph4/2rh1_1_ligph4.mol2 -int data/ints/2rh1_1_ints.mol2 -o data/labels/2rh1_1_labels.npy```

Training labels of cavity points for classification:
0: non-interacting
1: interacting
```python scripts/cavity_labels.py -c data/cavities/2rh1_1_cavityALL.mol2 -lph data/ligph4/2rh1_1_ligph4.mol2 -int data/ints/2rh1_1_ints.mol2 -o data/labels/2rh1_1_labels.npy```


### 3. Preparation of training test data

After generation of descriptors and labels for the entire database, split into application, balanced training and test sets:
```python scripts/prepare_data.py -f data/scpdb_2022.list -d data/desc -l data/labels -a data/scpdb_2022_annotation.tsv --split random > prepare_data.log```
outputs:
- `features_split.json`: training and test sets of points
- `external_test_pdb.json`: leave-out entire cavities, for application
- `prepare_data.log`: statistics of data splitting into positive, negative and application set

Several splitting schemes are implemented for left-out cavities: random, GPCR, time, kinase, protease, nuclear receptor.

### 4. Training models
Example for random left-out. Random Forest (rf) or XGBoost (xgb) classifiers were tested. Each of the seven ph4 types (CA hydrophobic, CZ aromatic, O h-bond acceptor, OG h-bond acceptor and donor, OD1 negative ionizable, N h-bond donor, NZ positive ionizable) is trained individually.
```python ../sharing/scripts/train.py -t data/split_train_test/random/features_split.json -d data/descriptors -clf rf```
outputs:
- `<ph4>.report`: statistics of cross-validation, training and external tests
- `<ph4>.model`: trained model binary, can be loaded with pickle

Further analyses are needed to assess the models.

### 5. Application of trained models on left-out cavities and generate pruned cavities
Example for random left-out.
```python scripts/application.py -v data/split_train_test/random/external_test_pdb.json -c data/scpdb_2022/cavities --models models/rf/$f/CA.model models/rf/random/CZ.model models/rf/random/O.model models/rf/random/OG.model models/rf/random/OD1.model models/rf/random/N.model models/rf/random/NZ.model --modelnames CA CZ O OG OD1 N NZ -d data/descriptors -o data/pred_pharm/rf/random/``` <br>

For all left-out.
```for f in random kinase gpcr nuclear_receptor protease time; do echo predicting for $f............; python scripts/application.py -v data/split_train_test /$f/external_test_pdb.json -c data/scpdb_2022/cavities --models models/rf/$f/CA.model models/rf/$f/CZ.model models/rf/$f/O.model models/rf/$f/OG.model models/rf/$f/OD1.model models/rf/$f/N.model models/rf/$f/NZ.model --modelnames CA CZ O OG OD1 N NZ -d data/descriptors -o data/pred_pharm/rf/$f/ done```

outputs:
- pruned cavities: only points precticted to be interacting are kept in the pruned cavities.
- figures: balanced accuracy, specificity, sensitivity, degree of prunning


### 6. Analysis of predictions: pruned cavities
Example for the prediction of GPCR left-out cavities 
```python scripts/statistics_application.py -c data/cavities -pr data/pred_pharm/rf/gpcr/ -ap data/split_train_test/gpcr/external_test_pdb.json -d ../data/descriptors -l ../data/labels -an scpdb_2022_annotation.tsv```




## Graph-based comparison and alignment of ligands to protein cavities

Transformation of ligands/molecules into pharmacophoric (ph4) features:
```python scripts/ligand_to_ph4.py --database --duplicate -l data/ligands.list -idir data/ligands -odir data/ligph4``` <br>

Align by searching common substructures in the ph4 features:
* *Using IChem VolSite grid-sampled cavities*
`python scripts/ligcare_graph.py -c data/cavities/2rh1_1_cavityALL.mol2 -p data/proteins/2rh1_1_protein.mol2 -l data/ligph4/2rh1_1_ligph4.mol2 -lig data/ligands/2rh1_1_ligand.mol2`

* *Using projected irregular protein features*
Compute protph4 features...
`python scripts/protein_to_ph4.py -p data/proteins/2rh1_1_protein.mol2 -c data/cavities/2rh1_1_cavityALL.mol2 -o data/protph4/2rh1_1_photph4.mol2`
... and compare
`python scripts/ligcare_graph.py -c data/protph4/2rh1_1_photph4.mol2 -p data/proteins/2rh1_1_protein.mol2 -l data/ligph4/2rh1_1_ligph4.mol2 -lig data/ligands/2rh1_1_ligand.mol2` <br>

outputs:
- `rot_<ligand_name>_x.mol2`: ligand poses, x = {1, 2, 3, etc.}
- `rot_<ligph4_name>_x.mol2`: ligand ph4 poses, x = {1, 2, 3, etc.}

By default, the top 20 solutions ranked by best RMSE of the aligned feature points are output.


## Point cloud-based comparison and alignment of ligands to protein cavities
Compute the ligvoxelplus representation of the ligand/molecule:
`python scripts/ligvoxelplus.py -i data/ligands/2rh1_1_ligand.mol2 -o data/ligvoxelplus/2rh1_1_ligvoxelplus.mol2` <br>

Align by point cloud registration with [ProCare](https://github.com/kimeguida/ProCare) (see requirements and usage) and apply roto/translation to the ligand atoms:
`(procare) $ python procare_launcher.py -t  data/cavities/2rh1_1_cavityALL.mol2 -s data/ligvoxelplus/2rh1_1_ligvoxelplus.mol2 --transform --ligandtransform data/ligands/2rh1_1_ligand.mol2`


## References

- Eguida M, 2022. Comparison of protein cavities by point cloud processing: principles and applications in drug design. [PhD thesis](https://www.theses.fr/s269955)
- Desaphy, J.; Azdimousa, K.; Kellenberger, E.; Rognan, D. Comparison and Druggability Prediction of Protein–Ligand Binding Sites from Pharmacophore-Annotated Cavity Shapes. J. Chem. Inf. Model. 2012, 52 (8), 2287–2299. [https://doi.org/10.1021/ci300184x](https://doi.org/10.1021/ci300184x)
- Da Silva, F.; Desaphy, J.; Rognan, D. IChem: A Versatile Toolkit for Detecting, Comparing, and Predicting Protein–Ligand Interactions. ChemMedChem 2018, 13 (6), 507–510. [https://doi.org/10.1002/cmdc.201700505](https://doi.org/10.1002/cmdc.201700505)
- Eguida, M., Rognan, D. A Computer Vision Approach to Align and Compare Protein Cavities: Application to Fragment-Based Drug Design. J. Med. Chem. 2020. [https://doi.org/10.1021/acs.jmedchem.0c00422](https://doi.org/10.1021/acs.jmedchem.0c00422)
