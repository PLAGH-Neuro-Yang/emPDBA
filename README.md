# **emPDBA: an ensemble model for Protein-DNA Binding Affinity prediction**



The emPDBA is a method for protein-DNA binding affinity prediction. This program can be run on Linux system. 

## Preparation
Before running the program, a PDB file of a protein-DNA complex with hydrogen atoms is required as input, which can be obtained through AutoDock or other softwares. Additionally, other softwares and packages are needed to be installed to ensure the program runs properly. 

emPDBA uses the following dependencies:
- python 3.8
- numpy
- sklearn
- xgboost
- joblib
- vmd

## Example
A protein-DNA complex with PDB file ‘5EGB.pdb’, is used as an example to show the prediction process. Hydrogen atoms are added by using AutoDock.

## How to run
**Step 1:** Put 5EGB.pdb, emPDBA.py, and code and feature folders in the same directory.

**Step 2:** Use the following bash command to run emPDBA.
		
```
python emPDBA.py [pdbfile] [type]
```

for example:
		
```
python emPDBA.py 5EGB.pdb D
```

The parameter **[pdbfile]** is the PDB file name of the protein-DNA complex structure.

The parameter **[type]** is the complex structure type. ‘S’ stands for the complexes with single-stranded DNAs, ‘D’ represents the complexes with double-stranded ones (emPDBA.py will classify this type into Double I, Double II and Double III based on the percentage of interface residues in protein) and ‘M’ represents the miscellaneous complexes.

Then, emPDBA.py program will perform feature extraction and prediction process, which will take a few seconds.

The prediction result of the input protein-DNA complex will be shown in the output of emPDBA.py program, i.e. on bash command interface.

## Help
For any questions, please contact us by yangshuang@emails.bjut.edu.cn or chunhuali@bjut.edu.cn.
