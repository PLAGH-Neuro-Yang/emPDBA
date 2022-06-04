# **emPDBA: an ensemble model for Protein-DNA Binding Affinity prediction**



The emPDBA program is a method for protein-DNA binding affinity prediction. This program can be run on Linux system. 

---


## Preparation
Before running the program, a PDB file of a protein-DNA complex with hydrogen atoms is required as input. The PDB file with hydrogen atoms added can be obtained through AutoDock or other software and algorithms. Also, several software and packages need to be installed to ensure the program runs properly. 

emPDBA uses the following dependencies:
- python 3.8
- numpy
- joblib
- vmd


## Example
A protein-DNA complex PDB file, ‘5EGB.pdb’, is used as an example to show the process. This PDB file is the experimental structure of the protein-DNA complex with PDB ID: 5EGB. Hydrogen atoms are added by using AutoDock.


## How to run
**Step 1:** Put the PDB file and emPDBA.py in the same directory with code and feature folders.

**Step 2:** Use the following bash command to run emPDBA.
		
```
python emPDBA.py [pdbfile] [type]
```

> for example:
```
python emPDBA.py 5EGB.pdb D
```


> The parameter [pdbfile] should be the PDB file name of the protein-DNA complex structure.

> The parameter [type] should be the DNA structure type. ‘S’ for single-stranded DNA-protein complex, with ‘D’ for double-stranded ones and ‘M’ for the miscellaneous complex.

Then, emPDBA program will perform feature extraction and prediction process, which will take a few seconds.

The prediction result of the input complex will be shown as the output of emPDBA.py program, or it can be seen on bash command interface.


---

## Help
For any questions, please contact us by yangshuang@emails.bjut.edu.cn or chunhuali@bjut.edu.cn.
