# **emPDBA: an ensemble model for Protein-DNA Binding Affinity prediction**



The emPDBA program is a method for protein-DNA binding affinity prediction. This program can be run on Linux system. 

---


## Preparation
Before running the program, a PDB file of a protein-DNA complex with hydrogen atoms is required as input. The PDB file with hydrogen atoms added can be obtained through AutoDock or other software and algorithms. Also, several software and packages need to be installed to ensure the program runs properly. 

emPDBA uses the following dependencies:
- python 3.8
- numpy
- sklearn
- xgboost
- joblib
- networkx
- vmd


## Example
A protein-DNA complex PDB file, ‘5EGB.pdb’, is used as an example to show the process. Hydrogen atoms are added by using AutoDock.


## How to run
**Step 1:** Add hydrogen atoms to the PDB file provided through AutoDock or other software and agorithms. This step can be skipped when using the example complex ‘5EGB.pdb’.

**Step 2:** Put the PDB file with hydrogen atoms (e.g. 5EGB.pdb), emPDBA.py, and code and feature folders in the same directory.

**Step 3:** Use the **MIBPB** tool (https://weilab.math.msu.edu/MIBPB/) [*Theoretical Chemistry Accounts*. 2017; 136(4):55.] to calculate the volume and reaction field energy of complex, protein and DNA. The data are stored in the 'volume.data' and 'reaction_field_energy_file.data' files in the order of complex, protein and DNA, separeted by spaces. Put the 'volume.data' and 'reaction_field_energy_file.data' files in the feature folder.

**Step 4:** Use the following bash command to run **emPDBA**.
		
```
python emPDBA.py [pdbfile] [type]
```

> for example:
```
python emPDBA.py 5EGB.pdb D
```


> The parameter [pdbfile] should be the PDB file name of the protein-DNA complex structure.

> The parameter [type] should be the DNA structure type. ‘D’ for double-stranded DNA-protein complex (emPDBA.py will classify this type into Double I, Double II and Double III based on the percentage of interface residues in protein) and ‘M’ for the miscellaneous complex. Please enter 'D' if the DNA in the user-supplied complex is double-stranded, otherwise enter 'M'.

Then, emPDBA program will perform feature extraction and prediction process, which will take a few seconds.

The prediction result of the input protein-DNA complex will be shown in the output of emPDBA.py program, i.e. on bash command interface.


---

## Help
For any questions, please contact us by yangshuang@emails.bjut.edu.cn or chunhuali@bjut.edu.cn.
