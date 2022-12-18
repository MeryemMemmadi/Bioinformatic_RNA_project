# Bioinformatic of RNA Project - RNA folding problem
## Compute an estimation of the Gibbs free energy


### Execute Code

Clone the git in your computer

```bash
$ git clone https://github.com/MeryemMemmadi/Bioinformatic_RNA_project
```

Execute the program for the 4GXY.pdb
```bash
$  python3 RNA_script3.py -f 4gxy.pdb
```

### Guidelines
For a given ribonucleotide chain, the RNA folding problem consists in finding the native fold
among the astronomically large number of possible conformations. The native fold being the
one with the lowest Gibbs free energy, the objective function should be an estimator of this
energy.


To answer this problematic, 3 scripts were developped :

#### Trainning script 

``RNA_script1.py`` is a python script that compute interatomic distances from the given dataset of PDB files. only C3’ atoms, “intrachain” distances and residues separated 
by at least 3 positions on the sequence are taked into account.

The observed frequencies of two residues i and j separated by a distance bin r is calculated as follows : 

$$ f_{i,j} ^{OBS}(r) = { N_{i,j}(r) / N_{i,j} } $$

where \N_{i,j}(r) is the count of i and j within the distance bin r, and *N_(i,j)* is the count of i andj for all distance bins.

The reference frequency is the same formula, except that the different residue types (A, U, C,G) are indistinct (“X”):

$$ f_{X,X} ^{REF}(r) = { N_{X,X}(r) / N_{X,X}} $$

Finally, the score (pseudo-energy)ūi,j(r) is computed as follows:

$$ u_{i,j}(r) = { -log \left( f _{i,j} ^{OBS}(r) / f_{i,j} ^{REF}(r) \right) } $$

The training script generate 10 files of 20 lines (1 line = 1 scoring value).The maximum scoring value is arbitrarily set to 10.


#### Plotting script

``RNA_script2.r`` is a  script R that plot the scoring profiles;


#### Compute the estimated Gibbs free energy

``RNA_script3.py`` is a python script that use the objective function to evaluate predicted structures from the RNA-Puzzles dataset. This script is partially 
similar to the first one, as it will compute all the distances for a given structure (same thresholds: 20 Å and i, i+4). For each distance, a scoring value will
be computed, using a linear interpolation :

$$ G_e = y_1 + (d - x_1) * {(y_2-y_1) \over (x_2-x_1)}$$

By summing all these scores, the script will calculate the estimated Gibbs free energy of the evaluated RNA conformation.

---


### Data
The directory liste_pdb contains all the pdb used as a trainning set for the function that estimate the Gibbs free energy.
The Downloaded PDB for which you want to compute the Gibbs free energy have to be in the same directory as the 3 scripts.
You can get the PDB files on the [Protein Databank Website](https://www.rcsb.org/).


---
### Code

#### Script 1 
Execute the first script with python and give the directory where the PDB used as a trainnning set are stored.

```bash
$ python3 RNA_script1.py -path liste_pdb/
```
generate as an output 10 files (one for each pair) of 20 lines (1 line = 1 scoring value).

#### Script 2 
Execute the second script on R to get the different plot of the scores

```bash
$ python3 RNA_script2.py
```
Generate as output 10 plot of the scores of each pair.

#### Script 3 
Execute the last script with python with the pdb file for which you want to estimate the Gibbs free energy, as an argument. 

```bash
$ python3 RNA_script3.py -f XXXX.pdb
```



