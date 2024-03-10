# Spatial-Division

This repository contains an implementation of a cell-specific spatial lineage tracing phylogeographic reconstruction method called spalin (placeholder name). This code uses sc-mail (https://github.com/raphael-group/sc-mail/) to solve the likelihood of the sequence mutation portion of spatial lineage tracing data, and extends that codebase to solve the spatial likelihood and combine them together. 

In order to run this program, there are a few prerequisites:

## Prerequisite Libraries/Packages
1. This implementation requires MOSEK (a package which handles optimization problems): Please refer to the installation page for MOSEK at: https://www.mosek.com/products/academic-licenses/
2. Besides the above, the following Python packages are required: cmake, cvxpy, numpy, scipy, setuptools, treeswift.

   
## Example Command
The following command runs our method. Inside the problin/ directory:
```
$ python run_problin.py -c character_matrix.csv -t true_tree.nwk -S leaf_locations.txt --delimiter comma -p k10_priors.csv -o test_output.txt --nInitials 10 --ultrametric --divide --radius 5
```
which will output to test_output.txt the quantities of interest (branch lengths, spatial sigma, mutation rate lambda, and runtime). true_tree.nwk is only used for the tree topology, and can thus have any branch lengths. 

Examples of the input files character_matrix.csv, true_tree.nwk, leaf_locations.txt, k10_priors.csv are given in the examples.zip file inside problin/. We briefly describe each of these files below.

## true_tree.nwk
In the example data files provided, the header contains first the name of the cell, followed by an enumeration of the sites of the lineage tracing experiment (s1,s2,...)
each row of the matrix should have its first column be the name of the cell, and then the state at each of the sites for that cell. For example:

| cell  | s1 | s2 | s3 |
| ------------- | ------------- | ------------- | ------------- |
| cell_1  | 2  | 0  | 1  |
| cell_2  | 2  | 2  | 1  |
| cell_3  | 1 | 2  | 0  |

The states should be made into integers. 

## true_tree.nwk
A newick tree file containing the tree topology. The leaves of this tree should match what was provided in the character matrix. 

## leaf_locations.txt
Each line in this .txt file is corresponds to a cell, followed by its x and y coordinates, for example
```
cell_1,3.3,4.5
cell_2,5.3,2.1
cell_3,1.0,3.2
```
could be a valid leaf_locations.txt file for the character matrix given mentioned above. Omitting the spatial locations will instead have the software optimize only the sequence mutation likelihood. 

## prior.csv
Describes a prior for the transition rate matrix. Call be omitted if we assume a uniform prior. Otherwise, it takes the form of a headerless csv where each row is:
```
cell_id, state, probability 
```

## Additional Flags
```
--divide --radius #
```

This flag causes the software to optimize using the symmetric displacement model with a average radius amount equal to the value of #. Omitting these commands while still providing the spatial locations will instead have the model optimize a traditional Brownian motion likelihood. 
