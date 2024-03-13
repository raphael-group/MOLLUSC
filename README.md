# Spatial-Division

This repository contains an implementation of a cell-specific spatial lineage tracing phylogeographic reconstruction method called spalin (placeholder name). This code uses sc-mail (https://github.com/raphael-group/sc-mail/) to solve the likelihood of the sequence mutation portion of spatial lineage tracing data, and extends that codebase to solve the spatial likelihood and combine them together. Please see that repository for more information on the possible input flags and options.

In order to run this program, there are a few prerequisites:

## Prerequisite Libraries/Packages
1. This implementation requires MOSEK (a package which handles optimization problems): Please refer to the installation page for MOSEK at: https://www.mosek.com/products/academic-licenses/
2. Besides the above, the following Python packages are required: cmake, cvxpy, numpy, scipy, setuptools, treeswift.

   
## Example Command 
The following command runs our method with the symmetric displacement model. Inside the problin/ directory:
```
$ python run_problin.py -c character_matrix.csv -t true_tree.nwk -S leaf_locations.txt --delimiter comma -p k10_priors.csv -o test_output.txt --nInitials 10 --ultrametric --divide --radius 5
```
which will output to test_output.txt the quantities of interest (branch lengths, spatial sigma, mutation rate lambda, and runtime). true_tree.nwk is only used for the tree topology, and can thus have any branch lengths. 

## Data Modalities Utilized 
There are three main modes that this method can perform: (1) The sequence only model, (2) Sequence + Location models jointly, and (3) Location only 


### Sequence Only Model
To run this method using only the sequence mutation model, simply omit the --spatial (-S) input flag. For example: 
```
$ python run_problin.py -c character_matrix.csv -t true_tree.nwk  --delimiter comma -p k10_priors.csv -o test_output.txt --nInitials 10 --ultrametric
```

### Sequence + Location
To run the method using both forms of data, include the spatial input flag. Furthermore, if you include the flags --divide & --radius, then the symmetric displacement model will be used with radius amount equal to the number after the --radius flag. For example, to run symmetric displacement with radius = 5, use
```
$ python run_problin.py -c character_matrix.csv -t true_tree.nwk -S leaf_locations.txt --delimiter comma -p k10_priors.csv -o test_output.txt --nInitials 10 --ultrametric --divide --radius 5
```

to run with Brownian motion, you can either set --radius to be 0, or omit the --divide flag entirely: 
```
$ python run_problin.py -c character_matrix.csv -t true_tree.nwk -S leaf_locations.txt --delimiter comma -p k10_priors.csv -o test_output.txt --nInitials 10 --ultrametric  --divide --radius 0
```

### Location Only
You can use the location only model by including the flag --spatial_only in conjunction on top of the flags mentioned above, for example:
```
$ python run_problin.py -c character_matrix.csv -t true_tree.nwk -S leaf_locations.txt --delimiter comma -p k10_priors.csv -o test_output.txt --nInitials 10 --ultrametric --divide --radius 5 --spatial_only
```


## Input File Descriptions
The following are the types of files one needs to run this method

- priors.csv (optional)
- character matrix file
- leaf_locations.txt
- true_tree.nwk

### priors.csv (optional)
(1) The file contains an example of a prior file that can be provided to our model for the Q matrix (see PMM paper). Note this file is not required if one assumes a uniform distribution for this prior. Each site and each mutated state should have a row in this matrix, structured as:

```
site_number, state, probability 
```

### character matrix file
(2) The first row is a header describing the contents of the subsequent matrix. The first column of the matrix is the name of a cell (in this case, 216_* refers to the fact that the leaves of the tree correspond to the 216th time sample of the video frame data from the intMEMOIR experiment). The remaining columns correspond to a given mutation site across the cells. In other words, for a given row of this matrix, the first entry is the name of the cell, the remaining  entries correspond to state the mutation sites. 

| cell_name  | site_1 | site_2 |
| ------------- | ------------- | ------------- |
| cell_1  | 1  | 2  |
| cell_2  | 0  | 2  |
| cell_3  | 1  | 0  |

### leaf locations file 
(3) This file contains the X-Y coordinates for each of the cells at the leaves of the tree. These should correspond to the cells in the character matrix in (1). Each line of this file should be of the form:

```
cell_name, x_coordinate, y_coordinate
```

### true tree (newick file)
(4) This file contains the true tree topology for this experiment, in newick format. the labels on the nodes should match the aforementioned files. 

You can see examples of these files in the included zip file, as well as in the data repository. 
