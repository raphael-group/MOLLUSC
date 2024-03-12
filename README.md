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

