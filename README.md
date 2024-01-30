# Spatial-Division

This repository contains an implementation of a cell-specific spatial lineage tracing phylogeographic reconstruction method called spalin (placeholder name). This code uses sc-mail (https://github.com/raphael-group/sc-mail/) to solve the likelihood of the sequence mutation portion of spatial lineage tracing data, and extends that codebase to solve the spatial likelihood and combine them together. 

In order to run this program, there are a few prerequisites:

## Prerequisite Libraries/Packages
1. This implementation requires MOSEK (a package which handles optimization problems): Please refer to the installation page for MOSEK at: https://www.mosek.com/products/academic-licenses/
2. Besides the above, the following Python packages are required: cmake, cvxpy, numpy, scipy, setuptools, treeswift.

