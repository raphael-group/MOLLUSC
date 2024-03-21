# MOLLUSC (Maximum Likelihood Estimation Of Lineage and Location Using Single Cell Spatial Lineage tracing Data)

This repository contains an implementation of MOLLUSC, a cell-specific spatial lineage tracing phylogeographic reconstruction method.

## Prerequisite Libraries/Packages
The following Python packages are required: cmake, cvxpy, numpy, scipy, setuptools, treeswift.

## Usage
A description of all the other input flags (as well as the ones described above for the different data input modalities) can be seen by running:
```
python run_mollusc.py --help
```

## Input File Descriptions
The following are the types of files one needs to run this method

- priors.csv (optional)
- character matrix file
- leaf_locations.txt
- true_tree.nwk

You can see examples of these files in the `example` folder.

### character matrix file
The first row is a header describing the contents of the subsequent matrix. The first column of the matrix is the name of a cell (in this case, 216_* refers to the fact that the leaves of the tree correspond to the 216th time sample of the video frame data from the intMEMOIR experiment). The remaining columns correspond to a given mutation site across the cells. In other words, for a given row of this matrix, the first entry is the name of the cell, the remaining  entries correspond to state the mutation sites. 

| cell_name  | site_1 | site_2 |
| ------------- | ------------- | ------------- |
| cell_1  | 1  | 2  |
| cell_2  | 0  | 2  |
| cell_3  | 1  | 0  |

### tree topology (newick file)
This file contains the tree topology for this experiment, in newick format. The labels on the leaf nodes must match the cell names specified in the character matrix.

### leaf locations file
This file contains the X-Y coordinates for each of the cells at the leaves of the tree. These should correspond to the cells in the character matrix. Each line of this file should be of the form:

```
cell_name, x_coordinate, y_coordinate
```

### priors.csv (optional)
The file contains prior mutation probabilities that can be provided to our model for the Q matrix of the PMM model (see the original paper for more details). Each site and each mutated state should have a row in this matrix, structured as:

```
site_number, state, probability 
```
The prior is not required if one assumes a uniform distribution for mutation states.

## Examples
The software can run in one of the following three main modes: (1) The sequence only model, (2) Sequence + Location models jointly, and (3) Location only. 
For each of these modalities we show an example command with the input files in the `example` directory.

### Sequence Only Model

To run this method using only the sequence mutation model, input the character matrix via `-c` and the tree topology via `-t`. For example: 
```
python run_mollusc.py -c example/character_matrix.csv -t example/true_tree.nwk --delimiter comma -p example/k10_priors.csv --nInitials 1 --randseeds 3103 -o example/sequence_only_example.txt --timescale 215 -v
```

### Sequence + Location
To run the method using both forms of data, include the spatial locations via `-S` in addition to `-c` and `-t` as described above. Furthermore, if you include the flags --divide & --radius, then the symmetric displacement model will be used with radius amount equal to the number after the --radius flag. For example, to run symmetric displacement with radius = 5, use
```
python run_mollusc.py -c example/character_matrix.csv -t example/true_tree.nwk --delimiter comma -p example/k10_priors.csv --nInitials 1 --randseeds 3103 -o example/sym_displacement_example.txt --timescale 215 -S example/leaf_locations.txt --divide --radius 5 -v
```

To run with Brownian motion, you can either set --radius to be 0, or omit the --divide flag entirely: 
```
python run_mollusc.py -c example/character_matrix.csv -t example/true_tree.nwk --delimiter comma -p example/k10_priors.csv --nInitials 1 --randseeds 3103 -o example/brownian_example.txt --timescale 215 -S example/leaf_locations.txt -v
```

### Location Only
You can use the location only model by including the flag --spatial_only in conjunction on top of the flags mentioned above. 
For example, to only use spatial data with symmetric displacement and radius = 5, use
```
python run_mollusc.py -c example/character_matrix.csv -t example/true_tree.nwk --delimiter comma -p example/k10_priors.csv --nInitials 1 --randseeds 3103 -o example/sym_displacement_spatialonly_example.txt --timescale 215 -S example/leaf_locations.txt --divide --radius 5 --spatial_only -v
```
If you instead want to use the Brownian motion model, use
```
python run_mollusc.py -c example/character_matrix.csv -t example/true_tree.nwk --delimiter comma -p example/k10_priors.csv --nInitials 1 --randseeds 3103 -o example/brownian_spatialonly_example.txt --timescale 215 -S example/leaf_locations.txt --spatial_only -v
```

