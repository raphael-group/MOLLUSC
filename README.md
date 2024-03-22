# MOLLUSC 
(Maximum Likelihood Estimation Of Lineage and Location Using Single Cell Spatial Lineage tracing Data)

This repository contains an implementation of MOLLUSC, a cell-specific spatial lineage tracing phylogeographic reconstruction method.

## Prerequisite Libraries/Packages
The following Python packages are required: cmake, cvxpy, numpy, scipy, setuptools, treeswift.

## Usage
A description of all the other input flags (as well as the ones described above for the different data input modalities) can be seen by running:
```
python run_mollusc.py --help
```

## Input File Descriptions
The following are the required inputs to MOLLUSC
- Character matrix
- Tree topology
- Leaf locations
- [optional] Prior mutation probabilities (associated with the input character matrix)

These input files must follow the formats described below.

### Character matrix
The first row is a header describing the contents of the subsequent matrix. The first column of the matrix is the name of a cell (in this case, 216_* refers to the fact that the leaves of the tree correspond to the 216th time sample of the video frame data from the intMEMOIR experiment). The remaining columns correspond to a given mutation site across the cells. In other words, for a given row of this matrix, the first entry is the name of the cell, the remaining  entries correspond to state the mutation sites. 

| cell_name  | site_1 | site_2 |
| ------------- | ------------- | ------------- |
| cell_1  | 1  | 2  |
| cell_2  | 0  | 2  |
| cell_3  | 1  | 0  |

See [example/character_matrix.csv](example/character_matrix.csv) for an example.

### Tree topology
The input tree topology must be given in [newick format](https://en.wikipedia.org/wiki/Newick_format#:~:text=In%20mathematics%2C%20Newick%20tree%20format,Maddison%2C%20Christopher%20Meacham%2C%20F.). 
In addition, the labels on the leaf nodes must match the cell names specified in the character matrix.

See [example/input_tree.nwk](example/input_tree.nwk) for an example.

### Leaf locations
This file contains the X-Y coordinates for each of the cells at the leaves of the tree. These should correspond to the cells in the character matrix. Each line of this file should be of the form:

```
cell_name, x_coordinate, y_coordinate
```
See [example/input_locations.txt](example/input_locations.txt) for an example.

### Prior mutation probabilities (optional)
The file contains prior mutation probabilities that can be provided to our model for the Q matrix of the PMM model (see the original paper for more details). Each site and each mutated state should have a row in this matrix, structured as:

```
site_number, state, probability 
```
The prior is not required if one assumes a uniform distribution for mutation states.
See [example/mutation_priors.csv](example/mutation_priors.csv) for an example.

## Examples
The software can run in one of the following three main modes: (1) The sequence only model, (2) Sequence + Location models jointly, and (3) Location only. 
For each of these modalities we show an example command with the input files in the `example` directory.

### Sequence Only Model

To run this method using only the sequence mutation model, input the character matrix via `-c` and the tree topology via `-t`. For example: 
```
python run_mollusc.py -c example/character_matrix.csv -t example/input_tree.nwk --delimiter comma -p example/mutation_priors.csv --nInitials 1 --randseeds 3103 -o sequence_only_example.txt --timescale 215 -v
```
The above command will produce an output file `sequence_only_example.txt`. An example output is provided in `example/use_case_1/sample_output.txt`.

### Sequence + Location
To run the method using both forms of data, include the spatial locations via `-S` in addition to `-c` and `-t` as described above. Furthermore, if you include the flags --divide & --radius, then the symmetric displacement model will be used with radius amount equal to the number after the --radius flag. For example, to run symmetric displacement with radius = 5, use
```
python run_mollusc.py -c example/character_matrix.csv -t example/input_tree.nwk --delimiter comma -p example/mutation_priors.csv --nInitials 1 --randseeds 3103 -o sym_displacement_example.txt --timescale 215 -S example/input_locations.txt --divide --radius 5 -v
```
The above command will produce an output file `sym_displacement_example.txt`. An example output is provided in `example/use_case_2/sample_output.txt`.

To run with Brownian motion, you can either set --radius to be 0, or omit the --divide flag entirely: 
```
python run_mollusc.py -c example/character_matrix.csv -t example/input_tree.nwk --delimiter comma -p example/mutation_priors.csv --nInitials 1 --randseeds 3103 -o browninan_example.txt --timescale 215 -S example/input_locations.txt --divide --radius 5 -v
```
The above command will produce an output file `brownian_example.txt`. An example output is provided in `example/use_case_3/sample_output.txt`.

### Location Only
You can use the location only model by including the flag --spatial_only in conjunction on top of the flags mentioned above. 
For example, to only use spatial data with symmetric displacement and radius = 5, use
```
python run_mollusc.py -c example/character_matrix.csv -t example/input_tree.nwk --delimiter comma -p example/mutation_priors.csv --nInitials 1 --randseeds 3103 -o sym_displacement_spatialonly_example.txt --timescale 215 -S input_locations.txt --divide --radius 5 --spatial_only -v
```
The above command will produce an output file `sym_displacement_spatialonly_example.txt`. An example output is provided in `example/use_case_4/sample_output.txt`.

If you instead want to use the Brownian motion model, use
```
python run_mollusc.py -c example/character_matrix.csv -t example/input_tree.nwk --delimiter comma -p example/mutation_priors.csv --nInitials 1 --randseeds 3103 -o sym_displacement_spatialonly_example.txt --timescale 215 -S input_locations.txt --spatial_only -v
```
The above command will produce an output file `brownian_spatialonly_example.txt`. An example output is provided in `example/use_case_5/sample_output.txt`.
