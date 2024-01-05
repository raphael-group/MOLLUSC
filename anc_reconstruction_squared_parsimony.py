import sys
import treeswift
import numpy as np
import os
import math
import warnings
from scipy.stats import norm
from scipy import optimize
from sys import argv
from math import sin,cos, sqrt

def get_tree(file_name):
	#input: file name

	f = open(file_name, 'r')
	data = f.read()
	newick_tree_start = data.find('[&R] ') + 5
	new_tree_end = data.find(';', newick_tree_start)
	tree_string = data[newick_tree_start:new_tree_end+1]
	# print(tree_string)

	newick_tree = treeswift.read_tree(tree_string, 'newick')
	f.close()
	return newick_tree

def normalize_tree(tree):
	height = tree.height()
	tree.scale_edges(1/height)
	return tree

def initialize_internals(length):
	x_avg = sum(value[0] for value in leaf_location_dict.values()) / len(leaf_location_dict)
	y_avg = sum(value[1] for value in leaf_location_dict.values()) / len(leaf_location_dict)

	x_vec = np.random.normal(x_avg, 5, length)
	y_vec = np.random.normal(y_avg, 5, length)
	initial = []
	for i in range(length):
		initial.append(x_vec[i])
		initial.append(y_vec[i])
	return initial


def optimize_internal_locations():
	# initialize the internal locations
	count = 0
	success = False
	best_result = []
	best_nllh = np.inf
	while (count < 30):
		print("iteration",count)
		x0 = initialize_internals(num_internal) + [5] #last entry is the starting value of sigma
		out = optimize.minimize(get_nllh, x0)
		count+=1
		if out.success == True:
			best_nllh = out.fun
			best_result = out.x
	if len(best_result) == 0:
		print("error, hit max initials")
		return x0
	return best_result

def get_nllh(internal_locations_sigma):
	llh = 0
	try:
		current_sigma = internal_locations_sigma[-1]

		for node in current_tree.traverse_postorder():
			if node.is_leaf():
				x,y = leaf_location_dict[node.label]
			else:
				x = internal_locations_sigma[internal_mapping[node.label]]
				y = internal_locations_sigma[internal_mapping[node.label] + 1]

			if node.is_root(): # maximum likelihood of the root location is the root location
				llh += math.log(norm.pdf(x, loc=x, scale=sqrt(current_sigma**2 * node.edge_length)))
				llh += math.log(norm.pdf(y, loc=y, scale=sqrt(current_sigma**2 * node.edge_length)))
			else:
				parent_x = internal_locations_sigma[internal_mapping[node.parent.label]]
				parent_y = internal_locations_sigma[internal_mapping[node.parent.label] + 1]

				if displacement_dict != None:
					parent_x += displacement_dict[node.label][0]
					parent_y += displacement_dict[node.label][0]

				llh += math.log(norm.pdf(x, loc=parent_x, scale=sqrt(current_sigma**2 * node.edge_length)))
				llh += math.log(norm.pdf(y, loc=parent_y, scale=sqrt(current_sigma**2 * node.edge_length)))
	except:
		return np.inf

	return -llh


def get_sigma(file_name):
	f = open(file_name, 'r')
	data = f.read()
	start = data.find('Sigma: ') + 7
	end = data[start:].find('\n')

	return float(data[start:start+end-1])

def get_locations(file_name):
	locations = {}
	with open(file_name,'r') as fin: 
		for line in fin:
			cellID,x,y = line.strip().split(",")
			locations[cellID] = (float(x),float(y))
	return locations

def get_displacement_amounts(file_name, tree):
	displacements = {}
	inThetas = False

	thetas = {}
	with open(file_name,'r') as fin: 
		for line in fin:
			if inThetas == True:
				cellID,theta = line.strip().split(" ")
				# fix up formatting (Getting rid of parenthesis and commas)
				thetas[cellID] = float(theta)
			if "Thetas" in line:
				inThetas = True

	if inThetas == False:
		return None

	for node in tree.traverse_preorder():
		sum_of_displacement = 0
		if node.parent == None:
			continue
		left_node = node.parent.child_nodes()[0]
		multiplier = -1
		if node.label == left_node: # was a left child, do addition
			multiplier = 1

		x_disp = multiplier * 5 * math.cos(thetas[node.parent.label])
		y_disp = multiplier * 5 * math.sin(thetas[node.parent.label])
		displacements[node.label] = (x_disp, y_disp)

	return displacements

def save_internal_locations(internal_locations_sigma, name_to_save):
	print("saving results")
	with open(name_to_save,'w') as fout: 
		for node in current_tree.traverse_postorder():
			if node.is_leaf():
				pass
			else:
				x_index = internal_mapping[node.label]
				y_index = x_index + 1
				x_loc = internal_locations_sigma[x_index]
				y_loc = internal_locations_sigma[y_index]
				fout.write(node.label + "," + str(x_loc) + "," + str(y_loc) + '\n')
		fout.write("sigma," + str(internal_locations_sigma[-1])+'\n')
	return


np.seterr(invalid='ignore')
current_tree = None
leaf_location_dict = {}
displacement_dict = None
internal_mapping = {}
num_internal = 0

leaf_location_file = argv[1]
brlen_sigma_results = argv[2]
output = argv[3]

#python anc_reconstruction_squared_parsimony.py /Users/gary/Documents/Projects/analyzing_spalin_data/real_location_data/s10c1/leaf_locations.txt /Users/gary/Documents/Projects/analyzing_spalin_data/real_location_data/s10c1/k10_r4_spalin-divide-beta5.txt /Users/gary/Documents/Projects/analyzing_spalin_data/real_location_data/s10c1/k10_r4_spalin-divide-beta5_locations.txt   

leaf_location_dict = get_locations(leaf_location_file)

file_name = brlen_sigma_results
current_tree = get_tree(file_name)

for node in current_tree.traverse_postorder():
	node.edge_length = 1

# only if we're using spalin-divide
displacement_dict = get_displacement_amounts(file_name, current_tree)


# construct mapping dictionary
if num_internal == 0: # if we haven't done this for the current tree
	index = 0
	for node in current_tree.traverse_postorder():
		if node.is_leaf() == False:
			internal_mapping[node.label] = index
			index+=2
			num_internal += 1

results = optimize_internal_locations()
save_internal_locations(results, output)


