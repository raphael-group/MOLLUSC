import sys
import treeswift
import numpy as np
import os
import math
import warnings
from scipy.stats import multivariate_normal
from scipy import optimize
from sys import argv
from math import sin,cos

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
	while (count < 100):
		print("iteration",count)
		x0 = initialize_internals(num_internal) + [5] #last entry is the starting value of sigma
		out = optimize.minimize(get_multivariate_normal_nllh, x0)
		count+=1
		if out.success == True:
			best_nllh = out.fun
			best_result = out.x
			break
	if len(best_result) == 0:
		print("error, hit max initials")
		return x0
	return best_result

def get_multivariate_normal_nllh(internal_locations_sigma):
	x_mean_vec, y_mean_vec = construct_mean_vector(internal_locations_sigma[0:-1])
	cov_matrix = construct_covariance_matrix(internal_locations_sigma[-1])

	all_locations_x = []
	all_locations_y = []
	for node in current_tree.traverse_postorder():
		if node.is_leaf():
			all_locations_x.append(leaf_location_dict[node.label][0])
			all_locations_y.append(leaf_location_dict[node.label][1])
		else:
			x_index = internal_mapping[node.label]
			y_index = x_index + 1
			all_locations_x.append(internal_locations_sigma[x_index])
			all_locations_y.append(internal_locations_sigma[y_index])

	x_llh = multivariate_normal.pdf(all_locations_x, mean=x_mean_vec, cov=cov_matrix)
	y_llh = multivariate_normal.pdf(all_locations_y, mean=y_mean_vec, cov=cov_matrix)
	if x_llh == 0 or y_llh == 0:
		return np.inf
	return - (math.log(x_llh) + math.log(y_llh))

def dist_from_root(node):
	# for some reason the distance_between() function of treeswift is throwing some error..
	dist = 0
	for p_node in node.traverse_ancestors():
		dist += p_node.edge_length
	return dist

def construct_covariance_matrix(sigma):
	# return an n x n matrix, where n is the number of nodes in the tree
	# entry ij is sigma * d(root, lca(i,j))

	# find the root of the tree:
	root = None
	num_nodes = 0
	cov_matrix = []
	for node in current_tree.traverse_preorder():
		num_nodes += 1
		if node.is_root():
			root = node

	for node_i in current_tree.traverse_postorder():
		entry = [0] * num_nodes
		j_index = 0
		for node_j in current_tree.traverse_postorder():

			lca = current_tree.mrca( {node_i.label, node_j.label} )
			entry[j_index] = sigma**2 * dist_from_root(lca)
			j_index += 1
		cov_matrix.append(entry)
	return cov_matrix

def construct_mean_vector(internal_locations):
	mean_vector_x = []
	mean_vector_y = []
	for node in current_tree.traverse_postorder():
		if node.is_root():
			x_index = internal_mapping[node.label]
			y_index = x_index + 1
			# should put x0 here (maybe assume x0 location = given tree root location?)
			x_mean_for_current_node = internal_locations[x_index]
			y_mean_for_current_node = internal_locations[y_index]
		else:
			x_index = internal_mapping[node.parent.label]
			y_index = x_index + 1
			x_mean_for_current_node = internal_locations[x_index]
			y_mean_for_current_node = internal_locations[y_index]

			if displacement_dict != None:
				x_mean_for_current_node += displacement_dict[node.label][0]
				y_mean_for_current_node += displacement_dict[node.label][1]

		mean_vector_x.append(x_mean_for_current_node)
		mean_vector_y.append(y_mean_for_current_node)

	return mean_vector_x, mean_vector_y

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


