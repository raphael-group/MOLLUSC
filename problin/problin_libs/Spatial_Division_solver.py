from treeswift import *
import math
from math import log,exp,sqrt,pi
from random import random, seed
from scipy import optimize
from scipy.sparse import csr_matrix
import warnings
import numpy as np
from problin_libs import min_llh, eps
from problin_libs.ML_solver import ML_solver
from problin_libs.SpaLin_solver import SpaLin_solver
from copy import deepcopy
from scipy.stats import norm 


class Spatial_Division_solver(SpaLin_solver):
	def __init__(self,treeTopo,data,prior,params={'nu':0,'phi':0,'sigma':0}):
		super(Spatial_Division_solver,self).__init__(treeTopo,data,prior,params)
		self.beta = params['beta']
		if self.beta == None:
			# ignore the nonlinear constraint
			pass
		self.leaf_locations = data['locations']
		self.params.sigma = params['sigma']
		self.params.displacement_amounts = {}

		self.num_internal = 0
		for tree in self.trees:
			for node in tree.traverse_postorder():
				if node.num_children() == 2:
					self.num_internal += 1

	def llh_of_separation_force(self, locations, displacement_amounts):
		llh = 0
		mutation_rate = 0.006
		sigma = self.params.sigma

		for i in [0,1]: # do this for the x and y coordinates
			for one_tree in self.trees:
				this_tree_locations = deepcopy(locations)
				tree = deepcopy(one_tree)
				#tree.scale_edges(1/mutation_rate)
				tree.scale_edges(sigma**2)


				stored_pruned_nodes = dict() #key = label, value = (value from data, mean)
				stored_means_nodes_in_tree = dict()

				# compute initial expected values for all the leaf nodes
				for node in tree.traverse_leaves():
					sum_of_displacement = 0
					for node_p in node.traverse_ancestors():
						if (node_p.label == 'virtual'):
							pass
						else:
							sum_of_displacement += displacement_amounts[node_p.label][i]
					stored_means_nodes_in_tree[node.label] = sum_of_displacement

		
				
				while tree.num_nodes() > 1:
					reset = 0
					for node_i in tree.traverse_leaves():
						for node_j in tree.traverse_leaves():
							if (node_i.label != node_j.label) and (node_i.parent == node_j.parent):
								node_k = node_i.parent

								#updated pruned nodes with (value_u_k, mu_u_k,brlen_u_k)

								pruned_node_value = this_tree_locations[node_i.label][i] - this_tree_locations[node_j.label][i]
								pruned_node_mean = stored_means_nodes_in_tree[node_i.label] - stored_means_nodes_in_tree[node_j.label]
								pruned_node_brlen = node_i.edge_length + node_j.edge_length
								stored_pruned_nodes[node_k.label] = (pruned_node_value, pruned_node_mean, pruned_node_brlen)

								#set and update values for the parent
								new_brlen_k = node_k.edge_length + (node_i.edge_length * node_j.edge_length)/(node_i.edge_length + node_j.edge_length)
								new_mean_k = (node_j.edge_length / (node_i.edge_length + node_j.edge_length)) * stored_means_nodes_in_tree[node_i.label] + (node_i.edge_length / (node_i.edge_length + node_j.edge_length)) * stored_means_nodes_in_tree[node_j.label]

								# initialize the internal node if it doesn't exist:
								if node_k.label not in this_tree_locations:
									this_tree_locations[node_k.label] = (0,0)
								if i == 0:
									this_tree_locations[node_k.label] = ((node_j.edge_length / (node_i.edge_length + node_j.edge_length)) * this_tree_locations[node_i.label][0] + (node_i.edge_length / (node_i.edge_length + node_j.edge_length)) * this_tree_locations[node_j.label][0],this_tree_locations[node_k.label][1])
								else:
									this_tree_locations[node_k.label] = (this_tree_locations[node_k.label][0],(node_j.edge_length / (node_i.edge_length + node_j.edge_length)) * this_tree_locations[node_i.label][1] + (node_i.edge_length / (node_i.edge_length + node_j.edge_length)) * this_tree_locations[node_j.label][1])
								stored_means_nodes_in_tree[node_k.label] = new_mean_k
								node_k.set_edge_length(new_brlen_k)

								#update the tree
								node_k.remove_child(node_i)
								node_k.remove_child(node_j)

								reset = 1
								break
						if reset == 1: # do this so that we recheck the num_nodes() > 2 condition after each prune
							break

				llh_one_tree = 0
				for pruned_node_label in stored_pruned_nodes:
					(pruned_value, pruned_mean, pruned_brlen) = stored_pruned_nodes[pruned_node_label]
					likelihood = norm.pdf(pruned_value, pruned_mean, sigma * sqrt(pruned_brlen))

					if likelihood == 0: # should set this tolerance thing a bit more intelligently
						llh_one_tree += -np.inf
					else:
						llh_one_tree += log(likelihood)
				llh += llh_one_tree
		return llh

	def __llh__(self):
		return self.lineage_llh() + self.llh_of_separation_force(self.leaf_locations, self.params.displacement_amounts)
		return final_llh

	def ini_sep(self,fixed_seps = None):
		if fixed_seps is not None:
			return fixed_seps
		else:
			forces = []     
			idx = 0
			if self.beta != None:
				sep_force = self.beta / sqrt(2)
			else:
				sep_force = 1 / sqrt(2)
			for tree in self.trees:
				for node in tree.traverse_postorder():
					if node.num_children() == 2:
					# initialize as symmetric

						# coordinate displacements for left child
						forces.append(sep_force) # x coor
						forces.append(sep_force) # y coor

						# coordinate displacements for right child
						forces.append(-sep_force) # x coor
						forces.append(-sep_force) # y coor

			return forces


	def ini_all(self,fixed_phi=None,fixed_nu=None):
		x0_brlens = self.ini_brlens()
		x0_nu = self.ini_nu(fixed_nu=fixed_nu)
		x0_phi = self.ini_phi(fixed_phi=fixed_phi)
		x0_sigma = 22 # hard code for now     
		x0_sep_forces = self.ini_sep()


		ini_params = (['brlens','nu','phi','sep_forces','sigma'],{'brlens':x0_brlens,'nu':[x0_nu],'phi':[x0_phi], 'sep_forces': x0_sep_forces,'sigma':[x0_sigma]})
		return ini_params 

	def optimize_one(self,randseed,fixed_phi=None,fixed_nu=None,optimize_brlens=True,verbose=1,ultra_constr=False):
		# optimize using a specific initial point identified by the input randseed
		# verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent
		warnings.filterwarnings("ignore")
		def nllh(x): 
			self.x2params(x,fixed_nu=fixed_nu,fixed_phi=fixed_phi,include_brlens=optimize_brlens)            
			return -self.__llh__()
		
		seed(a=randseed)
		#x0 = self.ini_all(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
		param_cats,ini_params = self.ini_all(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
		x0 = []
		for p in param_cats:
			if p != 'brlens' or optimize_brlens:
				x0 = x0 + ini_params[p]
		self.az_partition()

		bounds = self.get_bound(fixed_phi=fixed_phi,fixed_nu=fixed_nu,include_brlens=optimize_brlens)
		constraints = []    
		param_length = len(x0)
		if optimize_brlens:
			A = []
			b = []
			idx = 0
			# matrix A is for the fixed nodes
			for tree in self.trees:
				for node in tree.traverse_postorder():
					if node.mark_fixed:
						a = [0]*len(x0)
						a[idx] = 1
						A.append(a)
						b.append(node.edge_length)
					idx += 1   
			# include padding for matrix A to account for the other params
			# fixed_constraint_matrix = np.zeros((param_length, param_length))
			# fixed_constraint_matrix[:A.shape[0],:A.shape[1]] = A
			if len(A) > 0:    
				# fixed_constraint_bound = np.zeros(param_length)
				# fixed_constraint_bound[:len(b)] = b
				constraints.append(optimize.LinearConstraint(csr_matrix(A),b,b,keep_feasible=False))
			if ultra_constr:
				M = self.ultrametric_constr()
				# print(len(x0))
				# ultametric_constraint_matrix = np.zeros((param_length, param_length))
				# ultametric_constraint_matrix[:M.shape[0], M.shape[1]] = M
				constraints.append(optimize.LinearConstraint(csr_matrix(M),[0]*len(M),[0]*len(M),keep_feasible=False))

			# create the constraints for the separation amount parameters
			# find the index where the values start
			index_of_displacement_start = self.num_edges + 2

			# x,y displacements of sister nodes are symmetric
			curr_internal = 0
			symmetric_x_constraint_matrix = []
			for tree in self.trees:
				for node in tree.traverse_postorder():
					if node.num_children() == 2:
						a = [0] * len(x0)
						a[index_of_displacement_start + curr_internal] = 1
						a[index_of_displacement_start + curr_internal+2] = 1
						symmetric_x_constraint_matrix.append(a)
						curr_internal += 4
			constraints.append(optimize.LinearConstraint(csr_matrix(symmetric_x_constraint_matrix),[0]*len(symmetric_x_constraint_matrix),[0]*len(symmetric_x_constraint_matrix),keep_feasible=False))

			curr_internal = 0
			symmetric_y_constraint_matrix = []
			for tree in self.trees:
				for node in tree.traverse_postorder():
					if node.num_children() == 2:
						a = [0] * len(x0)
						a[index_of_displacement_start + curr_internal+1] = 1
						a[index_of_displacement_start + curr_internal+3] = 1
						symmetric_y_constraint_matrix.append(a)
						curr_internal += 4
			constraints.append(optimize.LinearConstraint(csr_matrix(symmetric_y_constraint_matrix),[0]*len(symmetric_y_constraint_matrix),[0]*len(symmetric_y_constraint_matrix),keep_feasible=False))

			# sum of squares of x,y displacements are equal to a constant
			def sum_of_squares(x):
				# return x^2 + y^2 for all internal nodes
				curr_internal = 0
				array_of_sums = []
				for tree in self.trees:
					for node in tree.traverse_postorder():
						if node.num_children() == 2:
							x_squared = x[index_of_displacement_start + curr_internal]**2
							y_squared = x[index_of_displacement_start + curr_internal+1]**2
							sum1 = x_squared + y_squared
							array_of_sums.append(sum1)

							x_squared = x[index_of_displacement_start + curr_internal+2]**2
							y_squared = x[index_of_displacement_start + curr_internal+3]**2
							sum2 = x_squared + y_squared
							array_of_sums.append(sum2)
							curr_internal += 4
				return array_of_sums
			if self.beta != None:
				constraints.append(optimize.NonlinearConstraint(sum_of_squares, [self.beta**2 * 0.9] * self.num_internal * 2, [self.beta**2 * 1.1] * self.num_internal * 2))

		disp = (verbose > 0)
		out = optimize.minimize(nllh, x0, method="SLSQP", options={'disp':disp,'iprint':3,'maxiter':1000}, bounds=bounds,constraints=constraints)

		if out.success:

			self.x2params(out.x,fixed_phi=fixed_phi,fixed_nu=fixed_nu,include_brlens=optimize_brlens)
			params = self.params
			f = out.fun
		else:
			f,params = None,None
		status = "optimal" if out.success else out.message
		return f,status


	def x2params(self,x,fixed_nu=None,fixed_phi=None,include_brlens=True):
		if include_brlens:
			self.x2brlen(x)
		self.x2nu(x,fixed_nu=fixed_nu,include_brlens=include_brlens)
		self.x2phi(x,fixed_phi=fixed_phi,include_brlens=include_brlens)
		i = self.num_edges + 2
		self.x2displacement(x,i)
		self.params.sigma = x[-1] 

	def x2displacement(self,x,idx):
		x_displacement_amounts = x[idx:-1]
		i = 0
		for tree in self.trees:
			for node in tree.traverse_postorder():
				if node.num_children() == 2:
				# initialize as symmetric
					children = node.child_nodes()
					self.params.displacement_amounts[children[0].label] = (x_displacement_amounts[i],x_displacement_amounts[i+1])
					self.params.displacement_amounts[children[1].label] = (x_displacement_amounts[i+2],x_displacement_amounts[i+3])

					i += 4

	def bound_sep_force(self):
		num_nodes = 0
		for tree in self.trees:
			for node in tree.traverse_postorder():
				# number of separation forces equal to the 
				if node.is_leaf() == False:
					num_nodes += 2
		if self.beta != None:
			lower_bound = [-self.beta] * (num_nodes * 2) # times 2 for both x and y coordinate
			upper_bound = [self.beta] * (num_nodes * 2)
		else:
			lower_bound = [-100] * (num_nodes * 2) # times 2 for both x and y coordinate
			upper_bound = [-100] * (num_nodes * 2)			
		return lower_bound, upper_bound

	def get_bound(self,keep_feasible=False,fixed_phi=None,fixed_nu=None,include_brlens=True):
		br_lower,br_upper = self.bound_brlen() if include_brlens else ([],[])
		phi_lower,phi_upper = self.bound_phi(fixed_phi=fixed_phi)
		nu_lower,nu_upper = self.bound_nu(fixed_nu=fixed_nu)
		sigma_lower,sigma_upper = self.bound_sigma()
		sep_force_lower, sep_force_upper = self.bound_sep_force()

		combined_lower = br_lower+[nu_lower,phi_lower]+sep_force_lower + [sigma_lower]
		combined_upper = br_upper+[nu_upper,phi_upper]+sep_force_upper +[sigma_upper]
		bounds = optimize.Bounds(combined_lower,combined_upper,keep_feasible=keep_feasible)
		return bounds   

			   




