from treeswift import *
import math
from math import log,exp,sqrt,pi,comb
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
	def __init__(self,treeTopo,data,prior,params={'nu':0,'phi':0,'sigma':0,'lambda_param': 1, 'radius': 0,'thetas': {}}):
		super(Spatial_Division_solver,self).__init__(treeTopo,data,prior,params)

		self.params.radius = params['radius']
		self.optimize_sigma = False
		self.optimize_radius = False
		if self.params.radius == None:
			print("We need to optimize the radius")
			self.optimize_radius = True
			self.params.radius = 1
		self.leaf_locations = data['locations']
		self.optimize_sigma = True

		#self.params.lambda_param = params['lambda_param']
		try:
			
			self.brlen_lower_bound = params['brlen_lower']
			self.brlen_upper_bound = params['brlen_upper']
		except:
			# since we make a mySolver again.
			pass

		if params['sigma'] == None:
			self.params.sigma = 10 # just some default
		else:
			self.params.sigma = params['sigma']
			self.optimize_sigma = False

		self.params.thetas= params['thetas']

		self.onlyspatial = False
		try:
			self.onlyspatial = params["spatialOnly"]
		except:
			pass
		self.num_internal = 0
		for tree in self.trees:
			for node in tree.traverse_postorder():
				if node.num_children() == 2:
					self.num_internal += 1

	def get_params(self):
		return {'lambda_param': self.params.lambda_param, 'phi':self.params.phi,'nu':self.params.nu,'sigma':self.params.sigma,'locations':self.leaf_locations, 'radius': self.params.radius, 'thetas': self.params.thetas}

	def llh_of_separation_force(self, locations, thetas):
		llh = 0
		sigma = self.params.sigma

		for i in [0,1]: # do this for the x and y coordinates
			for one_tree in self.trees:
				this_tree_locations = deepcopy(locations)
				tree = deepcopy(one_tree)
				tree.scale_edges(sigma**2)


				stored_pruned_nodes = dict() #key = label, value = (value from data, mean)
				stored_means_nodes_in_tree = dict()

				# compute initial expected values for all the leaf nodes
				for node in tree.traverse_leaves():
					sum_of_displacement = 0
					prev_node_label = node.label
					for node_p in node.traverse_ancestors():
						if node_p.label == node.label: # will always start with current node
							continue
						left_child = node_p.child_nodes()[0]

						multiplier = -1
						if left_child.label == prev_node_label: # was a left child, do addition
							multiplier = 1

						if i == 0: # doing x coordinate
							sum_of_displacement += multiplier * self.params.radius * math.cos(thetas[node_p.label])
						if i == 1: # doing y coordinate
							sum_of_displacement += multiplier * self.params.radius * math.sin(thetas[node_p.label])
						prev_node_label = node_p.label

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
		if self.onlyspatial:
			return self.llh_of_separation_force(self.leaf_locations, self.params.thetas)
		else:
			return self.lineage_llh() + self.llh_of_separation_force(self.leaf_locations, self.params.thetas)

	def ini_thetas(self):
		thetas = []     
		idx = 0
		for tree in self.trees:
			for node in tree.traverse_postorder():
				if node.num_children() == 2:
					rand_theta = np.random.uniform(0, 2*pi)
					thetas.append(rand_theta) # x coor
		return thetas


	def ini_all(self,fixed_phi=None,fixed_nu=None):
		x0_brlens = self.ini_brlens()
		x0_nu = self.ini_nu(fixed_nu=fixed_nu)
		x0_phi = self.ini_phi(fixed_phi=fixed_phi)
		x0_sigma = 22 # hard code for now     
		x0_thetas = self.ini_thetas()
		x0_radius = 1 #hard code for now to some unreasonable number

		if self.optimize_radius == True:
			ini_params = (['brlens','nu','phi','lambda_param','thetas','radius','sigma'],{'brlens':x0_brlens,'nu':[x0_nu],'phi':[x0_phi], 'thetas': x0_thetas,'lambda_param': [1], 'radius': [x0_radius], 'sigma':[x0_sigma]})
		else:
			ini_params = (['brlens','nu','phi','lambda_param','thetas','sigma'],{'brlens':x0_brlens,'nu':[x0_nu],'phi':[x0_phi], 'thetas': x0_thetas,'lambda_param': [1], 'sigma':[x0_sigma]})
		

		return ini_params 

	def optimize_one(self,randseed,fixed_phi=None,fixed_nu=None,optimize_brlens=True,verbose=1,ultra_constr=False):
		# optimize using a specific initial point identified by the input randseed
		# verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent
		warnings.filterwarnings("ignore")
		def nllh(x): 
			self.x2params(x,fixed_nu=fixed_nu,fixed_phi=fixed_phi,include_brlens=optimize_brlens)            
			return -self.__llh__()
		
		seed(a=randseed)
		np.random.seed(randseed)
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
		sum_constraints = []
		for tree in self.trees:
			for node in tree.traverse_leaves():
				one_constraint = [0]*len(x0)
				for node_p in node.traverse_ancestors(): 
					idx = self.node_label_to_index[node_p.label]
					one_constraint[idx] = 1
				sum_constraints.append(one_constraint)

		constraints.append(optimize.LinearConstraint(csr_matrix(sum_constraints),[self.brlen_upper_bound - 0.01]*len(sum_constraints),[self.brlen_upper_bound + 0.01]*len(sum_constraints),keep_feasible=False))


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
		self.x2lambda(x)
		self.x2thetas(x)

		if self.optimize_sigma:
			self.params.sigma = x[-1] 
		if self.optimize_radius:
			self.params.radius = x[-2]  

	def x2thetas(self,x):
		idx = self.num_edges + 3
		thetas = x[idx:-1]
		i = 0
		for tree in self.trees:
			for node in tree.traverse_postorder():
				if node.num_children() == 2:
					self.params.thetas[node.label] = thetas[i]
					i += 1

	def bound_thetas(self):
		num_nodes = 0		
		return [-np.inf]*self.num_internal,[np.inf]*self.num_internal

	def get_bound(self,keep_feasible=False,fixed_phi=None,fixed_nu=None,include_brlens=True):
		br_lower,br_upper = self.bound_brlen() if include_brlens else ([],[])
		phi_lower,phi_upper = self.bound_phi(fixed_phi=fixed_phi)
		nu_lower,nu_upper = self.bound_nu(fixed_nu=fixed_nu)
		sigma_lower,sigma_upper = self.bound_sigma()
		theta_lower, theta_upper = self.bound_thetas()
		lambda_lower, lambda_upper = (0,100000)
		radius_lower, radius_upper = (0.001, 30)

		if self.optimize_radius:
			combined_lower = br_lower+[nu_lower,phi_lower,lambda_lower]+theta_lower+[radius_lower]+[sigma_lower]
			combined_upper = br_upper+[nu_upper,phi_upper,lambda_upper]+theta_upper+[radius_upper]+[sigma_upper]
		else:
			combined_lower = br_lower+[nu_lower,phi_lower,lambda_lower]+theta_lower+[sigma_lower]
			combined_upper = br_upper+[nu_upper,phi_upper,lambda_upper]+theta_upper+[sigma_upper]

		bounds = optimize.Bounds(combined_lower,combined_upper,keep_feasible=keep_feasible)
		return bounds   