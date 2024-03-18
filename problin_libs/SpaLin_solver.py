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
from copy import deepcopy


class SpaLin_solver(ML_solver):
    def __init__(self,treeTopo,data,prior,params={'nu':0,'phi':0,'sigma':0,'lambda_param':1}):
        super(SpaLin_solver,self).__init__(treeTopo,data,prior,params)
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

        self.onlyspatial = False
        try:
            self.onlyspatial = params["spatialOnly"]
            if self.onlyspatial:
                print("only spatial")
        except:
            pass

    
    def get_params(self):
        return {'lambda_param': self.params.lambda_param,'phi':self.params.phi,'nu':self.params.nu,'sigma':self.params.sigma,'locations':self.leaf_locations}
    
    def show_params(self):                   
        nllh = self.negative_llh() 

        print("Negative llh: " + str(nllh))
        if self.onlyspatial == False: 
            print("Dropout rate: " + str(self.params.phi))
            print("Silencing rate: " + str(self.params.nu))
            print("Lambda: " + str(self.params.lambda_param))
        print("Sigma: " + str(self.params.sigma))
        

    def construct_distance_matrix(self,tree,locations):
        distances = dict()
        for node1 in tree.traverse_postorder():
            for node2 in tree.traverse_postorder():
                if node1.is_leaf() and node2.is_leaf():
                    x1,y1 = locations[node1.label]
                    x2,y2 = locations[node2.label]
                    distances[(node1.label,node2.label)] = (x1 - x2)**2 + (y1 - y2)**2
        return distances

    def spatial_llh_marginalized(self,locations):
        # given:
        # matrix D_ij as squared generalized distances
        # number of characters p
        # times of the tip populations (calculated from branch lengths) t_i
        # ancestors for each node
        llh = 0
        for one_tree in self.trees:
            tree = deepcopy(one_tree)

            tree.scale_edges(self.params.sigma**2)


            # construct distance matrix as dictionary of leaf names
            distance_matrix = self.construct_distance_matrix(tree, locations)

            # following the algorithm of Felsenstein 1973: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1762641/
            S = 0
            T = 1

            # while the tree has more than one leaf population
            while tree.num_nodes() > 2:
                for node in tree.traverse_preorder():
                    children = node.child_nodes()
                    if node.num_children() == 2 and children[0].is_leaf() and children[1].is_leaf():
                        v1 = children[0].get_edge_length()
                        v2 = children[1].get_edge_length()
                        
                        f = v1 / (v1+v2)

                        S += distance_matrix[(children[0].label,children[1].label)] / (v1 + v2)
                        T = T * (v1 + v2)

                        # update the distance matrix
                        for other_leaf in tree.traverse_postorder():
                            if other_leaf.is_leaf():
                                if (other_leaf.label != children[0].label) and (other_leaf.label != children[1].label):
                                    new_value = 0
                                    new_value += (1-f) * distance_matrix[(children[0].label, other_leaf.label)]
                                    new_value += f * distance_matrix[(children[1].label, other_leaf.label)]
                                    new_value -= f * (1-f) * distance_matrix[(children[0].label,children[1].label)]

                                    distance_matrix[(children[0].label, other_leaf.label)] = new_value
                                    distance_matrix[(other_leaf.label, children[0].label)] = new_value

                        # remove one of the children
                        node.remove_child(children[1])

                        # update the branch length of the other child
                        children[0].set_edge_length(node.get_edge_length() + (v1 * v2)/(v1 + v2))

                        # if this was a bifurcation, then we also remove this node and set the child to this node's parent
                        if node.num_children() == 1:
                            if tree.num_nodes() > 2:
                                parent = node.get_parent()

                                parent.remove_child(node)
                                children[0].set_parent(parent)
                                parent.add_child(children[0])
            if T == 0:
                llh += -np.inf
            else:
                llh += -log(T) - S/2
        return llh

    
    def spatial_llh(self,locations):
        llh = 0
        for tree in self.trees:
            for node in tree.traverse_preorder():
                if node.is_root() or not node.label in locations or not node.parent.label in locations:
                    continue
                d = node.edge_length
                curr_sigma = self.params.sigma*sqrt(d)
                x,y = locations[node.label]
                x_par,y_par = locations[node.parent.label]
                llh -= (0.5*((x-x_par)/curr_sigma)**2 + log(curr_sigma))
                llh -= (0.5*((y-y_par)/curr_sigma)**2 + log(curr_sigma))
        return llh 

    def __llh__(self):
        if self.onlyspatial:
            return self.spatial_llh_marginalized(self.leaf_locations)
        else:
            return self.lineage_llh() + self.spatial_llh_marginalized(self.leaf_locations)

    def optimize_brlen(self,x0,fixed_phi=None,fixed_nu=None,verbose=1,ultra_constr=False):
        # optimize branch lengths, phi, and nu using a specific initial point identified by the input randseed
        # verbose level: 1 --> show all messages; 0 --> show minimal messages; -1 --> completely silent
        warnings.filterwarnings("ignore")
        def nllh(x): 
            self.x2brlen(x)
            self.x2nu(x,fixed_nu=fixed_nu,include_brlens=True)
            self.x2phi(x,fixed_phi=fixed_phi,include_brlens=True)
            return -self.__llh__()        
        #seed(a=randseed)
        #x0 = self.ini_brlens() + [self.ini_nu(fixed_nu=fixed_nu),self.ini_phi(fixed_phi=fixed_phi)]
        self.az_partition()
        br_lower,br_upper = self.bound_brlen()  
        phi_lower,phi_upper = self.bound_phi(fixed_phi=fixed_phi)
        nu_lower,nu_upper = self.bound_nu(fixed_nu=fixed_nu)
        bounds = optimize.Bounds(br_lower+[nu_lower,phi_lower],br_upper+[nu_upper,phi_upper],keep_feasible=False)
        constraints = []    

        A = []
        b = []
        idx = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if node.mark_fixed:
                    a = [0]*len(x0)
                    a[idx] = 1
                    A.append(a)
                    b.append(node.edge_length)
                idx += 1   
        if len(A) > 0:     
            constraints.append(optimize.LinearConstraint(csr_matrix(A),b,b,keep_feasible=False))
        if ultra_constr:
            M = self.ultrametric_constr()
            constraints.append(optimize.LinearConstraint(csr_matrix(M),[0]*len(M),[0]*len(M),keep_feasible=False))
        disp = (verbose > 0)
        out = optimize.minimize(nllh, x0, method="SLSQP", options={'disp':disp,'iprint':1,'maxiter':1000}, bounds=bounds,constraints=constraints)
        if out.success:
            self.x2brlen(out.x)
            self.x2nu(out.x,fixed_nu=fixed_nu,include_brlens=True)
            self.x2phi(out.x,fixed_phi=fixed_phi,include_brlens=True)
            params = self.params
            f = out.fun
        else:
            f,params = None,None
        status = "optimal" if out.success else out.message
        return f,status        
    
    def ini_all(self,fixed_phi=None,fixed_nu=None):
        x0_brlens = self.ini_brlens()
        x0_nu = self.ini_nu(fixed_nu=fixed_nu)
        x0_phi = self.ini_phi(fixed_phi=fixed_phi)
        #x_lin = self.ini_brlens() + [self.ini_nu(fixed_nu=fixed_nu),self.ini_phi(fixed_phi=fixed_phi)]
        x0_sigma = 22 # hard code for now        
        ini_params = (['brlens','nu','phi','lambda_param', 'sigma'],{'brlens':x0_brlens,'nu':[x0_nu],'phi':[x0_phi],'lambda_param':[1], 'sigma':[x0_sigma]})
        #return x_lin, x_spa + [x_sigma]
        return ini_params 
    
    def x2params(self,x,fixed_nu=None,fixed_phi=None,include_brlens=True):
        if include_brlens:
            self.x2brlen(x)
        self.x2nu(x,fixed_nu=fixed_nu,include_brlens=include_brlens)
        self.x2phi(x,fixed_phi=fixed_phi,include_brlens=include_brlens)
        self.x2lambda(x)
        if self.optimize_sigma:
            self.params.sigma = x[-1]  
   
               
    def bound_sigma(self):
        return (eps,np.inf)    

    def get_bound(self,keep_feasible=False,fixed_phi=None,fixed_nu=None,include_brlens=True):
        br_lower,br_upper = self.bound_brlen() if include_brlens else ([],[])
        phi_lower,phi_upper = self.bound_phi(fixed_phi=fixed_phi)
        nu_lower,nu_upper = self.bound_nu(fixed_nu=fixed_nu)
        lambda_lower, lambda_upper = (0,100000)
        sigma_lower,sigma_upper = self.bound_sigma()
        bounds = optimize.Bounds(br_lower+[nu_lower,phi_lower, lambda_lower]+[sigma_lower],br_upper+[nu_upper,phi_upper, lambda_upper]+[sigma_upper],keep_feasible=keep_feasible)
        return bounds   
