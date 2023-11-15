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
    def __init__(self,treeTopo,data,prior,params={'nu':0,'phi':0,'sigma':0}):
        super(SpaLin_solver,self).__init__(treeTopo,data,prior,params)
        self.given_locations = data['locations']
        self.params.sigma = params['sigma']
        self.inferred_locations = {} if 'locations' not in params else params['locations']
        for x in self.given_locations:
            self.inferred_locations[x] = self.given_locations[x]
    
    def get_params(self):
        return {'phi':self.params.phi,'nu':self.params.nu,'sigma':self.params.sigma,'locations':self.inferred_locations}
    
    def show_params(self):                   
        nllh = self.negative_llh() 
        print("Negative llh: " + str(nllh))
        print("Dropout rate: " + str(self.params.phi))
        print("Silencing rate: " + str(self.params.nu))
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
        mutation_rate = 0.006 # hardcoding this because i'm lazy

        llh = 0
        for one_tree in self.trees:
            tree = deepcopy(one_tree)
            # need to correct the brlens with the mutation rate param as well as multiply by variance
            tree.scale_edges(1/mutation_rate)
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
        return self.lineage_llh() + self.spatial_llh_marginalized(self.inferred_locations)
        return final_llh
    
    '''
    def optimize_one(self,randseed,fixed_phi=None,fixed_nu=None,verbose=1,ultra_constr=False,optimize_brlens=True):
        param_cats,ini_params = self.ini_all(fixed_phi=fixed_phi,fixed_nu=fixed_nu)
        x0_brlens = ini_params['brlens'] if optimize_brlens else []
        x0_phi = ini_params['phi']
        x0_nu = ini_params['nu']
        x_lin_ini = x0_brlens + x0_phi + x0_nu
        x_spa_ini = ini_params['spatial']           
        for i in range(30):
            if optimize_brlens:
                f,status = self.optimize_brlen(x_lin_ini,fixed_phi=fixed_phi,fixed_nu=fixed_nu,verbose=verbose,ultra_constr=ultra_constr)
            f,status = self.optimize_locations()
        return f,status '''

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
            M = self.ultrametric_constr(padding={'phi','nu'})
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

    def optimize_locations(self):
        idx = 0
        for tree in self.trees:
            for node in tree.traverse_preorder():
                node.idx = idx
                idx += 1
        N = idx
        A_x = []
        b_x = []
        A_y = []
        b_y = []
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if node.label in self.given_locations:
                    a = [0]*N
                    a[node.idx] = 1
                    A_x.append(a)
                    A_y.append(a)
                    b_x.append(self.given_locations[node.label][0])
                    b_y.append(self.given_locations[node.label][1])
                elif node.is_root():
                    a = [0]*N
                    a[node.idx] = 1
                    a[node.children[0].idx] = -1
                    A_x.append(a)
                    A_y.append(a)
                    b_x.append(0)
                    b_y.append(0)
                else:
                    p = node.parent
                    c1,c2 = node.children
                    a = [0]*N
                    a[node.idx] = 1/node.edge_length + 1/c1.edge_length + 1/c2.edge_length
                    a[c1.idx] = -1/c1.edge_length
                    a[c2.idx] = -1/c2.edge_length
                    a[p.idx] = -1/node.edge_length
                    A_x.append(a)
                    A_y.append(a)
                    b_x.append(0)
                    b_y.append(0)        
            x = np.dot(np.linalg.pinv(A_x),b_x)
            y = np.dot(np.linalg.pinv(A_y),b_y)
            for tree in self.trees:
                for node in tree.traverse_postorder():
                    self.inferred_locations[node.label] = (x[node.idx],y[node.idx])
        f = self.negative_llh()
        status = "optimal"
        return f,status            
    
    def ini_all(self,fixed_phi=None,fixed_nu=None):
        x0_brlens = self.ini_brlens()
        x0_nu = self.ini_nu(fixed_nu=fixed_nu)
        x0_phi = self.ini_phi(fixed_phi=fixed_phi)
        #x_lin = self.ini_brlens() + [self.ini_nu(fixed_nu=fixed_nu),self.ini_phi(fixed_phi=fixed_phi)]
        min_x = min([x[0] for x in self.given_locations.values()])
        min_y = min([x[1] for x in self.given_locations.values()])
        max_x = max([x[0] for x in self.given_locations.values()])
        max_y = max([x[1] for x in self.given_locations.values()])
        x0_spa = []
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if not node.label in self.given_locations:
                    x0_spa += [min_x+(max_x-min_x)*random(),min_y+(max_y-min_y)*random()]
        x0_sigma = 22 # hard code for now        
        ini_params = (['brlens','nu','phi','spatial','sigma'],{'brlens':x0_brlens,'nu':[x0_nu],'phi':[x0_phi],'spatial':x0_spa,'sigma':[x0_sigma]})
        #return x_lin, x_spa + [x_sigma]
        return ini_params 
    
    def x2params(self,x,fixed_nu=None,fixed_phi=None,include_brlens=True):
        if include_brlens:
            self.x2brlen(x)
        self.x2nu(x,fixed_nu=fixed_nu,include_brlens=include_brlens)
        self.x2phi(x,fixed_phi=fixed_phi,include_brlens=include_brlens)
        i = self.num_edges + 2
        for tree in self.trees:
            for node in tree.traverse_postorder():
                if not node.label in self.given_locations:
                    self.inferred_locations[node.label] = (x[i],x[i+1])
                    i += 2
        self.params.sigma = x[-1]        
               
    def bound_locations(self,lower=-np.inf,upper=np.inf):
        N = 0
        for tree in self.trees:
            for node in tree.traverse_postorder():
                N += 2*int(not node.label in self.given_locations)
        return [lower]*N,[upper]*N

    def bound_sigma(self):
        return (eps,np.inf)    

    def get_bound(self,keep_feasible=False,fixed_phi=None,fixed_nu=None,include_brlens=True):
        br_lower,br_upper = self.bound_brlen() if include_brlens else ([],[])
        phi_lower,phi_upper = self.bound_phi(fixed_phi=fixed_phi)
        nu_lower,nu_upper = self.bound_nu(fixed_nu=fixed_nu)
        spa_lower,spa_upper = self.bound_locations()
        sigma_lower,sigma_upper = self.bound_sigma()
        bounds = optimize.Bounds(br_lower+[nu_lower,phi_lower]+spa_lower+[sigma_lower],br_upper+[nu_upper,phi_upper]+spa_upper+[sigma_upper],keep_feasible=keep_feasible)
        return bounds   
