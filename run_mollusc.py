#! /usr/bin/env python
import os
import pickle
import problin_libs as problin
from problin_libs.sequence_lib import read_sequences, read_priors
from problin_libs.ML_solver import ML_solver
from problin_libs.EM_solver import EM_solver
from problin_libs.SpaLin_solver import SpaLin_solver
from problin_libs.Spatial_Division_solver import Spatial_Division_solver
from problin_libs.Topology_search_parallel import Topology_search_parallel as Topology_search_parallel
from problin_libs.Topology_search import Topology_search as Topology_search_sequential
#from problin_libs.Topology_search import Topology_search
from treeswift import *
import random
import argparse
import timeit
from sys import argv,exit,stdout
from copy import deepcopy

def main():
    parser = argparse.ArgumentParser()

    # input arguments
    parser.add_argument("-t","--topology",required=True,help="Input trees in newick format.") 
    parser.add_argument("-c","--characters",required=True,help="The input character matrix. Must have header.")
    parser.add_argument("-p","--priors",required=False, default="uniform", help="The input prior matrix Q. Default: if not specified, use a uniform prior.")
    parser.add_argument("--delimiter",required=False,default="tab",help="The delimiter of the input character matrix. Can be one of {'comma','tab','whitespace'} .Default: 'tab'.")
    parser.add_argument("-m","--maskedchar",required=False,default="-",help="Masked character. Default: if not specified, assumes '-'.")
    parser.add_argument("-o","--output",required=True,help="The output file.")
   
    # which problem are you solving? 
    # parser.add_argument("--solver",required=False,default="EM",help="Specify a solver. Options are 'Spalin', 'Scipy' or 'EM'. Default: 'Spalin' if -S is given else 'EM'")
    parser.add_argument("--topology_search",action='store_true', required=False,help="Perform topology search using NNI operations. Always return fully resolved (i.e. binary) tree.")
    parser.add_argument("--resolve_search",action='store_true', required=False,help="Resolve polytomies by performing topology search ONLY on branches with polytomies. This option has higher priority than --topoloy_search.")
    parser.add_argument("-L","--compute_llh",required=False,help="Compute likelihood of the input tree using the input (phi,nu). Will NOT optimize branch lengths, phi, or nu. The input tree MUST have branch lengths. This option has higher priority than --topoloy_search and --resolve_search.")
    parser.add_argument("-S","--spatial",required=False,help="Spatial locations of a subset of cells.")
    parser.add_argument("--spatial_only", required=False,action='store_true',default=False,help="Ignore the sequence information and only use spatial information")
    # problem formulation
    parser.add_argument("--noSilence",action='store_true',help="Assume there is no gene silencing, but allow missing data by dropout in sc-sequencing.")
    parser.add_argument("--noDropout",action='store_true',help="Assume there is no sc-sequencing dropout, but allow missing data by gene silencing.")

    # constraints
    parser.add_argument("--timescale", type=int,required=False,default=1,help="Total time")

    # miscellaneous
    parser.add_argument("-v","--verbose",required=False,action='store_true',help="Show verbose messages.")
    parser.add_argument("--nInitials",type=int,required=False,default=1,help="The number of initial points. Default: 1.")
    parser.add_argument("--randseeds",required=False,help="Random seeds. Can be a single interger number or a list of intergers whose length is equal to the number of initial points (see --nInitials).")
    parser.add_argument("--randomreps", required=False, default=1, type=int, help="Number of replicates to run for the random strategy of topology search.")
    parser.add_argument("--maxIters", required=False, default=500, type=int, help="Maximum number of iterations to run topology search.")
    parser.add_argument("--parallel", required=False,action='store_true', help="Turn on parallel version of topology search.")

    parser.add_argument("--divide", required=False,action='store_true', help="Model division of cells")
    parser.add_argument("--radius", required=False,default=None,type=float,help="Fixed displacement of cell division")


    if len(argv) == 1:
        parser.print_help()
        exit(0)
    
    print("Launching " + problin.PROGRAM_NAME + " version " + problin.PROGRAM_VERSION)
    print(problin.PROGRAM_NAME + " was called as follow: " + " ".join(argv))
    start_time = timeit.default_timer()
    
    args = vars(parser.parse_args())
    
    # preprocessing: read and analyze input
    delim_map = {'tab':'\t','comma':',','whitespace':' '}
    delimiter = delim_map[args["delimiter"]]
    msa, site_names = read_sequences(args["characters"],filetype="charMtrx",delimiter=delimiter,masked_symbol=args["maskedchar"])
    prefix = '.'.join(args["output"].split('.')[:-1])

    with open(args["topology"],'r') as f:
        input_trees = []
        for line in f:
            input_trees.append(line.strip())

    k = len(msa[next(iter(msa.keys()))])
    if args["compute_llh"]:
        fixed_phi,fixed_nu = [float(x) for x in args["compute_llh"].strip().split()]
    else:    
        fixed_phi = 0 if args["noDropout"] else None
        fixed_nu = 0 if args["noSilence"] else None

    if args["randseeds"] is None:
        random_seeds = None
    else:
        random_seeds = [int(x) for x in args["randseeds"].strip().split()]
        if args["nInitials"] != 1 and len(random_seeds) == 1:
            random_seeds = random_seeds[0]

    if args["priors"] == "uniform":
        print("No prior file detected, using uniform prior probabilities for each alphabet on each site.")
        # use the uniform Q matrix
        Q = []
        for i in range(k):
            M_i = set(msa[x][i] for x in msa if msa[x][i] not in [0,"?"])
            m_i = len(M_i)
            q = {x:1/m_i for x in M_i}
            q[0] = 0
            Q.append(q)
    else:
        Q = read_priors(args["priors"], site_names)



    known_locations = {}
    if args["spatial"] is not None:
        with open(args["spatial"],'r') as fin: 
            for line in fin:
                cellID,x,y = line.strip().split(",")
                known_locations[cellID] = (float(x),float(y))

        if args["divide"]:
            print("Optimization by MOLLUSC-Divide: Scipy-SLSQP using spatial information")
            selected_solver = Spatial_Division_solver
        else:
            print("Optimization by MOLLUSC: Scipy-SLSQP using spatial information") 
            selected_solver = SpaLin_solver
    else:            
        print("Optimization by generic solver (Scipy-SLSQP)")       
        selected_solver = ML_solver   
    if args["divide"]:
        if args["spatial"] is None:
            print("Can't use divide without spatial information")




    # main tasks        
    data = {'charMatrices':[msa],'locations':known_locations} 
    prior = {'Q':[Q]} 
    

    params = {'brlen_lower': 0.005 * args["timescale"], 'brlen_upper': args["timescale"], 'nu':fixed_nu if fixed_nu is not None else problin.eps,'phi':fixed_phi if fixed_phi is not None else problin.eps, 'radius': args['radius'],  'thetas': {}, 'spatialOnly': args["spatial_only"], "sigma": None, 'lambda_param': 1}  
    Topology_search = Topology_search_sequential if not args["parallel"] else Topology_search_parallel

    myTopoSearch = Topology_search(input_trees, selected_solver, data=data, prior=prior, params=params)

    if args["compute_llh"]:
        print("Compute likelihood of the input tree and specified parameters without any optimization")
        mySolver = myTopoSearch.get_solver()
        nllh = mySolver.negative_llh()
        opt_trees = myTopoSearch.treeTopoList
        opt_params = myTopoSearch.params
        print("Tree neagtive log-likelihood: " + str(nllh))
        print("Tree log-likelihood: " + str(-nllh))
    else:
        # setup the strategy
        my_strategy = deepcopy(problin.DEFAULT_STRATEGY)
        # enforce ultrametric or not?
        my_strategy['ultra_constr'] = True
        # resolve polytomies or not?
        my_strategy['resolve_search_only'] = args["resolve_search"]
        # how many initials for optimization of numerical parameters?        
        my_strategy['initials'] = args["nInitials"]

        # full search or local search to only resolve polytomies? 
        if not args["resolve_search"] and not args["topology_search"]:
            print("Optimizing branch lengths, phi, and nu without topology search")
            mySolver = myTopoSearch.get_solver()
            nllh = mySolver.optimize(initials=args["nInitials"],fixed_phi=fixed_phi,fixed_nu=fixed_nu,verbose=args["verbose"],random_seeds=random_seeds,ultra_constr=True,optimize_brlens=True)      
            myTopoSearch.update_from_solver(mySolver)
            opt_trees = myTopoSearch.treeTopoList
            opt_params = myTopoSearch.params
        else:
            if args["resolve_search"]:
                print("Starting local topology search to resolve polytomies")
            else:
                print("Starting topology search")                 
            if args["parallel"]:
                print("Running topology search in parallel...")
            else:
                print("Running topology search sequentially...")
            randval = int(random.random() * 1000)
            checkpoint_file = f"{prefix}._ckpt.{randval}.txt"


            opt_trees,max_score,opt_params = myTopoSearch.search(maxiter=args["maxIters"], verbose=args["verbose"], strategy=my_strategy, nreps=args['randomreps'],checkpoint_file=checkpoint_file) 
            nllh = -max_score        
    
    # post-processing: analyze results and output 
    outfile = args["output"]   
    with open(outfile,'w') as fout:
        fout.write("Newick tree (s):\n") 
        for tree in opt_trees:
            fout.write(tree + "\n")

        fout.write("Negative-llh: " +  str(nllh) + "\n")
        if args['spatial_only'] == False:
            fout.write("Dropout rate: " + str(opt_params['phi']) + "\n")
            fout.write("Silencing rate: " + str(opt_params['nu']) + "\n")
            fout.write("Lambda: " + str(opt_params['lambda_param']) + "\n")
            
        if args["spatial"] is not None:
            fout.write("Sigma: " + str(opt_params['sigma']) + "\n")
        if args["spatial"] is not None:
            has_newly_inferred = False
            inferred_locations = {}
            mySolver = myTopoSearch.get_solver()

            if args["divide"] == True:
                fout.write("Radius " + str(opt_params['radius']) + '\n')
                fout.write("Thetas \n")
                for cell in mySolver.params.thetas:
                    theta = mySolver.params.thetas[cell]
                    fout.write(str(cell) + " " +  str(theta) + "\n")


    stop_time = timeit.default_timer()
    print("Runtime (s):", stop_time - start_time)

if __name__ == "__main__":
    main()
