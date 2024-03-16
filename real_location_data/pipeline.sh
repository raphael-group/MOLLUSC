#! /bin/bash

#bash copy_data.sh
#for x in s*; do
#    echo $x
#    python ~/my_gits/problin/simulate_sequence.py -k 10 --alphabetsize 2 -t $x/true_tree.nwk -r 5 --mu 0.006 -p $x/k10
#done    
for r in {1..5}; do
    for x in s*; do
        echo $x $r
        echo spalin leaves
        # python /Users/gary/Documents/Projects/spalin/problin/run_problin.py -c $x/k10_r$r\_character_matrix.csv -t $x/true_tree.nwk -S $x/leaf_locations.txt --delimiter comma -p $x/k10_priors.csv -o $x/k10_r$r\_spalin_leaves.txt --nInitials 10 --ultrametric --randomreps 3 > $x/k10_r$r\_spalin_leaves.log 2>&1
        
        echo spalin all
        # python /Users/gary/Documents/Projects/spalin/problin/run_problin.py -c $x/k10_r$r\_character_matrix.csv -t $x/true_tree.nwk -S $x/true_locations.txt --delimiter comma -p $x/k10_priors.csv -o $x/k10_r$r\_spalin_all.txt --nInitials 10 --ultrametric --randomreps 3 > $x/k10_r$r\_spalin_all.log 2>&1
        
        echo problin
        # python /Users/gary/Documents/Projects/spalin/problin/run_problin.py -c $x/k10_r$r\_character_matrix.csv -t $x/true_tree.nwk --delimiter comma -p $x/k10_priors.csv -o $x/k10_r$r\_problin.txt --nInitials 10 --ultrametric --randomreps 3 > $x/k10_r$r\_problin.log 2>&1
    
        # topology search
        echo spalin leaves topology search
        python /Users/gary/Documents/Projects/spalin/problin/run_problin.py -c $x/k10_r$r\_character_matrix.csv -t $x/true_tree.nwk -S $x/leaf_locations.txt --delimiter comma -p $x/k10_priors.csv -o $x/k10_r$r\_spalin_leaves_topology_search.txt --nInitials 10 --ultrametric --randomreps 3 --topology_search --verbose> $x/k10_r$r\_spalin_leaves_topology_search.log 2>&1

        echo spalin all topology search
        # python /Users/gary/Documents/Projects/spalin/problin/run_problin.py -c $x/k10_r$r\_character_matrix.csv -t $x/true_tree.nwk -S $x/true_locations.txt --delimiter comma -p $x/k10_priors.csv -o $x/k10_r$r\_spalin_all_topology_search.txt --nInitials 10 --ultrametric --randomreps 3 --topology_search > $x/k10_r$r\_spalin_all_topology_search.log 2>&1

        echo problin topology search
        # python /Users/gary/Documents/Projects/spalin/problin/run_problin.py -c $x/k10_r$r\_character_matrix.csv -t $x/true_tree.nwk --delimiter comma -p $x/k10_priors.csv -o $x/k10_r$r\_problin_topology_search.txt --nInitials 10 --ultrametric --randomreps 3 --topology_search > $x/k10_r$r\_problin_topology_search.log 2>&1
    done 
done    
