data_dir=/Users/uym2/my_gits/spalin_intMemoir/processed_data/brlen_estimation

for d in $data_dir/s*c*; do
    n=`wc -l $d/characters.txt | awk '{print $1;}'`
    if [ $n -gt 10 ]; then
        s=`basename $d`
        echo $s
        mkdir $s
        cp $d/{true_tree.nwk,true_locations.txt,leaf_locations.txt} $s
    fi    
done    
