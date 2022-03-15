#!/bin/bash

REPO_DIR="$(cd .. && pwd)"
MAX_JOBS=1
EMAILS="/home/rafael/wiki_dump.txt"

waitforjobs() {
    while test $(jobs -p | wc -w) -ge "$1"; do wait -n; done
}

experiment() {
    nkeys=$1 && shift
    nqrys=$1 && shift
    frac_full=$1 && shift
    membudg=$1 && shift
    filter=$1 && shift
    expdir=$1 && shift

    echo "$nkeys $nqrys $frac_full $membudg $filter $expdir"

    cd "$expdir"
    touch ./experiment_result
    echo -e "### BEGIN EXPERIMENT ###" >> ./experiment_result

    echo -e "\n\n### BEGIN EXPERIMENT DESCRIPTION ###" >> ./experiment_result
    echo -e "\tFilter name:\t$filter" >> ./experiment_result
    echo -e "\tNumber of keys:\t$nkeys" >> ./experiment_result
    echo -e "\tNumber of queries:\t$nqrys" >> ./experiment_result
    echo -e "\tApproximate fraction of full queries:\t$frac_full" >> ./experiment_result
    echo -e "\tFilter bits-per-key:\t$membudg" >> ./experiment_result

    echo python "$REPO_DIR"/workload_gen/gen_strings.py "$EMAILS" "./my_data" "$nkeys" "$nqrys" "$frac_full"
    python "$REPO_DIR"/workload_gen/gen_strings.py "$EMAILS" "./my_data" "$nkeys" "$nqrys" "$frac_full"
    if [ $filter = "SuRF" ]; then
        surfhlen=$1 && shift
        surfrlen=$1 && shift
        echo -e "\tSurf hash suffix length:\t$surfhlen" >> ./experiment_result
        echo -e "\tSurf real suffix length:\t$surfrlen" >> ./experiment_result
        echo -e "### END EXPERIMENT DESCRIPTION ###\n\n" >> ./experiment_result
        "$REPO_DIR"/standalone/run_workload_string ./my_data/data.txt ./my_data/txn.txt ./my_data/upper_bound.txt "$membudg" "$filter" "$surfhlen" "$surfrlen" >> ./experiment_result
    else
        dfsdiff=$1 && shift
        bfsdiff=$1 && shift
        cutoff=$1 && shift
        echo -e "\tBfs diffidence:\t$bfsdiff" >> ./experiment_result
        echo -e "\tDfs diffidence:\t$dfsdiff" >> ./experiment_result
        echo -e "\tCutoff:\t$cutoff" >> ./experiment_result
        echo -e "### END EXPERIMENT DESCRIPTION ###\n\n" >> ./experiment_result
        echo "$REPO_DIR"/standalone/run_workload_string ./my_data/data.txt ./my_data/txn.txt ./my_data/upper_bound.txt "$membudg" "$filter" "$dfsdiff" "$bfsdiff" "$cutoff" 
        "$REPO_DIR"/standalone/run_workload_string ./my_data/data.txt ./my_data/txn.txt ./my_data/upper_bound.txt "$membudg" "$filter" "$dfsdiff" "$bfsdiff" "$cutoff" >> ./experiment_result
    fi
    echo -e "### END EXPERIMENT ###" >> ./experiment_result
}

nkeys_arr=(1000)
nqrys_arr=(1000000)
dfsdiff_arr=(10000)
bfsdiff_arr=(64)
membudg_arr=(14)
surfhlen_arr=(0)
surfrlen_arr=(0)
expdir_arr=()
cutoff_arr=(0)
fracfull_arr=(0 0.5)

#DST experiments
for nkeys in "${nkeys_arr[@]}"; do
    for nqrys in "${nqrys_arr[@]}"; do
	for dfsdiff in "${dfsdiff_arr[@]}"; do
	    for bfsdiff in "${bfsdiff_arr[@]}"; do
		for membudg in "${membudg_arr[@]}"; do
		    for cutoff in "${cutoff_arr[@]}"; do
			for fracfull in "${fracfull_arr[@]}"; do
			    expdir=$(mktemp -d /tmp/dst_experiment.XXXXXXX)
			    expdir_arr+=("$expdir")
			    waitforjobs "$MAX_JOBS"
			    experiment $nkeys $nqrys $fracfull $membudg "DST" "$expdir" $dfsdiff $bfsdiff $cutoff
			done
		    done
		done
	    done
	done
    done
done

#SuRF experiments
for nkeys in "${nkeys_arr[@]}"; do
    for nqrys in "${nqrys_arr[@]}"; do
	for membudg in 10; do
	    for surfhlen in "${surfhlen_arr[@]}"; do
		for surfrlen in "${surfrlen_arr[@]}"; do
		    for fracfull in "${fracfull_arr[@]}"; do
			expdir=$(mktemp -d /tmp/dst_experiment.XXXXXXX)
			expdir_arr+=("$expdir")
			waitforjobs "$MAX_JOBS"
			experiment $nkeys $nqrys $fracfull $membudg "SuRF" "$expdir" $surfhlen $surfrlen
		    done
		done
	    done
	done
    done
done

wait
for expdir in "${expdir_arr[@]}"; do
    cat "$expdir/experiment_result"
done
