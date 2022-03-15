#!/bin/bash

REPO_DIR="$(cd .. && pwd)"
MAX_JOBS=1

waitforjobs() {
    while test $(jobs -p | wc -w) -ge "$1"; do wait -n; done
}

experiment() {
    nkeys=$1
    keylen=$2
    nqrys=$3
    kdist=$4
    qdist=$5
    mrange=$6
    dfsdiff=$7
    bfsdiff=$8
    pqratio=$9
    membudg="${10}"
    filter="${11}"
    expdir="${12}"
    corrd="${13}"

    cd "$expdir"
    touch ./experiment_result

    echo -e "Running $filter experiment with the following parameters:\n" >> ./experiment_result
    echo -e "\tNumber of keys:\t$nkeys" >> ./experiment_result
    echo -e "\tNumber of bits in a key:\t$keylen" >> ./experiment_result
    echo -e "\tNumber of queries:\t$nqrys" >> ./experiment_result
    echo -e "\tKey distribution:\t$kdist" >> ./experiment_result
    echo -e "\tQuery distribution:\t$qdist" >> ./experiment_result
    echo -e "\tMaximum range:\t$mrange" >> ./experiment_result
    echo -e "\tPoint-Query ratio:\t$pqratio" >> ./experiment_result
    echo -e "\tCorrelation degree:\t$corrd" >> ./experiment_result


    "$REPO_DIR"/workload_gen/workload "$nkeys" "$keylen" "$nqrys" "$mrange" "$pqratio" "$kdist" "$qdist" "$corrd"
    if [ $filter = "SuRF" ]; then
	surfhlen="${14}"
	surfrlen="${15}"
	echo -e "\tSurf hash suffix length:\t$surfhlen" >> ./experiment_result
	echo -e "\tSurf real suffix length:\t$surfrlen" >> ./experiment_result
	"$REPO_DIR"/standalone/disk_experiment ./my_data/data.txt ./my_data/txn.txt ./my_data/upper_bound.txt "$mrange" "$dfsdiff" "$bfsdiff" "$membudg" "$filter" "$surfhlen" "$surfrlen" >> ./experiment_result
    else
	echo -e "\tBfs diffidence:\t$bfsdiff" >> ./experiment_result
	echo -e "\tDfs diffidence:\t$dfsdiff" >> ./experiment_result
	echo -e "\tFilter bits-per-key:\t$membudg" >> ./experiment_result
	"$REPO_DIR"/standalone/disk_experiment ./my_data/data.txt ./my_data/txn.txt ./my_data/upper_bound.txt "$mrange" "$dfsdiff" "$bfsdiff" "$membudg" "$filter" >> ./experiment_result
    fi

    rm -rf ./my_data/
}

nkeys_arr=(100000)
keylen_arr=(64) # 64)
nqrys_arr=(1000000)
kdist_arr=("uniform") # "normal" "zipfian")
qdist_arr=("uniform") # "normal" "zipfian" "adversary")
#mrange_arr=(2048)
mrange_arr=(4096)
dfsdiff_arr=(1000)
bfsdiff_arr=(32)
pqratio_arr=(0.0 0.2 0.4 0.6 0.8 1.0)
membudg_arr=(8)
surfhlen_arr=(0)
surfrlen_arr=(0)
expdir_arr=()
corrd_arr=(1)

#DST experiments
for nkeys in "${nkeys_arr[@]}"; do
    for keylen in "${keylen_arr[@]}"; do
	for nqrys in "${nqrys_arr[@]}"; do
	    for kdist in "${kdist_arr[@]}"; do
		for qdist in "${qdist_arr[@]}"; do
		    for mrange in "${mrange_arr[@]}"; do
			for dfsdiff in "${dfsdiff_arr[@]}"; do
			    for bfsdiff in "${bfsdiff_arr[@]}"; do
				for pqratio in "${pqratio_arr[@]}"; do
				    for membudg in "${membudg_arr[@]}"; do
					for corrd in "${corrd_arr[@]}"; do
					    expdir=$(mktemp -d /tmp/dst_experiment.XXXXXXX)
					    expdir_arr+=("$expdir")
					    waitforjobs "$MAX_JOBS"
					    experiment $nkeys $keylen $nqrys $kdist $qdist $mrange $dfsdiff $bfsdiff $pqratio $membudg "DST" "$expdir" "$corrd_arr" &
					done
				    done
				done
			    done
			done
		    done
		done
	    done
	done
    done
done

#SuRF experiments
for nkeys in "${nkeys_arr[@]}"; do
    for keylen in "${keylen_arr[@]}"; do
	for nqrys in "${nqrys_arr[@]}"; do
	    for kdist in "${kdist_arr[@]}"; do
		for qdist in "${qdist_arr[@]}"; do
		    for mrange in "${mrange_arr[@]}"; do
			for dfsdiff in 1000; do
			    for bfsdiff in 32; do
				for pqratio in "${pqratio_arr[@]}"; do
				    for membudg in 10; do
					for surfhlen in "${surfhlen_arr[@]}"; do
					    for surfrlen in "${surfrlen_arr[@]}"; do
						for corrd in "${corrd_arr[@]}"; do
						    expdir=$(mktemp -d /tmp/dst_experiment.XXXXXXX)
						    expdir_arr+=("$expdir")
						    waitforjobs "$MAX_JOBS"
						    experiment $nkeys $keylen $nqrys $kdist $qdist $mrange $dfsdiff $bfsdiff $pqratio $membudg "SuRF" "$expdir" $corrd $surfhlen $surfrlen &
						done
					    done
					done
				    done
				done
			    done
			done
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
