#!/usr/local/bin/bash

REPO_DIR="$(cd .. && pwd)"
MAX_JOBS=1

waitforjobs() {
    while test $(jobs -p | wc -w) -ge "$1"; do wait -n; done
}

experiment() {
    nkeys=$1 && shift
    keylen=$1 && shift
    nqrys=$1 && shift
    kdist=$1 && shift
    qdist=$1 && shift
    pqratio=$1 && shift
    membudg=$1 && shift
    filter=$1 && shift
    expdir=$1 && shift
    corrd=$1 && shift

    cd "$expdir"
    touch ./experiment_result
    echo -e "### BEGIN EXPERIMENT ###" >> ./experiment_result

    echo -e "\n\n### BEGIN EXPERIMENT DESCRIPTION ###" >> ./experiment_result
    echo -e "\tFilter name:\t$filter" >> ./experiment_result
    echo -e "\tNumber of keys:\t$nkeys" >> ./experiment_result
    echo -e "\tNumber of bits in a key:\t$keylen" >> ./experiment_result
    echo -e "\tNumber of queries:\t$nqrys" >> ./experiment_result
    echo -e "\tKey distribution:\t$kdist" >> ./experiment_result
    echo -e "\tQuery distribution:\t$qdist" >> ./experiment_result
    echo -e "\tMaximum range:\t$mrange" >> ./experiment_result
    echo -e "\tPoint-Query ratio:\t$pqratio" >> ./experiment_result
    echo -e "\tFilter bits-per-key:\t$membudg" >> ./experiment_result
    echo -e "\tCorrelation degree:\t$corrd" >> ./experiment_result


    echo "$REPO_DIR"/workload_gen/workload "$nkeys" "$keylen" "$nqrys" "$mrange" "$pqratio" "$kdist" "$qdist" "$corrd"
    "$REPO_DIR"/workload_gen/workload "$nkeys" "$keylen" "$nqrys" "$mrange" "$pqratio" "$kdist" "$qdist" "$corrd"
    if [ $filter = "SuRF" ]; then
        surfhlen=$1 && shift
        surfrlen=$1 && shift
        echo -e "\tSurf hash suffix length:\t$surfhlen" >> ./experiment_result
        echo -e "\tSurf real suffix length:\t$surfrlen" >> ./experiment_result
        echo -e "### END EXPERIMENT DESCRIPTION ###\n\n" >> ./experiment_result
        "$REPO_DIR"/standalone/run_workload_disk ./my_data/data.txt ./my_data/txn.txt ./my_data/upper_bound.txt "$membudg" "$filter" "$surfhlen" "$surfrlen" >> ./experiment_result
    else
        dfsdiff=$1 && shift
        bfsdiff=$1 && shift
        cutoff=$1 && shift
        echo -e "\tBfs diffidence:\t$bfsdiff" >> ./experiment_result
        echo -e "\tDfs diffidence:\t$dfsdiff" >> ./experiment_result
        echo -e "\tCutoff:\t$cutoff" >> ./experiment_result
        echo -e "### END EXPERIMENT DESCRIPTION ###\n\n" >> ./experiment_result
        "$REPO_DIR"/standalone/run_workload_disk ./my_data/data.txt ./my_data/txn.txt ./my_data/upper_bound.txt "$membudg" "$filter" "$dfsdiff" "$bfsdiff" "$cutoff" >> ./experiment_result
    fi
    echo -e "### END EXPERIMENT ###" >> ./experiment_result

    rm -rf ./my_data/
}

nkeys_arr=(10000000)
keylen_arr=(64)
nqrys_arr=(10000)
kdist_arr=("uniform") # "normal" "zipfian")
qdist_arr=("adversary") # "normal" "zipfian" "adversary")
mrange_arr=(16)
dfsdiff_arr=(10000)
bfsdiff_arr=(64)
pqratio_arr=(0.0)
membudg_arr=(9.2)
surfhlen_arr=(0)
surfrlen_arr=(1)
expdir_arr=()
corrd_arr=(1)
cutoff_arr=(0)

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
                                            for cutoff in "${cutoff_arr[@]}"; do
                                                expdir=$(mktemp -d ~/tmp/dst_experiment.XXXXXXX)
                                                expdir_arr+=("$expdir")
                                                waitforjobs "$MAX_JOBS"
                                                experiment $nkeys $keylen $nqrys $kdist $qdist $pqratio $membudg "DST" $expdir $corrd $dfsdiff $bfsdiff $cutoff&
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

#SuRF experiments
for nkeys in "${nkeys_arr[@]}"; do
    for keylen in "${keylen_arr[@]}"; do
        for nqrys in "${nqrys_arr[@]}"; do
            for kdist in "${kdist_arr[@]}"; do
                for qdist in "${qdist_arr[@]}"; do
                    for mrange in "${mrange_arr[@]}"; do
                        for pqratio in "${pqratio_arr[@]}"; do
                            for membudg in 10; do
                                for surfhlen in "${surfhlen_arr[@]}"; do
                                    for surfrlen in "${surfrlen_arr[@]}"; do
                                        for corrd in "${corrd_arr[@]}"; do
                                            expdir=$(mktemp -d ~/tmp/dst_experiment.XXXXXXX)
                                            expdir_arr+=("$expdir")
                                            waitforjobs "$MAX_JOBS"
                                            experiment $nkeys $keylen $nqrys $kdist $qdist $pqratio $membudg "SuRF" "$expdir" $corrd $surfhlen $surfrlen &
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
