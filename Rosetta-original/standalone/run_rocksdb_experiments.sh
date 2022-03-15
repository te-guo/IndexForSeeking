#!/bin/bash

#
# These are the rocksdb experiments
#
# pre-req: range_query
#

if [ "$#" != "5" ]; then
    echo "Usage: $0 [path to experiment binary] [path to workload generator binary] [scratch directory to perform experiments in] [scratch directory to store experiment data in] [bulk]"
    exit
fi

EXPER_BIN="$(readlink -f $1)" && shift
WORKL_BIN="$(readlink -f $1)" && shift
SCRATCH_DIR="$(readlink -f $1)" && shift
DATA_DIR="$(readlink -f $1)" && shift
BULK_DIR="$(readlink -f $1)" && shift
MAX_JOBS=1

waitforjobs() {
    while test $(jobs -p | wc -w) -ge "$1"; do wait -n; done
}

fetch_data() {
    nkeys=$1 && shift
    keylen=$1 && shift
    nqrys=$1 && shift
    mrange=$1 && shift
    pqratio=$1 && shift
    kdist=$1 && shift
    qdist=$1 && shift
    corrd=$1 && shift

    path="$DATA_DIR/$nkeys/$keylen/$nqrys/$mrange/$pqratio/$kdist/$qdist/$corrd"

    if [ ! -e "$path/my_data" ]; then
	(
	    mkdir -p "$path" && cd "$path"
	    "$WORKL_BIN" "$nkeys" "$keylen" "$nqrys" "$mrange" "$pqratio" "$kdist" "$qdist" "$corrd"
	    shuf ./my_data/data.txt -o ./my_data/data.txt
	    mkfifo rand1 rand2
	    tee rand1 rand2 < /dev/urandom > /dev/null &
	    shuf --random-source=rand1 ./my_data/txn.txt -o ./my_data/txn.txt &
	    shuf --random-source=rand2 ./my_data/upper_bound.txt -o ./my_data/upper_bound.txt &
	    wait
	)
    fi
    cp -r "$path/my_data" ./
}

experiment() {
    nkeys=$1 && shift
    keylen=$1 && shift
    nqrys=$1 && shift
    kdist=$1 && shift
    qdist=$1 && shift
    mrange=$1 && shift
    dfsdiff=$1 && shift
    bfsdiff=$1 && shift
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
    echo -e "\tCorrelation degree:\t$corrd" >> ./experiment_result
    echo -e "\tFilter bits-per-key:\t$membudg" >> ./experiment_result

    fetch_data "$nkeys" "$keylen" "$nqrys" "$mrange" "$pqratio" "$kdist" "$qdist" "$corrd"

    if [ $filter = "SuRF" ]; then
	surftype=$1 && shift
        echo -e "\tSuRF type:\t$surftype" >> ./experiment_result
        echo -e "### END EXPERIMENT DESCRIPTION ###\n\n" >> ./experiment_result
        "$EXPER_BIN"_bulk_write "surf" "write" "$nkeys" "$nqrys" ./db/ "$BULK_DIR" ./my_data/ "$membudg" 0 "$surftype" >> ./experiment_result && "$EXPER_BIN" "surf" "read" "$nkeys" "$nqrys" ./db/ ./my_data/ "$membudg" 1 "$surftype" >> ./experiment_result
    elif [ $filter = "DST" ]; then
        echo -e "\tBfs diffidence:\t$bfsdiff" >> ./experiment_result
        echo -e "### END EXPERIMENT DESCRIPTION ###\n\n" >> ./experiment_result
        "$EXPER_BIN"_bulk_write "dst" "write" "$nkeys" "$nqrys" ./db/ "$BULK_DIR" ./my_data/ "$membudg" 0 "$bfsdiff" >> ./experiment_result && "$EXPER_BIN" "dst" "read" "$nkeys" "$nqrys" ./db/ ./my_data/ "$membudg" 1 "$bfsdiff" >> ./experiment_result
    elif [ $filter = "None" ]; then
        echo -e "### END EXPERIMENT DESCRIPTION ###\n\n" >> ./experiment_result
        "$EXPER_BIN"_bulk_write "none" "write" "$nkeys" "$nqrys" ./db/ "$BULK_DIR" ./my_data/ "$membudg" 0 >> ./experiment_result && "$EXPER_BIN" "none" "read" "$nkeys" "$nqrys" ./db/ ./my_data/ "$membudg" 1 >> ./experiment_result
    elif [ $filter = "Prefix" ]; then
	prefix=$1 && shift
        echo -e "### END EXPERIMENT DESCRIPTION ###\n\n" >> ./experiment_result
        "$EXPER_BIN"_bulk_write "prefix" "write" "$nkeys" "$nqrys" ./db/ "$BULK_DIR" ./my_data/ "$membudg" 0 "$prefix" >> ./experiment_result && "$EXPER_BIN" "prefix" "read" "$nkeys" "$nqrys" ./db/ ./my_data/ "$membudg" 1 "$prefix" >> ./experiment_result
    fi
    echo -e "### END EXPERIMENT ###" >> ./experiment_result

    rm -rf ./db/
    rm -rf ./my_data/
    rm -rf ./bulk
}

nkeys_arr=(100000000)
keylen_arr=(64)
nqrys_arr=(100000)
kdist_arr=("uniform")
qdist_arr=("uniform")
mrange_arr=(32)
dfsdiff_arr=(10000)
bfsdiff_arr=(64)
pqratio_arr=(0.9375 0.5)
membudg_arr=(10)
#membudg_arr=(0 10)
surfhlen_arr=(0)
surfrlen_arr=(0)
surftype_arr=(2)
expdir_arr=()
corrd_arr=(1)
prefix_arr=(3)

# #DST experiments
# for nkeys in "${nkeys_arr[@]}"; do
#     for keylen in "${keylen_arr[@]}"; do
#         for nqrys in "${nqrys_arr[@]}"; do
#             for kdist in "${kdist_arr[@]}"; do
#                 for qdist in "${qdist_arr[@]}"; do
#                     for mrange in "${mrange_arr[@]}"; do
#                         for dfsdiff in "${dfsdiff_arr[@]}"; do
#                             for bfsdiff in "${bfsdiff_arr[@]}"; do
#                                 for pqratio in "${pqratio_arr[@]}"; do
#                                     for membudg in "${membudg_arr[@]}"; do
#                                         for corrd in "${corrd_arr[@]}"; do
#                                             expdir=$(mktemp -d "$SCRATCH_DIR"/dst_experiment.XXXXXXX)
#                                             expdir_arr+=("$expdir")
#                                             waitforjobs "$MAX_JOBS"
#                                             experiment $nkeys $keylen $nqrys $kdist $qdist $mrange $dfsdiff $bfsdiff $pqratio $membudg "DST" "$expdir" "$corrd" &
#                                         done
#                                     done
#                                 done
#                             done
#                         done
#                     done
#                 done
#             done
#         done
#     done
# done

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
				    for membudg in "${membudg_arr[@]}"; do
					for corrd in "${corrd_arr[@]}"; do
					    for surftype in "${surftype_arr[@]}"; do
						expdir=$(mktemp -d "$SCRATCH_DIR"/dst_experiment.XXXXXXX)
						expdir_arr+=("$expdir")
						waitforjobs "$MAX_JOBS"
						experiment $nkeys $keylen $nqrys $kdist $qdist $mrange $dfsdiff $bfsdiff $pqratio $membudg "SuRF" "$expdir" $corrd $surftype &
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

# #Prefix-Bloom-Filter experiments
# for nkeys in "${nkeys_arr[@]}"; do
#     for keylen in "${keylen_arr[@]}"; do
#         for nqrys in "${nqrys_arr[@]}"; do
#             for kdist in "${kdist_arr[@]}"; do
#                 for qdist in "${qdist_arr[@]}"; do
#                     for mrange in "${mrange_arr[@]}"; do
#                         for dfsdiff in 1000; do
# 			    for bfsdiff in 32; do
# 				for pqratio in "${pqratio_arr[@]}"; do
# 				    for membudg in "${membudg_arr[@]}"; do
# 					for corrd in "${corrd_arr[@]}"; do
# 					    for prefix in "${prefix_arr[@]}"; do
# 						expdir=$(mktemp -d "$SCRATCH_DIR"/dst_experiment.XXXXXXX)
# 						expdir_arr+=("$expdir")
# 						waitforjobs "$MAX_JOBS"
# 						experiment $nkeys $keylen $nqrys $kdist $qdist $mrange $dfsdiff $bfsdiff $pqratio $membudg "Prefix" "$expdir" $corrd $prefix &
# 					    done
# 					done
# 				    done
# 				done
# 			    done
# 			done
# 		    done
# 		done
# 	    done
# 	done
#     done
# done

# #No Filter experiments
# for nkeys in "${nkeys_arr[@]}"; do
#     for keylen in "${keylen_arr[@]}"; do
#         for nqrys in "${nqrys_arr[@]}"; do
#             for kdist in "${kdist_arr[@]}"; do
#                 for qdist in "${qdist_arr[@]}"; do
#                     for mrange in "${mrange_arr[@]}"; do
#                         for dfsdiff in 1000; do
# 			    for bfsdiff in 32; do
# 				for pqratio in "${pqratio_arr[@]}"; do
# 				    for membudg in "${membudg_arr[@]}"; do
# 					for corrd in "${corrd_arr[@]}"; do
# 					    expdir=$(mktemp -d "$SCRATCH_DIR"/dst_experiment.XXXXXXX)
# 					    expdir_arr+=("$expdir")
# 					    waitforjobs "$MAX_JOBS"
# 					    experiment $nkeys $keylen $nqrys $kdist $qdist $mrange $dfsdiff $bfsdiff $pqratio $membudg "None" "$expdir" $corrd $prefix &
# 					done
# 				    done
# 				done
# 			    done
# 			done
# 		    done
# 		done
# 	    done
# 	done
#     done
# done

wait
for expdir in "${expdir_arr[@]}"; do
    cat "$expdir/experiment_result"
done
