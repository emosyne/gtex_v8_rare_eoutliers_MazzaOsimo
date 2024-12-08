#!/bin/bash

set -o nounset -o errexit -o pipefail

## Calculate PEER factors for each tissue.
## The nubmer of PEER factors is determined by the number of samples in the tissue.
## 15 factors for < 150 samples; 30 factors for between 150 and 250 samples; 35 factors for > 250 samples

echo "Calculate PEER factors for each tissue."


runPeer() {
    echo $1
    traitsFileName=$1
    prefix=${traitsFileName%.tpm.log2.ztrans.txt}
    tissue=$(basename ${prefix})  # Extract the filename
    tissue=${tissue%.log2.ztrans.txt}  # Remove the .log2.ztrans.txt part to get the tissue name
    nsamples=$(cat $traitsFileName | wc -l) # this is actually n samples + 1
    if [ $nsamples -le 150 ]; then
        maxFactorsN=15
    elif [ $nsamples -le 249 ]; then
        maxFactorsN=30
    elif [ $nsamples -le 349 ]; then
        maxFactorsN=45
    else
        maxFactorsN=60
    fi
    maxIterations=10000
    boundTol=0.001
    varTol=0.00001
    e_pa=0.1
    e_pb=10
    a_pa=0.001
    a_pb=0.1
    outdir=${WorkDir}/preprocessing_v8/PEER_v8/${tissue}_Factors"$maxFactorsN"
    indir=${WorkDir}/preprocessing_v8/PEER_v8/${tissue}_Factors"$maxFactorsN"
    echo $outdir

    mkdir -p $outdir

    ## actual calculation of peer factors
    echo "Rscript ${scriptdir}/correction/calculate_PEER_factors.R $traitsFileName $maxFactorsN $maxIterations \
            $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue" 
    Rscript ${scriptdir}/correction/calculate_PEER_factors.R $traitsFileName $maxFactorsN $maxIterations \
            $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue 2>&1 | tee -a ${outdir}/log_calculate_PEER_factors.txt
    
   
}

export -f runPeer




parallel --jobs 10 runPeer ::: ${peerdir}/*.log2.ztrans.txt

# Process the First 10 Files
# parallel --jobs 10 runPeer ::: $(ls ${peerdir}/*.log2.ztrans.txt | head -n 10)
# parallel --jobs 10 runPeer ::: $(ls ${peerdir}/*.log2.ztrans.txt | tail -n +11 | head -n 10)
# parallel --jobs 10 runPeer ::: $(ls ${peerdir}/*.log2.ztrans.txt | tail -n +21 | head -n 10)
# parallel --jobs 10 runPeer ::: $(ls ${peerdir}/*.log2.ztrans.txt | tail -n +31 | head -n 10)
# parallel --jobs 10 runPeer ::: $(ls ${peerdir}/*.log2.ztrans.txt | tail -n +41 | head -n 10)


echo "DONE!"
