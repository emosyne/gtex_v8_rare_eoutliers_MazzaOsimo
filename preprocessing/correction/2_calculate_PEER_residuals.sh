#!/bin/bash

set -o nounset -o errexit -o pipefail

## Calculate PEER factors for each tissue.
## The nubmer of PEER factors is determined by the number of samples in the tissue.
## 15 factors for < 150 samples; 30 factors for between 150 and 250 samples; 35 factors for > 250 samples

echo "Calculate PEER residuals for each tissue."


runPeer() {
    echo $1
    traitsFileName=$1
    # export traitsFileName="${peerdir}/Whole_Blood.log2.ztrans.txt"
    prefix=${traitsFileName%.tpm.log2.ztrans.txt}
    tissue=$(basename ${prefix})  # Extract the filename
    export tissue=${tissue%.log2.ztrans.txt}  # Remove the .log2.ztrans.txt part to get the tissue name
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
    export indir=${WorkDir}/preprocessing_v8/PEER_v8/${tissue}_Factors"$maxFactorsN"
    echo $outdir

    if [ ! -d "$outdir" ]; then
        echo "Directory $outdir does not exist. Skipping..."
        continue
    fi

    # computing residuals
    echo "computing residuals Rscript ${scriptdir}/correction/calculate_PEER_residuals.R ${traitsFileName} ${peerdir}/covariates.txt \
            ${indir}/factors.tsv ${gtex_eqtl_dir}/${tissue}.v8.egenes.txt.gz \
            ${WorkDir}/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed_withdups.txt \
            ${tissue}.peer.v8ciseQTLs.ztrans.txt"
    Rscript ${scriptdir}/correction/calculate_PEER_residuals.R ${traitsFileName} ${peerdir}/covariates.txt \
            ${indir}/factors.tsv ${gtex_eqtl_dir}/${tissue}.v8.egenes.txt.gz \
            ${WorkDir}/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed_withdups.txt \
            ${WorkDir}/preprocessing_v8/PEER_v8/${tissue}.peer.v8ciseQTLs.ztrans.txt 2>&1 | tee -a ${outdir}/${tissue}_calculate_PEER_residuals_log.txt  
}

export -f runPeer




parallel --jobs 10 runPeer ::: ${peerdir}/*.log2.ztrans.txt
# parallel --jobs 10 runPeer ::: ${peerdir}/Whole_Blood.log2.ztrans.txt ${peerdir}/Adipose_Subcutaneous ${peerdir}/ ${peerdir}/



echo "DONE!"
