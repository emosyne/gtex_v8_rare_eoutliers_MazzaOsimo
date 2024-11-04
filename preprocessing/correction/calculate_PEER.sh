#!/bin/bash

set -o nounset -o errexit -o pipefail

## Calculate PEER factors for each tissue.
## The nubmer of PEER factors is determined by the number of samples in the tissue.
## 15 factors for < 150 samples; 30 factors for between 150 and 250 samples; 35 factors for > 250 samples

echo "Calculate PEER factors for each tissue."

peerdir=${TEMPDIR}/preprocessing_v8/PEER_v8
scriptdir=`dirname \$(readlink -f "\$0")`
gtex_eqtl_dir=${GTEX_base}/eqtl/GTEx_Analysis_v8_eQTL

runPeer() {
    traitsFileName=$1
    prefix=${traitsFileName%.tpm.log2.ztrans.txt}
    tissue=`basename $prefix`
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
    outdir=${prefix}_Factors"$maxFactorsN"
    indir=${prefix}_Factors"$maxFactorsN"
    echo $outdir

    mkdir -p $outdir

    ## actual calculation of peer factors
    conda deactivate
    source activate eoutliers_calc_R_env
    
    echo "Rscript calculate_PEER_factors.R $traitsFileName $maxFactorsN $maxIterations $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue" > ${outdir}/log.txt
    Rscript ${scriptdir}/calculate_PEER_factors.R $traitsFileName $maxFactorsN $maxIterations \
            $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue >> ${outdir}/log.txt 2>&1
    
    # computing residuals
    Rscript ${scriptdir}/calculate_PEER_residuals.R $traitsFileName ${peerdir}/covariates.txt \
            ${indir}/factors.tsv ${gtex_eqtl_dir}/${tissue}.v8.egenes.txt.gz \
        $TEMPDIR/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt \
        ${prefix}.peer.v8ciseQTLs.ztrans.txt &> ${outdir}/log.residuals.txt  
    
    conda deactivate
    source activate expr_preprocessing_bash_py2_env
}

export -f runPeer
export scriptdir
export peerdir
export gtex_eqtl_dir

parallel --jobs 10 runPeer ::: ${peerdir}/*.log2.ztrans.txt

echo "DONE!"
