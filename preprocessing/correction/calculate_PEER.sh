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
            $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue" | tee -a ${outdir}/log.txt
    Rscript ${scriptdir}/correction/calculate_PEER_factors.R $traitsFileName $maxFactorsN $maxIterations \
            $boundTol $varTol $e_pa $e_pb $a_pa $a_pb $outdir $tissue >> ${outdir}/log.txt 2>&1
    
    # computing residuals
    echo "computing residuals Rscript ${scriptdir}/correction/calculate_PEER_residuals.R ${traitsFileName} ${peerdir}/covariates.txt \
            ${indir}/factors.tsv ${gtex_eqtl_dir}/${tissue}.v8.egenes.txt.gz \
            ${WorkDir}/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed_withdups.txt \
            ${tissue}.peer.v8ciseQTLs.ztrans.txt"
    Rscript ${scriptdir}/correction/calculate_PEER_residuals.R ${traitsFileName} ${peerdir}/covariates.txt \
            ${indir}/factors.tsv ${gtex_eqtl_dir}/${tissue}.v8.egenes.txt.gz \
            ${WorkDir}/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed_withdups.txt \
            ${tissue}.peer.v8ciseQTLs.ztrans.txt | tee -a ${outdir}/log.residuals.txt  
}

export -f runPeer

# # parallel --jobs 10 runPeer ::: ${peerdir}/*.log2.ztrans.txt
# Process the First 10 Files
parallel --jobs 10 runPeer ::: $(ls ${peerdir}/*.log2.ztrans.txt | head -n 1)
# parallel --jobs 10 runPeer ::: $(ls ${peerdir}/*.log2.ztrans.txt | tail -n +11 | head -n 10)
# parallel --jobs 10 runPeer ::: $(ls ${peerdir}/*.log2.ztrans.txt | tail -n +21 | head -n 10)
# parallel --jobs 10 runPeer ::: $(ls ${peerdir}/*.log2.ztrans.txt | tail -n +31 | head -n 10)
# parallel --jobs 10 runPeer ::: $(ls ${peerdir}/*.log2.ztrans.txt | tail -n +41 | head -n 10)

# # Get the list of files
# export files=(${peerdir}/*.log2.ztrans.txt)

# # Loop through the files in chunks of 10
# for ((i=0; i<${#files[@]}; i+=10)); do
#   sbatch <<EOT
# #!/bin/bash
# #SBATCH -J peer_$i
# #SBATCH -A MURRAY-SL3-CPU
# #SBATCH -p cclake
# #SBATCH --nodes=1
# #SBATCH --ntasks=16
# #SBATCH --mem=32G
# #SBATCH --time=12:00:00
# #SBATCH --mail-user=efo22@cam.ac.uk
# #SBATCH --mail-type=BEGIN,END,FAIL
# #SBATCH --output=${BASEDIR}/peer_job_$i.out
# #SBATCH --error=${BASEDIR}/peer_job_$i.err

# source /home/efo22/miniconda3/etc/profile.d/conda.sh
# source activate eoutliers_calc_R_env2

# parallel --jobs 10 runPeer ::: "${files[@]:i:10}"
# EOT
# done

echo "DONE!"
