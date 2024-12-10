#!/bin/bash
#SBATCH --job-name=filter_tissues_individuals
#SBATCH --output=slurm_%x_%j.out
#SBATCH -A MURRAY-SL3-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=1:00:00
#SBATCH --mail-user=efo22@cam.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL


#! sbatch directives end here (put any additional directives above this line)




###SCRIPT STARTS FROM HERE



export BASEDIR=/rds/project/rds-qBQA9s264aY/share/eosimo_fmazzarotto/gtex_v8_rare_eoutliers_MazzaOsimo
## Results paths
export WorkDir=${BASEDIR}/temp_workdir
#create folders
cd ${WorkDir}
mkdir -p preprocessing_v8
mkdir -p preprocessing_v8/PEER_v8
mkdir -p data_v8
mkdir -p data_v8/outliers
#define variables
export peerdir=${WorkDir}/preprocessing_v8/PEER_v8
export SAMPLE_TISSUES=${WorkDir}/preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt #created below

# GTEx public data:
export GTEX_base=/home/efo22/murray/share/eosimo_fmazzarotto/resources/DB/GTEx
export GTEX_expr=${GTEX_base}/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
export GTEX_SAMPLES=${GTEX_base}/sample_attrib/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
export GTEX_SUBJECTSv8=${GTEX_base}/sample_attrib/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
export gtex_eqtl_dir=${GTEX_base}/eqtl
export white_ids=${GTEX_base}/sample_attrib/gtex_v10_whiteRace_ids.txt
# this file comes from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz converted into txt with agat:
export gencode_gene_genetypes=/home/efo22/murray/share/eosimo_fmazzarotto/resources/DB/gencode/gencode.v26.GRCh38.genes_genetypes_autosomal.txt
## GTEx restricted data:
export GTEX_WGS=${GTEX_base}/WGS/phg001796.v1.GTEx_v9_WGS_phased.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_944Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz
export GTEX_PCs=${GTEX_base}/sample_attrib/phg001796.v1.GTEx_v9.genotype-qc.MULTI/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_support_files/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.20_Genotype_PCs.eigenvec.txt




source /home/efo22/miniconda3/etc/profile.d/conda.sh



# ## PART1: perform steps up and including PEER correction
### taken from multiple files, mostly based on https://github.com/nmferraro5/gtex_v8_rare_eoutliers/blob/716d43303dd31fae20869d4957b631393fb409c5/correction.md


# source activate expr_preprocessing_bash_py2_env


# #create sample tissues file
# # The generated mapping file excludes flagged individuals and samples.
# # Creates `preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt` and some files under `preprocessing_v8/PEER_v8/`.
# cat $GTEX_SAMPLES | tail -n+2 | cut -f1,7,17 | sed 's/ - /_/' | sed 's/ /_/g' | sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/' | awk '$3=="RNASEQ" {print $1"\t"$2}' > ${SAMPLE_TISSUES}

# head ${SAMPLE_TISSUES}


# #split expression by tissue
# echo "running split_expr_by_tissues.py"
# END='.tpm.txt'
# python2 ${BASEDIR}/preprocessing/split_expr_by_tissues.py --gtex $GTEX_expr --out $peerdir --sample $SAMPLE_TISSUES --end $END

# END='.reads.txt'
# python2 ${BASEDIR}/preprocessing/split_expr_by_tissues.py --gtex $GTEX_expr --out $peerdir --sample $SAMPLE_TISSUES --end $END

# ls ${peerdir} 


# # Creates one file with read counts and one with tpm per tissue in `preprocessing_v8` folder

# ### Transforming data prior to PEER correction

# conda deactivate
# source activate eoutliers_calc_R_env2
# echo "running MazzaOsimo_preprocess_expr.R"
# Rscript ${BASEDIR}/preprocessing/MazzaOsimo_preprocess_expr.R


# ### Generate list of top eQTLs for each gene in each tissue, extract from VCF, convert to number alternative alleles

# conda deactivate
# source activate expr_preprocessing_bash_py2_env
# echo "running get_eqtl_genotypes.sh"
# bash ${BASEDIR}/preprocessing/get_eqtl_genotypes.sh

# ls $WorkDir/preprocessing_v8/

# conda deactivate
# source activate eoutliers_calc_R_env2
# echo "running process_gtex_v8_cis_eqtl_genotypes.R"
# Rscript ${BASEDIR}/preprocessing/process_gtex_v8_cis_eqtl_genotypes.R

# conda deactivate

# # Generates several intermediate files in `preprocessing_v8` and relies on `process_gtex_v8_cis_eqtl_genotypes.R` to generate final `gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt` in `preprocessing_v8`



### Actually run PEER correction and compute residuals
# Creates a file for each tissue  under `preprocessing/PEER_v8/` with scaled and corrected log2(tpm) values.

# source activate eoutliers_calc_R_env2
# echo "calculating peer factors"
# # this script takes more than 36 hours running in parallel 10 files at a time on 32 nodes:
# bash ${BASEDIR}/preprocessing/correction/1_calculate_PEER_factors.sh
# # Relies on `preprocessing/correction/calculate_PEER_factors.R` 

# echo "calculating peer residuals"
# bash ${BASEDIR}/preprocessing/correction/2_calculate_PEER_residuals.sh
# # Relies on `preprocessing/correction/calculate_PEER_residuals.R`.



# # ### Generate files with data on what tissues are available per individual
# 
# bash ${BASEDIR}/preprocessing/get_tissue_by_individual.sh
# # Generates `preprocessing_v8/gtex_tissues_all_normalized_samples.txt` and `preprocessing_v8/gtex_individuals_all_normalized_samples.txt`


# ### Combine PEER-corrected data into a single flat file (and compress output file)
# conda deactivate
# source activate expr_preprocessing_bash_py2_env
# python2 ${BASEDIR}/preprocessing/gather_filter_normalized_expression.py 
# gzip ${WorkDir}/preprocessing_v8/gtex_normalized_expression.txt
# # Creates and compresses `preprocessing_v8/gtex_normalized_expression.txt.gz`.



# in order to now create normalized_expression_v8ciseQTLs_removed.txt by removing tissue and global expression outliers from 
# `preprocessing_v8/gtex_normalized_expression.txt.gz`, one needs to calculate /call outliers, originally done in a separate md file:

## PART 2: Outlier calling 
##  taken from https://github.com/nmferraro5/gtex_v8_rare_eoutliers/blob/716d43303dd31fae20869d4957b631393fb409c5/outlier_calling.md


### Run outlier calling
# conda deactivate
# source activate eoutliers_calc_R_env2
# Rscript ${BASEDIR}/outlier_calling/call_outliers.R \
#         --Z.SCORES=${WorkDir}/preprocessing_v8/gtex_normalized_expression.txt.gz \
#         --OUT.PREFIX=${WorkDir}/data_v8/outliers/gtex.outlier.controls.v8ciseQTLs \
#         --N.PHEN=5 

# # Generates one file (${WorkDir}/data_v8/outliers/gtex.outlier.controls.v8ciseQTLs.medz.txt) 
# # with columns gene ind N Df MedZ Y, where N refers to the number of individuals tested for a given gene, Df refers to the number of measurements 
# # for that gene in that individual and Y indicates outlier or control. 
# # The `N.PHEN` argument specifies the minimum number of measurements (i.e. tissues) available for a gene-individual required for inclusion. We use 5. 
# # `GLOBAL` specifies a file with a list of individual IDs to remove as global outliers. If NA, no individuals are removed.




### Identify global outlier individuals 
# # OULIERS: path to the Z-score data: 
# Rscript ${BASEDIR}/outlier_calling/identify_global_outliers.R \
#         --OUTLIERS=${WorkDir}/data_v8/outliers/gtex.outlier.controls.v8ciseQTLs.medz.txt \
#         --METHOD=proportion

# # Removes global outlier individuals, either defined based on the proportion of genes called as outliers per individual relative to the population or 
# the number of genes called as outliers, determined by setting 'proportion' or 'number' in `METHOD`. 
# Writes out a new outlier file:
        #  data_v8/outliers/gtex.outlier.controls.v8ciseQTLs.medz_globalOutliersRemoved.txt.
# I've also added two outputs:
        # preprocessing_v8/gtex_normalized_expression_global_outliers_removed.txt.gz",
        # data_v8/outliers/globalOutliers.txt

# ### From the multi-tissue outliers, determine which tissues have extreme effects and outlier sharing across tissues
# Rscript ${BASEDIR}/outlier_calling/extract_extreme_tissues.R \
#     --Z.SCORES=${WorkDir}/data_v8/outliers/gtex.outlier.controls.v8ciseQTLs.medz.txt \
#     --EXP.DATA=${WorkDir}/preprocessing_v8/gtex_normalized_expression_global_outliers_removed.txt.gz
# # Generates figures in `figures/GTEXv8_pair_jaccard.pdf`.


### Select tissues and individuals for downstream analyses (still from correction.md)
conda deactivate
source activate R_dplyr
Rscript ${BASEDIR}/preprocessing/filter_tissues_individuals.R
# Must be run from the upper level directory of the repo (e.g., the location of this readme).
# Generates `preprocessing_v8/gtex_2017-06-05_v8_design_passed.txt`, `preprocessing_v8/gtex_2017-06-05_v8_individuals_passed.txt`, `preprocessing_v8/gtex_2017-06-05_v8_tissues_passed.txt`, 
# and `preprocessing_v8/gtex_2017-06-05_v8_normalized_expression.subset.txt.gz`. 
# Also produces summary figures `figures/gtex_v8_design.pdf`. 
# The subset file filtered for missingness is used in correlation-outlier calling. No missingness filter is applied for multi-tissue eOutlier calling.
