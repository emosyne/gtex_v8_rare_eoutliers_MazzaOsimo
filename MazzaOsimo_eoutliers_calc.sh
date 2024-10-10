#!/bin/bash 

#SBATCH -J eoutliers_calc
#SBATCH -A MURRAY-SL2-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#! The Cascade Lake (cclake) nodes have 56 CPUs (cores) each and
#! 3420 MiB of memory per CPU.
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=1:00:00

#! sbatch directives end here (put any additional directives above this line)




###SCRIPT STARTS FROM HERE


module load miniconda/3
source activate expr_preprocessing

#define variables
BASEDIR=/rds/project/rds-qBQA9s264aY/share/eosimo_fmazzarotto/gtex_v8_rare_eoutliers_MazzaOsimo
TEMPDIR=${BASEDIR}/temp_workdir

#create folders
cd ${TEMPDIR}
mkdir -p data_v8
mkdir -p features_v8
mkdir -p figures
mkdir -p paper_figures
mkdir -p preprocessing_v8
mkdir -p preprocessing_v8/PEER_v8

GTEX_RNAv8=/rds/project/rds-qBQA9s264aY/share/eosimo_fmazzarotto/resources/DB/GTEx_expr
OUT=${TEMPDIR}/preprocessing_v8/PEER_v8
GTEX=${GTEX_RNAv8}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
GTEX_SAMPLES=${GTEX_RNAv8}/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
SAMPLE_TISSUES=${TEMPDIR}/preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt
END='.tpm.txt'



#create sample tissues file
cat $GTEX_SAMPLES | tail -n+2 | cut -f1,7,17 | sed 's/ - /_/' | sed 's/ /_/g' | sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/' | awk '$3=="RNASEQ" {print $1"\t"$2}' > ${SAMPLE_TISSUES}

#split expression by tissue
python2 $BASEDIR/preprocessing/split_expr_by_tissues.py --gtex $GTEX --out $OUT --sample $SAMPLE_TISSUES --end $END

# ```
# Creates one file with read counts and one with tpm per tissue in `preprocessing_v8` folder

# ### Transforming data prior to PEER correction
# ```
# Rscript preprocessing/preprocess_expr.R
# ```
# The generated mapping file excludes flagged individuals and samples.
# Creates `preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt` and some files under `preprocessing_v8/PEER_v8/`.

# ### Generate list of top eQTLs for each gene in each tissue, extract from VCF, convert to number alternative alleles
# ```
# bash get_eqtl_genotypes.sh
# ```
# Generates several intermediate files in `preprocessing_v8` and relies on `process_gtex_v8_cis_eqtl_genotypes.R` to generate final `gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt` in `preprocessing_v8`

# ### Actually run PEER correction and compute residuals
# ```
# bash preprocessing/correction/calculate_PEER.sh
# ```
# Relies on `preprocessing/correction/calculate_PEER_factors.R` and `preprocessing/correction/calculate_PEER_residuals.R`.

# Creates a file for each tissue  under `preprocessing/PEER_v8/` with scaled and corrected log2(tpm) values.

# ### Combine PEER-corrected data into a single flat file (and compress output file)
# ```
# python gather_filter_normalized_expression.py 
# gzip ${RAREDIR}/preprocessing_v8/gtex_2017-06-05_normalized_expression.txt
# ```

# Creates and compresses `preprocessing_v8/gtex_2017-06-05_normalized_expression.txt.gz`.

# ### Generate files with data on what tissues are available per individual
# ```
# bash get_tissue_by_individual.sh
# ```
# Generates `gtex_2017-06-05_tissues_all_normalized_samples.txt` and `gtex_2017-06-05_tissues_all_normalized_samples.txt` in `preprocessing_v8`

# ### Select tissues and individuals for downstream analyses
# ```
# Rscript preprocessing/filter_tissues_individuals.R
# ```
# Must be run from the upper level directory of the repo (e.g., the location of this readme).

# Generates `preprocessing_v8/gtex_2017-06-05_v8_design_passed.txt`, `preprocessing_v8/gtex_2017-06-05_v8_individuals_passed.txt`, `preprocessing_v8/gtex_2017-06-05_v8_tissues_passed.txt`, and `preprocessing_v8/gtex_2017-06-05_v8_normalized_expression.subset.txt.gz`. Also produces summary figures `figures/gtex_v8_design.pdf`. The subset file filtered for missingness is used in correlation-outlier calling. No missingness filter is applied for multi-tissue eOutlier calling.

