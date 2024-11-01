#!/bin/bash 

#SBATCH -J eoutliers_calc_scz_mazzaosimo
#SBATCH -A MURRAY-SL2-CPU
#SBATCH -p cclake
#SBATCH --nodes=1
#! The Cascade Lake (cclake) nodes have 56 CPUs (cores) each and
#! 3420 MiB of memory per CPU.
#SBATCH --ntasks=16
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --mail-user=efo22@cam.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

#! sbatch directives end here (put any additional directives above this line)




###SCRIPT STARTS FROM HERE


source activate expr_preprocessing2

#define variables
export BASEDIR=/rds/project/rds-qBQA9s264aY/share/eosimo_fmazzarotto/gtex_v8_rare_eoutliers_MazzaOsimo
## Results paths
export TEMPDIR=${BASEDIR}/temp_workdir

#create folders
cd ${TEMPDIR}
mkdir -p data_v8
mkdir -p features_v8
mkdir -p figures
mkdir -p paper_figures
mkdir -p preprocessing_v8
mkdir -p preprocessing_v8/PEER_v8

export GTEX_base=/home/efo22/murray/share/eosimo_fmazzarotto/resources/DB/GTEx
export OUT=${TEMPDIR}/preprocessing_v8/PEER_v8/
export GTEX_expr=${GTEX_base}/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
export GTEX_SAMPLES=${GTEX_base}/sample_attrib/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
export SAMPLE_TISSUES=${TEMPDIR}/preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt




#create sample tissues file
# The generated mapping file excludes flagged individuals and samples.
# Creates `preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt` and some files under `preprocessing_v8/PEER_v8/`.
cat $GTEX_SAMPLES | tail -n+2 | cut -f1,7,17 | sed 's/ - /_/' | sed 's/ /_/g' | sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/' | awk '$3=="RNASEQ" {print $1"\t"$2}' > ${SAMPLE_TISSUES}

head ${SAMPLE_TISSUES}

#split expression by tissue
python2 $BASEDIR/preprocessing/split_expr_by_tissues.py --gtex $GTEX_expr --out $OUT --sample $SAMPLE_TISSUES --end '.reads.txt'

ll $OUT


# Creates one file with read counts and one with tpm per tissue in `preprocessing_v8` folder

### Transforming data prior to PEER correction
export GTEX_SUBJECTSv8=${GTEX_base}/sample_attrib/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
export GTEX_PCs=${GTEX_base}/sample_attrib/phg001796.v1.GTEx_v9.genotype-qc.MULTI/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_support_files/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.20_Genotype_PCs.eigenvec.txt

Rscript ../preprocessing/MazzaOsimo_preprocess_expr.R


### Generate list of top eQTLs for each gene in each tissue, extract from VCF, convert to number alternative alleles
## restricted data:
export GTEX_WGS=${GTEX_base}/WGS/phg001796.v1.GTEx_v9_WGS_phased.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_944Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz



bash ../preprocessing/get_eqtl_genotypes.sh
# Rscript process_gtex_v8_cis_eqtl_genotypes.R

# Generates several intermediate files in `preprocessing_v8` and relies on `process_gtex_v8_cis_eqtl_genotypes.R` to generate final `gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt` in `preprocessing_v8`

# ### Actually run PEER correction and compute residuals

# bash preprocessing/correction/calculate_PEER.sh

# Relies on `preprocessing/correction/calculate_PEER_factors.R` and `preprocessing/correction/calculate_PEER_residuals.R`.

# Creates a file for each tissue  under `preprocessing/PEER_v8/` with scaled and corrected log2(tpm) values.

# ### Combine PEER-corrected data into a single flat file (and compress output file)

# python gather_filter_normalized_expression.py 
# gzip ${TEMPDIR}/preprocessing_v8/gtex_2017-06-05_normalized_expression.txt


# Creates and compresses `preprocessing_v8/gtex_2017-06-05_normalized_expression.txt.gz`.

# ### Generate files with data on what tissues are available per individual

# bash get_tissue_by_individual.sh

# Generates `gtex_2017-06-05_tissues_all_normalized_samples.txt` and `gtex_2017-06-05_tissues_all_normalized_samples.txt` in `preprocessing_v8`

# ### Select tissues and individuals for downstream analyses

# Rscript preprocessing/filter_tissues_individuals.R

# Must be run from the upper level directory of the repo (e.g., the location of this readme).

# Generates `preprocessing_v8/gtex_2017-06-05_v8_design_passed.txt`, `preprocessing_v8/gtex_2017-06-05_v8_individuals_passed.txt`, `preprocessing_v8/gtex_2017-06-05_v8_tissues_passed.txt`, and `preprocessing_v8/gtex_2017-06-05_v8_normalized_expression.subset.txt.gz`. Also produces summary figures `figures/gtex_v8_design.pdf`. The subset file filtered for missingness is used in correlation-outlier calling. No missingness filter is applied for multi-tissue eOutlier calling.

