source activate expr_preprocessing_bash_py2_env

#create sample tissues file
# The generated mapping file excludes flagged individuals and samples.
# Creates `preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt` and some files under `preprocessing_v8/PEER_v8/`.
cat $GTEX_SAMPLES | tail -n+2 | cut -f1,7,17 | sed 's/ - /_/' | sed 's/ /_/g' | sed 's/(//' | sed 's/)//' | sed 's/c-1/c1/' | awk '$3=="RNASEQ" {print $1"\t"$2}' > ${SAMPLE_TISSUES}

head ${SAMPLE_TISSUES}



#split expression by tissue
echo "running split_expr_by_tissues.py"
END='.tpm.txt'
python2 ${scriptdir}/split_expr_by_tissues.py --gtex $GTEX_expr --out $peerdir --sample $SAMPLE_TISSUES --end $END

END='.reads.txt'
python2 ${scriptdir}/split_expr_by_tissues.py --gtex $GTEX_expr --out $peerdir --sample $SAMPLE_TISSUES --end $END

ls ${peerdir} 


# Creates one file with read counts and one with tpm per tissue in `preprocessing_v8` folder

### Transforming data prior to PEER correction

conda deactivate
source activate eoutliers_calc_R_env
echo "running MazzaOsimo_preprocess_expr.R"
Rscript ${scriptdir}/MazzaOsimo_preprocess_expr.R


### Generate list of top eQTLs for each gene in each tissue, extract from VCF, convert to number alternative alleles

conda deactivate
source activate expr_preprocessing_bash_py2_env
echo "running get_eqtl_genotypes.sh"
bash ${scriptdir}/get_eqtl_genotypes.sh

ls $WorkDir/preprocessing_v8/

conda deactivate
source activate eoutliers_calc_R_env
echo "running process_gtex_v8_cis_eqtl_genotypes.R"
Rscript ${scriptdir}/process_gtex_v8_cis_eqtl_genotypes.R

# Generates several intermediate files in `preprocessing_v8` and relies on `process_gtex_v8_cis_eqtl_genotypes.R` to generate final `gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt` in `preprocessing_v8`

