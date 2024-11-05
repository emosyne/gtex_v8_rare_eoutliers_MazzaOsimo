#!/bin/bash

### Generate list of top eQTLs for each gene in each tissue
#public data
zcat $GTEX_base/GTEx_Analysis_v8_eQTL/*.v8.egenes.txt.gz | \
        cut -f14-17 | grep -v variant_pos | sed 's/^X/23/g' | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3,$4}' | sort -k1,1 -k2,2n | \
        sed 's/^23/X/g' | uniq > $WorkDir/preprocessing_v8/gtex_2017-06-05_v8_cis_eQTLs.bed         #change sorting so it matches vcf i.e. chr1, chr2, etc.

wait

### Extract these sites from the GTEx v8 VCF using bedtools
zcat $GTEX_WGS | head -4000 | grep '#' > $WorkDir/preprocessing_v8/gtex_v9_genotypes_cis_eQTLs.vcf

wait

bedtools intersect -a $WorkDir/preprocessing_v8/gtex_2017-06-05_v8_cis_eQTLs.bed -b $GTEX_WGS -loj -sorted | \
        grep -v '\-1' | cut -f6- | sort -k1,1 -k2,2n | uniq >> $WorkDir/preprocessing_v8/gtex_v9_genotypes_cis_eQTLs.vcf

wait

bgzip -f $WorkDir/preprocessing_v8/gtex_v9_genotypes_cis_eQTLs.vcf

wait

tabix -f -p vcf $WorkDir/preprocessing_v8/gtex_v9_genotypes_cis_eQTLs.vcf.gz

wait
### Convert the cis-eQTL genotypes in VCF format to the number of alternate alleles using VCFTools

vcftools --gzvcf $WorkDir/preprocessing_v8/gtex_v9_genotypes_cis_eQTLs.vcf.gz --out $WorkDir/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs --012 --maf 0.01


