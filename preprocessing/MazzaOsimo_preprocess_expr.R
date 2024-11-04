#!/usr/bin/env Rscript

## R script to split the entire GTEx matrices into files for each tissue
## then process RPKM and read count matrices so that columns match covariate files.
## Prepare matrices for input to PEER.


# A set of covariates identified using the Probabilistic Estimation of Expression Residuals (PEER) method                       (Stegle et al., PLoS Comp. Biol., 2010 ), calculated for the normalized expression matrices (described below). For eQTL analyses, the number of PEER factors was determined as function of sample size (N): 15 factors for N<150, 30 factors for 150≤ N<250, 45 factors for 250≤ N<350, and 60 factors for N≥350, as a result of optimizing for the number of eGenes discovered. For sQTL analyses, 15 PEER factors were computed for each tissue.

## Load required packages
require(data.table)
require(ggplot2)
require(stringr)

##------------- FUNCTIONS



## For each tissue, read in the TPM and read files from the given directory.
## Subset and reorder columns of TPM file to match the given covariates.
## Filter for genes with > 20% individuals with TPM > 0.1 and read count > 6.
## Log2(tpm + 2) transform the data, then z-transform.
## Finally, output the transformed matrix to a new file.
ztrans.tissue = function(tissue, dir, covs, read.filt = 6, tpm.filt = 0.1) {
    tpm = fread(paste0(dir, tissue, '.tpm.txt'))
    reads = fread(paste0(dir, tissue, '.reads.txt'))
    ## sanity checks
    stopifnot(sum(colnames(tpm) != colnames(reads)) == 0)
    stopifnot(sum(tpm$Gene != reads$Gene) == 0)

    covariates.subset = covs$SUBJID[covs$SUBJID %in% colnames(tpm)]
    genes = tpm$Gene
    tpm = tpm[, c(covariates.subset), with = F]
    reads = reads[, c(covariates.subset), with = F]
    ind.filt = round(0.2*ncol(tpm))
    zero.filt = round(0.05*ncol(tpm))
    indices.keep = rowSums(tpm > tpm.filt & reads > read.filt) >= ind.filt
    zero.keep = rowSums(tpm > 0) >= zero.filt
    indices.keep = indices.keep[which(indices.keep %in% zero.keep)]
    tpm = tpm[indices.keep, ]
    tpm.out = scale(t(log2(tpm + 2))) 
    colnames(tpm.out) = genes[indices.keep]
    write.table(tpm.out, paste0(dir, tissue, '.log2.ztrans.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
    return()
}

##------------- MAIN

dir = Sys.getenv('TEMPDIR')
peer.dir = Sys.getenv('peerdir')
pc.file = Sys.getenv('GTEX_PCs')
subject.file = Sys.getenv('GTEX_SUBJECTSv8')


## Import file mapping sample identifiers to tissues -- Restricted to samples that pass RNA-seq QC (marked as RNASEQ in column 28 of the sample file)
map.file = Sys.getenv('SAMPLE_TISSUES')


## Read in list of tissues and keep those with more than 50 samples
tissue.list = read.table(map.file, header = F, stringsAsFactors = F)[,2]
tissues = names(table(tissue.list)[table(tissue.list) > 50]) 
tissues = tissues[tissues != 'Cells_Leukemia_cell_line_CML'] # exclude K-652 samples

# Read in data for covariates and build the covariate matrix
pcs = read.table(pc.file, header = T, stringsAsFactors = F)
## shorten names: e.g., only keep GTEX-1117F of GTEX-1117F-0003-SM-6WBT7
pcs$SUBJID = apply(str_split_fixed(pcs$IID, '-', 5)[, c(1,2)], 1, paste, collapse = '-')
pcs = pcs[, -c(1,2)]
sex = read.csv(subject.file, header = T, stringsAsFactors = F, sep = '\t')[, c(1,3)]
covariates = merge(pcs, sex, by = 'SUBJID')

# Output the covariates matrix
write.table(covariates, paste0(peer.dir, '/covariates.txt'),
            col.names = T, row.names = F, quote = F, sep = '\t')

sapply(tissues, ztrans.tissue, dir = peer.dir, covs = covariates)
