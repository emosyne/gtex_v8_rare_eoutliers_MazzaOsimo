#!/bin/bash

# Go though GTEx data to figure out which individuals have data for which tissues.
# Create a file that has 2 columns:
# one column with the individual id and one column with the tissue name.
# Also generate files with the lists of unique tissue names and individual IDs.

set -o nounset -o errexit -o pipefail

dir=${RAREDIR}/preprocessing_v8
out=${dir}/gtex_2017-06-05_tissue_by_ind.txt 
tissues=${dir}/gtex_2017-06-05_tissues_all_normalized_samples.txt
inds=${dir}/gtex_2017-06-05_individuals_all_normalized_samples.txt

# put header
echo -e "Tissue\tId" > $out

for f in ${dir}/PEER_v8/*.peer.ztrans.txt
do
    fname=`basename $f`
    tissue=${fname%.peer.ztrans.txt}
    head -n1 $f | awk -v tissue=$tissue '{for(i=2; i<=NF; i++){print tissue"\t"$i}}' >> $out
done

# get the lists of individuals and tissues separately
cut -f1 $out | sort | uniq | grep -P -v '^Tissue' > $tissues
cut -f2 $out | sort | uniq | grep -P -v '^Id' > $inds
