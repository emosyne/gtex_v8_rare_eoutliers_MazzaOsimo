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

#! Notes:
#! Charging is determined by cpu number*walltime.
#! The --ntasks value refers to the number of tasks to be launched by SLURM only. This
#! usually equates to the number of MPI tasks launched. Reduce this from nodes*56 if
#! demanded by memory requirements, or if OMP_NUM_THREADS>1.
#! Each task is allocated 1 CPU by default, and each CPU is allocated 3420 MiB
#! of memory. If this is insufficient, also specify
#! --cpus-per-task and/or --mem (the latter specifies MiB per node).

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-ccl              # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:
module load miniconda/3

#! Full path to application executable: 
application=""

#! Run options for the application:
options=""

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 56:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$[${numnodes}*${mpi_tasks_per_node}]

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.


#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"


###############################################################
### You should not have to change anything below this line ####
###############################################################

cd $workdir
echo -e "Changed directory to `pwd`.\n"

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"

if [ "$SLURM_JOB_NODELIST" ]; then
        #! Create a machine file:
        export NODEFILE=`generate_pbs_nodefile`
        cat $NODEFILE | uniq > machine.file.$JOBID
        echo -e "\nNodes allocated:\n================"
        echo `cat machine.file.$JOBID | sed -e 's/\..*$//g'`
fi

echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD


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
mkdir -p preprocessing_v8/PEER_v8/

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

