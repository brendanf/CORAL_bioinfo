#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name OptimOTU
#SBATCH --account project_2005718
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --nodes 1
#SBATCH --mem 64G
#SBATCH --time 24:00:00
#SBATCH --mail-type ALL
#SBATCH --gres=nvme:10

export OMP_STACKSIZE=8096
if [ -v SLURM_CPUS_ON_NODE ] ; then
  export OMP_THREAD_LIMIT=$SLURM_CPUS_ON_NODE
fi
if [ -d "$LOCAL_SCRATCH" ] ; then
  export TMPDIR=$LOCAL_SCRATCH
  export SINGULARITY_BIND="${LOCAL_SCRATCH}:$(pwd)/userdir"
  echo "bind paths: $SINGULARITY_BIND"
fi
export PATH="/projappl/project_2005718/OptimOTU_v2/bin:$PATH"
if [[ $1 == "test" ]] ; then
if [[ $2 == "" ]] ; then
echo "Testing outdated targets..."
R --vanilla --quiet --no-echo -e 'targets::tar_outdated(callr_function=NULL)'
else
echo "Testing outdated targets leading to $2"
R --vanilla --quiet --no-echo -e "targets::tar_outdated(any_of(strsplit('$2', '[ ,;]')[[1]]), callr_function=NULL)"
fi
elif [[ $1 == "" ]] ; then
echo "Building plan on local machine"
R --vanilla --quiet --no-echo -e 'targets::tar_make(callr_function=NULL, reporter="timestamp")'
else
echo "Building target(s) '$1' on local machine"
R --vanilla --quiet --no-echo -e "targets::tar_make(any_of(strsplit('$1', '[ ,;]')[[1]]), callr_function=NULL, reporter='timestamp')"
fi
