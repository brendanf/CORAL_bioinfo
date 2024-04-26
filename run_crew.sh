#!/usr/bin/env bash

#filter and assign taxonomy to demultiplexed Illumina reads

#SBATCH --job-name OptimOTU
#SBATCH --account project_2005718
#SBATCH --partition small
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 16G
#SBATCH --time 3-00:00:00
#SBATCH --mail-type ALL

export PATH="/projappl/project_2005718/OptimOTU_v2/bin:$PATH"

if [[ $1 == "test" ]] ; then
if [[ $2 == "" ]] ; then
echo "testing outdated targets..."
echo "NOTE: crew is not used for testing, you could have used 'run_node.sh'"
R --vanilla --quiet --no-echo -e 'targets::tar_outdated(callr_function=NULL)'
else
echo "testing outdated targets leading to $2"
echo "NOTE: crew is not used for testing, you could have used 'run_node.sh'"
R --vanilla --quiet --no-echo -e "targets::tar_outdated($2, callr_function=NULL)"
fi
elif [[ $1 == "" ]] ; then
echo "Building plan using crew"
R --vanilla --quiet --no-echo -f run_crew.R
else
echo "Building target '$1' using crew"
OPTIMOTU_TARGET="$1" R --vanilla --quiet --no-echo -f run_crew.R
fi
