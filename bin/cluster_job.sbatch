#!/bin/bash
#SBATCH -J nuxs_spl_mkr
#SBATCH -c 16 # Number of cores
#SBATCH -p shared
#SBATCH --mem 64000
#SBATCH -t 0-20:00:00 # Maximum execution time (D-HH:MM:SS)
#SBATCH -o /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/nkamp/LIV2/sources/NuXSSplMkr/bin/logging/output/nuxssplmkr_%A_%a.out # Standard output
#SBATCH -e /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/nkamp/LIV2/sources/NuXSSplMkr/bin/logging/error/nuxssplmkr_%A_%a.err # Standard error$

source /n/holylfs05/LABS/arguelles_delgado_lab/Lab/common_software/setup.sh
cd /n/holylfs05/LABS/arguelles_delgado_lab/Everyone/nkamp/LIV2/sources/NuXSSplMkr/bin
source calc_xsecs.sh
