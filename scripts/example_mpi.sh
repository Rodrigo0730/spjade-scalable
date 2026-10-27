#!/bin/bash -l
#SBATCH --job-name=name
#SBATCH --account=project_2012081
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=hugemem
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=30000

# Load r-env
module load r-env

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2012081/TMP" >> ~/.Renviron

#config
SETTING="setting_3"
DS_VALUES="40,80,120"

# Run the R script
srun apptainer_wrapper exec Rscript --no-save --slave analysisMPI.R --setting=${SETTING} --ds=${DS_VALUES}
