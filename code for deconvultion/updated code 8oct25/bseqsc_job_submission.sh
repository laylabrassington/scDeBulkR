#!/bin/bash
#SBATCH --mail-user=your_emaill@university.edu
#SBATCH --mail-type=ALL
#SBATCH -J bseqsc_job_date     	             
#SBATCH --time=80:00:00              
#SBATCH --mem=50GB
#SBATCH -o bseqsc_job_date_out_%j.out     

#reference script that has absolute path to the Rscript
#need to submit job in directory where example_r_script.R and example_submission.sh is 
/filepath/my_rscript.sh bseqsc_job_7oct25.R
