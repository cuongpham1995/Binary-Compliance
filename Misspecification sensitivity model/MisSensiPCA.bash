#!/bin/sh
#SBATCH --partition standard        # Partition [debug, preempt, standard, reserved]
#SBATCH --time 0-10:00:00           # wall time [D-HH:MM:SS]

##SBATCH --job-name myjob           # Job name
##SBATCH —output my_output%j        # File for standard out - here the %j will be replaced by the JobID
##SBATCH —error my_error%j          # File for standard error.  If not specified will go to same file as standard out.

##SBATCH --mem-per-cpu 5G
#SBATCH --nodes 1                   # Number of nodes [must be power of 2 within queue limitations]
##SBATCH --ntasks-per-node 8         # Number of tasks per node (There are 16 cores per node)
##SBATCH --cpus-per-task 1           # Number of cores per task (If not 1, then you should use openmp or a compiler that can auto-parallelize)

#SBATCH --mail-type ALL             # Send e-mail when... [BEGIN, END, FAIL, REQUEUE, or ALL]
#SBATCH --mail-user Cuong_Pham@URMC.Rochester.edu
#SBATCH --mem=5G

#SBATCH -a 1-500
#SBATCH -o log_%a.txt 
#SBATCH -e err_%a.txt
##SBATCH -o log_last.txt
##SBATCH -e err_last.txt 

echo This is job $SLURM_ARRAY_TASK_ID

module load r/3.6.1/b1
R --no-save < MisSensiPCA.R