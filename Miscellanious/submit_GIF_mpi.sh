#!/bin/bash

#SBATCH -J GIFs
#SBATCH -N 2 -n 4 -c 24 --hint=memory_bound
#SBATCH --partition=skl_fua_prod
#SBATCH --account=FUA35_FELTOR
#SBATCH --mem=182000
#SBATCH --time=24:00:00


###please specify a info file (for humans)
#SBATCH -o "info/GIF_DIR_128_ly_1e-07_nu_x2_mpi_4.info"

date
Directory="/marconi_scratch/userexternal/crodrigu/Large_output/hdiff_outputs/"
direc="DIR_128_ly_1e-07_nu_x2_mpi_4/"
module load ffmpeg/
source /marconi/home/userexternal/crodrigu/Plasma/analysispy/bin/activate
    Folders=$Directory$direc
	echo $Folders
	srun --mpi=pmi2 python "Python_GIF/GIF_main_mpi.py" $Folders $direc
    srun --mpi=pmi2 python "Python_GIF/GIF_contour_mpi.py" $Folders $direc

# done
date
