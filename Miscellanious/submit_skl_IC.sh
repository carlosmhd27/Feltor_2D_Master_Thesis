#!/bin/bash

#SBATCH -J IC
#SBATCH -N 1 -n 1 -c 48 --hint=memory_bound
#SBATCH --partition=skl_fua_prod 
#SBATCH --account=FUA35_FELTOR
#SBATCH --mem=182000
#SBATCH --time=24:00:00


###please specify a info file (for humans)
#SBATCH -o "info/IC_DIR_128_ly_64_lx_1e-07_nu_0.001_dt_x2.info"

NAME='Marconi Skylake'

echo "$NAME"

date
model=IC
Folders="/marconi_scratch/userexternal/crodrigu/hdiff_outputs/DIR_128_ly_64_lx_1e-07_nu_0.001_dt_x2/"
output=$Folders"output_"$model".nc"
cp inputs/input_data_IC_DIR_128_ly_64_lx_1e-07_nu_0.001_dt_x2.json $Folders

srun ./../Feltor_2D_Master_Thesis/FELTOR-model/convection_hpc       "inputs/input_data_"$model"_DIR_128_ly_64_lx_1e-07_nu_0.001_dt_x2.json" $output 

# wait
# source /marconi/home/userexternal/crodrigu/Plasma/analysispy/bin/activate
# module load ffmpeg/

# python "Python_GIF/GIF_main_"$model".py" $Folders
# python "Python_GIF/GIF_contour_"$model".py" $Folders
