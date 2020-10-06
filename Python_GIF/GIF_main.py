import os
import matplotlib
matplotlib.use("Agg")

from GIF_modules import Analyzed, animate, init
from datetime    import datetime

import numpy                as np
import matplotlib.animation as animation
import matplotlib.pyplot    as plt
import time

## General variables
fps = 5
i   = 1   ## to reduce the number of points we use

## We define the directory where the outputs are located and confirm it exists
dir_name  = '../marconi_outputs/hdiff_outputs/Per_Per/'
exists    = os.path.isdir("./" + dir_name)

## We define the file of the output and confirm it exists
if not exists:
    raise Exception('The directory does not exist.')
    
models = ['IC', 'HW_mod', 'HW_ord'] ## , 'IC', 'HW_mod', 'HW_ord' 'HW_ord_dt_0.01', 'HW_ord_dt_0.001', 'HW_ord_dt_0.0001', 'HW_ord_dt_1e-05'
extra  = '_Neumann'

for model in models:
    
    GIF_name = dir_name + f'video_{model + extra}.gif'
    if os.path.exists(GIF_name):
        overwrite = input("The GIF file already exists, do you wanna overwrite it?([Y]/n)\n")
        if overwrite not in ['y', 'ye', 'yes', 'Y', 'YE', 'YES', None]:
            GIF_name = input('Please, give a new name')
            if GIF_name == None:
                raise Exception('You must change the name of the GIF file')
            

    info_file = dir_name + f"info_{model + extra}.txt"
    if os.path.exists(info_file):
        os.remove(info_file)

    start = time.time()
    File_name = os.path.join(dir_name, f'output_{model}.nc')
    existance = os.path.exists(File_name)

    ## Initiation of the figure
    fig, ax = plt.subplots(2, 3, figsize=(24, 12))

    information = open(info_file, 'a')

    if not existance:
        information.write(f'{File_name} does not exist.')

    else:

        ## We analysis the file and obtain the values we wanna measure
        Analytics = Analyzed(File_name)


        information.write('+' * 50 + "\n")
        information.write(f'model: {model + extra}, with {len(Analytics.ions) // i} out of {len(Analytics.ions)} points in the GIF' + "\n")
        information.write(f'The times goes from {Analytics.time[0]} to {Analytics.time[-1]}' + "\n")
        information.write('+' * 50 + "\n")
        information.write(f'The file started at {datetime.now()}'+ '\n')
        information.write('+' * 50 + "\n")
        information.close()

        ## Generate the init function for the GIF, it needs the model the fig and the parameters
        init_ = lambda : init(model = model + extra, Analytics = Analytics, ax = ax, fig = fig)

        ## We define the format for writing the mp4 file
#         Writer = animation.writers['ffmpeg']
        writer = animation.PillowWriter(fps=fps, metadata=dict(artist='Me'), bitrate = 600)

        ## We define the animation class and we save the animation, or show it 
        ani = animation.FuncAnimation(fig, animate, np.arange(0, len(Analytics.ions), i),
                                      fargs = (Analytics, ax),
                                      init_func=init_, interval = 100)

        ani.save(GIF_name)
        # plt.show()
        plt.clf()

        stop = time.time()
        information = open(info_file, 'a')
        information.write('After {:.1f}h, I am done with model: {}'.format((stop - start) / 3600, model + extra) + 2 * "\n")
    information.close()
    
print('Done')
