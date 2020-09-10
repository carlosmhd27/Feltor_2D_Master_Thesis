import os
import matplotlib
matplotlib.use("Agg")

from GIF_modules import Analyzed, animate, init

import numpy                as np
import matplotlib.animation as animation
import matplotlib.pyplot    as plt

## General variables
fps = 10
model = 'IC'

## Initiation of the figure
fig, ax = plt.subplots(2, 3, figsize=(24, 12))

## We define the directory where the outputs are located and confirm it exists
dir_name  = 'path/to/outputs/'
directory = './' + dir_name
exists    = os.path.isdir("./" + dir_name)

if not exists:
    raise Exception('The directory does not exist.')
    
## We define the file of the output and confirm it exists
# for model in ['IC', 'HW_ord', 'HW_mod']:
File_name = os.path.join(dir_name, 'output_{}.nc'.format(model))
existance = os.path.exists(File_name)

if not existance:
    raise Exception('The File does not exist.')

print('+' * 50)
print(r'model: {}'.format(model))
print('+' * 50)

## We analysis the file and obtain the values we wanna measure
Analytics = Analyzed(File_name)

## Generate the init function for the GIF, it needs the model the fig and the parameters
init_ = lambda : init(model = model, Analytics = Analytics, ax = ax, fig = fig)

## We define the format for writing the mp4 file
Writer = animation.writers['ffmpeg']
writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=600)

## We define the animation class and we save the animation, or show it 
ani = animation.FuncAnimation(fig, animate, np.arange(1, len(Analytics.ions), 3), fargs = (Analytics, ax),
                              init_func=init_, interval = 100)

ani.save('video_{}.mp4'.format(model), writer = writer)
# plt.show()

print('Done')
