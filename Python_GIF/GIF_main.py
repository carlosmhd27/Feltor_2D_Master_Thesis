from GIF_modules import Analyzed, animate, init
from datetime    import datetime
from os.path     import join, exists, isdir
from os          import remove
# from sys         import argv
from time        import time
from numpy       import arange
from mpi4py      import MPI

from matplotlib.pyplot    import clf, show, subplots
from matplotlib.animation import writers, FuncAnimation, PillowWriter
from matplotlib           import use

use("Agg")

## General variables
fps = 5
i   = 1   ## to reduce the number of points we use

## We define the directory where the outputs are located and confirm it exists
dir_name  = '../marconi_outputs/hdiff_outputs/Per_Per/'

assert isdir("./" + dir_name), 'The directory does not exist.'
    
models = ['IC', 'HW_mod', 'HW_ord'] ## , 'IC', 'HW_mod', 'HW_ord' 'HW_ord_dt_0.01', 'HW_ord_dt_0.001', 'HW_ord_dt_0.0001', 'HW_ord_dt_1e-05'
extra  = '_Neumann'

for model in models:
    
    GIF_name = join(dir_name, f'video_{model + extra}.gif')
    if exists(GIF_name):
        overwrite = input("The GIF file already exists, do you wanna overwrite it?([Y]/n)\n")
        if overwrite not in ['y', 'ye', 'yes', 'Y', 'YE', 'YES', '']:
            GIF_name = input('Please, give a new name')
            if GIF_name == '':
                continue
            

    info_file = join(dir_name, f"info_{model + extra}.txt")
    if exists(info_file):
        remove(info_file)

    start     = time()
    File_name = join(dir_name, f'output_{model}.nc')
    existance = exists(File_name)

    ## Initiation of the figure
    fig, ax = subplots(2, 3, figsize=(24, 12))

    if not existance:
        with open(info_file, 'a') as information:
            information.write(f'{File_name} does not exist.')
        continue
        
    ## We analysis the file and obtain the values we wanna measure
    Analytics = Analyzed(File_name)

    with open(info_file, 'a') as information:
        information.write('+' * 50 + "\n")
        information.write(f'model: {model + extra}, with {len(Analytics.ions) // i} out of {len(Analytics.ions)} points in the GIF' + "\n")
        information.write(f'The times goes from {Analytics.time[0]} to {Analytics.time[-1]}' + "\n")
        information.write('+' * 50 + "\n")
        information.write(f'The file started at {datetime.now()}'+ '\n')
        information.write('+' * 50 + "\n")


    ## Generate the init function for the GIF, it needs the model the fig and the parameters
    init_    = lambda : init(model = model + extra, Analytics = Analytics, ax = ax, fig = fig)
    animate_ = lambda pos, An, a: animate(pos, An, a, fig = fig)

    ## We define the format for writing the mp4 file
#         Writer = writers['ffmpeg']
    writer = PillowWriter(fps=fps, metadata=dict(artist='Me'), bitrate = 600)

    ## We define the animation class and we save the animation, or show it 
    ani = FuncAnimation(fig, animate_, np.arange(0, len(Analytics.ions), i),
                                  fargs = (Analytics, ax),
                                  init_func=init_, interval = 100) ## np.arange(1, len(Analytics.ions))

    ani.save(GIF_name)
    # show()
    clf()

    stop   = time()
    needed = stop - start
    nd_hr  = needed // 3600
    nd_mn  = needed % 60
    with open(info_file, 'a') as information:
        information.write('After {} h and {:.1f} min, I am done with model: {}'.format(nd_hr, nd_mn, model + extra) + 2 * "\n")
    
print('Done')
