from GIF_modules import Analyzed, animate, init
from datetime    import datetime
from os.path     import join, exists, isdir
from os          import remove
from sys         import argv
from time        import time
from numpy       import arange
from mpi4py      import MPI

from matplotlib.pyplot    import clf, show, subplots
from matplotlib.animation import writers, FuncAnimation
from matplotlib           import use

use("Agg")

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()
wrtr = 'pillow'
## General variables
fps = 5
i   = 1   ## to reduce the number of points we use

## We define the directory where the outputs are located and confirm it exists
# try:
#    dir_out = sys.argv[1]
# except:
#    dir_out = 'DIR_NEU_3'

#dir_name  = '/marconi_scratch/userexternal/crodrigu/Large_output/hdiff_outputs/'
dir_name = argv[1]
if rank == 0:
    print(dir_name)

assert isdir(dir_name), 'The directory does not exist.'

extra  = '_Marconi'
if len(argv) > 1:
    if  argv[2][-1] == '/':
        extra += '_' + argv[2][:-1]
    else:
        extra += '_' + argv[2]

if 'Mix' in dir_name:
    models = ['HW_ord', 'HW_mod']
elif 'tanh' in dir_name:
    models = ["IC_HW_mod", "IC_HW_ord", 'HW_mod_IC', 'HW_ord_IC']
elif 'Complete' in dir_name and size ==2:
    models = ['HW_mod_IC', 'HW_ord_IC']
elif 'complete' in dir_name and size==5:
    models = ['HW_mod_IC', 'HW_ord_IC', 'HW_ord', 'HW_mod', 'IC']
else:
    models = ['IC', 'HW_mod', 'HW_ord']

if size < len(models):
    print(f'Not all the videos will be created, size = {size} out of {len(models)}')

if rank < len(models):
    model = models[rank]
    print(f'Hello, I am processor {name} and I will take care of model: {model}.\n I hope I do it properly')
    GIF_name = join(dir_name, f'video_{model + extra}')
    GIF_name += '.mp4' if wrtr =='ffmpeg' else '.gif'
    #if exists(GIF_name):
        # overwrite = input("The GIF file already exists, do you wanna overwrite it?([Y]/n)\n")
    #   overwrite = 'Y'
    #   if overwrite not in ['y', 'ye', 'yes', 'Y', 'YE', 'YES', '']:
        # GIF_name = input('Please, give a new name')
    #       if GIF_name == '':
#           raise Exception('We need a name for the GIF')

    info_file = join(dir_name, f"info_{model + extra}.txt")
    if exists(info_file):
        remove(info_file)

    start     = time()
    File_name = join(dir_name, f'output_{model}.nc')
    existance = exists(File_name)

    ## Initiation of the figure
    fig, ax = subplots(2, 3, figsize=(24, 15))

    if not existance:
        with open(info_file, 'a') as information:
            information.write(f'{File_name} does not exist.')
        raise Exception('The output file does not exit.')

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
    init_    = lambda : init(model = model, Analytics = Analytics, ax = ax, fig = fig, extra = extra)
    animate_ = lambda pos, An, a: animate(pos, An, a, fig = fig)

    ## We define the format for writing the mp4 file
    Writer = writers[wrtr]
    writer = Writer(fps=fps, metadata=dict(artist='Me')) #, bitrate = 600

    ## We define the animation class and we save the animation, or show it 
    ani = FuncAnimation(fig, animate_, arange(0, len(Analytics.ions), i),
                        fargs = (Analytics, ax),
                        init_func=init_, interval = 100) ## arange(1, len(Analytics.ions))

    ani.save(GIF_name)
    # show()
    clf()

    stop   = time()
    needed = stop - start
    nd_hr  = needed // 3600
    nd_mn  = (needed / 3600 - nd_hr) * 60
    with open(info_file, 'a') as information:
        information.write('After {} h and {:.1f} min, I am done with model: {}'.format(nd_hr, nd_mn, model + extra) + 2 * "\n")
    print('After {} h and {:.1f} min, I am done with model: {}'.format(nd_hr, nd_mn, model + extra) + 2 * "\n")

    #except:
    #   print('Something went wrong')
