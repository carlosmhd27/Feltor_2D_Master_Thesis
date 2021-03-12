from sys import path
path.insert(1, '/marconi/home/userexternal/crodrigu/Plasma/prod/2D_FELTOR_Analysis/')

from Analysis          import Analyse
from os.path           import join, exists, isdir
from os                import chdir
from numpy             import average, shape, amax, absolute
from time              import time
# from scipy.io          import savemat, loadmat
from matplotlib.pyplot import show, subplots, savefig, plot, clf, subplot, suptitle

home_dir     = '/marconi_scratch/userexternal/crodrigu/'
# chdir(join(home_dir, 'Statistics/'))
output_dir_1 = 'hdiff_outputs/DIR_128_ly_0.0001_nu_x2_mpi/'
output_dir_2 = 'hdiff_outputs/DIR_128_ly_0.0001_nu_x2_new/'
Dir_1        = join(home_dir, output_dir_1)
Dir_2        = join(home_dir, output_dir_2)

if 'Mix' in output_dir_1:
    models = ['HW_ord', 'HW_mod']
if 'tanh' in output_dir_1:
    models = ['IC_HW_mod', 'IC_HW_ord', 'HW_mod_IC', 'HW_ord_IC']
else:
    models = ['IC', 'HW_mod', 'HW_ord']

nt_md = {'IC':1086, 'HW_mod':1247, 'HW_ord':330}

for model in models:

    File_1 = join(Dir_1, f'output_{model}.nc')
    File_2 = join(Dir_2, f'output_{model}.nc')
    # tic    = time()

    Analitics_1 = Analyse(File_1)
    nt1, nx1, ny1 = shape(Analitics_1.ions)
    print(average(Analitics_1.Mass), Analitics_1.Mass[0])
    Analitics_2 = Analyse(File_2)
    nt2, nx2, ny2 = shape(Analitics_2.ions)

    assert (nx2 == nx1 and ny2 == ny1), 'They do not even have the same spacial shape, bro!'

    print(f'The shape of the output {output_dir_1} is {(nt1, nx1, ny1)} and the last mesaured time is {Analitics_1.time[-1]}')
    print(f'The shape of the output {output_dir_2} is {(nt2, nx2, ny2)} and the last mesaured time is {Analitics_2.time[-1]}')

    nt1 = nt_md[model]

    if nt1 == nt2:
        dens_diff = average((Analitics_1.ions        - Analitics_2.ions)        ** 2)
        last_diff = average((Analitics_1.ions[-1]    - Analitics_2.ions[-1])    ** 2)
    elif nt1 < nt2:
        dens_diff = average((Analitics_1.ions[:nt1]  - Analitics_2.ions[:nt1])  ** 2)
        last_diff = average((Analitics_1.ions[nt1-1] - Analitics_2.ions[nt1-1]) ** 2)
    elif nt1 > nt2:
        dens_diff = average((Analitics_1.ions[:nt2]  - Analitics_2.ions[:nt2])  ** 2)
        last_diff = average((Analitics_1.ions[nt2-1] - Analitics_2.ions[nt2-1]) ** 2)

    print(f'The average error in density is {dens_diff}')
    print(f'The error at the position {nt1} point is {last_diff}')
    print(f'out of {max(amax(absolute(Analitics_1.ions)), amax(absolute(Analitics_2.ions)))}')


    plot(Analitics_1.time, Analitics_1.Mass)
    show()
