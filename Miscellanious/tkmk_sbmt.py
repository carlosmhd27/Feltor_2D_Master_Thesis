from sys import path
path.insert(1, '/marconi/home/userexternal/crodrigu/Plasma/prod/2D_FELTOR_Analysis/')

from Analysis        import units
from json            import dump, load
from os              import mkdir, replace, remove, system
from os.path         import join, isdir, exists
from scipy.constants import e, m_e, m_p, m_n
from pandas          import read_csv
from numpy           import ceil
csv_file = "characteristics/Tokamaks.csv"
tokamaks = read_csv(csv_file, index_col = 0)
orig_dir = "/marconi_scratch/userexternal/crodrigu/hdiff_outputs"
sbmt_dir = "submition/"
inpt_dir = "inputs/"

coment  = ''
comment = True

old_models = ["HW_mod", "HW_ord", "IC"]
Mix_models = ["IC_HW_mod", "IC_HW_ord", 'HW_mod_IC', 'HW_ord_IC']
#simply = {"HW_mod": "HWM_lynu", "HW_ord": "HWO_lynu", "IC": "IC_lynu"}
#models = ["IC_HW_mod", "IC_HW_ord", 'HW_mod_IC', 'HW_ord_IC']
models = ['HW_mod_IC', 'HW_ord_IC']
#models = ["HW_mod", "HW_ord", "IC"]
simply = {"IC_HW_mod": "ICHWM", "IC_HW_ord": "ICHWO",
          "HW_mod_IC": "HWMIC", "HW_ord_IC": "HWOIC",
          "HW_mod": "HWM", "HW_ord": "HWO", "IC": "IC"}

in_model = {"IC_HW_mod": "IC_HW", "IC_HW_ord": "IC_HW",
          "HW_mod_IC": "HW_IC", "HW_ord_IC": "HW_IC",
          "HW_mod": "HW", "HW_ord": "HW", "IC": "IC"}

in_val =  {}

in_val['lx']         = 200     #64  #Now 128
in_val['Nx']         = 500     #128 #Now 256
in_val['ly']         = 128     #64
in_val['Ny']         = 256     #128
in_val["nu_perp"]    = 0.001   #0.005 and 0.0001
in_val['maxout']     = 2500
in_val['dt']         = 0.0025
in_val['itstp']      = 200
in_val['n_out']      = 3
in_val['Nx_out']     = 32
in_val['Ny_out']     = 32
in_val['tanh_width'] = 0.01
in_val['x_b']        = 50      #in_val['lx'] / 2.

T_e = 25

for tkmk in tokamaks.index:
    if tkmk != 'EAST':
        continue
    new_dir   = f"tanh_Decristoforo_{tkmk}_T_e_{T_e}_nu_{in_val['nu_perp']}"
    directory = join(orig_dir, new_dir)
	
    if 'MPI' in new_dir or 'mpi' in new_dir:
        mk_t_mpi = True
    else:
        mk_t_mpi, np = False, 1

    if not isdir(directory):
        mkdir(directory)
    else:
        print('The directory already exists')

    tk_var = tokamaks.loc[tkmk]
    unit   = units(T_e = T_e, **tk_var)

    ## Parameters to submit
    input_prog = f"inputs/input_data_\"$model\"_{new_dir}.json"
    program = 'convection_mpi' if mk_t_mpi else 'convection_hpc'
    srun = 'srun --mpi=pmi2' if mk_t_mpi else 'srun'
    
    for model in models:
        if mk_t_mpi and 'MPI' not in simply[model]::
            simply[model] += 'MPI'		

        if model == 'IC' and 'Mix' in new_dir:
            continue

        File = join(inpt_dir, f"input_data_{model}.json")
        nw_npt = join(inpt_dir, f"input_data_{model}_{new_dir}.json")
        with open(File, 'r') as json_file:
            data = load(json_file)
        json_file.close()

        for var in in_val.keys():
            data[var] = in_val[var]

        data['model']     = in_model[model]
        data['modified']  = 1 if 'mod' in model else 0
                            #unit.kappa and unit.alpha
        data['curvature'] = unit.kappa if "IC" in model or 'Mix' in new_dir else 0
        data['adiabatic'] = 6e-4
        print(f"alpha = {data['adiabatic']} and kappa = {data['curvature']}")
        
        with open(nw_npt, 'w') as outfile:
            dump(data, outfile, indent=1)

        bash_name = f"submit_skl_{model}.sh"
        new_bash  = "new_" + bash_name
        bash_name = join(sbmt_dir, bash_name)
        new_bah   = join(sbmt_dir, new_bash)
        if exists(new_bash):
            remove(new_bash)

        if coment != 'A' and coment != '' :
            coment = input("Should I comment the GIF part?(yes/[A]ll/no)\n")
            if coment in ['y', 'ye', 'yes', 'Y', 'YE', 'YES', '', 'A', 'All', 'all', 'a']:
                comment = True
            else:
                comment = False

        with open(bash_name) as bash_file, open(new_bash, "a") as output:
            for line in bash_file:
                if "model=" in line:
                    line = "model=" + model + "\n"
                elif "Folders=" in line:
                    line = "Folders=\""+directory+"/\"\n"
                elif "cp" in line:
                    line = f"cp {nw_npt} $Folders\n"
                elif "#SBATCH -o" in line:
                    line = f"#SBATCH -o \"info/{model}_{new_dir}.info\"\n"
                elif "#SBATCH -J" in line:
                    line = f"#SBATCH -J {tkmk}_{simply[model]}\n"
                elif '#SBATCH -N' in line and mk_t_mpi:
                    line = f'#SBATCH -N {int(ceil(np / 2))} -n {np} --ntasks-per-socket=1 -c 24 --hint=memory_bound\n'
                elif '#SBATCH -N' in line and not mk_t_mpi:
                    line = f'#SBATCH -N 1 -n 1 -c 48 --hint=memory_bound\n'
                elif "run" in line and mk_t_mpi:
                    line = f'{srun} {program} "{input_prog}" $output {np0} {np1}\n'
                elif "run" in line and not mk_t_mpi:
                    line = f'{srun} {program} "{input_prog}" $output \n'
                elif "GIF" in line or 'source' in line or 'ffmpeg' in line or 'wait' in line:
                    if comment and '#' not in line:
                        line = '# ' + line
                    if not comment and '# ' in line:
                        line = line[3:]

                output.write(line)
        replace(new_bash, bash_name)

    GIF_file = "submit_GIF_mpi.sh"
    new_GIF  = f"submit_GIF_mpi_{new_dir}.sh"
    GIF_file = join(sbmt_dir, GIF_file)
    new_GIF  = join(sbmt_dir, new_GIF)
    if exists(new_GIF):
        remove(new_GIF)

    n = len(models) if 'Mix' not in new_dir else 2
    with open(GIF_file) as bash_file, open(new_GIF, "a") as output:
        for line in bash_file:
            if "direc=" in line:
                line = f"direc=\"{new_dir}/\"\n"
            elif "Directory=" in line:
                line = f"Directory=\"{orig_dir}/\"\n"
            elif "#SBATCH -o" in line:
                line = f"#SBATCH -o \"info/GIF_{new_dir}.info\"\n"
            elif '#SBATCH -N' in line:
                line = f'#SBATCH -N {int(ceil(n / 2))} -n {n} -c 24 --hint=memory_bound\n'
            elif 'run -n' in line and 'main' in line:
                line = f'\tsrun --mpi=pmi2 python \"~/Plasma/prod/Python_GIF/GIF_main_mpi.py\" $Folders $direc\n'
            elif 'run -n' in line and 'contour' in line:
                line = f'\tsrun --mpi=pmi2 python \"~/Plasma/prod/Python_GIF/GIF_contour_mpi.py\" $Folders $direc\n'

            output.write(line)
    # replace(new_GIF, GIF_file)

    submit_file = "submit.sh"
    new_submit  = "new_" + submit_file
    if exists(new_submit):
        remove(new_submit)

    with open(submit_file) as bash_file, open(new_submit, "a") as output:
        for line in bash_file:
            if "for model" in line:
                line = "for model in"
                for model in models:
                    line += " \"" +  model + "\""
                line += ";do\n"
            output.write(line)
    replace(new_submit, submit_file)

    submit = input("Should I submit the job?([Y]/n)\n")
    if submit in ['y', 'ye', 'yes', 'Y', 'YE', 'YES', '']:
    	system(". ./submit.sh")
    	print(f'Reactor {tkmk} submited')
    else:
        print('Done without submiting')
