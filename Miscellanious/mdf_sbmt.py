from json              import dump, load
from os                import mkdir, replace, remove, system
from os.path           import join, isdir, exists
from subprocess        import call
from scipy.constants   import e, m_e, m_p, m_n

## CONSTANTS
T_e       = 20 
B         = 1
kappa     = 0.015
m_i       = m_p + m_n
T_old     = 25; R_old = 1.65; B_old = 1.9
alpha_old = 0.000274

# omega_ci = e * B / (m_i)
## The charge of the electron is gone as the Temperature is in eV
rho_s    = (T_e * m_i) ** 0.5 / (B)
R_0      = rho_s / kappa

def calculate_alpha(R, T = 20, B = 1):
	alpha = alpha_old * (R_old **  2) * (T_old ** (-5/2)) * B_old
	alpha = alpha     * (R     ** -2) * (T     ** ( 5/2)) / B
	return alpha

## Directories
orig_dir = "/marconi_scratch/userexternal/crodrigu/hdiff_outputs"
sbmt_dir = "submition/"
inpt_dir = "inputs/"
in_val   = {}

coment  = ''
comment = True

smpl_models = ["HW_mod", "HW_ord", "IC"]
Mix_models  = ["IC_HW_mod", "IC_HW_ord", 'HW_mod_IC', 'HW_ord_IC']
#simply = {"HW_mod": "HWM_lynu", "HW_ord": "HWO_lynu", "IC": "IC_lynu"}
#models = ["IC_HW_mod", "IC_HW_ord", 'HW_mod_IC', 'HW_ord_IC']
simply = {"IC_HW_mod": "ICHWM", "IC_HW_ord": "ICHWO", 
		  "HW_mod_IC": "HWMIC", "HW_ord_IC": "HWOIC",
		  "IC": "IC", "HW_ord": "HWO", "HW_mod":"HWM"}


in_val['lx']         = 64       #64
in_val['Nx']         = 256      #128
in_val['ly']         = 128      #64
in_val['Ny']         = 512      #128
in_val['nu_perp']    = 0.0001   #0.005
in_val['maxout']     = 6250     #2500 for 0.0025 | 6250 for 0.001
in_val['dt']         = 0.001    #0.0025
in_val['itstp']      = 500      #200  for 0.0025 | 500  for 0.001
in_val['n_out']      = 3
in_val['Nx_out']     = 32
in_val['Ny_out']     = 32
in_val['tanh_width'] = 0.01
in_val["x_b"]        = in_val['lx'] / 2.

np, np0, np1 = 4, 2, 2

assert (np == np0 * np1), 'The number of processors is wrongly selected'

new_dir = f"DIR_{in_val['ly']}_ly_{in_val['nu_perp']}_nu_x2_mpi_{np}_1"
mk_t_mpi = True if ('MPI' in new_dir or 'mpi' in new_dir) else False

directory  = join(orig_dir, new_dir)

if not isdir(directory):
    mkdir(directory)
else:
    print('The directory already exists')

if 'Mix' in new_dir:
    models = ["HW_mod", "HW_ord"]
if 'tanh' in new_dir:
    models = ["IC_HW_mod", "IC_HW_ord", 'HW_mod_IC', 'HW_ord_IC']
else:
    models = ["HW_mod", "HW_ord", "IC"]

## Parameters to submit
input_prog = f"inputs/input_data_\"$model\"_{new_dir}.json"
program = 'convection_mpi' if mk_t_mpi else 'convection_hpc'
srun = 'srun --mpi=pmi2' if mk_t_mpi else 'srun'

for model in models:
    if mk_t_mpi:
        simply[model] += 'MPI'  
    if model == 'IC' and 'Mix' in new_dir:
        continue

    File = join(inpt_dir, f"input_data_{model}.json")
    nw_npt = join(inpt_dir, f"input_data_{model}_{new_dir}.json")
    with open(File, 'r') as json_file:
        data = load(json_file)
    
    if model in Mix_models:
    	data['model'] = "IC_HW" if "IC_HW" in model else "HW_IC"
    if model in smpl_models:
        data['model'] = 'IC' if 'IC' in model else 'HW'

    data['modified'] = 1 if 'mod' in model else 0

    for var in in_val.keys():
        data[var] = in_val[var]

    data['adiabatic'] = 0.05 # calculate_alpha(R_0)
    cond = "IC" in model or 'Mix' in new_dir
    data['curvature'] = kappa if cond else 0
    with open(nw_npt, 'w') as outfile:
        dump(data, outfile, indent=1)
       
    bash_name = f"submit_skl_{model}.sh"
    new_bash  = "new_" + bash_name
    bash_name = join(sbmt_dir, bash_name)
    new_bah   = join(sbmt_dir, new_bash)

    if exists(new_bash):
        remove(new_bash)
	
    if coment != 'Y' and coment != '' :
        coment = input("Should I comment the GIF part?(y/[Y]/n)\n")
        if coment in ['y', 'ye', 'yes', 'Y', 'YE', 'YES', '']:
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
            elif '#SBATCH -N' in line and mk_t_mpi:
                line = f'#SBATCH -N {int(ceil(np / 2))} -n {np} --ntasks-per-socket=1 -c 24 --hint=memory_bound\n'
            elif '#SBATCH -N' in line and not mk_t_mpi:
                line = f'#SBATCH -N 1 -n 1 -c 48 --hint=memory_bound\n'
            elif "#SBATCH -o" in line:
                line = f"#SBATCH -o \"info/{model}_{new_dir}.info\"\n"
            elif "#SBATCH -J" in line:
                line = f"#SBATCH -J {simply[model]}\n"
            elif "run" in line and mk_t_mpi:
                line = f'{srun} {program} "{input_prog}" $output {np0} {np1}\n'
            elif "run" in line and not mk_t_mpi:
                line = f'{srun} {program} "{input_prog}" $output \n'
            elif "GIF" in line or "source" in line or "ffmpeg" in line or "wait" in line:
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
with open(GIF_file) as fl, open(new_GIF, "a") as output:
    for line in fl:
        if "direc=" in line:
            line = f"direc=\"{new_dir}/\"\n"
        if "Directory=" in line:
            line = f"Directory=\"{orig_dir}/\"\n"
        if "#SBATCH -o" in line:
            line = f"#SBATCH -o \"info/GIF_{new_dir}.info\"\n"
        if '#SBATCH -N' in line:
            line = f'#SBATCH -N {int(ceil(n / 2))} -n {n} --ntasks-per-socket=1 -c 24 --hint=memory_bound\n'
        if 'run -n' in line and 'main' in line:
            line = f'\tsrun --mpi=pmi2 python \"~/Plasma/prod/Python_GIF/GIF_main_mpi.py\" $Folders $direc\n'
        if 'mpirun -n' in line and 'contour' in line:
            line = f'\tsrun --mpi=pmi2 python \"~/Plasma/prod/Python_GIF/GIF_contour_mpi.py\" $Folders $direc'

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

# Subprocess can not run source as it is not an executable command
submit = input("Should I submit the job?([Y]/n)\n")
if submit in ['y', 'ye', 'yes', 'Y', 'YE', 'YES', '']:
	system(". ./submit.sh")
print('Done')
