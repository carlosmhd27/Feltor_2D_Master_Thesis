from json              import dump, load
from os                import mkdir, replace, remove, system, path
from os.path           import join, isfile
from subprocess        import call
from numpy             import ceil

sbmt_dir = 'submition/'
# to_sbmt  = 'Folders_for_GIFs.txt'
to_sbmt  = 'GIFs_submitions.txt'
with open(to_sbmt, 'r') as file:
    for line in file:

        if line == '\n':
            break
        elif '\n' in line:
            line = line[:-1]
        if '/' not in line:
            GIF_file = join(sbmt_dir, line)
        else:
            GIF_file = line
        if isfile(GIF_file):
            system(f'sbatch {GIF_file}')
            print(f'{GIF_file} submited')
        else:
            print(f'{GIF_file} does not exist')

