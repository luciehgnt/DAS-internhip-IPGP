# This code stacks all correlograms given by "cross_corr.py".


import os
import multiprocessing as mp

path1 = "/gpfs/scratch/hugonnet/ccfiles/" # storage of correlograms for cross-correlation
path2 = "/gpfs/users/hugonnet/DAS/Sergi-ts-PWS/src/" # location of ts_pws file

#######################################

# define which part of fiber is used
start_channel = 56
end_channel = 430
ref_channel = 250		# reference point for correlation. currently used: 100, 200, 300, 400, 500, ... (depends on the part of the fiber that is being used)

#######################################

numbers = list(range(57, 182))
numbers

os.chdir(path1)

def task_stack(i_channel) :
    
    folder = path1 + f'all_correlograms_{ref_channel}_{i_channel}/'
    print(folder)
    os.system(f'cp {path2}ts_pws {folder}')
    os.chdir(folder)
    os.system('find . -type f -empty -print -delete')
    print(f'tf-pws for correlogram {ref_channel} and {i_channel}')
    os.system(f'ls --color=never *.sac > tf-pws-{ref_channel}_{i_channel}.in')
    os.system(f'./ts_pws tf-pws-{ref_channel}_{i_channel}.in osac="tspws_pcc_{ref_channel}_{i_channel}"')


if __name__ == '__main__' :

    n_processes = 32    # number of CPUs used for the task

    with mp.Pool(processes = n_processes) as pool :
        pool.map(task_stack, numbers)