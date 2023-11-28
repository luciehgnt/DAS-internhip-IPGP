# This code copies correlograms into specific folders that are then used
# for stacking in stack_cc.py

import os
import multiprocessing as mp

path1 = "/gpfs/scratch/hugonnet/ccfiles/" # storage of correlograms for cross-correlation

#######################################

# define which part of fiber is used - not used when start_position is not 'mid_all'
start_channel = 56
end_channel = 430

ref_channel = 250		# reference point for correlation. currently used: 100, 200, 300, 400, 500, ... (depends on the part of the fiber that is being used)

#######################################

#numbers = list(range(56, 431))
numbers = list(range(301, 431))
numbers

os.chdir(path1)


def task_cp(i_channel) :
    
    destination = path1 + f'all_correlograms_{ref_channel}_{i_channel}/'
    print(destination)
    
    if not os.path.exists(destination) : 
        os.mkdir(destination)

    os.system(f'find . -type f -wholename "*/Correlogram_{ref_channel}_{i_channel}_pcc/{ref_channel}_{i_channel}_*.sac" -exec cp {{}} {destination} \;')


if __name__ == '__main__' :

    n_processes = 32     # number of CPUs used for the task

    with mp.Pool(processes = n_processes) as pool :
        pool.map(task_cp, numbers)
