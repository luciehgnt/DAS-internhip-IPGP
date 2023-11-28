# This code completes the cross-correlation between the ref channel
# and all the other channels between start_channel and _end_channel,
# and stores the results (correlograms) into folders.


import os
import glob
import numpy as np
import xarray as xr
from obspy.signal.filter import bandpass
import multiprocessing as mp

from tools_pcc import *
from stack_functions_new import *

path1 = "/gpfs/sratch/hugonnet/ncfiles/"
path2 = "/gpfs/sratch/hugonnet/ccfiles/" # storage of correlograms for cross-correlation
path3 = "/gpfs/sratch/hugonnet/ccfiles/" # storage of correlograms for auto-correlation

files = glob.glob(path1 + '*.nc')
size_files = len(files)
files.sort()


#################################################

corr_type = 'phase-cross'   # can be 'phase-cross' (cross-correlation), 'phase-auto' (auto-correlation), 'normal', or 'onebit'

# define variable for correlation
window = 300 			# duration of time window in second
lag = window/2 			# shift of starting time on the next iteration
time = 0 			# time initialization

# filter boundary (Hz)
lowcut = 1
highcut = 10

start_position = 'mid_all' 

# define which part of fiber is used - not used when start_position is not 'mid_all'
start_channel = 56
end_channel = 430

ref_channel = 250		# reference point for correlation. currently used: 100, 200, 300, 400, 500, ... (depends on the part of the fiber that is being used)
nb_channel = 150		# number of channel used - not used when start_position is not 'mid_all'


#################################################


def task_corr(file) :
	
	if corr_type == 'phase-cross' :
		folder_destination = file.replace(path1, path2).replace(".nc", "/")

	elif corr_type == 'phase-auto' :
		folder_destination = file.replace(path1, path3).replace(".nc", "/")
		
		
	if not os.path.exists(folder_destination):
		os.mkdir(folder_destination)


	ds = xr.open_dataset(file)
	channels_coord = ds['channels']

	iDas_nc = xr.DataArray( 
		data=ds['strain_rates'],
		dims=("times", "channels"), 
		coords={"times": ds['strain_rates'].coords["times"], "channels": channels_coord}
		)
    
	dt = ( iDas_nc['times'][1].item() - iDas_nc['times'][0].item() ) * 1e-9 # s
	fs = 1/dt # Hz

	f = open('error.txt','w') # files with few channels are not used and are stored into the file 'f'

	if iDas_nc.shape[1] < 463: 
		print(f'This file only has {iDas_nc.shape[1]} channels. Skipping this file')
		f.write(f'{iDas_nc.shape[1]} channels\n')
		# skipping file if the channel is not complefind . -type f -path "*/C*/" -exec mv {} DossierDestination \;te

	else :

		print(f'Filtering using bandpass filter between {lowcut} and {highcut} Hz ...')
		traces_filt = iDas_nc.copy()
		print('no gauge length change')
	
    	# for j in [56,3249]:
		for j in range (start_channel,end_channel+1):
			try:
				traces_filt[:,j] = bandpass(traces_filt[:,j],lowcut,highcut,fs, zerophase=True) # obspy butterworth filter (faster)
			except:
				print('filtrage pas abouti pour j = ', j)
				#None
		
		length = iDas_nc.shape[0]	# total length of data

		# cross correlate

		if corr_type == 'phase-cross' :
			print('\nStarting the phase cross-correlation  ...')
			pcc, t = cross_correlate(traces_filt,folder_destination,window,lag,length,time,start_channel,end_channel,nb_channel,
                	                      ref_channel,start_position,onebit=False,phase_cc=True,decimate=1,resample=1)
			
		elif corr_type == 'phase-auto' :
			print('\nStarting the phase auto-correlation  ...')
			pcc_auto, t = auto_correlate(traces_filt,folder_destination,window,lag,length,time,start_channel,end_channel,nb_channel,
            	                          ref_channel,start_position,onebit=False,phase_cc=True,decimate=1,resample=1)


	f.close()
	i += 1


if __name__ == '__main__' :

    n_processes = 32 	# number of CPUs used for the task

    with mp.Pool(processes = n_processes) as pool :
        pool.map(task_corr, files)


