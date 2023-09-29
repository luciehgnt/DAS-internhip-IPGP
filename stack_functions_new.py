import matplotlib.pylab as plt
import numpy as np
import xarray as xr
from scipy import signal
from tools_pcc import *
import os
from obspy import Stream,Trace,UTCDateTime
from obspy.core import Stats
from obspy import read


from xdas.io.febus import correct_gps_time
from xdas.io.febus import read as read_das

user = "mba"

if user == "ipgp" :
    path1 = "/home/invite/Documents/DAS-Lucie/fichiers-nc/"

elif user == "mba" :
    path1 = "/Volumes/DISKENOIR/fichiers-nc/"
    path2 = "/Volumes/DISKENOIR/Correlograms/Cross-correlation"


ds = xr.open_dataset(path1 + "SR_FEBUS_STROMBOLI_2022-09-18_17-35-44_UTC.nc")
array=ds['strain_rates']
channels_coord = ds['channels']
sr_array_with_channels = xr.DataArray(
    data=ds['strain_rates'],  # Remplacez par les données réelles de 'strain_rates'
    dims=("times", "channels"),  # Spécifiez les dimensions
    coords={"times": ds['strain_rates'].coords["times"], "channels": channels_coord},
)

iDas_nc = sr_array_with_channels
Das_nc = sr_array_with_channels.coords

dt = 0.02



def normalize(trace):	# normalize trace
		return trace / max(trace)

def cut_traces(trace, start, duration): # netcdf trace afrom tools import *s input
	dt = trace['times'][1] - trace['times'][0]
	dt = float(dt) * 1e-9
	traces_cut = trace[int(start/dt):int(start/dt)+int(duration/dt),:]
	return traces_cut

def filter_bandpass(trace, lowcut, highcut, fs, order=2):
	nyq = 0.5 * fs
	low = lowcut/nyq
	high = highcut/nyq
	#order = 2
	b,a = signal.butter(order, [low, high], 'bandpass', analog=False)
	print(b,a)
	trace_filtered = signal.filtfilt(b, a, trace, axis=0)
	return trace_filtered

def make_date(traces):
	date = UTCDateTime(str(traces[0].time.data))
	year = str(date.year)
	month = str(date.month).zfill(2)
	day = str(date.day).zfill(2)
	hour = str(date.hour).zfill(2)
	minute = str(date.minute).zfill(2)
	fdate = f'{year}{month}{day}_{hour}{minute}'
	return fdate

def make_dir(traces, pcc, parent_dir, folder, stats, start,i, time, window):
	try:
		os.mkdir(parent_dir + folder)
	except:
		None
		# print(f'Folder {parent_dir + folder} already exists. Continue writing ...')

	stats.npts = len(pcc)
	correlogram = Trace(data=pcc, header=stats)
	fdate = make_date(traces)
	fname = f'correlogram_{start}_{i}_{int(time*dt)}_{window}_{fdate}'
	correlogram.write(parent_dir + folder + fname + '.sac', format='sac')		

# define function for cross correlation: one-bit cc, normal cc, and phase cc
def cross_correlate(cc,traces,parent_dir,window,lag,length,time,channel,start,start_position,
			onebit,phase_cc,decimate=1,resample=1):
	dt = traces['times'][1].item() - traces['times'][0].item()
	dt = float(dt) * 1e-9
	#dt = dt*resample

	t = np.arange(-window, window, dt)

	# looping cross correlation over varying time window
	iteration = 1
	#pcc_all = np.zeros([channel,cc.shape[0],cc.shape[1]])
	stats = Stats()
	stats.delta = dt
	stats.sampling_rate = 1/stats.delta
	stats.network = 'Correlogram'
	stats.station = 'STRMBL'
	stats.location = 'Italy'

	while time < length:
		if iteration%decimate == 0:
			j = 0
			print(f'Start of the {iteration}th time window. Starttime = {time*dt} s') 

			# cross correlate all selected channel in time window
			if start_position == 'mid':
				ranges = np.arange(start-channel,start+channel+1,1)
			elif start_position == 'early':
				ranges = np.arange(start, start+channel+1, 1)
			elif start_position == 'end':
				ranges = np.arange(start, start-channel-1, -1)
			elif start_position == 'mid_all':
				ranges = np.arange(57, 430+1, 1)

			for i in ranges:
			#for i in range(57,60): # use this when only interested in some channels
				
                # if onebit == True:
				# 	x1 = one_bit(traces[time:time+int(window/dt),start])
				# 	x2 = one_bit(traces[time:time+int(window/dt),i])
				# else:
				# 	x1 = traces[time:time+int(window/dt),start]
				# 	x2 = traces[time:time+int(window/dt),i]
				
				if phase_cc == True:
					try:
						_t , pcc = pcc2(x1, x2, dt, -window, window)
						cc[:,j] = cc[:,j] + pcc
						
						folder = f'Correlogram_{start}_{i}_phase/'
						make_dir(traces, pcc, parent_dir, folder, stats, start, i, time, window)
							
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break
				else:
					try: # to account for cutted trace in last time window
						corr = signal.correlate(x1,x2)
						corr = np.append(corr,0) # pad with zero to similarize the array length
						cc[:,j] = cc[:,j] + corr  
						folder = f'Correlogram_{start}_{i}_normal/'
						make_dir(traces, corr, parent_dir, folder, stats, start,i, time, window)
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break

				j = j + 1
		time = int(time + lag/dt)
		iteration = iteration + 1
	return cc, t

# define function for auto correlation: one-bit cc, normal cc, and phase cc
def auto_correlate(cc,traces,parent_dir,window,lag,length,time,channel,start,start_position,
			onebit,phase_cc,decimate=1,resample=1):
	dt = traces['times'][1].item() - traces['times'][0].item()
	dt = float(dt) * 1e-9
	#dt = dt*resample

	t = np.arange(-window, window, dt)

	# looping auto correlation over varying time window
	iteration = 1
	#pcc_all = np.zeros([channel,cc.shape[0],cc.shape[1]])
	stats = Stats()
	stats.delta = dt
	stats.sampling_rate = 1/stats.delta
	stats.network = 'Correlogram'
	stats.station = 'STRMBL'
	stats.location = 'Italy'

	while time < length:
		if iteration%decimate == 0:
			j = 0
			print(f'Start of the {iteration}th time window. Starttime = {time*dt} s') 

			# cross correlate all selected channel in time window
			if start_position == 'mid':
				ranges = np.arange(start-channel,start+channel+1,1)
			elif start_position == 'early':
				ranges = np.arange(start, start+channel+1, 1)
			elif start_position == 'end':
				ranges = np.arange(start, start-channel-1, -1)
			elif start_position == 'mid_all':
				ranges = np.arange(57, 430+1, 1)

			for i in ranges:
			#for i in range(57,60):
				x1 = traces[time:time+int(window/dt),i]
				
				if phase_cc == True:
					try:
						_t , pcc = apcc2(x1, dt, -window, window)
						cc[:,j] = cc[:,j] + pcc # linear stack
						
						folder = f'Correlogram_{start}_{i}_phase/'
						make_dir(traces, pcc, parent_dir, folder, stats, start,i, time, window)
							
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break
				else:
					try: # to account for cutted trace in last time window
						corr = signal.correlate(x1,x2)
						corr = np.append(corr,0) # pad with zero to similarize the array length
						cc[:,j] = cc[:,j] + corr  
						folder = f'Correlogram_{start}_{i}_normal/'
						make_dir(traces, corr, parent_dir, folder, stats, start,i, time, window)
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break

				j = j + 1
		time = int(time + lag/dt)
		iteration = iteration + 1
	return cc, t
