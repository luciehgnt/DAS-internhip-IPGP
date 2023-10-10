import numpy as np
from scipy import signal
from tools_pcc import *
import os
from obspy import Trace
# from obspy import Stream, UTCDateTime
from obspy.core import Stats
# from obspy import read



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

def one_bit(trace):
	trace_onebit = np.where(trace < 0, -1, trace)
	trace_onebit = np.where(trace_onebit > 0, 1, trace_onebit)
	return trace_onebit

def make_date(traces):
	# date = UTCDateTime(traces['times'][0].values)
	# year = str(date.year)
	# month = str(date.month).zfill(2)
	# day = str(date.day).zfill(2)
	# hour = str(date.hour).zfill(2)
	# minute = str(date.minute).zfill(2)

	date = traces['times'][0].values
	# year = str(np.datetime64(date, 'Y'))
	# month = str(np.datetime64(date, 'M'))
	# day = str(np.datetime64(date, 'D'))
	# hour = str(np.datetime64(date, 'h'))
	minute = str(np.datetime64(date, 'm')) # returns (example) "2022-09-18T17:35"
	
	#fdate = f'{year}{month}{day}_{hour}{minute}'
	fdate = minute.replace("-","").replace("T","_").replace(":","")

	return fdate

def make_dir(traces, pcc, parent_dir, folder, stats, start_channel,i, time, window):
	dt = traces['times'][1].item() - traces['times'][0].item()
	dt = float(dt) * 1e-9
	try:
		os.mkdir(parent_dir + folder)
	except:
		None
		# print(f'Folder {parent_dir + folder} already exists. Continue writing ...')

	stats.npts = len(pcc)
	correlogram = Trace(data=pcc, header=stats)
	fdate = make_date(traces)
	fname = f'{start_channel}_{i}_{int(time*dt)}_{window}_{fdate}'
	correlogram.write(parent_dir + folder + fname + '.sac', format='sac')		

# define function for cross correlation: one-bit cc, normal cc, and phase cc
def cross_correlate(traces,parent_dir,window,lag,length,time,nb_channel,start_channel,start_position,
			onebit,phase_cc,decimate=1,resample=1):
	dt = traces['times'][1].item() - traces['times'][0].item()
	dt = float(dt) * 1e-9
	#dt = dt*resample

	t = np.arange(-window, window, dt)

	# make output array
	if start_position == 'mid':
		cc = np.zeros((2*int(window/dt),2*nb_channel+1))
	elif start_position == 'mid_all':
		cc = np.zeros((2*int(window/dt),430+1))
		# cc = np.zeros((2*int(window/dt),3249+1))
	else:
		cc = np.zeros((2*int(window/dt),nb_channel+1))

	# cross correlate all selected channel in time window
	if start_position == 'mid':
		ranges = np.arange(start_channel-nb_channel,start_channel+nb_channel+1,1)
	elif start_position == 'early':
		ranges = np.arange(start_channel, start_channel+nb_channel+1, 1)
	elif start_position == 'end':
		ranges = np.arange(start_channel, start_channel-nb_channel-1, -1)
	elif start_position == 'mid_all':
		ranges = np.arange(56, 430+1, 1)
		#ranges = np.arange(56, 3250, 1)

	# looping cross correlation over varying time window
	iteration = 1
	#pcc_all = np.zeros([nb_channel,cc.shape[0],cc.shape[1]])
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


			for i_channel in ranges:
			#for i in range(57,60): # use this when only interested in some channels
			
				if onebit == True:
					x1 = one_bit(traces[time:time+int(window/dt),start_channel])
					x2 = one_bit(traces[time:time+int(window/dt),i_channel])
				else:
					x1 = traces[time:time+int(window/dt),start_channel]
					x2 = traces[time:time+int(window/dt),i_channel]
				
				if phase_cc == True:

					try:
						_t , pcc = pcc2(x1, x2, dt, -window, window)
						cc[:,j] = cc[:,j] + pcc
						
						folder = f'Correlogram_{start_channel}_{i_channel}_pcc/'
						make_dir(traces, pcc, parent_dir, folder, stats, start_channel, i_channel, time, window)
					
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break

				else:
					try: # to account for cutted trace in last time window
						corr = signal.correlate(x1,x2)
						corr = np.append(corr,0) # pad with zero to similarize the array length
						cc[:,j] = cc[:,j] + corr  
						folder = f'Correlogram_{start_channel}_{i_channel}_normal/'
						make_dir(traces, corr, parent_dir, folder, stats, start_channel,i_channel, time, window)
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break

				j = j + 1
		time = int(time + lag/dt) # here, time is an index, not a time value
		iteration = iteration + 1
	return cc, t

# define function for auto correlation: one-bit cc, normal cc, and phase cc
def auto_correlate(traces,parent_dir,window,lag,length,time,nb_channel,start_channel,start_position,
			onebit,phase_cc,decimate=1,resample=1):
	
	dt = traces['times'][1].item() - traces['times'][0].item()
	dt = float(dt) * 1e-9
	#dt = dt*resample

	t = np.arange(-window, window, dt)

	# make output array
	if start_position == 'mid':
		cc = np.zeros((2*int(window/dt),2*nb_channel+1))
	elif start_position == 'mid_all':
		cc = np.zeros((2*int(window/dt),430+1))
		# cc = np.zeros((2*int(window/dt),3249+1))
	else:
		cc = np.zeros((2*int(window/dt),nb_channel+1))

	# cross correlate all selected channel in time window
	if start_position == 'mid':
		ranges = np.arange(start_channel-nb_channel,start_channel+nb_channel+1,1)
	elif start_position == 'early':
		ranges = np.arange(start_channel, start_channel+nb_channel+1, 1)
	elif start_position == 'end':
		ranges = np.arange(start_channel, start_channel-nb_channel-1, -1)
	elif start_position == 'mid_all':
		ranges = np.arange(56, 430+1, 1)
		#ranges = np.arange(56, 3250, 1)

	# looping auto correlation over varying time window
	iteration = 1
	#pcc_all = np.zeros([nb_channel,cc.shape[0],cc.shape[1]])
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

			for i_channel in ranges:
			#for i in range(57,60):
		  
				x1 = traces[time:time+int(window/dt),i_channel].values.ravel()
				# attention à la gestion des channels (ajouter 1 ? enlever 1 ? à réfléchir)
            
				if phase_cc == True:

					try:
						_t , pcc = apcc2(x1, dt, -window, window)
						cc[:,j] = cc[:,j] + pcc # linear stack
						
						folder = f'Correlogram_{start_channel}_{i_channel}_pcc/'
						make_dir(traces, pcc, parent_dir, folder, stats, start_channel,i_channel, time, window)
							
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break

				else:
					try: # to account for cutted trace in last time window
						corr = signal.correlate(x1,x1)
						corr = np.append(corr,0) # pad with zero to similarize the array length
						cc[:,j] = cc[:,j] + corr  
						folder = f'Correlogram_{start_channel}_{i_channel}_normal/'
						make_dir(traces, corr, parent_dir, folder, stats, start_channel,i_channel, time, window)
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break

				j = j + 1
		time = int(time + lag/dt) # here, time is an index, not a time value
		iteration = iteration + 1
	return cc, t
