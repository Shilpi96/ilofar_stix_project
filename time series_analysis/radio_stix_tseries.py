##### This code will uplot all the time series from all the instruments inclusing orfees, ilofar, wind and stix.

from scipy.ndimage import gaussian_filter1d
from sunpy.net import Fido, attrs as a
import sunpy.timeseries
from sunpy.timeseries import TimeSeries
from stixpy import timeseries
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LogNorm
import numpy as np
from datetime import datetime
from astropy.time import Time
import matplotlib.style
import matplotlib as mpl
plt.rcParams.update(plt.rcParamsDefault)
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
import pdb,pylab,datetime
from sunpy.time import TimeRange
from stixpy.science import ScienceData
from itertools import product
from matplotlib.dates import DateFormatter
from matplotlib.ticker import AutoMinorLocator
from matplotlib import dates
from astropy.io import fits
import astropy.units as u
matplotlib.rcParams.update({'font.size': 16})
# -

def stix_plot(file1, file2, axes ):
    
    # Get the header and data
    hdu1 = fits.open(file1)
    hdu2 = fits.open(file2)
    
    stix_data1 = ScienceData.from_fits(file1)
    stix_data2 = ScienceData.from_fits(file2)
    
    
    ####### To concatenate both science products
    
    ## Sum over energies
    counts1, errors1, times1, timedeltas1, energies1 =  stix_data1.get_data(detector_indices = [[0, 31]],pixel_indices = [[0, 11]],energy_indices=[[12,16],[17,22]])
    
    counts2, errors2, times2, timedeltas2, energies2 =  stix_data2.get_data(detector_indices = [[0, 31]],pixel_indices = [[0, 11]],energy_indices=[[12,16],[17,22]])
    
    ## now concat the counts array
    counts = np.concatenate((counts1.value, counts2.value), axis = 0)
    
   ## Get the concat time array
    times1 = times1.to_datetime()+ datetime.timedelta(seconds = hdu1[0].header['EAR_TDEL']) 
    times2 = times2.to_datetime() + datetime.timedelta(seconds = hdu2[0].header['EAR_TDEL'])
    
    times = np.concatenate((times1, times2))
    
    ####### Plot the timeseries
    nt, nd, npix, ne = counts1.shape
    labels = ['15-25 keV', '25-50 keV']
    color = 'Reds'
    cm = pylab.get_cmap(color)
    collist = cm(np.linspace(0, 255, 6).astype(int))
    
    for did, pid, eid in product(range(nd), range(npix), range(ne)):
    	
    	lines = axes.plot(times, counts[:, did, pid, eid], label = labels[eid],color = collist[eid*2+2])
    
    #axes[5].set_yscale('log')
    axes.legend(loc = 'upper right')
    axes.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    axes.set_ylabel('Counts')
    #axes.xaxis.set_major_formatter(DateFormatter("%H:%M"))
    axes.xaxis.set_minor_locator(AutoMinorLocator())
    

    return axes
    
def radio_timeseries(data,time,axes,freqs,color = ['red','green','blue'],freq = [80,150,260],j = 0,instr = 'stx'):
	
	idx = [np.abs(freqs - freq[i]).argmin() for i in range(len(freq))]
	
	times_mpl = [dates.date2num(t) for t in time]
	
	data1D = [data[idx[i], ::]/np.mean(data[idx[i], ::]) for i in range(len(freq))]
	#pdb.set_trace()
	labels = [str(freq[i])+' MHz' for i in range(len(freq))]
	
	if instr == 'wind':
		
		if j ==0:
			for i in range(len(freq)):
				axes[0].step(times_mpl, data1D[i],where ='mid', label = labels[i],color = color[i])
				
		else:
			for i in range(len(freq)):
				axes[0].step(times_mpl, data1D[i],where ='mid', label = labels[i],color = color[i])
				

	else:
		if j ==0:
			for i in range(len(freq)):
				axes[i].plot(times_mpl, data1D[i], label = labels[i],color = color[i])
		else:
			for i in range(len(freq)):
				axes[i+j].plot(times_mpl, data1D[i], label = labels[i],color = color[i])
	


path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/stix_data/'

file1 = path+'solo_L1_stix-sci-xray-cpd_20221111T112917-20221111T113422_V01_2211111572-50228.fits'
file2 = path+'solo_L1_stix-sci-xray-cpd_20221111T113420-20221111T114142_V01_2211112966-50227.fits'

# Creating the summary plot
fig, axes = plt.subplots(7, 1, sharex=True, figsize = (12,18))
plt.subplots_adjust(wspace=0, hspace=0.05)

### Creating colormap
color = 'Blues'
cm = pylab.get_cmap(color)
collist = cm(np.linspace(0, 255, 13).astype(int))

stix_plot(file1, file2, axes[6])
axes[6].set_yscale('log')

##### Plot wind/waves
waves_data = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/all_radio_data/waves_data.npy')
waves_time = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/all_radio_data/waves_time.npy',allow_pickle = True)
waves_freq = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/all_radio_data/waves_freq.npy')

radio_timeseries(waves_data.T, waves_time, axes,waves_freq, color = [collist[3]],freq = [9],j=0,instr = 'wind')
axes[0].set_ylim(0,100)

##### Plot I-LOFAR
path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/plot_high_res_ilofar/'
data = np.load(path+'realta_ilofar.npy')
time = np.load(path+'realta_time.npy',allow_pickle = True)
freqs = np.load(path+'freq.npy')

radio_timeseries(data,time,axes,freqs,color = [collist[4],collist[6],collist[8]], freq = [50,150,250],j=1)  ## j is the axis number

##### Plot orfees
orfees_data = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/orfees_data.npy')
orfees_time = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/orfees_time.npy', allow_pickle = True)
orfees_freq = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/orfees_freqs.npy')

radio_timeseries(orfees_data, orfees_time, axes, orfees_freq,color = [collist[10],collist[12]],freq = [400,700],j=4) 


[axes[i].legend() for i in range(6)]
[axes[i].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S')) for i in range(6)]
[axes[i].set_xlim(datetime.datetime(2022,11,11,11,32),datetime.datetime(2022,11,11,11,45)) for i in range(6)]
[axes[i].xaxis.set_major_locator(dates.MinuteLocator(interval=2)) for i in range(6)]
#[axes[i].set_yscale('log') for i in range(6)]
#[axes[i].set_ylim(10**-2,10**2) for i in range(6)]
axes[3].set_ylabel('Intensity')
#=pdb.set_trace()

plt.show()

