'''
########################################################################
   
    Description:
        Plots the dynamic spectra from ilofar, orfees, stix and also plot the stix timeseries on top of the spectra
    
     Author:
       Shilpi Bhunia, Dublin Institute for Advanced Studies
       shilpibhunia6@gmail.com
########################################################################
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import os, pdb
import pylab
from matplotlib import dates
import matplotlib.ticker as ticker
import argparse
from datetime import datetime
from matplotlib import dates
import matplotlib as mpl
from stixpy.product import Product
from astropy.io import fits
import astropy.units as u
from sunpy.net import Fido, attrs as a
from astropy.time import Time
from radiospectra.spectrogram import Spectrogram
from matplotlib.colors import LogNorm
from matplotlib import colors,dates
from stixpy.science import ScienceData
from matplotlib.ticker import AutoMinorLocator
from datetime import datetime, timedelta
from itertools import product
mpl.rcParams.update({'font.size': 16})

##### Function to plot stix spectra
def summary_plot(file1, file2, ax,Energy_indices = None  ):
    
    '''
    Returns a figure containing data from GOES XRS, STIX, STEREO WAVES (HFR and LFR),
    WIND WAVES (R1 and R2) and I-LOFAR instruments using a start and end date.
    ''' 
    
    # Get the header and data
    hdu1 = fits.open(file1)
    hdu2 = fits.open(file2)
    
    stix_data1 = Product(file1)
    stix_data2 = Product(file2)
    
    
    ####### Easy way to plot the stix timeseries 
    
    #stix_data2.plot_timeseries(energy_indices=[[1,7],[7,12],[12,17],[17,23],[23,28]])
    
    ####### To concatenate both science products
    
    ## Sum over energies
    counts1, errors1, times1, timedeltas1, energies1 =  stix_data1.get_data(detector_indices = [[0, 31]],pixel_indices = [[0, 11]],energy_indices = Energy_indices)
    
    counts2, errors2, times2, timedeltas2, energies2 =  stix_data2.get_data(detector_indices = [[0, 31]],pixel_indices = [[0, 11]],energy_indices = Energy_indices)
    
    ## now concat the counts array
    counts = np.concatenate((counts1.value, counts2.value), axis = 0)
    counts = counts *100
    ## Get the concat time array
    times1 = times1.to_datetime()+ timedelta(seconds = hdu1[0].header['EAR_TDEL']) ###### adding the time delay. Remember the value is for the average time of the entire science file and is from the center of the sun
    #pdb.set_trace()
    
    times2 = times2.to_datetime() + timedelta(seconds = hdu2[0].header['EAR_TDEL'])
    
    times = np.concatenate((times1, times2))
    
    ## Plot the spectrogram
    
    #fig,ax = plt.subplots()
    
    emin = energies1['e_low']
    emax = energies1['e_high'] 
    e_edges = np.hstack([energies1['e_low'], emax[-1]])
    yticks = emin.value.tolist()
    Yticks = [yticks[i*3] for i in range(10)]
    del Yticks[1]
    del Yticks[2]
    
    ylabels=[f"{n:.0f}-{x:.0f}" for n,x in zip(emin.value,emax.value)]
    #pdb.set_trace()
    emin = energies1['e_low']
    emax = energies1['e_high']
    ylabels=[f"{n:.0f}-{x:.0f}" for n,x in zip(emin.value,emax.value)]     ##### ticklabels
    Ylabels = [ylabels[i*3] for i in range(10)]
    del Ylabels[1]
    del Ylabels[2]
    
    #pdb.set_trace()
    
    ax.pcolormesh(times,e_edges[:-1],counts[:,0,0].T,norm = LogNorm())
    
    #ax.xaxis_date()
    ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    ax.xaxis.set_major_locator(dates.MinuteLocator(interval=2))
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    
    ax.set_yticks(Yticks)
    ax.set_yticklabels(Ylabels)
    ax.set_ylim((emin[0],emin[28]))
    ax.set_ylabel('Energy bins (keV)')
    return ax

def tseries_plot(counts, times, axes ):
    #pdb.set_trace()
    nt, nd, npix, ne = counts.shape
    #pdb.set_trace()
    labels = ['4-10 keV', '10-15 keV', '15-25 keV', '25-50 keV', '50-84 keV']
    color = 'Reds'
    cm = pylab.get_cmap(color)
    collist = cm(np.linspace(0, 255, 6).astype(int))
    counts = counts *100
    for did, pid, eid in product(range(nd), range(npix), range(1,4)):
    	
    	lines = axes.plot(times, counts[:, did, pid, eid], label = labels[eid], color = collist[eid+1])
    
    axes.set_yscale('log')
    axes.legend(loc = 'lower right',fontsize = 10)
    axes.set_ylabel(r"Counts s$^{-1}$ keV$^{-1}$")
    return lines, counts, times


##### Creating new freq axis with blank freq gaps of 87.6953125 - 110.546875 MHz and 180.28125 - 210.546875 MHz.

def freq_axis(freqs):
	gap1 = np.flipud(freqs[288]+(np.arange(59)*0.390625)) 
	gap2 = np.flipud(freqs[88]+(np.arange(57)*0.390625))
	ax_shape = 59+57-1
	new_freq = np.zeros(ax_shape+freqs.shape[0])
	#pdb.set_trace()
	new_freq[0:88] = freqs[0:88]
	new_freq[88:145]  = gap2[:57]
	new_freq[145:345] = freqs[88:288]
	new_freq[345:404] = gap1[:59]
	new_freq[404:] = freqs[289:]
	
	return new_freq
	

def plot_ilofarspectro(data,time,freqs,ax):
	
	new_freq = freq_axis(freqs)
	
	data = np.log10(data)
	data[np.where(np.isinf(data)==True)] = 0.0
	
	
	data2 = np.empty((new_freq.shape[0], data.shape[1]))    
	data2[:] = np.NaN
	data2[0:88] = data[0:88]
	data2[145:345] = data[88:288]
	data2[404:] = data[289:]
	
	times_mpl = [dates.date2num(t) for t in time]
	
	ax.pcolormesh(times_mpl,new_freq,data2,vmin=np.percentile(data, 1),vmax=np.percentile(data, 99.7), cmap = plt.get_cmap('viridis'))
	
	#ax = plt.gca()
	#ax.xaxis_date()
    
	#ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
	#ax.set_xlabel('Time (UT)')
	ax.set_ylabel('Frequency (MHz)')
	ax.set_ylim(new_freq[0],10.7421875)
	#plt.gca().invert_yaxis()
	#plt.show()

	return (data, freqs)

#### to replace the horizontal intensity values >12000 around 400-500 Mhz
def subt(data):
	data = data.T
	ind = np.where(data[0]>12000)
	#pdb.set_trace()
	for j in range(data.shape[0]):
		for i in range(ind[0].shape[0]):
			data[j][ind[0][i]] = data[j][670]
	#pdb.set_trace()
	return data.T


####### Plot ilofar spectra
path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/plot_high_res_ilofar/'
data = np.load(path+'high_res_ilofar.npy')
time = np.load(path+'time.npy',allow_pickle = True)
freqs = np.load(path+'freq.npy')
freq_axis(freqs)

# set up the plots 3 rows 1 column with shared common x-axis
fig, axes = plt.subplots(3, 1, sharex=True, figsize=(12, 12), dpi=100)
plt.subplots_adjust(hspace=0.05)

plot_ilofarspectro(data,time,freqs,axes[2])


####### Plot orfees spectra
###### Plot the spec
orfees_data = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/orfees_data.npy')
orfees_time = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/orfees_time.npy', allow_pickle = True)
orfees_freq = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/orfees_freqs.npy')

data = subt(orfees_data)
vmm = np.percentile(data, [10,99])
axes[1].pcolormesh(orfees_time, orfees_freq[0], data, norm = colors.LogNorm(vmin=vmm[0], vmax=vmm[1]))
axes[1].set_title(' ')
axes[1].set_ylabel(' ')
#axes[2].set_xlabel('Time (UT)')
axes[1].set_ylim(244.91,1000)

axes[2].invert_yaxis()
#axes[1].invert_yaxis()
#[axes[i].invert_yaxis() for i in range(3)]

#[axes[i].set_xlim(datetime.datetime(2022,11,11, 11,32,25), datetime.datetime(2022,11,11,13,44)) for i in range(3)]

axes[2].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
#axes[2].xaxis.set_major_locator(dates.MinuteLocator(interval=2))
axes[2].xaxis.set_minor_locator(AutoMinorLocator())

### Plot stix spectra
path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/stix_data/NEW_STX_DATA/'

file1 = path+'solo_L1_stix-sci-xray-cpd_20221111T112917-20221111T113422_V02_2211111572-50228.fits'
file2 = path+'solo_L1_stix-sci-xray-cpd_20221111T113420-20221111T114142_V02_2211112966-50227.fits'
summary_plot(file1,file2,axes[0])
axes[2].set_xlabel('Time (UT)')
#axes[3].invert_yaxis()
#[axes[i].set_xlim(datetime.datetime(2022,11,11, 11,32,25), datetime.datetime(2022,11,11,11,44)) for i in range(3)]

##### Plot stix timeseries
stx_counts = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/stix_tseries_data/counts.npy')
stx_time =  np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/stix_tseries_data/corr_times.npy', allow_pickle = True)

ax = axes[0].twinx()
lines, counts, times = tseries_plot( stx_counts,stx_time,ax)

#x_line = [datetime.datetime(2022,11,11, 11,33,15),datetime.datetime(2022,11,11, 11,34,7),datetime.datetime(2022,11,11, 11,35,2),datetime.datetime(2022,11,11, 11,38,49),datetime.datetime(2022,11,11, 11,39,16),datetime.datetime(2022,11,11, 11,40,39)]

x_line = [datetime(2022,11,11, 11,33,51, 470000),datetime(2022,11,11, 11,34,23, 470000),datetime(2022,11,11, 11,34,50),datetime(2022,11,11, 11,35,40), datetime(2022,11,11, 11,38,37,450000),datetime(2022,11,11, 11,40,45) ]

for j in range(len(axes)):
	[axes[j].axvline(x = x_line[i],linestyle = '--', c = 'white',alpha = 0.7) for i in range(len(x_line))]

[axes[i].set_xlim(datetime(2022,11,11, 11,32,30), datetime(2022,11,11,11,42,30)) for i in range(3)]

plt.show()
