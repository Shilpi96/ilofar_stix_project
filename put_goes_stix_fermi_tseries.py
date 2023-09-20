from sunpy.net import Fido, attrs as a
from stixpy.science import ScienceData
from itertools import product
import sunpy.timeseries
from sunpy.timeseries import TimeSeries
from stixpy import timeseries
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LogNorm
import datetime
import numpy as np

from astropy.time import Time

import radiospectra
import radiospectra.net
from astropy.io import fits
import matplotlib.style
import matplotlib as mpl
plt.rcParams.update(plt.rcParamsDefault)
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
import pdb
from sunpy.time import TimeRange
from matplotlib import dates
from matplotlib.ticker import AutoMinorLocator
import pylab

# -

def summary_plot(start, end, start_time, end_time, axes):
    
    '''
    Returns a figure containing data from GOES XRS and FERMI instruments using a start and end date.
    ''' 
    
    # Obtaining data...
    
    query = Fido.search(a.Time(start,end), a.Instrument.goes)
    query1 = Fido.search(a.Time(start,end),a.Instrument.gbm, a.Resolution("ctime"), a.Detector("n3"))
    print(query)
    
    
    # GOES
    
    goes_data = Fido.fetch(query['xrs'])
    goes = TimeSeries(goes_data[0])
    tr = TimeRange(start_time, end_time) #change this to what time range you want to plot
    goes_tr = goes.truncate(tr) #### need to truncate if goes have two timeseries of GOES-13 and GOES-15
    
    
    # FERMI
    
    fermi_data = Fido.fetch(query1['gbm'][0])
    #pdb.set_trace()
    fermi = TimeSeries(fermi_data[0]) 
    fermi = fermi.remove_column('100-300 keV')
    fermi = fermi.remove_column('300-800 keV')
    fermi = fermi.remove_column('800-2000 keV')
    fermi_tr = fermi.truncate(tr)
    
    
    axes[0].text(0.5,0.95, 'GOES',fontsize =20)
    
        
    goes_tr.plot(axes=axes[0], markersize=0.1)
    
    fermi_tr.plot(axes=axes[2], markersize=0.1)
    
def stix_plot(file1, file2, axes ):
    
    '''
    Returns a figure containing data from GOES XRS, STIX, STEREO WAVES (HFR and LFR),
    WIND WAVES (R1 and R2) and I-LOFAR instruments using a start and end date.
    ''' 
    
    # Get the header and data
    hdu1 = fits.open(file1)
    hdu2 = fits.open(file2)
    stix_data1 = ScienceData.from_fits(file1)
    stix_data2 = ScienceData.from_fits(file2)
    
    ####### Easy way to plot the stix timeseries 
    
    #stix_data2.plot_timeseries(energy_indices=[[1,7],[7,12],[12,17],[17,23],[23,28]])
    
    ####### To concatenate both science products
    
    ## Sum over energies
    counts1, errors1, times1, timedeltas1, energies1 =  stix_data1.get_data(detector_indices = [[0, 31]],pixel_indices = [[0, 11]],energy_indices=[[1,6],[7,11],[12,16],[17,22],[23,27]])
    
    counts2, errors2, times2, timedeltas2, energies2 =  stix_data2.get_data(detector_indices = [[0, 31]],pixel_indices = [[0, 11]],energy_indices=[[1,6],[7,11],[12,16],[17,22],[23,27]])
    
    ## now concat the counts array
    counts = np.concatenate((counts1.value, counts2.value), axis = 0)
    
    ## Get the concat time array
    times1 = times1.to_datetime()+ datetime.timedelta(seconds = hdu1[0].header['EAR_TDEL']) ###### adding the time delay. Remember the value is for the average time of the entire science file and is from the center of the sun
    
    times2 = times2.to_datetime() + datetime.timedelta(seconds = hdu2[0].header['EAR_TDEL'])
    
    times = np.concatenate((times1, times2))
    
    ####### Plot the timeseries
    #fig, axes = plt.subplots()
    nt, nd, npix, ne = counts1.shape
    labels = ['4-10 keV', '10-15 keV', '15-25 keV', '25-50 keV', '50-84 keV']
    color = 'Reds'
    cm = pylab.get_cmap(color)
    collist = cm(np.linspace(0, 255, 6).astype(int))
    
    for did, pid, eid in product(range(nd), range(npix), range(ne)):
    	
    	lines = axes[1].plot(times, counts[:, did, pid, eid], label = labels[eid],color = collist[eid+1])
    
    axes[1].set_yscale('log')
    axes[1].xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    axes[1].set_ylabel('Counts')
    #axes.xaxis.set_major_formatter(DateFormatter("%H:%M"))
    axes[1].xaxis.set_minor_locator(AutoMinorLocator())
    axes[1].tick_params(axis='x', rotation=90)
    fig.autofmt_xdate()
    fig.tight_layout()
    #plt.show()   
# +
start_day_sample = "2022/11/11 00:00:00"
end_day_sample = "2022/11/11 22:00:00"

start_time = "2022/11/11 11:32:00"
end_time = "2022/11/11 11:44:00"

# Creating the summary plot
fig, axes = plt.subplots(3, 1, sharex=True, figsize = (12,18))
plt.subplots_adjust(wspace=0, hspace=0.05)
summary_plot(start_day_sample, end_day_sample, start_time, end_time, axes)

##### plot stix data
path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/stix_data/'
file1 = path+'solo_L1_stix-sci-xray-cpd_20221111T112917-20221111T113422_V01_2211111572-50228.fits'
file2 = path+'solo_L1_stix-sci-xray-cpd_20221111T113420-20221111T114142_V01_2211112966-50227.fits'

stix_plot(file1,file2,axes)
[axes[i].legend(loc = 'upper right') for i in range(3)]

plt.show()
