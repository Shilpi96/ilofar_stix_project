## This will plot goes and stix timeseries and also polt doteed lines on top of stix to show the non thermal periods.

from sunpy.net import Fido, attrs as a
from stixpy.science import ScienceData
from itertools import product
from sunpy.timeseries import TimeSeries
from stixpy import timeseries
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LogNorm
import numpy as np
from stixpy.product import Product
from astropy.time import Time
import radiospectra
import radiospectra.net
from astropy.io import fits
import matplotlib.style
import matplotlib as mpl
plt.rcParams.update(plt.rcParamsDefault)
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
from sunpy.time import TimeRange
from matplotlib import dates
from matplotlib.ticker import AutoMinorLocator
import pylab,pdb
matplotlib.rcParams.update({'font.size': 16})
from datetime import datetime, timedelta
# -

def goes_plot(start, end, start_time, end_time, axes):
    
    '''
    Returns a figure containing data from GOES XRS using a start and end date.
    ''' 
    
    # Obtaining data...
    query = Fido.search(a.Time(start,end), a.Instrument.goes)
    print(query)
    # GOES
    
    goes_data = Fido.fetch(query['xrs'])
    goes = TimeSeries(goes_data[0])
    tr = TimeRange(start_time, end_time) #change this to what time range you want to plot
    goes_tr = goes.truncate(tr) #### need to truncate if goes have two timeseries of GOES-13 and GOES-15
    
    #axes.text(0.5,0.95, 'GOES',fontsize =20)
    goes_tr.plot(axes=axes, markersize=0.1)
    
    
def stix_plot(file1, file2, axes ):
    
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
    counts1, errors1, times1, timedeltas1, energies1 =  stix_data1.get_data(detector_indices = [[0, 31]],pixel_indices = [[0, 11]],energy_indices=[[1,6],[7,11],[12,16],[17,22],[23,27]])
    
    counts2, errors2, times2, timedeltas2, energies2 =  stix_data2.get_data(detector_indices = [[0, 31]],pixel_indices = [[0, 11]],energy_indices=[[1,6],[7,11],[12,16],[17,22],[23,27]])
    
    ## now concat the counts array
    counts = np.concatenate((counts1.value*100, counts2.value*100), axis = 0)
    
    ## Get the concat time array
    times1 = times1.to_datetime()+ timedelta(seconds = hdu1[0].header['EAR_TDEL']) ###### adding the time delay. Remember the value is for the average time of the entire science file and is from the center of the sun
    
    times2 = times2.to_datetime() + timedelta(seconds = hdu2[0].header['EAR_TDEL'])
    
    times = np.concatenate((times1, times2))
    
    ####### Plot the timeseries
    #fig, axes = plt.subplots()
    nt, nd, npix, ne = counts1.shape
    labels = ['4-10 keV', '10-15 keV', '15-25 keV', '25-50 keV', '50-84 keV']
    color = 'Reds'
    cm = pylab.get_cmap(color)
    collist = cm(np.linspace(0, 255, 6).astype(int))
    
    for did, pid, eid in product(range(nd), range(npix), range(ne)):
    	
    	lines = axes.plot(times, counts[:, did, pid, eid], label = labels[eid],color = collist[eid+1])
    
    axes.set_yscale('log')
    axes.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    axes.set_ylabel(r"Counts s$^{-1}$ keV$^{-1}$")
    axes.set_xlabel('Time (UT)')
    #axes.xaxis.set_major_formatter(DateFormatter("%H:%M"))
    axes.xaxis.set_major_locator(dates.MinuteLocator(interval=2))
    #axes[1].xaxis.set_minor_locator(dates.SecondLocator(interval=60))
    axes.xaxis.set_minor_locator(AutoMinorLocator())
    
    #fig.autofmt_xdate()
    #fig.tight_layout()
    #plt.show()   
# +
start_day_sample = "2022/11/11 00:00:00"
end_day_sample = "2022/11/11 22:00:00"

start_time = "2022/11/11 11:32:30"
end_time = "2022/11/11 11:44:45"

# Creating the summary plot
fig, axes = plt.subplots(2, 1, sharex=True, figsize = (12,10))
plt.subplots_adjust(wspace=0, hspace=0.065)
goes_plot(start_day_sample, end_day_sample, start_time, end_time, axes[0])

##### plot stix data
path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/stix_data/NEW_STX_DATA/'
file1 = path+'solo_L1_stix-sci-xray-cpd_20221111T112917-20221111T113422_V02_2211111572-50228.fits'
file2 = path+'solo_L1_stix-sci-xray-cpd_20221111T113420-20221111T114142_V02_2211112966-50227.fits'

stix_plot(file1,file2,axes[1])
[axes[i].legend(loc = 'upper right', fontsize = 'x-small') for i in range(2)]
[axes[i].set_xlim(datetime(2022,11,11,11,32,30), datetime(2022,11,11,11,44,45)) for i in range(2)]

x_line = [datetime(2022,11,11, 11,33,51, 470000),datetime(2022,11,11, 11,34,23, 470000),datetime(2022,11,11, 11,34,50),datetime(2022,11,11, 11,35,40), datetime(2022,11,11, 11,38,37,450000),datetime(2022,11,11, 11,40,45) ]

[axes[1].axvline(x = x_line[i],linestyle = '--', c = 'black',alpha = 0.3) for i in range(len(x_line))]


plt.show()
