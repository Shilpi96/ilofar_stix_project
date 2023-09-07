###### This script plots stix contours for all energy ranges overlaid on aia on a single plot. It also overlays the time stamps on the stix timeseries.

from sunpy.net import Fido, attrs as a
import sunpy.timeseries
from sunpy.timeseries import TimeSeries
from stixpy import timeseries
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LogNorm
import datetime
import numpy as np
from datetime import datetime
from astropy import units as u 
from astropy.time import Time
from astropy.io import fits
import matplotlib.style
import matplotlib as mpl
plt.rcParams.update(plt.rcParamsDefault)
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
matplotlib.rcParams.update({'font.size': 11.5})
import pdb
from sunpy.time import TimeRange
from stixpy.science import ScienceData
from itertools import product
from matplotlib.dates import DateFormatter
from datetime import datetime, timedelta
from matplotlib.ticker import AutoMinorLocator
from matplotlib import dates
import sunpy.visualization.colormaps as cm
import pylab, glob
import sunpy.map
from sunpy.coordinates.ephemeris import get_horizons_coord
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import matplotlib.gridspec as gridspec
sdoaia171 = matplotlib.colormaps['sdoaia171']
sdoaia1700 = matplotlib.colormaps['sdoaia1700']
# -
######## Function to Plot stix timeseries with marked time range
def summary_plot(counts, times, axes, tstart = datetime(2022,11,11,11,33,10), tend = datetime(2022,11,11,11,34,10), tdelta = 183 ):
    
    ####### Plot the timeseries
    nt, nd, npix, ne = counts.shape
    labels = ['4-10 keV', '10-15 keV', '15-25 keV', '25-50 keV', '50-84 keV']
    color = 'Reds'
    cm = pylab.get_cmap(color)
    collist = cm(np.linspace(0, 255, 6).astype(int))
    
    
    for did, pid, eid in product(range(nd), range(npix), range(ne)):
    	
    	lines = axes.plot(times, counts[:, did, pid, eid], label = labels[eid], color = collist[eid+1])
    
    axes.set_yscale('log')
    axes.legend(fontsize = 8)
    axes.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
    axes.xaxis.set_major_locator(dates.MinuteLocator(interval=2))
    axes.set_ylabel('Counts')
    #axes.xaxis.set_major_formatter(DateFormatter("%H:%M"))
    axes.xaxis.set_minor_locator(AutoMinorLocator())
    
    # adding time correction 
    tstart = tstart + timedelta(seconds = tdelta)  #### tdelta = 183.47 for 1st fits and 183.45 for 2nd fits file
    tend = tend + timedelta(seconds = tdelta)
    print(tstart,
    tend)
    axes.axvspan(tstart, tend, alpha = 0.1, color = 'black')
    

    return tstart

####### Create aia image with stix contours

##### Get a proper helioprojective header for stix and create a sunpy stix map

def stix_hpj_map(fitsfile, aia_map):
	
	hdu1 = fits.open(fitsfile)
	header = hdu1[0].header
	data = hdu1[0].data
	#pdb.set_trace()

	solo_hgs = get_horizons_coord('solo', time=header["DATE-OBS"])

	ref_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective(observer=solo_hgs, obstime=header["DATE-OBS"]))

	new_header = sunpy.map.make_fitswcs_header(data, 
                                       ref_coord,
                                       reference_pixel=[int(header["CRPIX1"])-1, int(header["CRPIX2"])-1]*u.pixel, 
                                       scale=[float(header["CDELT1"]), float(header["CDELT2"])]*u.arcsec/u.pixel)

	new_header['rsun_ref']= aia_map.meta['rsun_ref']

	stix_map = sunpy.map.Map(data, new_header)
	map_start = header['DATE-OBS']
	map_endt = header['DATE-END']
	#pdb.set_trace()
	
	return stix_map, map_start, map_endt	

###### Function to plot aia with stix contours

def img(fitsfile,ax, color = 'Blue'):

	stix_map,startt,endt = stix_hpj_map(fitsfile,aia_map)

	
	##### draw stix contours
	CS = stix_map.draw_contours(levels=np.arange(70, 100, 10)*u.percent, colors = color,axes = ax)
	h1,_ = CS.legend_elements()
	return h1, startt, endt
	
#### find the index of closest time of aia to stix time
def find_aia(aia_t,dtime):
	tidx = np.abs(np.array(aia_t) - dtime).argmin()
	return tidx
	
### read aia and stix fits file and plot 

path1 = '/home/shilpi/stix_idl_files/testing/'

fitsfile = sorted(glob.glob(path1+'*.fits'))
path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/'
counts = np.load(path+'stix_tseries_data/counts.npy')
times = np.load(path+'stix_tseries_data/corr_times.npy',allow_pickle = True)

#pdb.set_trace()
path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/aia_data/'

aia_map = sunpy.map.Map(path + 'data_1700/AIA.20221111_113940.1700.image_lev1.fits')
#aia_map.plot_settings['norm'] = plt.Normalize(100, 3000)

aia_list = sorted(glob.glob('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/aia_data/all_data_1700/*.fits'))

aia_t = []
for i in range(len(aia_list)):
		hdu = fits.open(aia_list[i])
		time = hdu[1].header['DATE-OBS'][:19]
		#pdb.set_trace()
		aia_t.append(datetime.strptime(time, '%Y-%m-%dT%H:%M:%S'))


#ax1.set_title(title, fontsize = 14)
colors = ['green','royalblue','Blue','black']
label = [ '10-15 keV','15-25 keV','25-50 keV ','4-10 keV ']
#color = 'coolwarm'
#cm = pylab.get_cmap(color)
#collist = cm(np.linspace(0, 255, 10).astype(int))

e_range = len(label)-1 ### this is the number of energy bins -1

for i in range(len(fitsfile)):
	h = []
	k = i+i*e_range
	print('making image')
	print('defining the axes')
	fig2 = plt.figure(figsize=(12, 10))
	outer = gridspec.GridSpec(3, 1, wspace=0.9, hspace=0.5)
	outer.update(left=0.1, right=0.9, top=0.95, bottom=0.05, wspace=0.3)
	inner0 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0:1,:3], wspace=0.1, hspace=0.4)
        
	######################
	####### plotting stix + aia image
	inner1 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[1:,:3], wspace=0.1, hspace=0.04)           ##### from 0th row and 1st column  

	
	ax1 = plt.subplot(inner1[0], projection = aia_map)
	aia_obs_time = aia_map.meta['date-obs'][11:19]

	##### Setting limit in both x and y axis
	xlims_world = [144, 234]*u.arcsec
	ylims_world = [122, 189]*u.arcsec
	world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_map.coordinate_frame)
	pixel_coords = aia_map.world_to_pixel(world_coords)

	# we can then pull out the x and y values of these limits.
	xlims_pixel = pixel_coords.x.value
	ylims_pixel = pixel_coords.y.value
	ax1.set_xlim(xlims_pixel)
	ax1.set_ylim(ylims_pixel)
	
	# x and y label and title
	ax1.set_xlabel('Helioprojective X coord',fontsize = 10 )
	ax1.set_ylabel('Helioprojective Y coord',fontsize = 10)
	#aia_map.plot(axes = ax1,cmap=sdoaia1700)
	for j in range(len(colors)):
		print(fitsfile[i+j])
		
		leg, startt, endt = img(fitsfile[k+j],ax1, color = colors[j])
		h.append(leg[0])
	
	
	date_format = '%Y-%m-%dT%H:%M:%S.%f'
	Start = datetime.strptime(startt, date_format)
	Endt = datetime.strptime(endt, date_format)
	#pdb.set_trace()
	#### Plotting stix timeseries
	ax = plt.subplot(inner0[0])
	
	tstart = summary_plot(counts, times, ax, tstart = Start, tend = Endt, tdelta = 183.45 )
	
	##### Plotting aia image
	index  =  find_aia(aia_t,tstart)
	aia_map1 = sunpy.map.Map(aia_list[index])
	aia_map1.plot_settings['norm'] = plt.Normalize(100, 5000)
	aia_map1.plot(axes = ax1,cmap=sdoaia1700)
	
	ax1.legend(h,label,fontsize = "9") 
	title = startt+" - "+endt[11:]
	ax1.set_title(title)
	
	plt.savefig('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/stix_images/testing_2/image_{:03d}.png'.format(i))

	#plt.show()  
	



