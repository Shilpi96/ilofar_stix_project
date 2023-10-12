#### This script will produce aia 171 & 1700 images with stix contours for 4-10 and 15-30 kev for 4 time intervals.

from astropy.io import fits
from astropy import units as u 
from astropy.coordinates import SkyCoord
import sunpy.map, matplotlib
import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
import numpy as np
from sunpy.time import parse_time
from astropy.visualization import ImageNormalize, SqrtStretch
from matplotlib import colors
sdoaia171 = matplotlib.colormaps['sdoaia171']

##### This code plots stix timeseries on top axis woth marked time range for which the image reconstruction has been done

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
matplotlib.rcParams.update({'font.size': 15})
import pdb
from sunpy.time import TimeRange
from stixpy.science import ScienceData
from itertools import product
from matplotlib.dates import DateFormatter
import datetime
from matplotlib.ticker import AutoMinorLocator
from matplotlib import dates
import sunpy.visualization.colormaps as cm
import pylab
import sunpy.map
from sunpy.coordinates.ephemeris import get_horizons_coord
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import matplotlib.gridspec as gridspec
import sunpy.visualization.colormaps as cm
sdoaia171 = matplotlib.colormaps['sdoaia171']
sdoaia1700 = matplotlib.colormaps['sdoaia1700']
# -



####### Create aia image with stix contours

##### Get a proper helioprojective header for stix and create a sunpy stix map

def stix_hpj_map(fitsfile, aia_map):
	
	hdu1 = fits.open(fitsfile)
	header = hdu1[0].header
	data = hdu1[0].data

	solo_hgs = get_horizons_coord('solo', time=header["DATE-OBS"])

	ref_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=frames.Helioprojective(observer=solo_hgs, obstime=header["DATE-OBS"]))

	new_header = sunpy.map.make_fitswcs_header(data, 
                                       ref_coord,
                                       reference_pixel=[int(header["CRPIX1"])-1, int(header["CRPIX2"])-1]*u.pixel, 
                                       scale=[float(header["CDELT1"]), float(header["CDELT2"])]*u.arcsec/u.pixel)

	new_header['rsun_ref']= aia_map.meta['rsun_ref']
	
	stix_map = sunpy.map.Map(data, new_header)
	
	return stix_map	

###### Function to plot aia with stix contours

def img(fitsfile,ax,aia_map,color = 'Blue'):

	stix_map = stix_hpj_map(fitsfile,aia_map)

	##### draw stix contours
	stix_map.draw_contours(levels=np.arange(70, 100, 10)*u.percent, colors = color,axes = ax)
	#pdb.set_trace()

	##### Setting limit in both x and y axis
	xlims_world = [144, 234]*u.arcsec
	ylims_world = [100, 189]*u.arcsec
	world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_map.coordinate_frame)
	pixel_coords = aia_map.world_to_pixel(world_coords)

	# we can then pull out the x and y values of these limits.
	xlims_pixel = pixel_coords.x.value
	ylims_pixel = pixel_coords.y.value
	ax.set_xlim(xlims_pixel)
	ax.set_ylim(ylims_pixel)
	

### read aia and stix fits file and plot 

path1 = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/stix_data/paper_map/'


fits_4 = [ path1 +'map_11:33:52-4-10.fits',path1 +'map_11:34:51-4-10.fits',path1 +'map_11:38:38-4-10.fits',path1 +'map_11:39:30-4-10.fits']
fits_15 = [ path1 +'map_11:33:52-15-30.fits',path1 +'map_11:34:51-15-30.fits',path1 +'map_11:38:38-15-30.fits',path1 +'map_11:39:30-15-30.fits']

path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/aia_data/'

data_171 = [path + 'all_data_171/AIA.20221111_113357.0171.image_lev1.fits',path + 'all_data_171/AIA.20221111_113457.0171.image_lev1.fits',path + 'all_data_171/AIA.20221111_113845.0171.image_lev1.fits',path + 'all_data_171/AIA.20221111_113933.0171.image_lev1.fits']

data_1700 = [path + 'all_data_1700/AIA.20221111_113404.1700.image_lev1.fits',path + 'all_data_1700/AIA.20221111_113452.1700.image_lev1.fits',path + 'all_data_1700/AIA.20221111_113852.1700.image_lev1.fits',path + 'all_data_1700/AIA.20221111_113940.1700.image_lev1.fits']

map_171 = [sunpy.map.Map(data_171[i]) for i in range(len(data_171))]
map_1700 = [sunpy.map.Map(data_1700[i]) for i in range(len(data_1700))]

#aia_map.plot_settings['norm'] = plt.Normalize(300, 3000)

############## Defining Axes
print('defining the axes')
fig2 = plt.figure(figsize=(9, 7))

####### Plot stix timeseries and stix contours on aia 1700

for i in range(4):

	map_1700[i].plot_settings['norm'] = plt.Normalize(300, 3000)
	#map_171[i].plot_settings['norm'] = plt.Normalize(300, 1000)
	
	ax1 = plt.subplot(2,4,i+1,projection = map_1700[i])
	ax2 = plt.subplot(2,4,i+5,projection = map_171[i])
	
	map_1700[i].plot(axes = ax1,cmap=sdoaia1700)
	map_171[i].plot(axes = ax2,cmap=sdoaia171)
	
	img(fits_4[i],ax1,map_1700[i],color = 'Black')
	img(fits_15[i],ax1,map_1700[i],color = 'Blue')
	img(fits_4[i],ax2,map_171[i],color = 'Black')
	img(fits_15[i],ax2,map_171[i],color = 'Blue')
	
	
	# x and y label and title
	ax1.set_title(fits_4[i][58:66])
	ax2.set_title('  ')
	axes = [ax1,ax2]
	[axes[i].set_xlabel('Solar X (arcsec)') for i in range(2)]
	[axes[i].set_ylabel('Solar Y (arcsec)') for i in range(2)]
	
	ax1.axes.coords[0].set_ticklabel_visible(False)
	[axes[i].axes.coords[1].set_ticklabel_visible(False) for i in range(2)]
	if i==0:ax1.axes.coords[1].set_ticklabel_visible(True)
	if i==0:ax2.axes.coords[1].set_ticklabel_visible(True)
	

plt.show()

