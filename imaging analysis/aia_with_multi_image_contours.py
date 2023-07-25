##### This code plots aia images with NRH radio contours for three timeslices

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
from stixpy.net.client import STIXClient
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
matplotlib.rcParams.update({'font.size': 14})
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
from scipy.io import readsav
# -



####### Create aia image with stix contours

	

###### Function to plot aia with stix contours

def img(fitsfile, aia_map,ax, fig,cmap = sdoaia171):

	aia_map.plot(axes = ax,cmap=sdoaia171)
	
	##### Setting limit in both x and y axis
	xlims_world = [100, 600]*u.arcsec
	ylims_world = [-350, 220]*u.arcsec
	world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_map.coordinate_frame)
	pixel_coords = aia_map.world_to_pixel(world_coords)

	# we can then pull out the x and y values of these limits.
	xlims_pixel = pixel_coords.x.value
	ylims_pixel = pixel_coords.y.value
	ax.set_xlim(xlims_pixel)
	ax.set_ylim(ylims_pixel)
	
	# x and y label and title
	ax.set_xlabel('Solar X coord',fontsize = 10 )
	ax.set_ylabel('Solar Y coord',fontsize = 10)
	
	
def nrh_sunpy_map(fname, dtime =  b'2022-11-11T11:39:42.570Z'):

	nrh = readsav(fname, python_dict=True)
	
	dtime = dtime + nrh["nrh_hdr"]['DATE_OBS'][0][19:]
	
	index = np.where(nrh["nrh_hdr"]['DATE_OBS'] ==dtime )
	
	#pdb.set_trace()

	b = {name:nrh["nrh_hdr"][name][index[0][0]] for name in nrh["nrh_hdr"].dtype.names}

	map_data = nrh["nrh_data"][index[0][0]]


	ref_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, 
                     frame=frames.Helioprojective(observer="earth", obstime=b["DATE_OBS"]), 
                     )

	header = sunpy.map.make_fitswcs_header(map_data, 
                                       ref_coord, 
                                       reference_pixel=[int(b["CRPIX1"])-1, int(b["CRPIX2"])-1]*u.pixel, 
                                       scale=[float(b["CDELT1"]), float(b["CDELT2"])]*u.arcsec/u.pixel, 
                                       wavelength=float(b["FREQ"])*u.MHz)

	nrh_map = sunpy.map.Map(map_data, header)
 
	return nrh_map

def plot_contours(fname, aia_map,ax, dtime =  b'2022-11-11T11:39:42', color = 'red',j=0 ):
	
	aia_map.plot(axes = ax,cmap=sdoaia171,title=False)
	nrh_map = nrh_sunpy_map(fname, dtime = dtime)
	
	with frames.Helioprojective.assume_spherical_screen(aia_map.observer_coordinate):
		nrh_map2 = nrh_map.reproject_to(aia_map.wcs)
		
	CS = nrh_map.draw_contours(levels=np.arange(85, 100, 5)*u.percent, colors = color, axes = ax)
	h,_ = CS.legend_elements()
	
	##### Setting limit in both x and y axis
	xlims_world = [0, 1000]*u.arcsec
	ylims_world = [-300, 400]*u.arcsec
	world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_map.coordinate_frame)
	pixel_coords = aia_map.world_to_pixel(world_coords)

	# we can then pull out the x and y values of these limits.
	xlims_pixel = pixel_coords.x.value
	ylims_pixel = pixel_coords.y.value
	ax.set_xlim(xlims_pixel)
	ax.set_ylim(ylims_pixel)
	
	if j ==0: 
		ax.set_ylabel(' ')
	else: 
		ax.set_ylabel('Solar Y')
	ax.set_xlabel('Solar X')
	ax.set_title(aia_map.meta['date-obs'][11:19])
	
	return h

#### nrh data
path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/all_radio_data/'
nrh_432 = path+'nrh_maps/nrh_maps_221111_113201_114259_432.sav'
nrh_408 = path+'nrh_maps/nrh_maps_221111_113201_114259_408.sav'
nrh_150 = path+'nrh_maps/nrh_maps_221111_113200_114259_150.sav'
nrh_173 = path+'nrh_maps/nrh_maps_221111_113200_114259_173.sav'
nrh_228 = path+'nrh_maps/nrh_maps_221111_113200_114300_228.sav'
nrh_270 = path+'nrh_maps/nrh_maps_221111_113200_114259_270.sav'
nrh_298 = path+'nrh_maps/nrh_maps_221111_113200_114259_298.sav'
nrh_327 = path+'nrh_maps/nrh_maps_221111_113200_114259_327.sav'


### read aia and stix fits file and plot 

path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/aia_data/'

aia_map1 = sunpy.map.Map(path + 'all_data_171/AIA.20221111_113845.0171.image_lev1.fits')
aia_map2 = sunpy.map.Map(path + 'all_data_171/AIA.20221111_113921.0171.image_lev1.fits')
aia_map3 = sunpy.map.Map(path + 'all_data_171/AIA.20221111_114045.0171.image_lev1.fits')
#aia_map.plot_settings['norm'] = plt.Normalize(100, 3000)


############## Defining Axes
print('defining the axes')

####### Plot stix timeseries and stix contours on aia 1700
fig2 = plt.figure(figsize=(12, 7))
outer = gridspec.GridSpec(3, 1, wspace=0.3, hspace=0.2)
outer.update(left=0.1, right=0.9, top=0.95, bottom=0.05, wspace=0.3)
inner0 = gridspec.GridSpecFromSubplotSpec(1, 1,
subplot_spec=outer[0:1,:3], wspace=0.2, hspace=0.4)
inner1 = gridspec.GridSpecFromSubplotSpec(1, 3,
	subplot_spec=outer[1:,:3], wspace=0.3, hspace=0.2)

##### Plot orfees
orfees_data = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/orfees_data.npy')
orfees_time = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/orfees_time.npy', allow_pickle = True)
orfees_freq = np.load('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/orfees_freqs.npy')
pdb.set_trace()
axes = plt.subplot(inner0[0])


axes.pcolormesh(orfees_time, orfees_freq[0], orfees_data,norm = colors.LogNorm(),cmap='viridis')


axes.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S')) 
axes.xaxis.set_major_locator(dates.MinuteLocator(interval=2))

axes.set_ylabel('Frequency (MHz)')
axes.set_ylim(top = 435)

axes.set_xlim(datetime.datetime(2022,11,11, 11,37), datetime.datetime(2022,11,11,11,43))

tmark = [datetime.datetime(2022,11,11, 11,38,46),datetime.datetime(2022,11,11, 11,39,22),datetime.datetime(2022,11,11, 11,40,45)]
for i in range(len(tmark)):
	axes.axvline(x=tmark[i], alpha = 0.7, color = 'white', linestyle = '--', linewidth = 0.8)


axes.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
axes.xaxis.set_major_locator(dates.MinuteLocator(interval=2))
axes.xaxis.set_minor_locator(AutoMinorLocator())

###### Plot aia with radio contours
ax = plt.subplot(inner1[0],projection = aia_map1)
ax1 = plt.subplot(inner1[1],projection = aia_map2)
ax2 = plt.subplot(inner1[2],projection = aia_map3)

color = 'coolwarm'
cm = pylab.get_cmap(color)
collist = cm(np.linspace(0, 255, 9).astype(int))
h1=[]
h2=[]
h3=[]

fname = [nrh_432, nrh_408, nrh_327, nrh_298, nrh_270, nrh_228, nrh_173, nrh_150]
for i in range(len(fname)):
	leg = plot_contours(fname[i], aia_map1,ax, dtime =  b'2022-11-11T11:38:46', color = collist[i],j=1)
	h1.append(leg[0])
for i in range(len(fname)):
	leg=plot_contours(fname[i], aia_map2,ax1, dtime =  b'2022-11-11T11:39:22', color = collist[i],j=0)
	h2.append(leg[0])
	
for i in range(len(fname)):
	leg=plot_contours(fname[i], aia_map3,ax2, dtime =  b'2022-11-11T11:40:45', color = collist[i],j=0)
	h3.append(leg[0])





plt.show()











