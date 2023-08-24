######### This plot makes aia with nrh radio contours

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy import units as u 
from astropy.coordinates import SkyCoord
import sunpy.map
from sunpy.net import Fido, attrs as a
from sunpy.time import parse_time
from astropy.visualization import ImageNormalize, SqrtStretch
from matplotlib import colors
from sunpy.coordinates import frames
import sunpy.coordinates
import pdb, matplotlib
from sunpy.coordinates.ephemeris import get_horizons_coord
from scipy.io import readsav
from astropy.visualization import ImageNormalize,LogStretch
import sunpy.visualization.colormaps as cm
sdoaia171 = matplotlib.colormaps['sdoaia171']
import pylab, glob
from datetime import datetime,timedelta
from sklearn.preprocessing import normalize
import matplotlib.gridspec as gridspec

def find_aia(aia_t,dtime):
	tidx = np.abs(np.array(aia_t) - dtime).argmin()
	return tidx

    
def nrh_sunpy_map(fname, dtime =  b'2022-11-11T11:39:42'):

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

def plot_contours(fname, aia_map, ax,dtime =  b'2022-11-11T11:39:42', color = 'Purple',label = '432' ):
	nrh_map = nrh_sunpy_map(fname, dtime = dtime)
	CS = nrh_map.draw_contours(axes=ax,levels=np.arange(85, 100, 5)*u.percent, colors = color)
	h1,_ = CS.legend_elements()
	return h1
		
def make_img(aia_list, tidx,dtime, ax):
	
	h = []

	for i in range(len(clist)):
		leg = plot_contours(fname[i], aia_map, ax,dtime =  dtime, color = clist[i] )
		h.append(leg[0])

	print('plotting the legend')
	ax.legend(h,frq,fontsize = "9")    
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
	
	
	
	
####### Read nrh data 

path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/all_radio_data/'

nrh_432 = path+'nrh_maps/nrh_maps_221111_113201_114259_432.sav'
nrh_408 = path+'nrh_maps/nrh_maps_221111_113201_114259_408.sav'
nrh_150 = path+'nrh_maps/nrh_maps_221111_113200_114259_150.sav'
nrh_173 = path+'nrh_maps/nrh_maps_221111_113200_114259_173.sav'
nrh_228 = path+'nrh_maps/nrh_maps_221111_113200_114300_228.sav'
nrh_270 = path+'nrh_maps/nrh_maps_221111_113200_114259_270.sav'
nrh_298 = path+'nrh_maps/nrh_maps_221111_113200_114259_298.sav'
nrh_327 = path+'nrh_maps/nrh_maps_221111_113200_114259_327.sav'

###### Plot radio contours on aia data

dtime = b'2022-11-11T11:39:19'
fname = [nrh_432, nrh_408, nrh_327, nrh_298, nrh_270, nrh_228, nrh_173, nrh_150]
frq = ['432 MHz','408 MHz', '327 MHz', '298 MHz','270 MHz','228 MHz', '173 MHz', '150 MHz']
clist = ['red','yellow','purple','blue','orange','magenta','cyan','yellowgreen']

times = [datetime(2022,11,11,11,33) +  timedelta(seconds = 1*i) for i in range(210)]

str_times = [datetime.strftime(times[i], '%Y-%m-%dT%H:%M:%S') for i in range(len(times))]
bytes_times = [bytes(str_times[i], 'utf-8') for i in range(len(times))]


##### load ilofar and orfees combined data
path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/plot_high_res_ilofar/'
data = np.load(path+'ilofar_orfees.npy')
freq = np.load(path+'ilofar_orfees_freq.npy',allow_pickle = True)
orfees_time = np.load(path+'time_1sec.npy',allow_pickle = True)

#data2 = normalize(data, axis=1, norm='l1')

####look for closest aia file to radio
aia_list = sorted(glob.glob('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/aia_data/all_data_171/*.fits'))

aia_t = []
for i in range(len(aia_list)):
		hdu = fits.open(aia_list[i])
		time = hdu[1].header['DATE-OBS'][:19]
		#pdb.set_trace()
		aia_t.append(datetime.strptime(time, '%Y-%m-%dT%H:%M:%S'))
		
for j in range(119,210):
	time_aia = times[j]
	tidx = find_aia(aia_t,time_aia)
	
	print('defining the axes')
	fig2 = plt.figure(figsize=(12, 7))
	outer = gridspec.GridSpec(2, 1, wspace=0.9, hspace=0.5)
	outer.update(left=0.1, right=0.9, top=0.95, bottom=0.05, wspace=0.3)
	inner0 = gridspec.GridSpecFromSubplotSpec(1, 1,
        subplot_spec=outer[0:1,:3], wspace=0.1, hspace=0.4)
        
        ##### Plotting the spectra
	ax1 = plt.subplot(inner0[0])
	#pdb.set_trace()
	ax1.pcolormesh(orfees_time,freq,data, cmap = plt.get_cmap('viridis'))
	ax1.set_xlim(datetime(2022,11,11,11,32),datetime(2022,11,11,11,37))
	ax1.axvline(times[j], alpha = 0.7, color = 'white', linestyle = '--', linewidth = 0.8)
        
        ##### Plotting aia map
	inner1 = gridspec.GridSpecFromSubplotSpec(1, 1,
	subplot_spec=outer[1:,:3], wspace=0.2, hspace=0.2) 
        
	aia_map = sunpy.map.Map(aia_list[tidx])
	ax = plt.subplot(inner1[0],projection = aia_map)
	aia_map.plot(axes=ax, title = False)
	make_img(aia_list, tidx, bytes_times[j],ax)
	
	ax.set_title(str_times[j])
	#plt.show()
	### save the figure
	print('Saving the figure')
	plt.savefig('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/nrh_movies/aia_all_radio_cntrs_images_1/image_{:03d}.png'.format(j))
	
   
#plt.show()


    
