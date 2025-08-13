##### this script look for the closet AIA file for a particular NRH imaging time and then plot NRH contours 
of multiple frequencies on the AIA image.

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
from sunpy.coordinates import frames, screens, sun
import sunpy.coordinates
import pdb, matplotlib
from sunpy.coordinates.ephemeris import get_horizons_coord
from scipy.io import readsav
from astropy.visualization import ImageNormalize,LogStretch
import sunpy.visualization.colormaps as cm
sdoaia171 = matplotlib.colormaps['sdoaia171']
import pylab, glob
from datetime import datetime,timedelta
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from astropy.time import Time
from sunpy.map import Map, make_fitswcs_header
from tqdm import tqdm
from matplotlib.colors import LogNorm

def find_aia(aia_t,dtime):
	tidx = np.abs(np.array(aia_t) - dtime).argmin()
	return tidx
def get_time(hdul):
    times = hdul[1].data['TIME']
    date_obs_str = hdul[1].header['DATE-OBS']
    reference_time = Time(f"{date_obs_str}T00:00:00", scale='utc')
    times_datetime = reference_time + times * u.s / 1000
    return np.array([datetime.strptime(times_datetime.isot[i], "%Y-%m-%dT%H:%M:%S.%f") for i in range(times_datetime.shape[0])])

def nrh_sunpy_map_from_cache(hdul, times, dtime):
    ind = np.abs(times - dtime).argmin()
    data = hdul[1].data['STOKESI'][ind]
    ref_coord = SkyCoord(0*u.arcsec, 0*u.arcsec,
                         frame=frames.Helioprojective(observer="earth", obstime=times[ind]))
    solar_pix = hdul[1].header['SOLAR_R']
    solar_r_arcsec = sun.angular_radius(times[ind]).value
    cdelta = solar_r_arcsec / solar_pix
    header = sunpy.map.make_fitswcs_header(
        data, ref_coord,
        reference_pixel=[hdul[1].header['CRPIX1'], hdul[1].header['CRPIX2']]*u.pixel,
        scale=[cdelta, cdelta]*u.arcsec/u.pixel,
        wavelength=hdul[1].header['FREQ']*u.MHz
    )
    return Map(data, header), times[ind]

    
def plot_contours(aia_map, dtime, fname, color, label, ax, legend_handles):
    global nrh_hdul, nrh_times
    nrh_map, imgtime = nrh_sunpy_map_from_cache(nrh_hdul[fname], nrh_times[fname], dtime)
    nrh_map.meta['rsun_ref'] = aia_map.meta['rsun_ref']
    with frames.Helioprojective.assume_spherical_screen(aia_map.observer_coordinate):
        nrh_map.draw_contours(axes=ax, levels=np.arange(92, 99, 3)*u.percent, colors=color)
    legend_handles.append(Line2D([0], [0], color=color, lw=2, label=label))
    return legend_handles, imgtime

def make_img(fnrh, aia_map, dtime, ax, labels,imgtime,collist):
    legend_handles = []
    for i in range(len(fnrh)):
        legend_handles, t = plot_contours(aia_map, dtime, fnrh[i], collist[i], labels[i], ax, legend_handles)
        imgtime.append(t)

    ax.legend(handles=legend_handles, fontsize='small')
    ax.patch.set_facecolor('black')
    xlims_world = [-1472, -805]*u.arcsec
    ylims_world = [0, 670]*u.arcsec
    world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_map.coordinate_frame)
    pixel_coords = aia_map.world_to_pixel(world_coords)
    ax.set_xlim(pixel_coords.x.value)
    ax.set_ylim(pixel_coords.y.value)
    return imgtime
	
	
####### Read nrh data 

root1 = '/home/shilpi/march_campaign/event_2025_03_28/data/'
fnrh = sorted(glob.glob(root1 + 'nrh2_*_b.fts'))
fnrh = fnrh[::-1]
freqs = [444,432,408,327,298.7,270.6,228,173.2,150.9]
labels = [f.split('_')[5] for f in fnrh]
nrh_hdul = {}
nrh_times = {}
for fname in tqdm(fnrh, desc="NRH"):
        hdul = fits.open(fname)
        times = get_time(hdul)
        nrh_hdul[fname] = hdul
        nrh_times[fname] = times


####look for closest aia file to radio
aia_list = sorted(glob.glob('/home/shilpi/march_campaign/event_2025_03_28/data/171/*.fits'))

aia_t = []
for i in range(len(aia_list)):
		hdu = fits.open(aia_list[i])
		time = hdu[1].header['DATE-OBS'][:19]
		#pdb.set_trace()
		aia_t.append(datetime.strptime(time, '%Y-%m-%dT%H:%M:%S'))
		
Times = datetime(2025,3,28,15,12,47,970000) ### time for NRH imaging
time_aia = Times
tidx = find_aia(aia_t,time_aia)
print('plotting aia with nrh contours')
aia_map = sunpy.map.Map(aia_list[tidx])
fig = plt.figure(figsize=(7, 7))
ax = plt.subplot(111,projection = aia_map)
aia_map.plot(axes=ax, clip_interval=(1, 99.95)*u.percent)
imgtime = []
collist = plt.cm.coolwarm(np.linspace(0, 1, len(fnrh)))
imgtimes = make_img(fnrh[:5], aia_map, time_aia, ax, labels,imgtime,collist)
str_times = datetime.strftime(imgtimes[0], '%Y-%m-%dT%H:%M:%S')
ax.set_title(str_times)
plt.show()
	




