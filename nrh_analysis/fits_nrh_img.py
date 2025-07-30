##### this script plot nrh image from the fits files
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from sunpy.map import Map, make_fitswcs_header
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames, sun
import matplotlib.pyplot as plt
import pdb
import sunpy.map
from datetime import datetime

def get_time(hdul):
	times = hdul[1].data['TIME']
	date_obs_str = hdul[1].header['DATE-OBS']
	reference_time = Time(f"{date_obs_str}T00:00:00", scale='utc')
	times_datetime = reference_time + times * u.s / 1000
	times = np.array([datetime.strptime(times_datetime.isot[i], "%Y-%m-%dT%H:%M:%S.%f") for i in range(times_datetime.shape[0])])
	return times
	
def calculate_cdelt_from_header(date,header):
    solar_r_pix = header.get('SOLAR_R')
    if solar_r_pix is None:
        raise ValueError("Le mot-clé SOLAR_R est manquant dans l'en-tête.")

    solar_r_pix = header.get('SOLAR_R')
    time_obj = Time(date)
    solar_r_arcsec = sun.angular_radius(time_obj)
    solar_r_arcsec = 960  # Rayon solaire en arcsec
    cdelt_arcsec = solar_r_arcsec / solar_r_pix
    return -cdelt_arcsec, cdelt_arcsec  # X négatif, Y positif

def nrh_sunpy_map(nrh, dtime =  datetime(2022, 11, 11, 11, 32, 1, 400000)):
    hdul = fits.open(fits_file)
    times = get_time(hdul)
    ind = np.abs(times-dtime).argmin()
    data = hdul[1].data['STOKESI'][ind]
    print(times[ind])
    #pdb.set_trace()
    ref_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, 
                     frame=frames.Helioprojective(observer="earth", obstime=times[ind]), 
                     )
    solar_pix = hdul[1].header['SOLAR_R']
    solar_r_arcsec = sun.angular_radius(times[ind]).value
    cdelta = solar_r_arcsec/solar_pix
    
    header = sunpy.map.make_fitswcs_header(data, 
                                       ref_coord, 
                                       reference_pixel=[int(hdul[1].header['CRPIX1']), int(hdul[1].header['CRPIX2'])]*u.pixel, 
                                       scale=[float(cdelta), float(cdelta)]*u.arcsec/u.pixel, 
                                       wavelength=float(hdul[1].header['FREQ'])*u.MHz)
    #pdb.set_trace()
    nrh_map = Map(data, header)
    freq = nrh_map.meta['wavelnth']
    
    return nrh_map

# read fits file
fits_file = '/home/shilpi/march_campaign/event_2025_03_28/data/nrh2_1509_h80_20250328_150700c04_b.fts'
nrh_map=nrh_sunpy_map(fits_file,dtime =  datetime(2025, 3, 28, 15,10,10,200000))
fig = plt.figure()
ax = fig.add_subplot(projection=nrh_map)
nrh_map.plot(axes=ax, cmap="viridis")
nrh_map.draw_limb(axes=ax, color='w')
#times = hdul[1].data['TIME']
#pdb.set_trace()
plt.show()
