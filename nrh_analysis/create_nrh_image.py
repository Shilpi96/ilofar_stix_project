from astropy.io import fits
from astropy import units as u 
from astropy.coordinates import SkyCoord
import sunpy.map
import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
from scipy.io import readsav
from sunpy.coordinates import frames
import numpy as np
import pandas as pd
from sunpy.time import parse_time
from astropy.visualization import ImageNormalize, SqrtStretch
from matplotlib import colors
import pdb

####### Read nrh data 

path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/all_radio_data/'
files = path+'nrh_fits/nrh2_1509_h70_20221111_113200c05_b.fts'
'''
aa = fits.open(files)
aa[1].header # header

tstart = aa[0].header["DATE_OBS"]
tend = aa[0].header["DATE_END"]

aa[1].data.dtype.names
aa[1].data["STOKESI"].shape
'''

nrh = readsav(path+'nrh_maps/nrh_maps_221111_113200_114259_298.sav', python_dict=True)
#pdb.set_trace()
b = {name:nrh["nrh_hdr"][name][0] for name in nrh["nrh_hdr"].dtype.names} # getting header info for the first time index
pdb.set_trace()
map_data = nrh["nrh_data"][0]


ref_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, 
                     frame=frames.Helioprojective(observer="earth", obstime=b["DATE_OBS"]), 
                     )

header = sunpy.map.make_fitswcs_header(map_data, 
                                       ref_coord, 
                                       reference_pixel=[int(b["CRPIX1"])-1, int(b["CRPIX2"])-1]*u.pixel, 
                                       scale=[float(b["CDELT1"]), float(b["CDELT2"])]*u.arcsec/u.pixel, 
                                       wavelength=float(b["FREQ"])*u.MHz)

nrh_map = sunpy.map.Map(map_data, header)
pdb.set_trace()
nrh_map.plot(cmap="viridis")
nrh_map.draw_limb(color='w')
nrh_map.draw_grid(color='w')
plt.show()
