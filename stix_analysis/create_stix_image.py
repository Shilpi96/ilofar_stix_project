from astropy.io import fits
from astropy import units as u 
from astropy.coordinates import SkyCoord
import sunpy.map
import matplotlib.pyplot as plt
from sunpy.net import Fido, attrs as a
import numpy as np
from sunpy.time import parse_time
from astropy.visualization import ImageNormalize, SqrtStretch
from matplotlib import colors
from sunpy.coordinates import frames
import pdb

path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/stix_data/stix_map/'
hdu = fits.open(path + 'mymap.fits')
header = hdu[0].header
data = hdu[0].data

ref_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, 
                     frame=frames.Helioprojective(observer="earth", obstime=header["DATE_OBS"]), 
                     )

new_header = sunpy.map.make_fitswcs_header(data, 
                                       ref_coord, 
                                       reference_pixel=[int(header["CRPIX1"])-1, int(header["CRPIX2"])-1]*u.pixel, 
                                       scale=[float(header["CDELT1"]), float(header["CDELT2"])]*u.arcsec/u.pixel)

stix_map = sunpy.map.Map(data, new_header)

pdb.set_trace()

fig = plt.figure()
ax = fig.add_subplot(projection=stix_map)
stix_map.plot()
stix_map.plot(cmap="viridis")
stix_map.draw_limb(color='w')
stix_map.draw_grid(color='w')

##### reprojecting stix to aia

path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/aia_data/data_171/'
map_aia = sunpy.map.Map(path + 'AIA.20221111_113009.0171.image_lev1.fits')


out_header = sunpy.map.make_fitswcs_header(stix_map.data,coordinate = map_aia.reference_coordinate.replicate(rsun=stix_map.reference_coordinate.rsun),scale=u.Quantity(map_aia.scale))

pdb.set_trace()

outmap = stix_map.reproject_to(out_header)
outmap.data[0:129] = stix_map.data[0:129]

fig = plt.figure()
ax1 = fig.add_subplot(121, projection=stix_map)
stix_map.plot(axes=ax1)
ax2 = fig.add_subplot(122, projection=outmap)
outmap.plot(axes=ax2, title='stix image as seen from SDO')

stix_map.draw_limb(color='blue')


plt.show()
