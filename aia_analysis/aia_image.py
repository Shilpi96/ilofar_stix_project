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

path = '/mnt/LOFAR-PSP/ilofar_STIX_shilpi/aia_data/data_171/'
aia_map = sunpy.map.Map(path + 'AIA.20221111_113009.0171.image_lev1.fits')
fig = plt.figure()
ax = fig.add_subplot(projection=aia_map)
aia_map.plot()

##### Setting limit in both x and y axis
xlims_world = [-300, 800]*u.arcsec
ylims_world = [-500, 700]*u.arcsec
world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_map.coordinate_frame)
pixel_coords = aia_map.world_to_pixel(world_coords)

# we can then pull out the x and y values of these limits.
xlims_pixel = pixel_coords.x.value
ylims_pixel = pixel_coords.y.value
ax.set_xlim(xlims_pixel)
ax.set_ylim(ylims_pixel)


plt.show()
