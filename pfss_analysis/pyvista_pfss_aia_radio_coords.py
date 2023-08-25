########## This script uses sunkit pyvista to plot aia in 3d with pfss and pos of radio coords.

import os
import astropy.units as u
import matplotlib.pyplot as plt
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

import pfsspy
import pfsspy.utils
import pfsspy.tracing as tracing
import numpy as np 
from astropy.coordinates import SkyCoord
import astropy.constants as const
from sunpy.coordinates import frames
from sunkit_pyvista import SunpyPlotter
from astropy.constants import R_sun
from matplotlib import colors


root = '/Users/shilpibhunia/Library/Mobile Documents/com~apple~CloudDocs/ilofar_stix_project/aia_hmi_data/'
aia = sunpy.map.Map(root+'AIA.20221111_113509.0171.image_lev1.fits')

hmi_file = [root+'hmi.synoptic_mr_polfil_720s.2264.Mr_polfil.fits']

hmi_map = sunpy.map.Map(hmi_file[0])
hmi_map = hmi_map.resample([360, 180] * u.pix)


nrho = 35
rss = 2.5
pfss_in = pfsspy.Input(hmi_map, nrho, rss)

hp_lon = np.linspace(-600, 600, 20) * u.arcsec
hp_lat = np.linspace(-600, 600, 20) * u.arcsec
# Make a 2D grid from these 1D points
lon, lat = np.meshgrid(hp_lon, hp_lat)
seeds = SkyCoord(lon.ravel(), lat.ravel(),
                 frame=aia.coordinate_frame)


pfss_out = pfsspy.pfss(pfss_in)
tracer = tracing.FortranTracer()
flines = tracer.trace(seeds, pfss_out)


####### Create 2d plot of aia with pfss

fig = plt.figure()
ax = plt.subplot(1, 1, 1, projection=aia)
aia.plot(ax)
# for fline in flines:
#     ax.plot_coord(fline.coords, alpha=0.8, linewidth=1, color='white')

for field_line in flines:
    color = {0: 'white', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity)
    coords = field_line.coords
    coords.representation_type = 'cartesian'
    ax.plot_coord(coords,
            	  color=color, linewidth=1)
# ax.set_xlim(500, 900)
# ax.set_ylim(400, 800)
plt.show()


######### Create 3d AIA with pfss

plotter = SunpyPlotter()
low_res_aia = aia.resample([512, 512] * u.pix)
# Plot a map
plotter.plot_map(low_res_aia, clip_interval=[1, 99] * u.percent)
# Add an arrow to show the solar rotation axis
plotter.plot_solar_axis()
'''
def my_fline_color_func(field_line):    ##### you can change the polarity color
    color = {0: "grey", -1: "tab:blue", 1: "tab:red"}.get(
                    field_line.polarity,
                )
    return colors.to_rgb(color)
'''
'''
norm = colors.LogNorm(vmin=1, vmax=1000) #### you can also chnage it to a colormap
cmap = plt.get_cmap("viridis")
return cmap(norm(np.abs(field_line.expansion_factor)))'''


# Plotting the field lines
print('Plot the field lines in 3d aia')
plotter.plot_field_lines(flines)
camera_coord = SkyCoord(
    0 * u.deg,
    0 * u.deg,
    6 * R_sun,
    frame=frames.HeliographicStonyhurst,
    obstime=low_res_aia.date,
)
plotter.set_camera_coordinate(camera_coord)

########### Plot the type III source positions for one soingle timestamp
x = np.array([0.2846517938334465, 0.316279727554574, 0.2846517938334465,
    0.2846517938334465, 0.31627972755455264, 0.3479076475794979,
    0.316279727554574, 0.3479076475794979])*R_sun
y = np.array([0.03162796761938121, 0.06325592085769638, 0.03162796761938121,
    0.03162796761938121, 0.03162796111366393, -0.06325590647663047,
    -0.06325592085767008, -0.06325590647663047])*R_sun
z = np.array([1.017510448666871, 1.019218634708425, 1.074262205710927,
    1.1050586662864654, 1.1239071197658113, 1.1543750974458407,
    1.2530858585022187, 1.3002348532319432])*R_sun

coords = SkyCoord(x, y, z, frame=frames.Heliocentric(observer="earth", obstime=aia.date))  #### heliocentric is typically used with cartesian coords
cools = plt.cm.viridis(np.linspace(0,1,  len(coords)))
########### Plot the coords of the sources
print('Plotting the positions of the radio sources')

for i, coord in enumerate(coords):
    plotter.plot_coordinates(coord, color=cools[i], opacity='geom_r')

plotter.show()
