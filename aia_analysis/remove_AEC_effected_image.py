import astropy.units as u
from sunpy.net import Fido, attrs
import sunpy.map
import os
import glob

list_171 = glob.glob('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/aia_data/data_171'+'/*.fits')
list_193 = glob.glob('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/aia_data/data_171'+'/*.fits')
list_211 = glob.glob('/mnt/LOFAR-PSP/ilofar_STIX_shilpi/aia_data/data_171'+'/*.fits')

###### Now we have to delete files for which the exposure time is less than the standard time. for 171 = 1.5 sec, 211 = 1.5 sec and 193 = 1.0 sec. because during solar eruptions the image can get saturated so the camera will try to close more quickly because it will try to receive less photons. That means we are not getting correct information. So need to delete these images.
######## Check the images to make sure we're not using AEC-affected images and delete these AEC-affected images from the list
min_exp_t_193 = 1.0
min_exp_t_211 = 1.5
min_exp_t_171 = 1.5


## 211 channel
for i in range(len(list_211)):
	aia = sunpy.map.Map(list_211[i])
	if aia.meta['exptime']<min_exp_t_211 : os.remove(list_211[i])

	else: 
		print('there is no AEC effected images in 211')

## 193 channel
for i in range(len(list_193)):
	aia = sunpy.map.Map(list_193[i])
	if aia.meta['exptime']<min_exp_t_193 : os.remove(list_193[i])

	else: 
		print('there is no AEC effected images in 193')

## 171 channel
for i in range(len(list_171)):
	aia = sunpy.map.Map(list_171[i])
	if aia.meta['exptime']<min_exp_t_171 : os.remove(list_171[i])

	else: 
		print('there is no AEC effected images in 171')

