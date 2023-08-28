##### Using serpentine this code plots stereo A 'SEPT' electron flux timeseries from (45-169 keV)

from seppy.loader.stereo import stereo_load
import datetime as dt
from matplotlib import pyplot as plt
import pdb
import numpy as np

instrument = "SEPT"
startdate = dt.datetime(2022, 11, 11)
enddate = dt.datetime(2022, 11, 12)
path = None
resample = "5min"
pos_timestamp = None

df, meta = stereo_load(instrument=instrument, startdate=startdate, enddate=enddate,
                       path=path, resample=resample, pos_timestamp=pos_timestamp)
#pdb.set_trace()
#df.Electron_Flux_1.plot(logy=True, ylabel=meta['Electron_Flux_UNITS'], label=meta['Electron_Bins_Text'][1][0], 
#                        title='STEREO/HET electrons')
##### save the data, e_bins and times
root = '/Users/shilpibhunia/Library/Mobile Documents/com~apple~CloudDocs/ilofar_stix_project/codes/epd_analysis/stereo_data'
np.save(root+'sept_data.npy',data)
np.save(root+'e_bins.npy',bins_array)
np.save(root+'times.npy',times)

fig, ax = plt.subplots()
channels = ["ch_"+str(i+2) for i in range(10)]
bins = meta["ch_strings"].tolist()
ax.set_prop_cycle('color', plt.cm.winter(np.linspace(0,1,len(channels))))

for i in range(len(channels)):
	df[channels[i]].plot(logy=True, ax = ax, label=bins[i])

ax.set_xlim(dt.datetime(2022,11,11, 10), dt.datetime(2022,11,11,14))
#ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., title=f'Electrons ({viewing})')
ax.legend(loc = 'lower right', fontsize=9)
ax.set_ylabel('Electron flux')
plt.show()
