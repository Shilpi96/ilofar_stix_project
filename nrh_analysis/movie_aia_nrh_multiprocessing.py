##### this code will use multiprocessing to plot NRH contours on AIA and create a movie

import warnings
from sunpy.util.exceptions import SunpyDeprecationWarning
warnings.filterwarnings("ignore", category=SunpyDeprecationWarning)
from astropy.io.fits.verify import VerifyWarning
warnings.filterwarnings("ignore", category=VerifyWarning)

import matplotlib
matplotlib.use('Agg')

import os,pdb
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames, sun
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from sunpy.map import Map
from datetime import datetime, timedelta
from astropy.time import Time
import numpy as np
import glob, sunpy.map
from tqdm import tqdm
import multiprocessing as mp

# ------------------ GLOBAL VARIABLES (for workers) ------------------
aia_maps = []
aia_times = []
nrh_hdul = {}
nrh_times = {}
fnrh = []
labels = []
data = None
freq = None
orfees_times = None
str_times = []
output_dir = ''

# ------------------ HELPERS ------------------
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

def find_aia(aia_times, dtime):
    return np.abs(np.array(aia_times) - dtime).argmin()

def plot_contours(aia_map, dtime, fname, color, label, ax, legend_handles):
    global nrh_hdul, nrh_times
    nrh_map, imgtime = nrh_sunpy_map_from_cache(nrh_hdul[fname], nrh_times[fname], dtime)
    nrh_map.meta['rsun_ref'] = aia_map.meta['rsun_ref']
    with frames.Helioprojective.assume_spherical_screen(aia_map.observer_coordinate):
        nrh_map.draw_contours(axes=ax, levels=np.arange(92, 99, 3)*u.percent, colors=color)
    legend_handles.append(Line2D([0], [0], color=color, lw=2, label=label))
    return legend_handles, imgtime

def make_img(fnrh, aia_map, dtime, ax, labels):
    collist = plt.cm.coolwarm(np.linspace(0, 1, len(fnrh)))
    legend_handles = []
    imgtime = []

    for i in range(len(fnrh)):
        legend_handles, t = plot_contours(aia_map, dtime, fnrh[i], collist[i], labels[i], ax, legend_handles)
        imgtime.append(t)

    ax.legend(handles=legend_handles, fontsize='x-small')
    ax.patch.set_facecolor('black')
    xlims_world = [-1950, -550]*u.arcsec
    ylims_world = [-100, 900]*u.arcsec
    world_coords = SkyCoord(Tx=xlims_world, Ty=ylims_world, frame=aia_map.coordinate_frame)
    pixel_coords = aia_map.world_to_pixel(world_coords)
    ax.set_xlim(pixel_coords.x.value)
    ax.set_ylim(pixel_coords.y.value)
    return imgtime

def process_frame(j):
    global fnrh, aia_maps, aia_times, labels, data, freq, orfees_times, str_times, output_dir

    time_aia = datetime(2025, 3, 28, 15, 12) + timedelta(seconds=0.25 * (j))
    print(time_aia)
    tidx = find_aia(aia_times, time_aia)
    aia_map = aia_maps[tidx]

    fig = plt.figure(figsize=(12, 7))
    outer = gridspec.GridSpec(2, 1, wspace=0.9, hspace=0.5)
    outer.update(left=0.1, right=0.9, top=0.95, bottom=0.05)

    inner0 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[0], wspace=0.1, hspace=0.4)
    ax1 = plt.subplot(inner0[0])
    vmm = np.percentile(data, [1, 99])
    ax1.imshow(data, cmap='Spectral_r', norm=LogNorm(vmin=vmm[0], vmax=6),
               aspect='auto', origin='lower',
               extent=(orfees_times[0], orfees_times[-1], freq[0][0], freq[0][-1]))
    ax1.set_xlim(datetime(2025, 3, 28, 15, 12), datetime(2025, 3, 28, 15, 18, 55))
    ax1.set_ylim(149, 500)

    inner1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer[1], wspace=0.2, hspace=0.2)
    ax = plt.subplot(inner1[0], projection=aia_map)
    aia_map.plot(axes=ax, clip_interval=(1, 99.95)*u.percent)

    imgtime = make_img(fnrh, aia_map, time_aia, ax, labels)
    str_times = datetime.strftime(imgtime[0], '%Y-%m-%dT%H:%M:%S.%fZ')
    #pdb.set_trace()
    ax.set_title(str_times)
    ax1.axvline(imgtime[0], alpha=0.7, color='black', linestyle='--', linewidth=0.8)

    outpath = os.path.join(output_dir, f"image_{j:03d}.png")
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)

# ------------------ MAIN FUNCTION ------------------
def main():
    global fnrh, labels, aia_maps, aia_times, nrh_hdul, nrh_times
    global data, freq, orfees_times, str_times, output_dir

    root = '/home/shilpi/march_campaign/event_2025_03_28/data/'
    output_dir = os.path.join(root, 'nrh_aia_imanges1')
    os.makedirs(output_dir, exist_ok=True)

    # Load spectrogram
    data = np.load(root + 'orfee_data.npy')
    freq = np.load(root + 'orfee_freqs.npy')
    orfee_time = np.load(root + 'orfee_times.npy', allow_pickle=True)
    orfees_times = np.array([t.datetime for t in orfee_time])

    # Load AIA maps
    print("Loading AIA maps...")
    aia_files = sorted(glob.glob(root + '171/*.fits'))
    for f in tqdm(aia_files, desc="AIA"):
        hdu = fits.open(f)
        time = datetime.strptime(hdu[1].header['DATE-OBS'][:19], '%Y-%m-%dT%H:%M:%S')
        aia_times.append(time)
        aia_maps.append(Map(f))
    #pdb.set_trace()
    # Load NRH maps
    print("Loading NRH maps...")
    fnrh = sorted(glob.glob(root + 'nrh2_*_b.fts'))
    fnrh = fnrh[::-1]
    #pdb.set_trace()
    labels = [f.split('_')[5] for f in fnrh]
    for fname in tqdm(fnrh, desc="NRH"):
        hdul = fits.open(fname)
        times = get_time(hdul)
        nrh_hdul[fname] = hdul
        nrh_times[fname] = times

    Times = [datetime(2025, 3, 28, 15, 12) + timedelta(seconds=0.25 * i) for i in range(1661)]
    str_times = [t.strftime('%Y-%m-%dT%H:%M:%S') for t in Times]


    indices = list(range(len(Times)))

    # Use multiprocessing Pool
    print("Starting multiprocessing...")
    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.map(process_frame, indices)

if __name__ == '__main__':
    main()


