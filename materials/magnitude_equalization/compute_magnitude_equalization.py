# %% Compute data driven magnitude equalization function
# for SUpDEq interpolation
#
# Fabian Brinkmann 06/2024
import warnings
import pyfar as pf                  # v0.6.6
import sofar as sf                  # v1.3.3
import numpy as np                  # v1.26.4
import matplotlib.pyplot as plt     # v3.9.0
import glob
import os
%matplotlib qt

warnings.filterwarnings(
    "ignore", category=pf.classes.warnings.PyfarDeprecationWarning)

# 'SOMICOM', 'HUTUBS_measured', and 'HUTUBS_simulated' implemented
database = 'SONICOM'

# specify which SOFA files are used
if database == 'SONICOM':
    # all data must in SONICOM folder in subfolder P*, e.g., SONICOM/P0001
    subjects = sorted(glob.glob('SONICOM/P*'))
    sampling_rate = 48
    kind = 'FreeFieldComp'
elif database == 'HUTUBS_measured':
    subjects = glob.glob('HUTUBS/HRIRs/*measured.sofa')
elif database == 'HUTUBS_simulated':
    subjects = glob.glob('HUTUBS/HRIRs/*simulated.sofa')
else:
    raise ValueError(f'Unknown database {database}')

if subjects:
    n_subjects = len(subjects)
else:
    raise ValueError(f'{database} database not found')

# Width of fractional octave smoothing
smoothing_width = 6

for nn, subject in enumerate(subjects):

    # load data
    if database == 'SONICOM':
        file = os.path.join(
            subject, 'HRTF', 'HRTF', f'{sampling_rate}kHz',
            f'{os.path.basename(subject)}_{kind}_{sampling_rate}kHz.sofa')
    elif database.startswith('HUTUBS'):
        file = subject

    sofa = sf.read_sofa(file, verbose=False)
    hrirs, sources, _ = pf.io.convert_sofa(sofa)

    # allocate array for magnitude reference
    if nn == 0:
        magnitude_equalization = np.zeros_like(hrirs.freq)

    # add current subject to magnitude reference
    if smoothing_width:
        hrirs = pf.dsp.smooth_fractional_octave(
        hrirs, smoothing_width)[0]

    magnitude_equalization += np.abs(hrirs.freq)

# %% finalize and plot reference

if database == 'SONICOM':
    offset = -15
    name = f'{database}_{sampling_rate}k_{kind}_smoothed_{smoothing_width}'
else:
    offset = 0
    name = f'{database}_smoothed_{smoothing_width}'

# average and make signal
magnitude_equalization_signal = pf.Signal(
    magnitude_equalization/n_subjects, hrirs.sampling_rate, hrirs.n_samples,
    'freq')

# write as sofa
mag_reference_sofa = sf.read_sofa(file, verbose=False)
mag_reference_sofa.Data_IR = magnitude_equalization_signal.time
sf.write_sofa(name + '.sofa', mag_reference_sofa)

# (ugly) plots as sanity check

# median plane
idx = np.where(sources.lateral == 0)[0]
idx = idx[np.argsort(sources.polar[idx])]
gains = np.arange(0, -5*len(idx), -5)
ax = pf.plot.freq(magnitude_equalization_signal[idx, 0] * 10**(gains/20))
ax.set_ylim(gains[-1] - 40 - offset, 20)
for nn, gain in enumerate(gains):
    plt.text(200, gain + offset,
             f'{np.round(sources.polar[idx[nn]] / np.pi * 180)} deg')
plt.savefig(f'{name}_median_plane.pdf')
plt.close()

# horizontal plane
idx = np.where(np.logical_and(
    sources.elevation == 0,
    (sources.azimuth / np.pi * 180) % 20 < .1))[0]
idx = idx[np.argsort(sources.azimuth[idx])]
gains = np.arange(0, -5*len(idx), -5)
ax = pf.plot.freq(magnitude_equalization_signal[idx, 0] * 10**(gains/20))
ax.set_ylim(gains[-1] - 40 - offset, 20)
for nn, gain in enumerate(gains):
    plt.text(200, gain + offset,
             f'{np.round(sources.azimuth[idx[nn]] / np.pi * 180)} deg')
plt.savefig(f'{name}_horizontal_plane.pdf')
plt.close()
