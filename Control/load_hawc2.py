# -*- coding: utf-8 -*-
"""How to load various HAWC2 formats to a NumPy array.
Also does a plot thing because why not.
"""
import h5py
import matplotlib.pyplot as plt
import numpy as np
import os


def load_hawc_ascii(dat_path):
    """Load HAWC ascii data to a NumPy array.
    Input path to .dat file.
    """
    data = np.loadtxt(dat_path)
    return data


def load_hawc_binary(dat_path):
    """Load HAWC binary file to a NumPy array.
    Input path to .dat file.
    """
    sel_path = dat_path.replace('.dat', '.sel')
    with open(sel_path, 'r') as f:
        content = f.readlines()
        # get number lines and channels
        nr_lines, nr_chan = [int(s) for s in content[8].split()[:2]]
        # get the scale factors
        i_scale = [i+1 for i in range(len(content)) if 'scale factors' in content[i].lower()][0]
        scale_factors = np.loadtxt(sel_path, skiprows=i_scale)
    with open(dat_path, 'rb') as fid:
        data = np.zeros((nr_lines, nr_chan))
        j = 0
        for i in range(nr_chan):
            fid.seek(i * nr_lines * 2, 0)
            data[:, j] = np.fromfile(fid, 'int16', nr_lines) * scale_factors[i]
            j += 1
    return data


cwd = os.getcwd()
cwd = cwd.replace('\\','/')
folder = cwd + '/res/'

itime = 0  # column index of time
irotspd = 9  # column index of rotor speed

# -------- ascii --------
case = 'C2'
path = folder + 'dtu_10mw_rwt_' + case + '.dat'
#%%
data = load_hawc_ascii(path)

# -------- hawc binary --------
# name = 'hawc_binary'
# dat_path = folder + 'dtu_10mw_rwt_' + name + '.dat'
# data = load_hawc_binary(dat_path)


# ====================== plot things for fun ======================
time = data[:, itime]
rotspd = data[:, irotspd]
fig, ax = plt.subplots(num=1, clear=True, figsize=(7, 3))
ax.plot(time, rotspd)
ax.set(xlabel='Time [s]', ylabel='Rotor speed [rad/s]')
fig.tight_layout()