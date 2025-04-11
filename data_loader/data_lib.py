# ============================================================================ #
# data_lib.py
#
# James Burgoyne jburgoyne@phas.ubc.ca
# Modified by Jonah Lee
# CCAT Prime 2024
#
# Map Maker Iterative data library.
# Copied & modified from https://github.com/ccat850/map_making_blasttng/
# ============================================================================ #

import os
import re
import gc
import numpy as np
import pandas as pd


# ============================================================================ #
def loadMasterData(roach, dir_master, dir_roach):
    '''Loads all the common data into memory.
    '''

    # load master fields
    # master_fields = ['time', 'time_usec', 'ra', 'dec', 'el', 'az', 'alt', 'lon', 'lat']
    master_fields = ['time', 'time_usec', 'el', 'az', 'alt', 'lon', 'lat']
    dat_raw = {
        field: np.load(dir_master + field + '.npy')
        for field in master_fields}

    # load common roach fields
    dat_raw['roach_time'] = np.load(dir_roach + f'ctime_built_roach{roach}.npy', mmap_mode='r')

    # combine time_usec into time
    dat_raw['time'] = dat_raw['time'].astype('float64') + dat_raw['time_usec']*1e-6
    del(dat_raw['time_usec'])

    # coord conversions (as per format file)
    # for field in ['ra', 'dec', 'el', 'az']:
    for field in ['el', 'az']:
        dat_raw[field]  = dat_raw[field].astype('float64')*8.38190317e-08
    for field in ['lon', 'lat']:
        dat_raw[field]  = dat_raw[field].astype('float64')*1.676380634e-07

    return dat_raw


# ============================================================================ #
def loadTargSweepsData(dir_targ):
    '''Loads and combines the target sweep files.

    dir_targ: (str or path) The absolute filename str or path.
    '''

    # load and combine targ files (If, Qf)
    pattern = r'^\d{9}\.dat$'
    files = os.listdir(dir_targ)
    matched_files = [f for f in files if re.match(pattern, f)]
    sorted_files = sorted(matched_files)
    dat_targs = np.array([
        np.fromfile(os.path.join(dir_targ, f), dtype = '<f')
        for f in sorted_files
    ])

    # load frequency file (Ff)
    Ff = np.loadtxt(dir_targ + 'sweep_freqs.dat')

    return dat_targs, Ff


# ============================================================================ #
def alignMasterAndRoachTods(dat_raw):
    '''Interpolate master arrays and align roach.
    '''

    # master_fields = ['time', 'ra', 'dec', 'el', 'az', 'alt', 'lon', 'lat']
    master_fields = ['time', 'el', 'az', 'alt', 'lon', 'lat']
    roach_time = dat_raw['roach_time']

    # interpolate master fields to roach_time length
    def interp(a, a_match):
        x_old = np.arange(len(a))
        x_new = np.linspace(0, len(a)-1, len(a_match))
        a_interp = np.interp(x_new, x_old, a)
        return a_interp

    dat_aligned = {}
    for field in master_fields:
        dat_aligned[field] = interp(dat_raw[field], roach_time)
        del(dat_raw[field])
        gc.collect()

    # dat_aligned = {
    #     field: interp(dat_raw[field], roach_time)
    #     for field in master_fields}

    # indices to align roach tods to master tods
    indices = np.searchsorted(roach_time, dat_aligned['time'], side='left')
    indices[indices == len(roach_time)] = len(roach_time) - 1 # max index bug fix

    # use aligned roach time as main time tod
    dat_aligned['time'] = roach_time[indices]

    return dat_aligned, indices


# ============================================================================ #
# samplingFrequency
def samplingFrequency(tod_time):
    '''Calculate fs assuming constant.'''

    # print(tod_time[:100])
    # print(np.diff(tod_time))

    # dt = tod_time[1] - tod_time[0]
    # dt = np.nanmean(np.diff(tod_time))

    i = 0
    dt = 0
    while dt == 0:
        dt = tod_time[i + 1] - tod_time[i]
        i += 1

    return 1/dt


# ============================================================================ #
# ds
def ds(X, Y):
    '''Spatial bin diff.'''

    i = 0
    ds = 0
    while ds == 0:
        ds = np.sqrt((X[i + 1] - X[i])**2 + (Y[i + 1] - Y[i])**2)
        i += 1

    return ds


# ============================================================================ #
def findAllKIDs(directory):
    '''Search given directory for KID files and return set of KID numbers.
    Note that the KID numbers are strings with leading zero, e.g. '0100'.
    '''

    files = os.listdir(directory)

    # Extract 'kid' values from filenames
    kid_values = set()
    for file in files:
        # Check if the file matches the expected format
        if file.startswith('i_kid') and file.endswith('.npy'):
            # Extract 'kid' value from the filename
            kid = int(file.split('_')[1][4:])
            kid_values.add(kid)

    # Sort 'kid' values and format them with leading zeros
    sorted_kid_values = sorted(kid_values)
    sorted_kid_values_strings = [f"{kid:04}" for kid in sorted_kid_values]

    return sorted_kid_values_strings


# ============================================================================ #
def loadKidRejects(file_rejects):
    '''
    '''

    # load rejects file
    dat = np.loadtxt(file_rejects, delimiter=' ', dtype=str)

    return dat


# ============================================================================ #
def loadKIDI(roach, kid, dir_roach) -> np.memmap:
    '''Preps KID I (in-phase) for on-demand loading.
    '''

    I: np.memmap = np.load(dir_roach + f'i_kid{kid}_roach{roach}.npy',
                allow_pickle=False, mmap_mode='r')

    return I

# ============================================================================ #
# ============================================================================ #
def loadKIDQ(roach, kid, dir_roach) -> np.memmap:
    '''Preps KID Q (quadrature) for on-demand loading.
    '''

    Q: np.memmap = np.load(dir_roach + f'q_kid{kid}_roach{roach}.npy',
                allow_pickle=False, mmap_mode='r')

    return Q


# ============================================================================ #
def loadKIDData(roach, kid, dir_roach) -> tuple[np.memmap, np.memmap]:
    '''Preps KID I and Q for on-demand loading.
    '''

    return loadKIDI(roach, kid, dir_roach), loadKIDQ(roach, kid, dir_roach)


# ============================================================================ #
def getTargSweepIQ(kid, dat_targs):
    '''Filter roach targs for target sweep for this kid.

    kid: (int) The KID number.
    dat_targs: (2D array; floats) Array of all target sweeps for this roach.
    '''

    I = dat_targs[:, 2*int(kid)]
    Q = dat_targs[:, 2*int(kid)+1]

    return I, Q


# ============================================================================ #
def abFromLayout(file_layout):
    '''Get the a,b coords from the detector layout file.
    file_layout: (str) Absolute file name of detector layout file.
    ab coordinates are micron offsets from geometric center of array.
    Use platescale to convert to on-sky coords.
    '''

    # Load layout file CSV
    data = np.loadtxt(file_layout, skiprows=1, delimiter=',')

    # ab = {kid_i: (a_i, b_i)}
    # ab = {
    #     f"{int(row[0]):04}": (row[1], row[2])
    #     for row in sorted(data, key=lambda x: int(x[0]))
    # }

    # prep the fields
    kids = [f"{kid:04}" for kid in sorted(data[:,0].astype(int))]
    a = data[:,1].astype(float)
    b = data[:,2].astype(float)
    # convert to dict:
    ab = {kid: (a[i],b[i]) for i,kid in enumerate(kids)}

    return ab


# ============================================================================ #
def load_kid_layout(file, rejects_file=None) -> dict[str, tuple[float, float]]:
    if file.endswith('.csv'):
        return load_kid_layout_csv(file, rejects_file)
    if file.endswith('.npy'):
        return load_kid_layout_npy(file, rejects_file)

    raise ValueError('Layout file must end with .csv or .npy')


# ============================================================================ #
def load_kid_layout_npy(file, rejects_file=None) -> dict[str, tuple[float, float]]:
    """Loads KID x/y coordinates on the image plane for a ROACH

    Returns a dictionary which maps KID IDs to coordinate pairs.
    """
    try:
        layouts_dict = np.load(file, allow_pickle=True).item()

        # parse keys in the form '0000' or 'roach1_0000'
        if isinstance(next(iter(layouts_dict.keys())), str):
            layouts_dict = {key[-4:]: val for key, val in layouts_dict.items()}
        # parse keys in the form of ints
        if isinstance(next(iter(layouts_dict.keys())), int):
            layouts_dict = {f'{key:04}': val for key, val in layouts_dict.items()}

        if rejects_file is not None:
            try:
                rejects = np.loadtxt(rejects_file, delimiter=' ', dtype=str)
                if isinstance(rejects[0], str): rejects = [key[-4:] for key in rejects]
                if isinstance(rejects[0], int): rejects = [f'{key:04}' for key in rejects]
                layouts_dict = {key: val for key, val in layouts_dict.items() if key in rejects}
            except FileNotFoundError as err:
                print(f"File {rejects_file} not found")
                raise err

        return layouts_dict

    except FileNotFoundError as err:
        print(f"File {file} not found")
        raise err


# ============================================================================ #
def load_kid_layout_csv(file, rejects_file=None) -> dict[str, tuple[float, float]]:
    """Loads KID x/y coordinates on the image plane for a ROACH

    Returns a dictionary which maps KID IDs to coordinate pairs.
    """
    try:
        df = pd.read_csv(file, index_col=0)

        if rejects_file is not None:
            try:
                rejects = np.loadtxt(rejects_file, delimiter=' ', dtype=str),
                df = df[~df.index.isin(rejects)]
            except FileNotFoundError as err:
                print(f"File {rejects_file} not found")
                raise err

        return {str(kid): (df['x'][kid], df['y'][kid]) for kid in df.index}

    except FileNotFoundError as err:
        print(f"File {file} not found")
        raise err