# ============================================================================ #
# config.py
#
# Jonah Lee
#
# Configuration file for loading and exporting of BLAST-TNG .g3 data.
# ============================================================================ #

g3_dir = ""

# base data directories
dir_root   = '/media/player1/blast2020fc1/fc1/'   # control computer
dir_conv   = dir_root + "converted/"              # control computer

# data directories and files
dir_master = dir_conv + 'master_2020-01-06-06-21-22/'

dir_targ_dict = {
    1: dir_root + 'roach_flight/roach1/targ/Tue_Jan__7_00_55_50_2020/',
    2: dir_root + 'roach_flight/roach2/targ/Tue_Jan__7_00_55_50_2020/',
    3: dir_root + 'roach_flight/roach3/targ/Tue_Jan__7_00_55_51_2020/',
    4: dir_root + 'roach_flight/roach4/targ/Tue_Jan__7_00_55_50_2020/',
    5: dir_root + 'roach_flight/roach5/targ/Tue_Jan__7_00_55_51_2020/',
}

dir_roach_dict = {
    1: dir_conv +'roach1_2020-01-06-06-22-01/',
    2: dir_conv +'roach2_2020-01-06-06-22-01/',
    3: dir_conv +'roach3_2020-01-06-06-21-56/',
    4: dir_conv +'roach4_2020-01-06-06-22-01/',
    5: dir_conv +'roach5_2020-01-06-06-22-01/',
}

# KID rejects list
file_rejects_dict = {
    1: dir_root + f'map_making/kid_rejects/kid_rejects_roach1.dat',
    2: dir_root + f'map_making/kid_rejects/kid_rejects_roach2.dat',
    3: dir_root + f'map_making/kid_rejects/kid_rejects_roach3.dat',
    4: dir_root + f'map_making/kid_rejects/kid_rejects_roach4.dat',
    5: dir_root + f'map_making/kid_rejects/kid_rejects_roach5.dat',
}

# common-mode file
file_commonmode = lambda roach: f'common_mode_roach{roach}.dat'

# log file
log_file = 'map_making.log'

# single KID maps output directory
dir_single = 'single_maps/'

# map aligning parameters output dir and file
dir_xform = 'align/'
# file_xform = dir_xform + f'align_roach{roach}.npy'
# file_source_coords = dir_xform + f'source_coords_roach{roach}.npy'