# ============================================================================ #
# main.py
#
# Jonah Lee
#
# Entry point for compilation of .g3 data
# ============================================================================ #

from frame_generators import BlastData, FrameGenManager, ScanFrameGenerator, CalFrameGenerator, ObservationFrameGenerator
from data_loader.roach import ScanPass
from data_loader import config
import os

from so3g import proj
import numpy as np

from spt3g import core
from spt3g import calibration

def add_bolometer_properties(frame, x_shifts="x_shifts", y_shifts="y_shifts", out_key="BolometerProperties"):
    if frame.type != core.G3FrameType.Calibration:
        return
    kids = sorted(frame[x_shifts].keys())
    bp_map_list = []
    for kid in kids:
        bp = calibration.BolometerProperties()
        bp.physical_name = kid
        bp.x_offset = frame[x_shifts][kid]
        bp.y_offset = frame[y_shifts][kid]
        bp_map_list.append((kid, bp))
    frame[out_key] = calibration.BolometerPropertiesMap(bp_map_list)

class GetFocalPlane:
    def __init__(self, x_shifts="x_shifts", y_shifts="y_shifts", out_key="BolometerProperties"):
        self.x_shifts = x_shifts
        self.y_shifts = y_shifts
        self.out_key = out_key
        self.focalplane = None

    def __call__(self, frame):
        if frame.type != core.G3FrameType.Calibration:
            return
        idxs = np.argsort(frame[self.x_shifts].keys())
        frame[self.x_shifts].values())
        proj.FocalPlane.from_xieta()
        frame[self.out_key] = calibration.BolometerPropertiesMap(bp_map_list)

if __name__ == '__main__':

    out_dir = os.path.join(config.g3_dir, config.version_dir)
    os.makedirs(out_dir, exist_ok=True)

    scan_pass = ScanPass.PASS_3

    # at the moment, the program runs out of memory with all 5 roaches
    data = BlastData(roach_ids=(1,), scan_pass=scan_pass)

    obs_generator = ObservationFrameGenerator(scan_pass=scan_pass)
    cal_generator = CalFrameGenerator(data, 1)
    scan_generator = ScanFrameGenerator(data, 1)

    # create a FrameGenManager which manages frame insertion between generator modules
    frame_gen_manager = FrameGenManager([obs_generator, cal_generator, scan_generator])

    pipe = core.G3Pipeline()
    pipe.Add(frame_gen_manager)
    pipe.Add(add_bolometer_properties)
    pipe.Add(core.G3Writer, filename=os.path.join(out_dir, 'roach1_pass3.g3'))

    pipe.Run()
