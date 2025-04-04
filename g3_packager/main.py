# ============================================================================ #
# main.py
#
# Jonah Lee
#
# Entry point for compilation of .g3 data
# ============================================================================ #

from frame_generators import BlastData, FrameGenManager, ScanFrameGenerator, CalFrameGenerator
from data_loader.roach import ScanPass
from data_loader import config
import os

from spt3g import core

if __name__ == '__main__':

    out_dir = os.path.join(config.g3_dir, config.version_dir)
    os.makedirs(out_dir, exist_ok=True)

    # at the moment, the program runs out of memory with all 5 roaches
    data = BlastData(roach_ids=(1,), scan_pass=ScanPass.PASS_3)
    scan_generator = ScanFrameGenerator(data, 1)
    cal_generator = CalFrameGenerator(data)
    frame_gen_manager = FrameGenManager(data, [cal_generator, scan_generator])

    pipe = core.G3Pipeline()
    pipe.Add(frame_gen_manager)
    pipe.Add(core.G3Writer, filename=os.path.join(out_dir, 'roach1_pass3.g3'))

    pipe.Run()
