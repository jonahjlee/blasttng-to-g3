# ============================================================================ #
# main.py
#
# Jonah Lee
#
# Entry point for compilation of .g3 data
# ============================================================================ #

import config
import mmi_data_lib as dlib

from spt3g import core
import numpy as np
import so3g
import os


if __name__ == '__main__':

    class ScanFrameGenerator(core.G3Module):
        def __init__(self, roaches=None):
            ...

        def Process(self, frame):
            ...

    out_dir = os.path.join(config.g3_dir, config.version_dir)
    os.makedirs(out_dir, exist_ok=True)

    pipe = core.G3Pipeline()
    pipe.Add(ScanFrameGenerator)
    pipe.Add(core.G3Writer, filename=os.path.join(out_dir, 'testing.g3'))
    pipe.Run()
    ...