# ============================================================================ #
# main.py
#
# Jonah Lee
#
# Entry point for compilation of .g3 data
# ============================================================================ #
from setuptools.command.bdist_egg import scan_module

from data_loader import config
from data_loader.config import roach_ids
from data_loader.roach import RoachPass, RoachID, ScanPass

from spt3g import core
import so3g
import numpy as np
import os


class BlastData:
    """Contains BLAST-TNG in the form of references RoachPass objects.

    Provides convenience methods for accessing the data contained in the RoachPass objects,
    for the purpose of G3 file packaging
    """

    def __init__(self, roach_ids: tuple[int] = None):
        self.roach_ids: tuple[int] = roach_ids if roach_ids is not None else (1, 2, 3, 4, 5)
        self.roaches: dict = self._load_roaches()

    def _load_roaches(self):
        """Loads a RoachPass objet for each Roach ID in self.roach_ids

        Loads all 3 passes by default
        """
        roaches = {roach_id: RoachPass(RoachID(roach_id), ScanPass.ALL, use_rejects_file=False)
                   for roach_id in self.roach_ids}

        return roaches

    def get_time(self, roach_id):
        return self.roaches[roach_id].dat_sliced['time']

    def get_azimuth(self, roach_id):
        return self.roaches[roach_id].dat_sliced['az']

    def get_elevation(self, roach_id):
        return self.roaches[roach_id].dat_sliced['el']

    def get_latitude(self, roach_id):
        return self.roaches[roach_id].dat_sliced['lat']

    def get_longitude(self, roach_id):
        return self.roaches[roach_id].dat_sliced['lon']

    def get_altitude(self, roach_id):
        return self.roaches[roach_id].dat_sliced['alt']

    def get_kid_i_q(self, roach_id, kid) -> tuple[np.ndarray, np.ndarray]:
        return self.roaches[roach_id].get_kid_i_q(kid)

class ScanFrameGenerator:
    def __init__(self, data: BlastData, ref_roach_id: int, scan_seconds: float=3, data_freq: float=476.5):
        self.data: BlastData = data
        self.ref_roach_id: id = ref_roach_id
        self.scan_idx: int = 0
        self.scan_seconds: float = scan_seconds
        self.data_freq: float = data_freq
        self.scan_len: int = int(self.scan_seconds * self.data_freq)
        self.done: bool = False

    def _get_scan_slice(self):
        """Determine the indices which define the start (inclusive) and stop (exclusive) for the current scan frame

        The resulting indices can be used equally in dat_sliced fields or in kid_i_q timestreams.
        """
        slice_i = self.scan_idx * self.scan_len
        slice_f = (self.scan_idx + 1) * self.scan_len - 1
        return slice_i, slice_f

    def __call__(self, frame):
        if self.done:
            return []

        out_frame = core.G3Frame(core.G3FrameType.Scan)

        start_i, stop_i = self._get_scan_slice()

        ts = so3g.G3SuperTimestream()
        ts.names = names
        ts.times = core.G3VectorTime(times * core.G3Units.s)
        ts.data = data

        return out_frame


if __name__ == '__main__':

    out_dir = os.path.join(config.g3_dir, config.version_dir)
    os.makedirs(out_dir, exist_ok=True)

    data = BlastData()
    generator = ScanFrameGenerator(data)

    pipe = core.G3Pipeline()

    pipe.Add(generator)
    pipe.Add(core.G3Writer, filename=os.path.join(out_dir, 'testfile.g3'))

    pipe.Run()

    breakpoint()
