# ============================================================================ #
# main.py
#
# Jonah Lee
#
# Entry point for compilation of .g3 data
# ============================================================================ #

from data_loader import config
from data_loader.config import roach_ids
from data_loader.roach import RoachPass, RoachID, ScanPass

# so3g must be imported before spt3g to avoid segmentation fault
import so3g
from spt3g import core
import numpy as np
import os


class BlastData:
    """Contains BLAST-TNG in the form of references RoachPass objects.

    Provides convenience methods for accessing the data contained in the RoachPass objects,
    for the purpose of G3 file packaging
    """

    def __init__(self, roach_ids: tuple[int] = None):
        self.roach_ids: tuple[int] = roach_ids if roach_ids is not None else (1, 2, 3, 4, 5)
        self.roaches: dict[int, RoachPass] = self._load_roaches()

    def _load_roaches(self):
        """Loads a RoachPass objet for each Roach ID in self.roach_ids

        Loads all 3 passes by default
        """
        roaches = {roach_id: RoachPass(RoachID(roach_id), ScanPass.ALL, use_rejects_file=False)
                   for roach_id in self.roach_ids}

        return roaches

    def get_time(self, roach_id) -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['time'] * core.G3Units.s

    def get_azimuth(self, roach_id) -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['az'] * core.G3Units.deg

    def get_elevation(self, roach_id) -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['el'] * core.G3Units.deg

    def get_latitude(self, roach_id)  -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['lat'] * core.G3Units.deg

    def get_longitude(self, roach_id)  -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['lon'] * core.G3Units.deg

    def get_altitude(self, roach_id)  -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['alt'] * core.G3Units.m

    def get_kid_i(self, roach_id, kid) -> np.ndarray:
        return self.roaches[roach_id].get_kid_i(kid)

    def get_kid_q(self, roach_id, kid) -> np.ndarray:
        return self.roaches[roach_id].get_kid_q(kid)

class ScanFrameGenerator:
    def __init__(self, data: BlastData, ref_roach_id: int,
                 scan_seconds: float=3, data_freq: float=476.5, max_scans: int=None):
        self.data: BlastData = data
        self.ref_roach_id: id = ref_roach_id
        self.scan_idx: int = 0
        self.scan_seconds: float = scan_seconds
        self.data_freq: float = data_freq
        self.scan_len: int = int(self.scan_seconds * self.data_freq)
        self.max_scans: int = max_scans
        self.done: bool = False

    def _get_scan_slice(self):
        """Determine the indices which define the start (inclusive) and stop (exclusive) for the current scan frame

        The resulting indices can be used equally in dat_sliced fields or in kid_i_q timestreams.
        """
        slice_i = self.scan_idx * self.scan_len
        slice_f = (self.scan_idx + 1) * self.scan_len - 1
        return slice_i, slice_f

    def _get_kid_data(self, start_i, stop_i) -> so3g.G3SuperTimestream:
        times = core.G3VectorTime(self.data.get_time(self.ref_roach_id)[start_i:stop_i])
        kid_i_q_data = None
        kid_i_q_names = None
        for id, roach in self.data.roaches.items():
            roach_i_names = [f'roach{id}_{kid}_I' for kid in roach.kids]
            roach_i = [roach.get_kid_i(kid)[start_i:stop_i] for kid in roach.kids]
            roach_q_names = [f'roach{id}_{kid}_Q' for kid in roach.kids]
            roach_q = [roach.get_kid_q(kid)[start_i:stop_i] for kid in roach.kids]
            kid_i_q_data = np.array(roach_i + roach_q)
            kid_i_q_names = roach_i_names + roach_q_names
        # see https://so3g.readthedocs.io/en/latest/cpp_objects.html#how-to-work-with-float-arrays
        quanta = 0.01 * np.ones(len(kid_i_q_names))
        ts = so3g.G3SuperTimestream(kid_i_q_names, times, kid_i_q_data, quanta)
        return ts

    def __call__(self, frame):
        if self.done:
            return []

        out_frame = core.G3Frame(core.G3FrameType.Scan)

        slice_i, slice_f = self._get_scan_slice()

        out_frame['data'] = self._get_kid_data(slice_i, slice_f)
        out_frame['time'] = core.G3Time.Now()

        t_i = out_frame['data'].times[0]
        t_f = out_frame['data'].times[-1]

        data_funcs = {
            "az": BlastData.get_azimuth,
            "el": BlastData.get_elevation,
            "lat": BlastData.get_latitude,
            "lon": BlastData.get_longitude,
            "alt": BlastData.get_altitude,
        }

        for key, func in data_funcs.items():
            for roach_id in self.data.roach_ids:
                ts = core.G3Timestream(func(data, roach_id)[slice_i:slice_f])
                ts.start = t_i
                ts.stop = t_f
                out_frame[key] = ts


        self.scan_idx += 1
        if self.max_scans is not None and self.scan_idx >= self.max_scans:
            self.done = True

        breakpoint()

        return out_frame


if __name__ == '__main__':

    out_dir = os.path.join(config.g3_dir, config.version_dir)
    os.makedirs(out_dir, exist_ok=True)

    # at the moment, the program runs out of memory with all 5 roaches
    data = BlastData(roach_ids=(1,))
    generator = ScanFrameGenerator(data, 1, max_scans=10)

    pipe = core.G3Pipeline()

    pipe.Add(generator)
    pipe.Add(core.G3Writer, filename=os.path.join(out_dir, 'testfile.g3'))

    pipe.Run()

    breakpoint()
