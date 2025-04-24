import os

from data_loader.roach import RoachPass, RoachID, ScanPass
from data_loader import data_lib as dlib
from data_loader import config
# so3g must be imported before spt3g to avoid segmentation fault
import so3g
from spt3g import core
import numpy as np
from astropy.wcs import WCS


class BlastData:
    """Contains BLAST-TNG in the form of references RoachPass objects.

    Provides convenience methods for accessing the data contained in the RoachPass objects,
    for the purpose of G3 file packaging
    """

    def __init__(self, roach_ids: tuple[int]=None, scan_pass: ScanPass=None):
        self.roach_ids: tuple[int] = roach_ids if roach_ids is not None else (1, 2, 3, 4, 5)
        self.scan_pass = scan_pass if scan_pass is not None else ScanPass.ALL
        self.roaches: dict[int, RoachPass] = self._load_roaches()

    def _load_roaches(self):
        """Loads a RoachPass object for each Roach ID in self.roach_ids

        Loads all 3 passes by default
        """
        roaches = {roach_id: RoachPass(RoachID(roach_id), self.scan_pass, use_rejects_file=False)
                   for roach_id in self.roach_ids}
        return roaches

    def get_time(self, roach_id) -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['time'] * core.G3Units.s

    def get_azimuth(self, roach_id) -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['az'] * core.G3Units.deg

    def get_elevation(self, roach_id) -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['el'] * core.G3Units.deg

    def get_latitude(self, roach_id) -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['lat'] * core.G3Units.deg

    def get_longitude(self, roach_id) -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['lon'] * core.G3Units.deg

    def get_altitude(self, roach_id) -> np.ndarray:
        return self.roaches[roach_id].dat_sliced['alt'] * core.G3Units.m

    def get_kid_i(self, roach_id, kid) -> np.ndarray:
        return self.roaches[roach_id].get_kid_i(kid)

    def get_kid_q(self, roach_id, kid) -> np.ndarray:
        return self.roaches[roach_id].get_kid_q(kid)

    def get_kid_target_sweeps(self, roach_id, kid) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Load target sweeps If/Qf/Ff data for KID

        see mmi_data_lib loadTargSweepsData/getTargSweepIQ for more"""
        dat_targs = self.roaches[roach_id].dat_targs
        Ff = self.roaches[roach_id].Ff
        kid_targ_i = dat_targs[:, 2 * int(kid)]
        kid_targ_q = dat_targs[:, 2 * int(kid) + 1]
        return kid_targ_i, kid_targ_q, Ff


class FrameGenManager(core.G3Module):
    """Contains a generator module for each type of frame produced and manages their frame output.

    FrameGenManager is intended to be the first module, so it will be passed an unlimited
    number of None frames until it returns [] and ends processing.
    """

    def __init__(self, generators=None):
        super().__init__()
        self.generators: list[core.G3Module] = generators if generators is not None else []
        self.generator_idx = 0
        self.done = False

    def Process(self, frame):
        if self.generator_idx >= len(self.generators):
            # end processing; all generators are done
            return []

        active_generator = self.generators[self.generator_idx]
        print(f"Adding a frame. {active_generator=}")
        out = active_generator.Process(frame)

        if active_generator.done:
            self.generator_idx += 1

        return out


class ObservationFrameGenerator(core.G3Module):
    def __init__(self, scan_pass: ScanPass=None):
        self.obs_num = -1 if scan_pass is None else scan_pass.value
        self.done = False
    def Process(self, _):
        assert not self.done, "Generator should not be called when it is done!"

        out_frame = core.G3Frame(core.G3FrameType.Observation)
        out_frame['SourceName'] = "RCW92"
        out_frame['ObservationNumber'] = self.obs_num

        self.done = True
        return out_frame


class ScanFrameGenerator(core.G3Module):
    def __init__(self, data: BlastData, ref_roach_id: int,
                 scan_seconds: float=3, data_freq: float=476.5, frame_limit: int=None):
        super().__init__()
        self.data: BlastData = data
        self.ref_roach_id: int = ref_roach_id
        self.frame_num: int = 0
        self.scan_seconds: float = scan_seconds
        self.data_freq: float = data_freq
        self.scan_len: int = int(self.scan_seconds * self.data_freq)
        self.max_frame_num: int = self._get_max_frame_num()
        self.frame_limit: int | None = frame_limit
        self.done: bool = False

    def _get_max_frame_num(self):
        """Determines the maximum ``frame_num`` given the length of the scan and the length of the data

        ``frame_num`` is zero for the first scan frame and increments for each frame.
        The last frame, with ``frame_num == max_frame_num`` may have fewer items than ``scan_len``,
        and will use up the remaining data in the ref_roach's main axis.
        """
        num_indices = len(self.data.roaches[self.ref_roach_id])
        # truncate towards zero since max_frame_num is one less than the total number of frames.
        # for example, if there are 4.6 frames worth of data, we would have 5 frames
        # (the fifth being a bit shorter) and max_frame_num would be 5
        max_frame_num = int(num_indices / self.scan_len)
        return max_frame_num

    def _get_scan_slice(self):
        """Determine the indices which define the start (inclusive) and stop (exclusive) for the current scan frame

        The resulting indices can be used equally in dat_sliced fields or in kid_i_q timestreams.
        """
        slice_i = self.frame_num * self.scan_len
        slice_f = (self.frame_num + 1) * self.scan_len - 1 if self.frame_num < self.max_frame_num else -1
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

    def Process(self, _):
        assert not self.done, "Generator should not be called when it is done!"

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
                ts = core.G3Timestream(func(self.data, roach_id)[slice_i:slice_f])
                ts.start = t_i
                ts.stop = t_f
                out_frame[key] = ts

        print(f"Generated Scan Frame {self.frame_num}/{self.max_frame_num}")

        if ((self.frame_limit is not None and self.frame_num >= self.frame_limit)
                or (self.frame_num >= self.max_frame_num)):
            # stop when we reach the manually imposed frame_limit or run out of data, whichever comes first
            self.done = True
        self.frame_num += 1

        return out_frame


class CalFrameGenerator(core.G3Module):
    """Generate Calibration Frame(s) for BLAST data.

    This means dat_targs and Ff, since these are the inputs required
    to compute DF given kid I/Q data.
    """
    def __init__(self, data: BlastData, ref_roach_id: int):
        super().__init__()
        self.data: BlastData = data
        self.done = False
        self.ref_roach_id: int = ref_roach_id

    def get_target_sweeps(self):
        iq_dict: dict[str, core.G3Timestream] = {}
        for roach_id in self.data.roach_ids:
            kids = self.data.roaches[roach_id].kids
            for kid in kids:
                If, Qf, Ff = self.data.get_kid_target_sweeps(roach_id, kid)
                iq_dict[f"roach{roach_id}_{kid}_I"] = core.G3Timestream(If.astype(np.float64))
                iq_dict[f"roach{roach_id}_{kid}_Q"] = core.G3Timestream(Qf.astype(np.float64))
                iq_dict[f"roach{roach_id}_{kid}_F"] = core.G3Timestream(Ff.astype(np.float64))
        # initialize the timestream map
        # start/stop are not set because they are not used for DF calculation
        ts = core.G3TimestreamMap(iq_dict)
        return ts

    @staticmethod
    def um_to_az_el_offsets(x_um, y_um, platescale):
        """Converts micron offsets on image plane to az/el on-sky offsets

        Returns values as angular G3 Units
        """
        w = WCS(naxis=2)
        w.wcs.crpix = [0, 0]  # center of the focal plane is tangent point
        w.wcs.crval = [0., 0.]  # source is at center in offsets map
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.cdelt = [platescale, platescale]

        x_deg, y_deg = w.wcs_pix2world(x_um, y_um, 0)

        # Normalize the angles to the range [-180, 180)
        x_deg_normalized = np.mod(x_deg + 180, 360) - 180
        y_deg_normalized = np.mod(y_deg + 180, 360) - 180

        return x_deg_normalized * core.G3Units.deg, y_deg_normalized * core.G3Units.deg

    def get_kid_shifts(self) -> tuple[core.G3MapDouble, core.G3MapDouble]:
        """Return source-determined KID ra_shifts and dec_shifts.

        Currently only supports roach 1
        """
        # note: file currently has only roach1
        shifts_file = os.path.join(os.path.dirname(__file__), "detector_layouts", f"roach1_shifts_radec.npy")
        kid_shifts: dict[str, tuple[float, float]] = np.load(shifts_file, allow_pickle=True).item()
        x_keyvals = [(k, v[0]) for k, v in kid_shifts.items()]
        y_keyvals = [(k, v[1]) for k, v in kid_shifts.items()]
        return core.G3MapDouble(x_keyvals), core.G3MapDouble(y_keyvals)

    def _get_cal_lamp_kid_data(self) -> so3g.G3SuperTimestream:
        kid_i_q_data = None
        kid_i_q_names = None
        for id, roach in self.data.roaches.items():
            roach_i_names = [f'roach{id}_{kid}_I' for kid in roach.kids]
            roach_i = [roach.get_kid_cal_lamp_i(kid) for kid in roach.kids]
            roach_q_names = [f'roach{id}_{kid}_Q' for kid in roach.kids]
            roach_q = [roach.get_kid_cal_lamp_q(kid) for kid in roach.kids]
            kid_i_q_data = np.array(roach_i + roach_q)
            kid_i_q_names = roach_i_names + roach_q_names
        # see https://so3g.readthedocs.io/en/latest/cpp_objects.html#how-to-work-with-float-arrays
        times = core.G3VectorTime(np.zeros(config.cal_f_offset - config.cal_i_offset))  # real times not needed, use placeholder for now
        quanta = 0.01 * np.ones(len(kid_i_q_names))
        ts = so3g.G3SuperTimestream(kid_i_q_names, times, kid_i_q_data, quanta)
        return ts

    def Process(self, _):
        assert not self.done, "Generator should not be called when it is done!"

        out_frame = core.G3Frame(core.G3FrameType.Calibration)

        out_frame["target_sweeps"] = self.get_target_sweeps()
        ra_shifts, dec_shifts = self.get_kid_shifts()
        out_frame["ra_shifts"] = ra_shifts
        out_frame["dec_shifts"] = dec_shifts
        out_frame["cal_lamp_data"] = self._get_cal_lamp_kid_data()

        self.done = True
        return out_frame  # insert the calframe into the pipeline
