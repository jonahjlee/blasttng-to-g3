# ============================================================================ #
# g3utils.py
#
# Jonah Lee
#
# Collection of miscellaneous tools for processing BLAST-TNG G3 files.
# Mostly contains G3 Modules for use in pipelines.
# Includes modules for debugging, signal processing, map-making etc.
# This file will probably be moved into a more organized package.
# ============================================================================ #

import so3g
from spt3g import core
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as au
from spt3g.core import G3Units as gu
from scipy.ndimage import gaussian_filter

def ordinal(n: int):
    """source: https://stackoverflow.com/questions/9647202/ordinal-numbers-replacement/20007730#20007730"""
    if 11 <= (n % 100) <= 13:
        suffix = 'th'
    else:
        suffix = ['th', 'st', 'nd', 'rd', 'th'][min(n % 10, 4)]
    return str(n) + suffix

class FrameCounter(core.G3Module):
    def __init__(self):
        super(FrameCounter, self).__init__()
        self.previous_type = None
        self.num_repeats = 0
    def Process(self, frame):
        type = frame.type
        if type == self.previous_type:
            self.num_repeats += 1
            print(f"{type} (x{self.num_repeats + 1})", end='\r')
        else:
            print()
            print(f"{type}", end='\r')
            if type == core.G3FrameType.EndProcessing: print()
            self.num_repeats = 0
        self.previous_type = type


class NthFrameGrabber:
    """Stores the nth frame of a given type"""
    def __init__(self, n=1, frame_type: core.G3FrameType=None, verbose=False):
        self.n = n
        self.frame_type = frame_type if frame_type is not None else core.G3FrameType.Scan
        self.num_seen = 0
        self.nth_frame = None
        self.verbose = verbose
    def __call__(self, frame):
        if self.nth_frame is not None:
            # already found the frame
            return
        if frame.type == self.frame_type:
            self.num_seen += 1
            if self.num_seen == self.n:
                self.nth_frame = frame
                if not self.verbose: return
                print(f"Found the {ordinal(self.n)} frame with type: {self.frame_type}!")
                print(f"The frame is now stored in {self}'s nth_frame attribute.")


class FirstFrameGrabber(NthFrameGrabber):
    """Stores the first frame of a given type"""
    def __init__(self, frame_type: core.G3FrameType=None, verbose=False):
        super().__init__(1, frame_type, verbose)
    @property
    def first_frame(self):
        return self.nth_frame


class LastFrameGrabber:
    """Grabs the last frame of a given type"""
    def __init__(self, frame_type: core.G3FrameType=None):
        self.frame_type = frame_type if frame_type is not None else core.G3FrameType.Scan
        self.last_frame = None
    def __call__(self, frame):
        if frame.type == self.frame_type:
            self.last_frame = frame
            return
        if frame.type == core.G3FrameType.EndProcessing:
            print(f"Found the last frame with type: {self.frame_type}!")
            print(f"The frame is now stored in {self}'s last_frame attribute.")


def add_radec_astropy(frame,
                      az: str="az", el: str="el",
                      lat: str="lat", lon: str="lon",
                      alt: str="alt",
                      data: str="data",
                      ra: str="ra", dec: str="dec"):
    """Use astropy coordinate transformations to convert az/el --> ra/dec.

    Keyword arguments indicate keys in the scan frame for inputs/outputs
    Requires lat/lon/alt to determine the telescope location as an input to astropy.coordinates.SkyCoord
    """
    if frame.type != core.G3FrameType.Scan:
        return

    az_deg = np.array(frame[az]) / gu.deg
    el_deg = np.array(frame[el]) / gu.deg
    lat_deg = np.array(frame[lat]) / gu.deg
    lon_deg = np.array(frame[lon]) / gu.deg
    alt_m = np.array(frame[alt]) / gu.m

    unix_times = np.array(frame[data].times) / gu.s
    times = Time(unix_times, format="unix")

    blasttng_loc = EarthLocation(lat=lat_deg, lon=lon_deg, height=alt_m)
    sky_coords = SkyCoord(alt=el_deg * au.deg,
                          az=az_deg * au.deg,
                          obstime=times,
                          frame='altaz',
                          location=blasttng_loc)
    t_i = frame[data].times[0]
    t_f = frame[data].times[-1]

    ra_ts = core.G3Timestream(sky_coords.icrs.ra.deg * gu.deg)
    ra_ts.start = t_i
    ra_ts.stop = t_f
    frame[ra] = ra_ts

    dec_ts = core.G3Timestream(sky_coords.icrs.dec.deg * gu.deg)
    dec_ts.start = t_i
    dec_ts.stop = t_f
    frame[dec] = dec_ts

def add_radec_so3g(frame,
                   az: str="az", el: str="el",
                   lat: str="lat", lon: str="lon",
                   alt: str="alt",
                   data: str="data",
                   ra: str="ra", dec: str="dec"):
    if frame.type != core.G3FrameType.Scan:
        return

    # Approx location for this time frame
    middle_idx = len(frame[lon]) // 2
    site = so3g.proj.EarthlySite(
        frame[lon][middle_idx] / gu.deg,
        frame[lat][middle_idx] / gu.deg,
        frame[alt][middle_idx] / gu.m
    )

    # coords() returns an array with shape (n_time, 4); each 4-tuple contains values (lon, lat, cos(gamma), sin(gamma))
    # Construct a CelestialSightLine to az_el to on-sky coords
    times = np.array(frame[data].times) / gu.s
    csl = so3g.proj.CelestialSightLine.naive_az_el(
        times,
        frame[az] / gu.rad,
        frame[el] / gu.rad,
        site=site
    )
    coords = csl.coords(so3g.proj.FocalPlane.boresight())

    x = np.mod(coords[0][:, 0], 2 * np.pi) * gu.rad
    y = coords[0][:, 1] * gu.rad

    frame[ra] = core.G3VectorDouble(x)
    frame[dec] = core.G3VectorDouble(y)


    def plot(self, ax=None):
        with np.errstate(invalid='ignore'):
            m = self.data / self.hits
        if ax is not None:
            ax.imshow(m, origin='lower')
            ax.set_xticks(range(self.nx + 1)[::10], [f"{ra:.2f}" for ra in self.ra_edges[::10] / core.G3Units.deg],
                          rotation=45)
            ax.set_yticks(range(self.ny + 1)[::10], [f"{dec:.2f}" for dec in self.dec_edges[::10] / core.G3Units.deg])
            ax.set_xlabel("RA (deg)")
            ax.set_ylabel("DEC (deg)")
            ax.set_title("Combined Map")
        else:
            plt.imshow(m, origin='lower')
            plt.xticks(range(self.nx + 1)[::10], [f"{ra:.2f}" for ra in self.ra_edges[::10] / core.G3Units.deg],
                       rotation=45)
            plt.yticks(range(self.ny + 1)[::10], [f"{dec:.2f}" for dec in self.dec_edges[::10] / core.G3Units.deg])
            plt.colorbar(label="DF")
            plt.xlabel("RA (deg)")
            plt.ylabel("DEC (deg)")
            plt.title("Combined Map")
            plt.show()

units_srings: dict[float, str] = {
    gu.deg: "deg",
    gu.degree: "degree",
    gu.degrees: "degrees",
    gu.arcmin: "arcmin",
    gu.arcsec: "arcsec",
    gu.rahour: "rahour",
    gu.raminute: "raminute",
    gu.rasecond: "rasecond",
    gu.rahr: "rahr",
}
class PlotRaDec:
    def __init__(self, ra_key="ra", dec_key="dec", ax=None, units=None):
        self.ra_key = ra_key
        self.dec_key = dec_key
        self.ax = ax
        self.units = units if units is not None else gu.deg
        self.units_str = units_srings[self.units]

        self.ra_data: list[np.ndarray] = []
        self.dec_data: list[np.ndarray] = []
    def __call__(self, frame):
        if frame.type == core.G3FrameType.Scan and self.ra_key in frame and self.dec_key in frame:
            self.ra_data.append(np.asarray(frame[self.ra_key]))
            self.dec_data.append(np.asarray(frame[self.dec_key]))
        if frame.type == core.G3FrameType.EndProcessing:
            ra: np.ndarray = np.concatenate(self.ra_data, axis=0).flatten()
            dec: np.ndarray = np.concatenate(self.dec_data, axis=0).flatten()
            if self.ax is not None:
                self.ax.plot(ra / self.units, dec / self.units)
                self.ax.set_xlabel(f"RA ({self.units_str})")
                self.ax.set_ylabel(f"DEC ({self.units_str})")
                plt.show()
            else:
                plt.plot(ra / self.units, dec / self.units)
                plt.xlabel(f"RA ({self.units_str})")
                plt.ylabel(f"DEC ({self.units_str})")
                plt.show()


# not a module, but used in DF modules below
def df_IQangle(I, Q, If, Qf, Ff, i_f0=None):
    '''Calculate df using IQ Angle Method.

    I: (1D array of floats) Timestream S21 real component.
    Q: (1D array of floats) Timestream S21 imaginary component.
    If: (1D array of floats) Target sweep S21 real component.
    Qf: (1D array of floats) Target sweep S21 imaginary component.
    Ff: (1D array of floats) Target sweep S21 frequency axis.
    '''

    if i_f0 is None:  # resonant frequency index
        i_f0 = np.argmin(np.abs(If + 1j * Qf))

    cI = (If.max() + If.min()) / 2  # centre of target IQ loop
    cQ = (Qf.max() + Qf.min()) / 2

    # target sweep
    If_c, Qf_c = If - cI, Qf - cQ  # shift center to origin
    θf = np.arctan2(Qf_c, If_c)  # find IQ angles

    # observations
    I_c, Q_c = I - cI, Q - cQ  # shift origin
    θ = np.arctan2(Q_c, I_c)  # find IQ angles

    # adjust frequencies for delta from f0
    Ff0 = Ff - Ff[i_f0]  # center Ff on f0

    # interpolate
    df = np.interp(θ, θf, Ff0, period=2 * np.pi)

    return df / Ff[i_f0]


def add_cal_lamp_df(frame, roach_id=1, iq_key="data"):
    if frame.type != core.G3FrameType.Calibration:
        return

    super_ts = frame[iq_key]

    kids = set([id_str[-6:-2] for id_str in super_ts.names])
    df_data = np.zeros((len(kids), super_ts.data.shape[1]))
    names = []
    for i, kid in enumerate(kids):
        i_idx = int(np.where(np.asarray(super_ts.names) == f"roach{roach_id}_{kid}_I")[0][0])
        q_idx = int(np.where(np.asarray(super_ts.names) == f"roach{roach_id}_{kid}_Q")[0][0])
        kid_i = super_ts.data[i_idx]
        kid_q = super_ts.data[q_idx]

        # load target sweeps
        If = np.array(frame["target_sweeps"][f"roach{roach_id}_{kid}_I"])
        Qf = np.array(frame["target_sweeps"][f"roach{roach_id}_{kid}_Q"])
        Ff = np.array(frame["target_sweeps"][f"roach{roach_id}_{kid}_F"])

        # build df tod
        names.append(f"roach{roach_id}_{kid}")
        df_data[i] = df_IQangle(kid_i, kid_q, If, Qf, Ff)

    times = super_ts.times
    quanta = np.ones(len(kids)) * np.std(df_data) / 10_000

    df_super_ts = so3g.G3SuperTimestream(names, times, df_data, quanta)

    frame["cal_lamp_df"] = df_super_ts


class AddSingleKidDF:
    def __init__(self, roach_id=1, kid="0000"):
        self.roach_id = roach_id
        self.kid = kid
        self.calframe = None

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration:
            self.calframe = frame
        if frame.type != core.G3FrameType.Scan:
            return
        assert self.calframe is not None, "failed to process scan frame: missing prior calibration frame!"

        # load I and Q
        ts: so3g.G3SuperTimestream = frame["data"]
        i_idx = int(np.where(np.asarray(ts.names) == f"roach{self.roach_id}_{self.kid}_I")[0][0])
        q_idx = int(np.where(np.asarray(ts.names) == f"roach{self.roach_id}_{self.kid}_Q")[0][0])
        kid_i = ts.data[i_idx]
        kid_q = ts.data[q_idx]

        # load target sweeps
        If = np.array(self.calframe["target_sweeps"][f"roach{self.roach_id}_{self.kid}_I"])
        Qf = np.array(self.calframe["target_sweeps"][f"roach{self.roach_id}_{self.kid}_Q"])
        Ff = np.array(self.calframe["target_sweeps"][f"roach{self.roach_id}_{self.kid}_F"])

        # build df tod
        df_tod = df_IQangle(kid_i, kid_q, If, Qf, Ff)

        t_i = frame["data"].times[0]
        t_f = frame["data"].times[-1]

        df_ts = core.G3Timestream(df_tod)
        df_ts.start = t_i
        df_ts.stop = t_f
        frame[f"roach{self.roach_id}_{self.kid}_DF"] = df_ts


class AddScanDF:
    def __init__(self, roach_id=1):
        self.roach_id = roach_id
        self.calframe = None

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration:
            self.calframe = frame
        if frame.type != core.G3FrameType.Scan:
            return
        assert self.calframe is not None, "failed to process scan frame: missing prior calibration frame!"

        # get an arbitrarily ordered list of unique kids from calframe keys
        kids: list[str] = list({key[7:11] for key in self.calframe["target_sweeps"].keys()})

        # inputs to G3SuperTimestream constructor
        times: core.G3VectorTime = frame["data"].times  # same timestamps as I/Q data
        names: list[str] = []
        df_data: np.ndarray = np.zeros(shape=(len(kids), len(times)))

        for i, kid in enumerate(kids):
            # load I and Q
            ts: so3g.G3SuperTimestream = frame["data"]
            i_idx = int(np.where(np.asarray(ts.names) == f"roach{self.roach_id}_{kid}_I")[0][0])
            q_idx = int(np.where(np.asarray(ts.names) == f"roach{self.roach_id}_{kid}_Q")[0][0])
            kid_i = ts.data[i_idx]
            kid_q = ts.data[q_idx]
            # load target sweeps
            If = np.array(self.calframe["target_sweeps"][f"roach{self.roach_id}_{kid}_I"])
            Qf = np.array(self.calframe["target_sweeps"][f"roach{self.roach_id}_{kid}_Q"])
            Ff = np.array(self.calframe["target_sweeps"][f"roach{self.roach_id}_{kid}_F"])
            # build df tod and update names/df_data
            df_data[i] = df_IQangle(kid_i, kid_q, If, Qf, Ff)
            names.append(f"roach{self.roach_id}_{kid}")

        compressed_resolution = 10000
        quanta: np.ndarray[float] = ((np.abs(df_data).max() / compressed_resolution)
                                     * np.ones(len(kids)))

        # add the G3SuperTimestream to the scan frame
        df_super_timestream = so3g.G3SuperTimestream(names, times, df_data, quanta)
        frame["df"] = df_super_timestream


class DetectorStats:
    """Determine the median and standard deviation for detectors over all scans

    Stores a `medians` and `stds` attribute, which have shape (n_scans, n_dets)
    """

    def __init__(self, data_key: str = "df"):
        self.data_key = data_key
        self.stds = []
        self.medians = []

    def __call__(self, frame):
        if frame.type != core.G3FrameType.Scan:
            return
        data = frame[self.data_key].data
        self.stds.append(np.std(data, axis=1))
        self.medians.append(np.median(data, axis=1))


def naive_normalize_df(frame, detector_medians=None, detector_stds=None):
    """Normalize tods by setting median to zero and stdev to 1"""
    if frame.type != core.G3FrameType.Scan:
        return

    # data has shape (n_dets, n_samps)
    data = frame["df"].data
    n_dets = data.shape[0]

    data_zeroed = data - detector_medians[:, None]
    norm_df = data_zeroed / detector_stds[:, None]

    out_super_ts = so3g.G3SuperTimestream(
        frame["df"].names,
        frame["df"].times,
        norm_df,
        np.ones(n_dets) * 0.00001,  # quanta - float resolution when compressed, current val is arbitrary
    )

    frame["norm_df"] = out_super_ts


class NormalizeDF():
    def __init__(self, detector_medians=None, in_key="df", out_key="norm_df", cal_df="cal_lamp_df"):
        self.in_key = in_key
        self.out_key = out_key
        self.cal_df = cal_df
        self.calframe = None
        self.detector_medians = detector_medians

    def __call__(self, frame):
        """Normalize tods by setting median to zero and cal lamp max to 1"""
        if frame.type == core.G3FrameType.Calibration:
            self.calframe = frame
        if frame.type != core.G3FrameType.Scan:
            return
        assert self.calframe is not None, "cannot normalize DF in scan without prior calibration data!"

        # data has shape (n_dets, n_samps)
        data = frame[self.in_key].data
        n_dets = data.shape[0]

        data_zeroed = data - self.detector_medians[:, None]
        norm_df = data_zeroed / np.max(self.calframe[self.cal_df].data, axis=1)[:, None]

        out_super_ts = so3g.G3SuperTimestream(
            frame[self.in_key].names,
            frame[self.in_key].times,
            norm_df,
            np.ones(n_dets) * 0.00001,  # quanta - float resolution when compressed, current val is arbitrary
        )

        frame[self.out_key] = out_super_ts


class GenericPlotter:
    def __init__(self, array_getter=None, label: str=None, subplots_args: dict=None, plot_args: dict=None):
        """array_getter is a callable which takes in a frame and returns an array-like object to plot"""
        assert array_getter is not None, "GenericPlotter was not given an array_getter"
        self.array_getter = array_getter
        self.timestreams: list[np.ndarray] = []
        self.label = label if label is not None else self.array_getter.__name__
        self.subplots_args = subplots_args if subplots_args is not None else {}
        self.plot_args = plot_args if plot_args is not None else {}
    def __call__(self, frame):
        if frame.type == core.G3FrameType.EndProcessing:
            flatts: np.ndarray = np.concatenate(self.timestreams, axis=0).flatten()
            fig, ax = plt.subplots(**self.subplots_args)
            ax.plot(flatts, **self.plot_args)
            ax.set_xlabel("index")
            ax.set_ylabel(self.label)
            ax.set_title(self.label + " vs. index")
            plt.show()
        if frame.type == core.G3FrameType.Scan:
            self.timestreams.append(self.array_getter(frame))


class TimeStreamPlotter(GenericPlotter):
    def __init__(self, ts_key=None, label: str=None, subplots_args: dict=None, plot_args: dict=None):
        assert ts_key is not None, "TimeStreamPlotter was not given a ts_key"
        self.ts_key = ts_key
        def timestream_array_getter(frame):
            return np.array(frame[self.ts_key])
        super().__init__(timestream_array_getter, label=label, subplots_args=subplots_args, plot_args=plot_args)


def remove_common_mode(frame):
    # skip frames that don't contain the input key
    if "df" not in frame:
        return

    # get the input timestream data
    ts_in: so3g.G3SuperTimestream = frame["df"]

    # use broadcast to remove common-mode from all timestreams
    tsarr: np.ndarray[np.f64] = ts_in.data
    common_mode = np.mean(tsarr, axis=0)
    out_arr = tsarr - common_mode

    # create output object with correct timestamps and units
    df_ctremoved = so3g.G3SuperTimestream(ts_in.names, ts_in.times, out_arr, ts_in.quanta)

    # save the common mode as well
    common_mode = core.G3Timestream(common_mode)
    common_mode.start = ts_in.times[0]
    common_mode.stop = ts_in.times[-1]

    # store the calibrated timestreams to the output key in the frame
    frame["df_ctremoved"] = df_ctremoved
    frame["common_mode"] = common_mode


class SingleMapBinner:
    def __init__(self, kid, timestreams="df", ra0=None, dec0=None, xlen=None, ylen=None, res=None):
        # center of the sky map
        assert ra0 is not None, "must set ra0!"
        assert dec0 is not None, "must set dec0!"
        self.ra0 = ra0
        self.dec0 = dec0

        self.xlen = xlen if xlen is not None else 1 * core.G3Units.deg
        self.ylen = ylen if ylen is not None else 1 * core.G3Units.deg

        self.res = res if res is not None else 1 * core.G3Units.arcmin

        # number of bins along each axis
        self.nx = int(self.xlen / self.res)
        self.ny = int(self.ylen / self.res)

        # bin edges
        self.ra_edges = np.linspace(-self.xlen / 2, self.xlen / 2, self.nx + 1) + self.ra0
        self.dec_edges = np.linspace(-self.ylen / 2, self.ylen / 2, self.ny + 1) + self.dec0

        self.kid = kid
        self.timestreams = timestreams

        # array for storing the binned timestream data
        self.data = np.zeros((self.ny, self.nx), dtype=float)

        # array for storing the number of times each pixel is "hit" in the timestreams
        self.hits = np.zeros((self.ny, self.nx), dtype=float)

    def source_coords(self):
        """Find source pixel coordinates in this KID's map

        Returns: (x, y)
        """
        zz = self.data / self.hits

        # gaussian smoothing
        smoothing_kernel = 2  # in pixel coords
        nonan_array = np.where(np.isnan(zz), np.nanmedian(zz), zz)
        smoothed_array = gaussian_filter(nonan_array, sigma=smoothing_kernel)

        # identify source in smoothed map by max value (pixel coords)
        max_coords = np.unravel_index(np.argmax(smoothed_array), smoothed_array.shape)
        return max_coords[::-1]

    def plot(self, ax=None):
        with np.errstate(invalid='ignore'):
            m = self.data / self.hits
        source_coords = self.source_coords()
        if ax is not None:
            ax.imshow(m, origin='lower')
            ax.set_xticks(range(self.nx + 1)[::10], [f"{ra:.2f}" for ra in self.ra_edges[::10] / core.G3Units.deg],
                          rotation=45)
            ax.set_yticks(range(self.ny + 1)[::10], [f"{dec:.2f}" for dec in self.dec_edges[::10] / core.G3Units.deg])
            ax.set_xlabel("RA (deg)")
            ax.set_ylabel("DEC (deg)")
            ax.set_title(f"{self.kid}")
            ax.plot(*source_coords, 'ro')
        else:
            plt.imshow(m, origin='lower')
            plt.xticks(range(self.nx + 1)[::10], [f"{ra:.2f}" for ra in self.ra_edges[::10] / core.G3Units.deg],
                       rotation=45)
            plt.yticks(range(self.ny + 1)[::10], [f"{dec:.2f}" for dec in self.dec_edges[::10] / core.G3Units.deg])
            plt.colorbar(label="DF")
            plt.xlabel("RA (deg)")
            plt.ylabel("DEC (deg)")
            plt.title(f"{self.kid}")
            plt.plot(*source_coords, 'ro')
            plt.show()

    def __call__(self, frame):
        if self.timestreams not in frame:
            return

        super_ts = frame[self.timestreams]

        kid_idx = int(np.where(np.array(super_ts.names) == self.kid)[0][0])
        kid_ts = super_ts.data[kid_idx]

        # naively use boresight pointing for now...
        y = np.asarray(frame["dec"])
        x = np.asarray(frame["ra"])

        # update data and hits, in-place
        self.data += np.histogram2d(y, x, bins=[self.dec_edges, self.ra_edges], weights=kid_ts)[0]
        self.hits += np.histogram2d(y, x, bins=[self.dec_edges, self.ra_edges])[0]


rcw_92_min_lat = -77.110
rcw_92_max_lat = -77.070
rcw_92_min_lon = 162.20
rcw_92_max_lon = 162.55
rcw_92_min_alt = 36030
rcw_92_max_alt = 36120
rcw_92_avg_lat = rcw_92_min_lat / 2. + rcw_92_max_lat / 2.
rcw_92_avg_lon = rcw_92_min_lon / 2. + rcw_92_max_lon / 2.
rcw_92_avg_alt = rcw_92_min_alt / 2. + rcw_92_max_alt / 2.
BLASTTNG_SITE = so3g.proj.EarthlySite(rcw_92_avg_lon, rcw_92_avg_lat, rcw_92_avg_alt)  # we could also add weather


class MapBinner:
    def __init__(self, timestreams="df", site=None, source_coords=None, ra0=None, dec0=None, xlen=None, ylen=None,
                 res=None,
                 select_kids: list[str] = None):
        self.timestreams = timestreams
        self.site = site if site is not None else BLASTTNG_SITE
        self.source_coords = source_coords
        assert source_coords is not None, "must set source_coords!"
        assert ra0 is not None, "must set ra0!"
        assert dec0 is not None, "must set dec0!"

        self.ra0 = ra0
        self.dec0 = dec0
        self.xlen = xlen if xlen is not None else 1 * core.G3Units.deg
        self.ylen = ylen if ylen is not None else 1 * core.G3Units.deg
        self.res = res if res is not None else 1 * core.G3Units.arcmin
        # number of bins along each axis
        self.nx = int(self.xlen / self.res)
        self.ny = int(self.ylen / self.res)
        # bin edges
        self.ra_edges = np.linspace(-self.xlen / 2, self.xlen / 2, self.nx + 1) + self.ra0
        self.dec_edges = np.linspace(-self.ylen / 2, self.ylen / 2, self.ny + 1) + self.dec0

        self.select_kids = select_kids

        # array for storing the binned timestream data
        self.data = np.zeros((self.ny, self.nx), dtype=float)
        # array for storing the number of times each pixel is "hit" in the timestreams
        self.hits = np.zeros((self.ny, self.nx), dtype=float)

    def _get_kids(self, super_ts) -> list[str]:
        """Determine the list of kids which should contribute to the map

        Returned KIDs exisit in both:
        - detector timestream data
        - select_kids, if provided
        """
        if self.select_kids is None: return super_ts.names
        return list(set(super_ts.names).intersection(set(self.select_kids)))

    def __call__(self, frame):
        if self.timestreams not in frame:
            return

        super_ts = frame[self.timestreams]

        common_kids = self._get_kids(super_ts)

        for kid in common_kids:
            kid_timestream_idx = int(np.where(np.asarray(super_ts.names) == kid)[0][0])
            kid_ts = super_ts.data[kid_timestream_idx]

            x = frame["ra"] + (self.source_coords[kid][0] - self.nx / 2) * self.res
            y = frame["dec"] + (self.source_coords[kid][1] - self.ny / 2) * self.res

            # update data and hits, in-place
            self.data += np.histogram2d(y, x, bins=[self.dec_edges, self.ra_edges], weights=kid_ts)[0]
            self.hits += np.histogram2d(y, x, bins=[self.dec_edges, self.ra_edges])[0]

    def plot(self, ax=None):
        with np.errstate(invalid='ignore'):
            m = self.data / self.hits
        if ax is not None:
            ax.imshow(m, origin='lower')
            ax.set_xticks(range(self.nx + 1)[::10], [f"{ra:.2f}" for ra in self.ra_edges[::10] / core.G3Units.deg],
                          rotation=45)
            ax.set_yticks(range(self.ny + 1)[::10], [f"{dec:.2f}" for dec in self.dec_edges[::10] / core.G3Units.deg])
            ax.set_xlabel("RA (deg)")
            ax.set_ylabel("DEC (deg)")
            ax.set_title("Combined Map")
        else:
            plt.imshow(m, origin='lower')
            plt.xticks(range(self.nx + 1)[::10], [f"{ra:.2f}" for ra in self.ra_edges[::10] / core.G3Units.deg],
                       rotation=45)
            plt.yticks(range(self.ny + 1)[::10], [f"{dec:.2f}" for dec in self.dec_edges[::10] / core.G3Units.deg])
            plt.colorbar(label="DF")
            plt.xlabel("RA (deg)")
            plt.ylabel("DEC (deg)")
            plt.title("Combined Map")
            plt.show()
