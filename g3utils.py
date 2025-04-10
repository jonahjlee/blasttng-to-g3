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
import spt3g.core.G3Units as gu

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


class FirstFrameGrabber:
    """Stores the first frame of a given type"""
    def __init__(self, frame_type: core.G3FrameType=None):
        self.frame_type = frame_type if frame_type is not None else core.G3FrameType.Scan
        self.first_frame = None
    def __call__(self, frame):
        if self.first_frame is not None:
            # already found the frame
            return
        if frame.type == self.frame_type:
            self.first_frame = frame
            print(f"Found the first frame with type: {self.frame_type}!")
            print(f"The frame is now stored in {self}'s first_frame attribute.")
            return


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


def plot_ra_dec(frame, ra_key="ra", dec_key="dec"):
    # skip any frame that doesn't contain the right key
    if ra_key not in frame or dec_key not in frame:
        return

    # plot coordinates in real units
    plt.plot(frame[ra_key] / gu.deg, frame[dec_key] / gu.deg)


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
        if frame.type != core.G3FrameType.Scan:
            return
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

        self.xlen = xlen if xlen is not None else 1 * gu.deg
        self.ylen = ylen if ylen is not None else 1 * gu.deg

        self.res = res if res is not None else 1 * gu.arcmin

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

    def plot(self, ax=None):
        with np.errstate(invalid='ignore'):
            m = self.data / self.hits
            if ax is not None:
                ax.imshow(m, origin='lower')
                ax.set_xticks(range(self.nx+1)[::10], [f"{ra:.2f}" for ra in self.ra_edges[::10] / gu.deg], rotation=45)
                ax.set_yticks(range(self.ny+1)[::10], [f"{dec:.2f}" for dec in self.dec_edges[::10] / gu.deg])
                ax.set_xlabel("RA (deg)")
                ax.set_ylabel("DEC (deg)")
                ax.set_title(f"{self.kid}")
            else:
                plt.imshow(m, origin='lower')
                plt.xticks(range(self.nx+1)[::10], [f"{ra:.2f}" for ra in self.ra_edges[::10] / gu.deg], rotation=45)
                plt.yticks(range(self.ny+1)[::10], [f"{dec:.2f}" for dec in self.dec_edges[::10] / gu.deg])
                plt.colorbar(label="DF")
                plt.xlabel("RA (deg)")
                plt.ylabel("DEC (deg)")
                plt.title(f"{self.kid}")
                plt.show()
