# ============================================================================ #
# signal.py
#
# Jonah Lee
#
# Various tools for G3 TOD processing
# ============================================================================ #

from spt3g import core
import so3g
import numpy as np


class DetectorStats:
    """
    G3 pipeline module.
    Determines the median and standard deviation for detectors over all scans

    Stores a `medians` and `stds` attribute, which have shape (n_scans, n_dets)
    """
    def __init__(self, data_key: str="df"):
        """
        Instantiate a DetectorStats object

        :param data_key: Key into a G3SuperTimestream for scan frames.
        """
        self.data_key = data_key
        self.stds = []
        self.medians = []

    def __call__(self, frame):
        if frame.type != core.G3FrameType.Scan:
            return
        data = frame[self.data_key].data
        self.stds.append(np.std(data, axis=1))
        self.medians.append(np.median(data, axis=1))


# used in DF modules below
def df_iqangle(I, Q, If, Qf, Ff, i_f_theta=None):
    """
    Calculate df using IQ Angle Method.

    :param I: (1D array of floats) Timestream S21 real component.
    :param Q: (1D array of floats) Timestream S21 imaginary component.
    :param If: (1D array of floats) Target sweep S21 real component.
    :param Qf: (1D array of floats) Target sweep S21 imaginary component.
    :param Ff: (1D array of floats) Target sweep S21 frequency axis.
    :param i_f_theta: Resonant frequency Index
    """

    if i_f_theta is None:  # resonant frequency index
        i_f_theta = np.argmin(np.abs(If + 1j * Qf))

    cI = (If.max() + If.min()) / 2  # centre of target IQ loop
    cQ = (Qf.max() + Qf.min()) / 2

    # target sweep
    If_c, Qf_c = If - cI, Qf - cQ  # shift center to origin
    theta_f = np.arctan2(Qf_c, If_c)  # find IQ angles

    # observations
    I_c, Q_c = I - cI, Q - cQ  # shift origin
    theta = np.arctan2(Q_c, I_c)  # find IQ angles

    # adjust frequencies for delta from f0
    Ff_theta = Ff - Ff[i_f_theta]  # center Ff on f0

    # interpolate
    df = np.interp(theta, theta_f, Ff_theta, period=2 * np.pi)

    return df / Ff[i_f_theta]


def compute_df_data(kids: list[str], super_ts, target_sweeps) -> np.ndarray:
    """
    Note: This is not a G3 pipeline module. See `add_cal_lamp_df` and `AddScanDF`, which call this function internally.
    Computes a DF (delta-frequency) ndarray with shape (n_dets, n_samps).

    :param kids: List of names of kids, e.g. "roach1_0000". Names should match `^roach[1-5]_\d{4}$`, and their order
                 will be preserved as the indices for axis 0 in the resulting 2d ndarray.
    :param super_ts: G3SuperTimestream containing KID I/Q data. Names are like above, except there are twice
                     as many because half have '_I' appended to indicate in-phase and half have '_Q' appended to indicate
                     quadrature
    :param target_sweeps: G3TimestreamMap containing calibration sweep data with indices appended by '_I', '_Q' or '_F'.
    :return: a DF ndarray with shape (n_dets, n_samps).
    """
    df_data = np.zeros((len(kids), super_ts.data.shape[1]))
    for i, kid in enumerate(kids):
        i_idx = int(np.where(np.asarray(super_ts.names) == f"{kid}_I")[0][0])
        q_idx = int(np.where(np.asarray(super_ts.names) == f"{kid}_Q")[0][0])
        kid_i = super_ts.data[i_idx]
        kid_q = super_ts.data[q_idx]

        If = np.array(target_sweeps[f"{kid}_I"])
        Qf = np.array(target_sweeps[f"{kid}_Q"])
        Ff = np.array(target_sweeps[f"{kid}_F"])

        df_data[i] = df_iqangle(kid_i, kid_q, If, Qf, Ff)
    return df_data


def add_cal_lamp_df(frame, iq_key: str="cal_lamp_data", target_sweeps_key="target_sweeps", out_key="cal_lamp_df"):
    """
    G3 pipeline module.
    Compute DF (delta-frequency) for the calibration lamp data stored in the calibration frame.

    :param frame: G3FrameObject passed in automatically by pipeline.
    :param iq_key: Key to I/Q G3SuperTimestream in calibration frame.
    :param out_key: Key to output DF G3SuperTimestream into in calibration frame.
    """
    if frame.type != core.G3FrameType.Calibration:
        return

    super_ts = frame[iq_key]

    kids = list({id_str[:-2] for id_str in super_ts.names})
    times = super_ts.times
    df_data = compute_df_data(kids, super_ts, frame[target_sweeps_key])
    quanta = np.ones(len(kids)) * np.std(df_data) / 10_000

    df_super_ts = so3g.G3SuperTimestream(kids, times, df_data, quanta)
    frame[out_key] = df_super_ts


class AddScanDF:
    """
    G3 pipeline module.
    Compute DF (delta-frequency) data for all scan frames.
    """
    def __init__(self, iq_key="data", target_sweeps_key="target_sweeps", out_key="df"):
        """
        Instantiate an AddScanDF object

        :param iq_key: Key to I/Q G3SuperTimestream in scan frame.
        :param target_sweeps_key: Key to calibration sweep G3SuperTimestream in the calibration frame.
        :param out_key: Key to output DF G3SuperTimestream into in scan frame.
        """
        self.calframe = None
        self.iq_key = iq_key
        self.target_sweeps_key = target_sweeps_key
        self.out_key = out_key

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration:
            self.calframe = frame
        if frame.type != core.G3FrameType.Scan:
            return

        assert self.calframe is not None, "failed to process scan frame: missing prior calibration frame!"
        target_sweeps = self.calframe[self.target_sweeps_key]

        kids = list({id_str[:-2] for id_str in target_sweeps.keys()})
        times = frame[self.iq_key].times
        df_data = compute_df_data(kids, frame[self.iq_key], target_sweeps)
        quanta = (np.abs(df_data).max() / 10_000) * np.ones(len(kids))

        df_super_timestream = so3g.G3SuperTimestream(kids, times, df_data, quanta)
        frame[self.out_key] = df_super_timestream


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
        cal_data = self.calframe[self.cal_df].data
        n_dets = data.shape[0]

        # set the median of the normalized array to zero
        data_zeroed = data - self.detector_medians[:, None]
        cal_data_zeroed = cal_data - self.detector_medians[:, None]

        # scale such that the maximum during calibration is set to 1
        norm_df = data_zeroed / np.max(cal_data_zeroed, axis=1)[:, None]

        out_super_ts = so3g.G3SuperTimestream(
            frame[self.in_key].names,
            frame[self.in_key].times,
            norm_df,
            np.ones(n_dets) * 0.00001,  # quanta - float resolution when compressed, current val is arbitrary
        )

        frame[self.out_key] = out_super_ts


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
