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


def add_cal_lamp_df(
    frame,
    iq_key: str="cal_lamp_data",
    target_sweeps_key: str="target_sweeps",
    out_key: str="cal_lamp_df",
    select_kids: list[str] | tuple[str] | set[str] = None,
    exclude_kids: list[str] | tuple[str] | set[str] = None
):
    """
    G3 pipeline module.
    Compute DF (delta-frequency) for the calibration lamp data stored in the calibration frame.

    :param frame: G3FrameObject passed in automatically by pipeline.
    :param target_sweeps_key: Key to target sweeps G3SuperTimestream in calibration frame.
    :param iq_key: Key to I/Q G3SuperTimestream in calibration frame.
    :param out_key: Key to output DF G3SuperTimestream into in calibration frame.
    :param select_kids: Optional, list of kids to include in combined map. If `None` (default), all kids are included.
    :param exclude_kids: Optional, list of kids to exclude form combined map. If `None` (default), no kids are excluded.
                         Note: If both select_kids and exclude_kids are provided, uses only KIDs that are in
                         select_kids and not in exclude_kids.
    """
    if frame.type != core.G3FrameType.Calibration:
        return

    super_ts = frame[iq_key]

    kids = {id_str[:-2] for id_str in super_ts.names}
    if select_kids is not None:
        kids = kids.intersection(set(select_kids))
    if exclude_kids is not None:
        kids = kids - set(exclude_kids)
    kid_list = sorted(list(kids))
    times = super_ts.times
    df_data = compute_df_data(kid_list, super_ts, frame[target_sweeps_key])
    quanta = np.ones(len(kid_list)) * np.std(df_data) / 10_000

    df_super_ts = so3g.G3SuperTimestream(kid_list, times, df_data, quanta)
    frame[out_key] = df_super_ts


class AddScanDF:
    """
    G3 pipeline module.
    Compute DF (delta-frequency) data for all scan frames.
    """
    def __init__(self,
                 iq_key="data", target_sweeps_key="target_sweeps", out_key="df",
                 select_kids: list[str] | tuple[str] | set[str] = None,
                 exclude_kids: list[str] | tuple[str] | set[str] = None):
        f"""
        Instantiate an AddScanDF object

        :param iq_key: Key to I/Q G3SuperTimestream in scan frame.
        :param target_sweeps_key: Key to calibration sweep G3SuperTimestream in the calibration frame.
        :param out_key: Key to output DF G3SuperTimestream into in scan frame.
        :param select_kids: If provided, only compute/store DF for kids in select_kids.
        :param exclude_kids: If provided, do not compute DF for kids in exclude_kids.
                             If neither select_kids nor exclude_kids are provided, use all detectors.
                             If both are provided, use set(select_kids) - set(exclude_kids).
        """
        self.iq_key = iq_key
        self.target_sweeps_key = target_sweeps_key
        self.out_key = out_key
        self.select_kids = select_kids
        self.exclude_kids = exclude_kids
        self.target_sweeps = None
        self.kid_list = None

    def _get_kid_list(self) -> list[str]:
        """
        Determine the ordered list of KIDs to output DF for.
        This determines the `.names` list for the resulting df G3SuperTimestream.
        """
        kids = {id_str[:-2] for id_str in self.target_sweeps.keys()}
        if self.select_kids is not None:
            kids = kids.intersection(set(self.select_kids))
        if self.exclude_kids is not None:
            kids = kids - set(self.exclude_kids)
        return sorted(list(kids))

    def __call__(self, frame):
        if frame.type == core.G3FrameType.Calibration:
            self.target_sweeps = frame[self.target_sweeps_key]
            self.kid_list = self._get_kid_list()
        if frame.type != core.G3FrameType.Scan:
            return

        assert self.target_sweeps is not None, ("failed to process scan frame: "
                                                "missing target sweeps from prior calibration frame!")

        times = frame[self.iq_key].times
        df_data = compute_df_data(self.kid_list, frame[self.iq_key], self.target_sweeps)
        quanta = (np.abs(df_data).max() / 10_000) * np.ones(len(self.kid_list))

        df_super_timestream = so3g.G3SuperTimestream(self.kid_list, times, df_data, quanta)
        frame[self.out_key] = df_super_timestream


class NormalizeDF():
    """
    G3 pipeline module.
    Normalizes tods by setting median to zero and calibration lamp max to 1
    """
    def __init__(self, detector_medians=None, in_key="df", out_key="norm_df", cal_df_key="cal_lamp_df"):
        """
        Instantiate an NormalizeDF module.

        Requires that cal lamp df and scan df have the same shape (n_dets, n_samps). If you are filtering detectors
        using select_kids and/or exclude_kids, make sure that this is consistent between calibration and scan df.

        :param detector_medians: Mapping from detector identifiers (^roach[1-5]_\d{4}$) to detector signal median value.
        :param in_key: Key to input G3SuperTimestream in scan frame.
        :param out_key: Key to output G3SuperTimestream into in scan frame.
        :param cal_df_key: Key to calibration lamp DF G3SuperTimestream in calibration frame.
        """
        self.in_key = in_key
        self.out_key = out_key
        self.cal_df = cal_df_key
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


def remove_common_mode(frame, in_key="df", out_key="df_ctremoved", ct_key=None):
    """
    G3 pipeline module.
    Removes the commmon-mode c(t) from detector timestreams.

    :param frame: G3FrameObject passed in automatically by pipeline.
    :param in_key: Key to input DF G3SuperTimestream in scan frame.
    :param out_key: Key to c(t) removed DF G3SuperTimestream into in scan frame.
    :param ct_key: Optional - key in which to store the common-mode as a G3Timestream.
                   If this is `None` (default), do not store the common-mode.
    """
    # skip frames that don't contain the input key
    if in_key not in frame:
        return

    # get the input timestream data
    ts_in: so3g.G3SuperTimestream = frame[in_key]

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
    frame[out_key] = df_ctremoved
    if ct_key is not None: frame[ct_key] = common_mode


def azelToMapPix(az, el, x_edges, y_edges):
    '''Convert az/el coords to map pix coords.

    az/el: (1D array of floats) Array of az/el coordinates.
    x_edges/y_edges: (1D array of floats) The map az/el bin edges.
    '''
    indices_x = np.searchsorted(x_edges, az, side='right') - 1
    indices_y = np.searchsorted(y_edges, el, side='right') - 1
    # Correct the indices if they go out of bounds
    indices_x = np.clip(indices_x, 0, len(x_edges) - 2)
    indices_y = np.clip(indices_y, 0, len(y_edges) - 2)

    return indices_x, indices_y


def common_mode_iter(frame, in_key="df_ctremoved", out_key="df_iterated", kid_shifts=None,
                     prev_map=None, ra0=None, dec0=None, xlen=None, ylen=None, res=None):
    if frame.type != core.G3FrameType.Scan:
        return

    assert kid_shifts is not None, ("kid_shifts is required for common_mode_iter.\n"
                                    "It can be obtained from the previous MapBinner's kid_shifts property")

    # figure out bins
    nx = int(xlen / res)
    ny = int(ylen / res)
    ra_edges = np.linspace(-xlen / 2, xlen / 2, nx + 1) + ra0
    dec_edges = np.linspace(-ylen / 2, ylen / 2, ny + 1) + dec0

    # get astronomical signal tod estimate from map for each kid
    ast_signal_estimate = np.zeros_like(frame[in_key].data)
    for i, kid in enumerate(frame[in_key].names):
        x = frame["ra"] + kid_shifts[kid][0]
        y = frame["dec"] + kid_shifts[kid][1]
        indices_x, indices_y = azelToMapPix(x, y, ra_edges, dec_edges)
        ast_signal_estimate[i] = prev_map[indices_y, indices_x]
    common_mode = np.nanmean(frame[in_key].data - ast_signal_estimate, axis=0)
    # plt.imshow(ast_signal_estimate); plt.title("est. ast.") ;plt.show();
    # plt.plot(common_mode); plt.title("ct it 2") ;plt.show();
    ct_removed = frame[in_key].data - common_mode[None, :]

    super_ts = so3g.G3SuperTimestream(frame[in_key].names,
                                      frame[in_key].times,
                                      ct_removed,
                                      frame[in_key].quanta)
    frame[out_key] = super_ts
