# ============================================================================ #
# maps.py
#
# Jonah Lee
#
# G3 pipeline modules for map-making (binning images, combining maps, etc.)
# ============================================================================ #

import so3g
from spt3g import core
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import tools


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


class SingleMapBinner:
    """
    G3 Pipeline Module.
    Bins one detector's TOD (time-ordered-data) into a flat sky map with Plate-Carree projection.
    """
    def __init__(self,
                 kid,
                 timestreams="df",
                 ra0: float=None, dec0: float=None,
                 xlen: float=None, ylen: float=None,
                 res: float=None,
                 roach: int=None):
        """
        Create a new SingleMapBinner.

        :param kid: Identifier for a detector
            If `kid` is a string, it must match the expression `^roach[1-5]_\d{4}$`, e.g. "roach1_0000"
            or be parsable as an int or float.
            If `kid` is (or is parsed as) an int or float, it will be parsed into a string matching the above pattern,
            and `roach` must be defined.
        :param timestreams: Key into detector G3SuperTimestream for scan frames
        :param ra0: G3Units angle - center of map in right ascension
        :param dec0: G3Units angle - center of map in declination
        :param xlen: G3Units angle - width of map in right ascension
        :param ylen: G3Units angle - height of map in declination
        :param res: G3Units angle - size of map square pixels
        :param roach: Roach number, used to find detector identifier string if `kid` is an int/float
        """
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

        self.kid = kid if roach is None else tools.kid_string(kid, roach)
        self.timestreams = timestreams

        # array for storing the binned timestream data
        self.data = np.zeros((self.ny, self.nx), dtype=float)

        # array for storing the number of times each pixel is "hit" in the timestreams
        self.hits = np.zeros((self.ny, self.nx), dtype=float)

    def source_coords(self):
        """
        Find source pixel coordinates in this KID's map.
        Applies a gaussian smoothing filter, then finds the location of the brightest pixel.

        :returns: (x, y)
        """
        with np.errstate(invalid='ignore'):
            zz = self.data / self.hits

        # gaussian smoothing
        smoothing_kernel = 2  # in pixel coords
        nonan_array = np.where(np.isnan(zz), np.nanmedian(zz), zz)
        smoothed_array = gaussian_filter(nonan_array, sigma=smoothing_kernel)

        # identify source in smoothed map by max value (pixel coords)
        max_coords = np.unravel_index(np.argmax(smoothed_array), smoothed_array.shape)
        return max_coords[::-1]

    def plot(self, ax=None, show=True):
        """
        Plot this SingleMapBinner's map.

        Image data for `matplotlib.pyplot.imshow` is given by `self.data / self.hits`.

        :param ax: Matplotlib axes instance. If None (default), will use plt.gca()
        :param show: If True (default), calls plt.show() automatically.
        """
        with np.errstate(invalid='ignore'):
            m = self.data / self.hits
        source_coords = self.source_coords()

        if ax is None: ax = plt.gca()

        ax.imshow(m, origin='lower')
        ax.set_xticks(range(self.nx + 1)[::10], [f"{ra:.2f}" for ra in self.ra_edges[::10] / core.G3Units.deg],
                      rotation=45)
        ax.set_yticks(range(self.ny + 1)[::10], [f"{dec:.2f}" for dec in self.dec_edges[::10] / core.G3Units.deg])
        ax.set_xlabel("Boresight RA (deg)")
        ax.set_ylabel("Boresight DEC (deg)")
        ax.set_title(f"{self.kid}")
        ax.plot(*source_coords, 'ro')

        if show: plt.show()

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


class MapBinner:
    def __init__(self, timestreams="df", source_coords: dict=None, stds: dict=None,
                 ra0: float=None, dec0: float=None,
                 xlen: float=None, ylen: float=None,
                 res: float=None,
                 select_kids: list[str] = None):
        """
        Create a new MapBinner.

        :param timestreams: Key into detector G3SuperTimestream for scan frames
        :param source_coords: Mapping from detector identifiers (^roach[1-5]_\d{4}$) to source (col, row) pixel coords
        :param stds: Mapping from detector identifiers (^roach[1-5]_\d{4}$) to detector signal standard deviation.
                     Detectors are weighted by 1/σ² when combining maps. If `None` (default), weights maps equally.
        :param ra0: G3Units angle - center of map in right ascension
        :param dec0: G3Units angle - center of map in declination
        :param xlen: G3Units angle - width of map in right ascension
        :param ylen: G3Units angle - height of map in declination
        :param res: G3Units angle - size of map square pixels
        :param select_kids: Optional, list of kids to include in combined map. If `None` (default), all kids are included.
        """
        self.timestreams = timestreams
        self.source_coords = source_coords
        self.stds = stds
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

    def _get_kids(self, ts_names: list[str]) -> list[str]:
        """Determine the list of kids which should contribute to the map

        Returned KIDs exisit in both:
        - detector timestream data
        - select_kids, if provided

        :param ts_names: List of detector names in G3SuperTimestream
        """
        if self.select_kids is None: return ts_names
        return list(set(ts_names).intersection(set(self.select_kids)))

    def __call__(self, frame):
        """Update MapBinner with a new frame. Called within a G3 pipeline."""
        if self.timestreams not in frame:
            return

        super_ts = frame[self.timestreams]

        common_kids = self._get_kids(super_ts)

        for kid in common_kids:
            kid_timestream_idx = int(np.where(np.asarray(super_ts.names) == kid)[0][0])
            kid_ts = super_ts.data[kid_timestream_idx]

            x = frame["ra"] + (self.nx / 2 - self.source_coords[kid][0]) * self.res
            y = frame["dec"] + (self.ny / 2 - self.source_coords[kid][1]) * self.res

            # update data and hits, in-place
            if self.stds is not None:
                kid_weight = self.stds[kid_timestream_idx] ** -2  # same index as super_ts
            else:
                kid_weight = 1
            kid_data = np.histogram2d(y, x, bins=[self.dec_edges, self.ra_edges], weights=kid_ts)[0] * kid_weight
            self.data += kid_data
            kid_hits = np.histogram2d(y, x, bins=[self.dec_edges, self.ra_edges])[0] * kid_weight
            self.hits += kid_hits

    def plot(self, ax=None, show=True):
        """
        Plot this MapBinner's map.

        Image data for `matplotlib.pyplot.imshow` is given by `self.data / self.hits`.

        :param ax: Matplotlib axes instance. If None (default), will use plt.gca()
        :param show: If True (default), calls plt.show() automatically.
        """
        with np.errstate(invalid='ignore'):
            m = self.data / self.hits

        if ax is None: ax = plt.gca()
        ax.imshow(m, origin='lower')
        ax.set_xticks(range(self.nx + 1)[::10], [f"{ra:.2f}" for ra in self.ra_edges[::10] / core.G3Units.deg],
                      rotation=45)
        ax.set_yticks(range(self.ny + 1)[::10], [f"{dec:.2f}" for dec in self.dec_edges[::10] / core.G3Units.deg])
        ax.set_xlabel("RA (deg)")
        ax.set_ylabel("DEC (deg)")
        ax.set_title("Combined Map")

        if show: plt.show()
