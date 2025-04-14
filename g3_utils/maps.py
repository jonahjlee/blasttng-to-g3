import so3g
from spt3g import core
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


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
        with np.errstate(invalid='ignore'):
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
            ax.set_xlabel("Boresight RA (deg)")
            ax.set_ylabel("Boresight DEC (deg)")
            ax.set_title(f"{self.kid}")
            ax.plot(*source_coords, 'ro')
        else:
            plt.imshow(m, origin='lower')
            plt.xticks(range(self.nx + 1)[::10], [f"{ra:.2f}" for ra in self.ra_edges[::10] / core.G3Units.deg],
                       rotation=45)
            plt.yticks(range(self.ny + 1)[::10], [f"{dec:.2f}" for dec in self.dec_edges[::10] / core.G3Units.deg])
            plt.colorbar(label="DF")
            plt.xlabel("Boresight RA (deg)")
            plt.ylabel("Boresight DEC (deg)")
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
