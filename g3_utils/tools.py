from spt3g import core
import numpy as np
import matplotlib.pyplot as plt
from spt3g.core import G3Units as gu


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


class GenericPlotter:
    def __init__(self, array_getter=None, label: str=None, ax=None, show=True):
        """Helper module to plot data in a G3Pipeline. Accumulates the data from array_getter into a 1-d array
        for the entire pipeline, then plots the data when an EndProcessing frame is encountered.

        :param array_getter: a callable which takes in a frame and returns an array-like object to plot
        :param label: a label for the plot y-axis and title
        :param ax: a matplotlib axes object. If `None`, uses `plt.gca()`
        :param show: If `True`, add labels and show when an EndProcessing frame is reached.
        """
        assert array_getter is not None, "GenericPlotter was not given an array_getter"
        self.array_getter = array_getter
        self.timestreams: list[np.ndarray] = []
        self.label = label if label is not None else self.array_getter.__name__
        self.ax = ax
        self.show = show

    def __call__(self, frame):
        if frame.type == core.G3FrameType.EndProcessing:
            ax = self.ax if self.ax is not None else plt.gca()
            flatts: np.ndarray = np.concatenate(self.timestreams, axis=0).flatten()
            ax.plot(flatts)
            if not self.show: return
            ax.set_xlabel("index")
            ax.set_ylabel(self.label)
            ax.set_title(self.label + " vs. index")
            plt.show()
        if frame.type == core.G3FrameType.Scan:
            self.timestreams.append(self.array_getter(frame))


class TimeStreamPlotter(GenericPlotter):
    def __init__(self, ts_key=None, label: str=None, ax=None, show=True):
        assert ts_key is not None, "TimeStreamPlotter was not given a ts_key"
        self.ts_key = ts_key
        def timestream_array_getter(frame):
            return np.array(frame[self.ts_key])
        super().__init__(timestream_array_getter, label=label, ax=ax, show=show)

