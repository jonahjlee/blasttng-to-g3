# ============================================================================ #
# coords.py
#
# Jonah Lee
#
# G3 pipeline modules to help with astronomical coordinate transformations,
# mainly converting Az/El to RA/Dec
# ============================================================================ #

import so3g
from spt3g import core
from spt3g import maps
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spt3g.core import G3Units as gu
import astropy.units as au
import numpy as np


def get_earthlocation(frame, lat: str, lon: str, alt: str, use_middle_loc: bool):
    """Get a vector of astropy EarthLocation objects from a frame

    :param frame: G3Frame
    :param lat: Key to G3Timestream of latitude
    :param lon: Key to G3Timestream of longitude
    :param alt: Key to G3Timestream of altitude
    :param use_middle_loc: Compute a single EarthLocation as the mean of the start and end times.
                           Improves performance significantly, but reduces precision.
    """
    if use_middle_loc:
        lat_deg = (frame[lat][-1] - frame[lat][0]) / gu.deg
        lon_deg = (frame[lat][-1] - frame[lat][0]) / gu.deg
        alt_m = (frame[lat][-1] - frame[lat][0]) / gu.deg
    else:
        lat_deg = np.array(frame[lat]) / gu.deg
        lon_deg = np.array(frame[lon]) / gu.deg
        alt_m = np.array(frame[alt]) / gu.m
    blasttng_loc = EarthLocation(lat=lat_deg, lon=lon_deg, height=alt_m)
    return blasttng_loc


def add_radec_astropy(
    frame,
    az: str="az", el: str="el",
    lat: str="lat", lon: str="lon",
    alt: str="alt",
    ra: str="ra", dec: str="dec",
    use_middle_loc = False
):
    """
    G3 pipeline module.
    Uses astropy coordinate transformations to convert az/el --> ra/dec.

    :param frame: G3Frame passed in automatically by pipeline
    :param az: Key for azimuth G3Timestream in scan frame
    :param el: Key for elevation angle G3Timestream in scan frame
    :param lat: Key for latitude G3Timestream in scan frame
    :param lon: Key for longitude G3Timestream in scan frame
    :param alt: Key for altitude G3Timestream in scan frame
    :param ra: Key in which to output right ascension G3Timestream in scan frame
    :param dec: Key in which to output declination G3Timestream in scan frame
    :param use_middle_loc: Compute a single EarthLocation as the mean of the start and end times.
                           Improves performance significantly, but reduces precision.
    """
    if frame.type != core.G3FrameType.Scan:
        return

    az_deg = np.array(frame[az]) / gu.deg
    el_deg = np.array(frame[el]) / gu.deg

    unix_times = np.array(frame[az].times) / gu.s
    times = Time(unix_times, format="unix")

    blasttng_loc = get_earthlocation(frame, lat, lon, alt, use_middle_loc)
    sky_coords = SkyCoord(alt=el_deg * au.deg,
                          az=az_deg * au.deg,
                          obstime=times,
                          frame='altaz',
                          location=blasttng_loc)
    t_i = frame[az].start
    t_f = frame[az].stop

    ra_ts = core.G3Timestream(sky_coords.icrs.ra.deg * gu.deg)
    ra_ts.start = t_i
    ra_ts.stop = t_f
    frame[ra] = ra_ts

    dec_ts = core.G3Timestream(sky_coords.icrs.dec.deg * gu.deg)
    dec_ts.start = t_i
    dec_ts.stop = t_f
    frame[dec] = dec_ts


def add_radec_spt3g(
    frame,
    az: str="az", el: str="el",
    lat: str="lat", lon: str="lon",
    alt: str="alt",
    ra: str="ra", dec: str="dec",
    use_middle_loc=False
):
    """
    G3 pipeline module.
    Uses spt3g coordinate transformations to convert az/el --> ra/dec.

    :param frame: G3Frame passed in automatically by pipeline
    :param az: Key for azimuth G3Timestream in scan frame
    :param el: Key for elevation angle G3Timestream in scan frame
    :param lat: Key for latitude G3Timestream in scan frame
    :param lon: Key for longitude G3Timestream in scan frame
    :param alt: Key for altitude G3Timestream in scan frame
    :param ra: Key in which to output right ascension G3Timestream in scan frame
    :param dec: Key in which to output declination G3Timestream in scan frame
    :param use_middle_loc: Compute a single EarthLocation as the mean of the start and end times.
                           Improves performance significantly, but reduces precision.
    """
    if frame.type != core.G3FrameType.Scan:
        return

    blasttng_loc = get_earthlocation(frame, lat, lon, alt, use_middle_loc)
    ra_ts, dec_ts = maps.convert_azel_to_radec(frame[az], frame[el], blasttng_loc)
    frame[ra] = ra_ts
    frame[dec] = dec_ts


def add_radec_so3g(
    frame,
    az: str="az", el: str="el",
    lat: str="lat", lon: str="lon",
    alt: str="alt",
    ra: str="ra", dec: str="dec"
):
    """
    G3 pipeline module.
    Use so3g coordinate transformations to convert az/el --> ra/dec.

    Telescope location passed into `so3g.proj.CelestialSightLine` is given by the average of
    the location at the start of frame and end of frame.

    Uses `so3g.proj.CelestialSightLine.naive_az_el`, which does not account for weather in calculations:

        "This will be off by several arcminutes… but less than a degree. The weather is ignored."

    :param frame: G3Frame passed in automatically by pipeline
    :param az: Key for azimuth G3Timestream in scan frame
    :param el: Key for elevation angle G3Timestream in scan frame
    :param lat: Key for latitude G3Timestream in scan frame
    :param lon: Key for longitude G3Timestream in scan frame
    :param alt: Key for altitude G3Timestream in scan frame
    :param ra: Key in which to output right ascension G3Timestream in scan frame
    :param dec: Key in which to output declination G3Timestream in scan frame
    """
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
    times = np.array(frame[az].times) / gu.s
    csl = so3g.proj.CelestialSightLine.naive_az_el(
        times,
        frame[az] / gu.rad,
        frame[el] / gu.rad,
        site=site
    )
    coords = csl.coords(so3g.proj.FocalPlane.boresight())

    ra_ts = core.G3Timestream(np.mod(coords[0][:, 0], 2 * np.pi) * gu.rad)
    dec_ts = core.G3Timestream(coords[0][:, 1] * gu.rad)

    ra_ts.start = frame[az].start
    ra_ts.stop = frame[az].stop
    dec_ts.start = frame[az].start
    dec_ts.stop = frame[az].stop

    frame[ra] = ra_ts
    frame[dec] = dec_ts
