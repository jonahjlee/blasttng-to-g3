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
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from spt3g.core import G3Units as gu
import astropy.units as au
import numpy as np

def add_radec_astropy(frame,
                      az: str="az", el: str="el",
                      lat: str="lat", lon: str="lon",
                      alt: str="alt",
                      ra: str="ra", dec: str="dec"):
    """
    G3 pipeline module.
    Uses astropy coordinate transformations to convert az/el --> ra/dec.

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

    az_deg = np.array(frame[az]) / gu.deg
    el_deg = np.array(frame[el]) / gu.deg
    lat_deg = np.array(frame[lat]) / gu.deg
    lon_deg = np.array(frame[lon]) / gu.deg
    alt_m = np.array(frame[alt]) / gu.m

    unix_times = np.array(frame[az].times) / gu.s
    times = Time(unix_times, format="unix")

    blasttng_loc = EarthLocation(lat=lat_deg, lon=lon_deg, height=alt_m)
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

def add_radec_so3g(frame,
                   az: str="az", el: str="el",
                   lat: str="lat", lon: str="lon",
                   alt: str="alt",
                   ra: str="ra", dec: str="dec"):
    """
    G3 pipeline module.
    Use so3g coordinate transformations to convert az/el --> ra/dec.

    Telescope location passed into `so3g.proj.CelestialSightLine` is given by the average of
    the location at the start of frame and end of frame.

    Uses `so3g.proj.CelestialSightLine.naive_az_el`, which does not account for weather in calculations:

        "This will be off by several arcminutes… but less than a degree. The weather is ignored."

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

    x = np.mod(coords[0][:, 0], 2 * np.pi) * gu.rad
    y = coords[0][:, 1] * gu.rad

    frame[ra] = core.G3VectorDouble(x)
    frame[dec] = core.G3VectorDouble(y)
