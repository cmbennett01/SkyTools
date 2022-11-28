#!/usr/bin/env python
"""
File: atlas.py
Author: Christopher Bennett
Email: chris@east20thst.net
Github: https://github.com/cmbennett01
Description: Script to lookup Uranometria 2000.0 chart number for a coordinate.
"""

from math import ceil
import numpy as np


# Structured array for U2K bands (lower cutoff and # panels
u2k_zones = np.array(
    [(84.5,  1), (73.5,  6), (62.0, 10), (51.0, 12), (40.0, 15),
     (29.0, 18), (17.0, 18), (5.5,  20), (0.0,  20), (0.0,   0)],
     dtype=[('lowDec', '>f4'), ('numZones', '>i4')])


def u2k_chart(coord_obj):
    """Uranometria 2000.0 Deep Sky Atlas chart number.
    u2k_chart(coord_obj)

    Parameters
    ----------
    coord_obj : astropy SkyCoord

    Returns
    -------
    tuple: vol, pg
        vol : int (Uranometria volume number)
        pg : int (Uranometria chart number)

    Notes
    -----
    Python implementation of `pyephem/libastro-3.7.7/atlas.c` by
    Brandon Rhodes <https://github.com/brandon-rhodes/pyephem>
    """

    ra = coord_obj.ra.hour
    dec = coord_obj.dec.deg

    if ra < 0.0 or ra >= 24.0 or dec < -90.0 or dec > 90.0:
        raise ValueError("Invalid ra and/or dec")

    if dec < 0.0:
        dec = -dec
        south = 1  # South is mirror of North
    else:
        south = 0

    panel, band = 1, 0

    # scan u2k_zones for the correct band:
    while u2k_zones[band]['numZones'] != 0 and \
            dec <= u2k_zones[band]['lowDec']:
        panel += u2k_zones[band]['numZones']  # accumulate total panels
        band += 1

    ra -= 12.0 / u2k_zones[band]['numZones']  # offset by half-width of panel
    if ra >= 24.0:  # reality check. shouldn't happen.
        ra -= 24.0

    if ra <= 0.0:  # offset could give negative ra
        ra += 24.0

    if south and u2k_zones[band]['numZones']:
        panel = 222 - panel - u2k_zones[band]['numZones']

    return south+1, panel+int(u2k_zones[band]['numZones']*(24.0 - ra)/24.0)


def sky_quad(coord_obj, lat=39.5):
    """Sky quadrant .
    u2k_chart(coord_obj)

    Parameters
    ----------
    coord_obj : astropy SkyCoord
    lat : float, default 39.5

    Returns
    -------
    quad : str (NP, NQ1, SQ1, NQ2 ... )

    Notes
    -----
    Returns 'NP' for North Circumpolar (determined by lat)
    NQ1, NQ2 ... for northern quads, SQ1, SQ2 ... for southern.
    """
    ra = coord_obj.ra.hour
    dec = coord_obj.dec.deg

    if ra < 0.0 or ra > 24.0 or dec < -90.0 or dec > 90.0:
        raise ValueError("Invalid ra and/or dec")

    if dec > 90-lat:
        return 'NP'

    if dec >= 0.0:
        return 'NQ'+str(ceil(ra/6))

    return 'SQ'+str(ceil(ra/6))
