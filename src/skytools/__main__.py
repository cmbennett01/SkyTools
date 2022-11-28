import argparse

import astropy.units as u
from astropy.coordinates import SkyCoord

from .atlas import u2k_chart, sky_quad

parser = argparse.ArgumentParser()
parser.add_argument("-ra", help="Right ascension in hms")
parser.add_argument("-dec", help="Declination in dms")
parser.add_argument('-n', '--name', help="Get object by name")
args = parser.parse_args()

if args.name:
    c = SkyCoord.from_name(args.name)
else:
    coords = args.ra + args.dec
    c = SkyCoord(coords, unit=(u.deg, u.hourangle))


ra = c.ra.to_string(unit='hour', precision=0, format='unicode') # type: ignore
dec = c.dec.to_string(unit='deg', precision=0, alwayssign=True, # type: ignore
                      format='unicode')
vol, pg = u2k_chart(c)
quad = sky_quad(c)

print(f"RA:{ra} Dec:{dec} Quad:{quad}")
print(f"Uranometria 2000.0: Vol. {vol}, p. {pg}")
