"""
reference impl to check 
"""

import numpy as np
import pandas as pd

import astropy.time  as at
import astropy.units as au
import astropy.coordinates as asc

SOURCES   = {
    '3C138': asc.SkyCoord (79.5687917*au.deg, 16.5907806*au.deg),
    '3C147': asc.SkyCoord (84.6812917*au.deg, 49.8285556*au.deg),
    '3C48' : asc.SkyCoord (24.4220417*au.deg, 33.1597417*au.deg)
}

SCHOICES  = list(SOURCES.keys())

######################
## taken from observatories.dat
## from tempo2
gmrt_ant_x         = 1656342.30
gmrt_ant_y         = 5797947.77
gmrt_ant_z         = 2073243.16
gmrt_el            = asc.EarthLocation.from_geocentric (gmrt_ant_x, gmrt_ant_y, gmrt_ant_z, unit="m")
######################
singles      = pd.read_csv ("casa_single_antpos.txt", names=[
  "idx", "name", "oname", "radius", "ru", "long", "lat", "bx", "by", "bz", "gx", "gy", "gz"
], sep='\\s+')
singles_el   = asc.EarthLocation.from_geocentric ( singles.gx, singles.gy, singles.gz, unit='m' )
singles['lat_deg'] = singles_el.lat.deg
singles['lon_deg'] = singles_el.lon.deg
######################

def write_code():
    """
    to write code
    """
    print (" GMRT POS")
    print ( "LONGITUDE_DEG", gmrt_el.lon.deg, sep='\t' )
    print ( "LATITUDE_DEG", gmrt_el.lat.deg, sep='\t' )

    print (" GMRT SINGLE ANTENNA POS")

    for name, lon, lat in zip ( singles['oname'], singles['lon_deg'], singles['lat_deg'] ):
        print (f"{{ antname_t{{\"{name.split(':')[0]}\"}}, antpos_t({lon:.6f}, {lat:.6f}) }},")

def get_args ():
    import argparse
    agp  = argparse.ArgumentParser("gmrt_coordinates", description='Computes parallactic angles and such')
    add  = agp.add_argument
    add ('-s','--source', help='Source', choices=SCHOICES, required=True, dest='source')
    add ('-t','--mjd', help='MJD', dest='mjd', type=float, required=True)
    return agp.parse_args()

def get_sla_lst (mjd, longitude):
    """
    returns Quantity

    from slalib
    """

    tu  =  ( mjd - 51544.5 ) / 36525.0

    lst =  ( np.mod ( mjd, 1.0 ) * 2.0 * np.pi + \
            ( 24110.54841 + ( 8640184.812866 +   \
                ( 0.093104 - 6.2e-6 * tu ) * tu ) * tu ) * ( np.pi/12.0/3600.0)
    )

    lst = np.mod ( lst, 2.0*np.pi )

    return asc.Angle(lst*au.radian) + longitude

def get_parallactic_angle ( sc, tobs, el ):
    """ source coordinates, mjd --> parallactic angle (degree) 
        Taken from :GMRT-FRB/get_par_angle.py:
    """

    lst   = tobs.sidereal_time ( 'mean', longitude=el.lon, model=None )
    lst2  = get_sla_lst ( tobs.mjd, el.lon )
    print ( lst.to_string(au.hourangle), lst2.to_string(au.hourangle), sep='\t' )
    h     = (lst - sc.ra).radian
    print ( lst.radian, h, sep='\t' )
    q     = np.arctan2 ( 
            np.sin ( h ), 
            np.tan ( el.lat.radian ) * np.cos ( sc.dec.radian ) - 
            np.sin ( sc.dec.radian ) * np.cos ( h )
    )
    return np.rad2deg(q)

if __name__ == "__main__":
    ### write_code only when we need to update code
    write_code()
    import sys
    sys.exit(0)
    ### write_code only when we need to update code

    args = get_args()

    src  = SOURCES [ args.source ]
    tobs = at.Time ( args.mjd, format='mjd' )

    ist  = tobs + 5.5*au.hour

    print ( f" At Time [UTC] = {tobs.isot}" )
    print ( f" At Time [IST] = {ist.isot}" )

    ###
    par  = get_parallactic_angle ( src, tobs, gmrt_el )
    print ( f"GMRT   {par:.3f} deg" )
    ###
    # for name, el in zip(singles['name'], singles_el):
        # par  = get_parallactic_angle ( src, tobs, el )
        # print ( f"{name}    {par:.3f} deg" )


