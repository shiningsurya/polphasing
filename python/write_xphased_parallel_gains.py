# coding: utf-8

import numpy as np
import pandas as pd

ORDER_ANTS = "C00 C01 C02 C03 C04 C05 C06 C08 C09 C10 C11 C12 C13 C14 E02 E03 E04 E05 E06 S01 S02 S03 S04 S06 W01 W02 W03 W04 W05 W06 C07 S05".split()

def get_args():
    import argparse
    agp = argparse.ArgumentParser('xphased_gains')
    add = agp.add_argument
    add('soldf', help='Solutions df (rr,ll)  we need normalized solutions')
    add ('-t','--tag', help='Output tag', dest='tag', required=True)
    add ('-x','--xphase', help='Xphase_df, output of :xphase_from_noisef:', required=True, dest='xphase')
    return agp.parse_args()

def write_gains ( ofile, df, corr, xphase=None, nchan=2048 ):
    """
    ofile and dataframe
    """
    with open(ofile,'w') as f:
        ## write header
        line = " ".join(ORDER_ANTS)
        line += "\n"
        f.write (line)
        ##
        for ichan in range(nchan):
            line = ""
            for ant in ORDER_ANTS:
                g  = df.gains.loc[ant, corr][ichan] 
                if xphase is not None:
                    _xp  = 1.0
                    try:
                        _xp = np.exp ( 1.0j * xphase['dphi'].loc[ant][ichan] )
                    except:
                        ## probably xphase does not have that antenna
                        _xp = 1.0
                    g = g * _xp
                line += f"{g:+.3f} "
            line += "\n"
            f.write (line)
        ##

if __name__ == "__main__":
    args = get_args()
    ###
    sf   = pd.read_pickle( args.soldf )
    otag = f"{args.tag}"
    ###
    xf   = pd.read_pickle ( args.xphase )
    ###
    write_gains ( otag+"_rr.gains", sf, 'rr', xphase=xf )
    write_gains ( otag+"_ll.gains", sf, 'll', xphase=None)




