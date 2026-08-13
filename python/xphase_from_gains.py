"""
complex gains we solve introduce cross hand phase

there is a phase ramp or a linear slope, which we can subtract out

whatever remains is technically not leakage but noise in the system 

we see this noise

this script just computes. plot using other script

this computation is pretty heavy
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm

from collections import defaultdict
###############################
def get_args():
    import argparse
    agp = argparse.ArgumentParser("xphase_from_gains", description="Computes crosshand phases and residuals and delays")
    add = agp.add_argument
    add ('tag', help="Stem/tag that was passed to polphase")
    return agp.parse_args()


if __name__ == "__main__":
    args    = get_args()
    stem    = args.tag
    odf     = stem + "_xphase_df.pkl"

    rgains  = pd.read_csv (f"{stem}_r.gains", sep='\\s+').map(complex)
    lgains  = pd.read_csv (f"{stem}_l.gains", sep='\\s+').map(complex)

    mm      = np.abs( rgains.sum(0) ) == 0.
    goodants =  list(mm[~mm].index)

    freqs   = np.linspace ( 550., 750., 2048, endpoint=True )
    freqs_ghz  = freqs * 1E-3

    # ants    = list(set(rgains.columns).intersection(lgains.columns))
    ants    = goodants

    def measure_d ( freqs_ghz, dphase, dmin=-300, dmax=300, dsize=1024 ):
        """
        measures delay bias and residuals
        """
        delays_grid = np.linspace ( dmin, dmax, dsize )
        dmags       = np.zeros_like ( delays_grid )
        for i, idelay in enumerate ( delays_grid ):
            dmags [ i ]  = np.abs ( np.nanmean ( np.exp ( 2.0j * ( dphase - ( 0.5 * np.pi * idelay * freqs_ghz ) ) ) ) )
        ####
        delay_ns    = delays_grid [ np.argmax ( dmags ) ]
        biaser      = np.exp ( 2.0j * ( dphase - ( 0.5 * np.pi * delay_ns * freqs_ghz ) ) )
        bias_rad    = 1.0 * np.angle ( np.nanmean ( biaser ) )
        model       = 0.5 * np.angle ( np.exp ( 2.0j * 0.5 * ( ( np.pi * delay_ns * freqs_ghz ) + bias_rad ) ) ) 
        residuals   = 0.5 * np.angle ( np.exp ( 2.0j * (dphase -  ( 0.5 * np.pi * delay_ns * freqs_ghz ) - 0.5 * bias_rad  ) ) )
        return {'delay':delay_ns, 'bias':bias_rad, 'res':residuals, 'model':model, 'delays':delays_grid, 'mag':dmags, 'dphase':dphase}

    dps   = defaultdict(list)
    for ant in tqdm (ants, desc='antennas', unit='ant'):
        dphi   = np.angle ( rgains[ant] / lgains[ant] )
        dphase = 0.5 * dphi
        sols   = measure_d ( freqs_ghz, dphase )
        #
        dps['ant'].append ( ant )
        dps['dphase'].append ( dphase )
        dps['dphi'].append ( dphi )
        dps['mag'].append ( sols['mag'] )
        dps['model'].append ( sols['model'] )
        dps['delay'].append ( sols['delay'] )
        dps['bias'].append ( sols['bias'] )
        dps['res'].append ( sols['res'] )

    dps   = pd.DataFrame ( dps ).set_index('ant').sort_index()
    dps.to_pickle(odf)

