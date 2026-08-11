"""

crosshand phase from noise diode

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm

from collections import defaultdict

def get_args():
    import argparse
    agp = argparse.ArgumentParser("xphase_from_noisef", description="Computes crosshand phases and residuals and delays")
    add = agp.add_argument
    add ('on', help=':write_baseline: of ON scan', )
    add ('off', help=':write_baseline: of OFF scan', )
    add ('-o', '--ofile', help='xphase dataframe', required=True, dest='ofile')
    return agp.parse_args()

def make_ngf ( f ):

    ff        = pd.read_csv (f,sep='\\s+')
    selfrow   = (ff['ant1'] == ff['ant2'])
    ff        = pd.DataFrame ( ff[['ant1','correlation','complex']][selfrow] )
    ff['complex'] = ff['complex'].apply(complex)
    gf        = ff.groupby(['ant1','correlation']).agg(list)
    gf['complex'] = gf['complex'].apply ( np.array )
    return gf

if __name__ == "__main__":
    args    = get_args()
    odf     = args.ofile
    ###################
    on      = make_ngf ( args.on )
    of      = make_ngf ( args.off )
    ### beauty of pandas
    oo      = on - of
    ###################
    ants    = sorted ( set(oo.index.get_level_values(0)) )

    freqs   = np.linspace ( 550., 750., 2048, endpoint=True )
    freqs_ghz  = freqs * 1E-3

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
        dphi   = np.angle ( oo.complex.loc[ant, 'rl'] )
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

