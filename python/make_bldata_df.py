# coding: utf-8
import numpy as np
import pandas as pd

freqs   = np.linspace ( 550., 750., 2048, endpoint=True )
freqs_ghz  = freqs * 1E-3

def get_args():
    import argparse
    agp = argparse.ArgumentParser("make_bldata_df", description="Makes a grouped dataframe out of bldata")
    add = agp.add_argument
    add ('bldata', help='Output of :write_baselines:')
    return agp.parse_args()


def measure_d ( freqs_ghz, dphi, dmin=-300, dmax=300, dsize=1024 ):
    """
    measures delay bias and residuals
    """
    delays_grid = np.linspace ( dmin, dmax, dsize )
    dmags       = np.zeros_like ( delays_grid )
    for i, idelay in enumerate ( delays_grid ):
        dmags [ i ]  = np.abs ( np.nanmean ( np.exp ( 1.0j * ( dphi - ( 2.0 * np.pi * idelay * freqs_ghz ) ) ) ) )
    ####
    delay_ns    = delays_grid [ np.argmax ( dmags ) ]
    biaser      = np.exp ( 1.0j * ( dphi - ( 2.0 * np.pi * delay_ns * freqs_ghz ) ) )
    bias_rad    = 1.0 * np.angle ( np.nanmean ( biaser ) )
    model       = np.angle ( np.exp ( 1.0j * ( ( 2.0 * np.pi * delay_ns * freqs_ghz ) + bias_rad ) ) ) 
    residuals   = np.angle ( np.exp ( 1.0j * (dphi -  ( 2.0 * np.pi * delay_ns * freqs_ghz ) - 1.0 * bias_rad  ) ) )
    return {'delay':delay_ns, 'bias':bias_rad, 'res':residuals, 'model':model, 'delays':delays_grid, 'mag':dmags, 'dphi':dphi}

if __name__ == "__main__":
    args = get_args ()
    bf   = pd.read_csv(args.bldata,sep='\\s+')
    bf['complex'] = bf['complex'].apply(complex)
    gf = bf[['ant1','ant2','correlation','complex']].groupby(['ant1','ant2','correlation']).agg(list)

    # gf['delay_ns'] = 0.
    # gf['bias_rad'] = 0.

    # for (ant1,ant2,corr), row in gf.iterrows():
        # _d = measure_d ( freqs_ghz, np.angle(row.complex) )
        # gf.loc[(ant1, ant2, corr), "delay_ns"] = _d['delay']
        # gf.loc[(ant1, ant2, corr), "bias_rad"] = _d['bias']

    gf.to_pickle (args.bldata + "_df.pkl")
