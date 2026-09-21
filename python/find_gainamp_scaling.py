# coding: utf-8

import numpy as np
import pandas as pd

def get_args():
    import argparse
    agp = argparse.ArgumentParser('find_gainamp_scaling')
    add = agp.add_argument
    add('selfcorr', help='Output of :write_self:')
    add ('-f','--full', help='Fullpolphase solution df', dest='sols', required=True)
    add ('-a','--atol', help='Absolute tolerance wrt 0.0', default=0.01, type=float, dest='atol')
    add ('--replace', help='Replace inverse of small gain amps with', default=100., type=float, dest='rval')
    add ('-r','--reference', help='Reference value of selfcorr', type=float, default=100., dest='ref')
    return agp.parse_args()

if __name__ == "__main__":
    args = get_args()

    SELF = args.selfcorr
    SOLS = args.sols
    ATOL = args.atol
    RVAL = args.rval
    BVAL = args.ref

    ## load the self baselines
    ss = pd.read_csv(SELF, sep='\\s+', usecols=['ant','corr','real']).groupby(['ant','corr']).agg(list)
    ss['real'] = ss['real'].apply(np.array)

    ## load the solved solutions
    sol = pd.read_pickle(SOLS)

    ## compute inverse amplitudes
    sol['amp']  = sol['gains'].apply(np.abs)
    sol['iamp'] = sol['amp'].apply(lambda f : np.where(np.isclose(f,0.,atol=ATOL),RVAL, 1.0/f))

    ## estimate the scaling
    fac    = ss.real * sol.iamp * sol.iamp
    val    = np.nanmedian(fac.loc[:,['rr','ll']].apply(np.nanmedian))

    sac    = np.sqrt ( BVAL / val )

    print (f" Amplitude scaling := {sac:.5f}")

