# coding: utf-8

from itertools import cycle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from numpy import conjugate as conj
###############################
"""
use output of :xphase_from_gains.py:
to plot delays against antenna
xphase against channel
residuals against channel
"""
def get_args ():
    import argparse
    agp  = argparse.ArgumentParser ("examine_xphase", description="xphase check")
    add  = agp.add_argument
    add ('xphase_df', help="output of :xphase_from_gains.py:")
    add ('-o', '--outfile', help='Output file', default=None, dest='ofile')
    add ('-d','--drop', help='Drop the following antennas', nargs='+', dest='badants', default=[])
    return agp.parse_args()


if __name__ == "__main__":
    args    = get_args()
    DFPKL   = args.xphase_df
    OFILE   = args.ofile

    xf    = pd.read_pickle(DFPKL)
    xf    = xf.drop(index=args.badants)

    ants  = sorted(set(xf.index))

    COLORS = cycle(['red','green','blue','black','orange','yellow', 'magenta'])

    fig     = plt.figure('xphase',figsize=(8,6))
    axde, axxp, axrs = fig.subplots ( 3,1, )

    axbi = axde.twinx()


    for iant,color in zip(ants,COLORS):
        axxp.plot ( xf.loc[iant].dphase, color=color, label=iant )
        axrs.plot ( xf.loc[iant].res, color=color, label=iant )

    axde.plot ( xf.delay, c='k', marker='.' )

    axbi.plot ( xf.bias, c='b', marker='.' )

    axde.grid(axis='x', which='major', ls=':', c='k', alpha=0.4)

    axde.set_ylabel ('Delay / ns')
    axbi.set_ylabel ('Bias / ns')

    axde.xaxis.tick_top()
    axde.xaxis.set_label_position('top')
    axde.set_xlabel ('Antenna')

    axbi.yaxis.label.set_color('blue')
    axbi.tick_params(axis='y', color='blue')

    axxp.set_ylabel ('Phase / rad')
    axrs.set_ylabel ('Residual / rad')
    axrs.set_xlabel ('Channel')

    axrs.sharex(axxp)

    axrs.axhline(0.,ls=':', c='k')

    axde.tick_params(axis='x', labelsize='x-small', labelrotation=45)

    if args.ofile:
        fig.savefig (args.ofile, dpi=300, bbox_inches='tight')
    else:
        plt.show ()
