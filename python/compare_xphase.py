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
    agp  = argparse.ArgumentParser ("compare_xphase", description="xphase check")
    add  = agp.add_argument
    add ('dfs', help="output of :xphase_from_gains.py:", nargs=2)
    add ('-l','--label', help='Labels', nargs=2, dest='labels', required=True)
    add ('-o', '--outfile', help='Output file', default=None, dest='ofile')
    add ('-d','--drop', help='Drop the following antennas', nargs='+', dest='badants', default=[])
    return agp.parse_args()


if __name__ == "__main__":
    args    = get_args()
    OFILE   = args.ofile

    xf    = [pd.read_pickle(ifs).drop(index=args.badants) for ifs in args.dfs]

    ants  = sorted(set(xf[0].index))

    COLORS = cycle(['red','green','blue','black','orange','yellow', 'magenta'])

    fig     = plt.figure('xphase',figsize=(8,6))
    axde, axdd = fig.subplots ( 2,1, sharex=True)

    axde.plot ( xf[0].delay, c='red', marker='.', label=args.labels[0] )
    axde.plot ( xf[1].delay, c='blue', marker='.', label=args.labels[1] )

    axde.set_title(f"red={args.labels[0]} blue={args.labels[1]}")

    dd    = xf[0].delay - xf[1].delay

    axdd.plot ( dd, c='k' )

    axde.grid(axis='x', which='major', ls=':', c='k', alpha=0.4)
    axdd.grid(axis='x', which='major', ls=':', c='k', alpha=0.4)

    axde.set_ylabel ('Delay / ns')
    axdd.set_ylabel ('Difference / ns')

    axde.xaxis.tick_top()
    axde.xaxis.set_label_position('top')
    axde.set_xlabel ('Antenna')

    axde.tick_params(axis='x', labelsize='x-small', labelrotation=45)
    axdd.tick_params(axis='x', labelsize='x-small', labelrotation=45)

    axdd.set_xlabel('Antenna')

    if args.ofile:
        fig.savefig (args.ofile, dpi=300, bbox_inches='tight')
    else:
        plt.show ()
