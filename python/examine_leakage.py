# coding: utf-8

from itertools import cycle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_args ():
    import argparse
    agp  = argparse.ArgumentParser ("examine_leakage", description="leakage check")
    add  = agp.add_argument
    add ('leakage_df', help="output of :make_solutions_df.py -l:")
    add ('-o', '--outfile', help='Output file', default=None, dest='ofile')
    add ('-a','--ants', help='Show the antenna', nargs='+', dest='ants', default=[])
    return agp.parse_args()

if __name__ == "__main__":
    args = get_args()

    det  = pd.read_pickle(args.leakage_df)

    ants = args.ants
    if len(ants) == 0:
        ants = sorted(set(det.index.get_level_values('ant')))
    ###################

    fig = plt.figure ('3345')

    axes = fig.subplots ( 2,2, sharex=True, sharey='col', gridspec_kw={'hspace':0., 'wspace':0.} )

    COLORS = cycle(['red','green','blue','black','orange','yellow', 'magenta'])
    for ant,color in zip(ants, COLORS):

        drl = det.gains.loc[ant, 'drl']
        dlr = det.gains.loc[ant, 'dlr']

        """
        rl: amp, phase
        lr: amp, phase
        """

        axes[0,0].semilogy ( np.abs(drl), c=color )
        axes[0,1].plot ( np.angle(drl), c=color )

        axes[1,0].semilogy ( np.abs(dlr), c=color )
        axes[1,1].plot ( np.angle(dlr), c=color )

        axes[1,0].set_xlabel ('Channel')
        axes[1,1].set_xlabel ('Channel')

        axes[0,0].set_ylabel('RL Amp')
        axes[0,1].set_ylabel('RL Phs / rad')

        axes[1,0].set_ylabel('LR Amp')
        axes[1,1].set_ylabel('LR Phs / rad')

        axes[0,1].yaxis.tick_right()
        axes[1,1].yaxis.tick_right()

        axes[0,1].yaxis.set_label_position('right')
        axes[1,1].yaxis.set_label_position('right')


    if args.ofile:
        fig.savefig (args.ofile, dpi=300, bbox_inches='tight')
    else:
        plt.show ()
