# coding: utf-8

from tqdm import tqdm

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#############################################

def get_args():
    import argparse
    agp = argparse.ArgumentParser("compare_polphase_quartical_diag", description="Comparing polphase and quartical diag_complex solutions")
    add = agp.add_argument
    add ('-p', '--polphase', help='tag of polphase solutions', required=True, dest='pptag')
    add ('-q', '--quartical', help='csv file of quartical solutions :make_quartical_df.py:', required=True, dest='qcsv')
    add ('-r', '--refant', help='Reference antenna', default='C00', dest='refant')
    add ('-O', '--outdir', help='put plots in this directory', default='compare_polphase_quartical_diag', dest='odir')
    return agp.parse_args()

# PPTAG = "full_jones_par/3C138_cpar_q"
# QSOLS = "quartical_3C138_full_gains.csv"
# ODIR  = "compare_polphase_quartical_3C138_cpar"
#############################################

def action ( drr, dll, mrr, mll, nrr, nll, ant, ofile ):
    """
    d?? and m?? are complex gains

    d?? are the quartical ones

    """
    axes = fig.subplots ( 2,2,sharex=True, sharey='col', gridspec_kw={'hspace':0., 'wspace':0.04},  )

    axes[0,0].semilogy ( np.abs(drr), c='k' )
    axes[1,0].semilogy ( np.abs(dll), c='k' )

    axes[0,0].semilogy ( np.abs(mrr), c='b' )
    axes[1,0].semilogy ( np.abs(mll), c='b' )

    axes[0,1].plot ( np.angle(drr), c='k', alpha=0.5 )
    axes[1,1].plot ( np.angle(dll), c='k', alpha=0.5 )

    axes[0,1].plot ( np.angle(mrr), c='b', alpha=0.5, )
    axes[1,1].plot ( np.angle(mll), c='b', alpha=0.5,)

    axes[0,1].plot ( np.angle(nrr), c='g', alpha=0.5 )
    axes[1,1].plot ( np.angle(nll), c='g', alpha=0.5 )

    axes[1,0].set_xlabel ('Channel')
    axes[1,1].set_xlabel ('Channel')

    axes[0,0].set_ylabel ('RR-Amp')
    axes[1,0].set_ylabel ('LL-Amp')

    axes[0,1].set_ylabel ('RR-Phs/rad')
    axes[1,1].set_ylabel ('LL-Phs/rad')

    axes[0,1].yaxis.tick_right()
    axes[0,1].yaxis.set_label_position('right')
    axes[1,1].yaxis.tick_right()
    axes[1,1].yaxis.set_label_position('right')

    fig.suptitle ( f"{ant}\nquartical=black\npolphase(ref={refant})=green polphase(noref)=blue" )
    # fig.suptitle ( f"{ant}\nquartical=black polphase=blue" )

if __name__ == "__main__":
    args    = get_args()
    PPTAG   = args.pptag
    QSOLS   = args.qcsv
    ODIR    = args.odir
    refant  = args.refant

    if not os.path.exists ( ODIR ): os.mkdir ( ODIR )

    ## read polphase
    urr  = pd.read_csv (f"{PPTAG}_rr.gains", sep='\\s+').map(complex)
    ull  = pd.read_csv (f"{PPTAG}_ll.gains", sep='\\s+').map(complex)
    ## reference polphase
    refrr  = np.exp(1.0j*np.angle(urr[refant]))
    refll  = np.exp(1.0j*np.angle(ull[refant]))
    ## quartical
    ff   = pd.read_csv ( QSOLS,).set_index(['antenna','correlation']).sort_index()
    ff['gains'] = ff['gains'].apply(complex)
    fflags = ff['gain_flags'] == 1.0
    ff.loc[fflags, 'gains']   = np.nan
    # ff.to_pickle(QSOLF)
    # adfdf
    ########################

    fig  = plt.figure ('polphase_quartical', figsize=(5,5))

    for ant in tqdm (urr.columns, desc='antennas'):
    # for ant in tqdm (['C02'], desc='antennas'):

        ofile = f"{ODIR}/{ant}.png"
        action ( 
            ff.loc[ant,'RR'].gains.values, 
            ff.loc[ant,'LL'].gains.values, 
            urr[ant],
            ull[ant],
            urr[ant] / refrr,
            ull[ant] / refll,
            ant,
            ofile
        )

        fig.savefig (ofile, bbox_inches='tight', dpi=300)
        fig.clf()

        # plt.show ()
        # break


