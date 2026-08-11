# coding: utf-8

from tqdm import tqdm

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_args():
    import argparse
    agp = argparse.ArgumentParser("compare_polphase_quartical", description="Comparing polphase and quartical complex solutions")
    add = agp.add_argument
    add ('-p', '--polphase', help='tag of polphase solutions', required=True, dest='pptag')
    add ('-q', '--quartical', help='csv file of quartical solutions :make_quartical_df.py:', required=True, dest='qcsv')
    add ('-O', '--outdir', help='put plots in this directory', default='compare_polphase_quartical_diag', dest='odir')
    return agp.parse_args()

def action ( drr, drl, dlr, dll, mrr, mrl, mlr, mll, ant, ofile ):
    """
    d?? and m?? are complex gains
    """
    mse  = lambda f : np.nanmean ( np.power ( np.abs(f), 2.0 ) )
    err  = mse ( drr - mrr )
    erl  = mse ( drl - mrl )
    elr  = mse ( dlr - mlr )
    ell  = mse ( dll - mll )

    axes = fig.subplots ( 4,2,sharex=True, sharey='col', gridspec_kw={'hspace':0., 'wspace':0.04},  )

    axes[0,0].semilogy ( np.abs(drr), c='k' )
    axes[1,0].semilogy ( np.abs(drl), c='k' )
    axes[2,0].semilogy ( np.abs(dlr), c='k' )
    axes[3,0].semilogy ( np.abs(dll), c='k' )

    axes[0,0].semilogy ( np.abs(mrr), c='b' )
    axes[1,0].semilogy ( np.abs(mrl), c='b' )
    axes[2,0].semilogy ( np.abs(mlr), c='b' )
    axes[3,0].semilogy ( np.abs(mll), c='b' )

    axes[0,1].plot ( np.angle(drr), c='k' )
    axes[1,1].plot ( np.angle(drl), c='k', alpha=0.5 )
    axes[2,1].plot ( np.angle(dlr), c='k', alpha=0.5 )
    axes[3,1].plot ( np.angle(dll), c='k' )

    axes[0,1].plot ( np.angle(mrr), c='b' )
    axes[1,1].plot ( np.angle(mrl), c='b', alpha=0.5 )
    axes[2,1].plot ( np.angle(mlr), c='b', alpha=0.5 )
    axes[3,1].plot ( np.angle(mll), c='b' )

    axes[3,0].set_xlabel ('Channel')
    axes[3,1].set_xlabel ('Channel')

    axes[0,0].set_ylabel ('RR-Amp')
    axes[1,0].set_ylabel ('RL-Amp')
    axes[2,0].set_ylabel ('LR-Amp')
    axes[3,0].set_ylabel ('LL-Amp')

    axes[0,1].set_ylabel ('RR-Phs/rad')
    axes[1,1].set_ylabel ('RL-Phs/rad')
    axes[2,1].set_ylabel ('LR-Phs/rad')
    axes[3,1].set_ylabel ('LL-Phs/rad')

    axes[0,1].yaxis.tick_right()
    axes[0,1].yaxis.set_label_position('right')
    axes[1,1].yaxis.tick_right()
    axes[1,1].yaxis.set_label_position('right')
    axes[2,1].yaxis.tick_right()
    axes[2,1].yaxis.set_label_position('right')
    axes[3,1].yaxis.tick_right()
    axes[3,1].yaxis.set_label_position('right')

    fig.suptitle ( f"{ant}\nquartical=black polphase=blue\nMSE: RR={err:.1e} RL={erl:.1e} LR={elr:.1e} LL={ell:.1e}" )
    # fig.savefig (ofile, bbox_inches='tight', dpi=300)
    # fig.clf()

if __name__ == "__main__":
    args    = get_args()
    PPTAG   = args.pptag
    QSOLS   = args.qcsv
    ODIR    = args.odir

    if not os.path.exists ( ODIR ): os.mkdir ( ODIR )

    ## read polphase
    urr  = pd.read_csv (f"{PPTAG}_rr.gains", sep='\\s+').map(complex)
    url  = pd.read_csv (f"{PPTAG}_rl.gains", sep='\\s+').map(complex)
    ulr  = pd.read_csv (f"{PPTAG}_lr.gains", sep='\\s+').map(complex)
    ull  = pd.read_csv (f"{PPTAG}_ll.gains", sep='\\s+').map(complex)
    ## quartical
    ff   = pd.read_csv ( QSOLS,).set_index(['antenna','correlation']).sort_index()
    ff['gains'] = ff['gains'].apply(complex)
    fflags = ff['gain_flags'] == 1.0
    ff.loc[fflags, 'gains']   = np.nan
    # ff.to_pickle(QSOLF)
    # adfdf
    ########################

    fig  = plt.figure ('polphase_quartical', figsize=(6,8))

    for ant in tqdm (urr.columns, desc='antennas'):
    # for ant in tqdm (['C00'], desc='antennas'):

        ofile = f"{ODIR}/{ant}.png"
        action ( 
            ff.loc[ant,'RR'].gains.values, 
            ff.loc[ant,'RL'].gains.values, 
            ff.loc[ant,'LR'].gains.values, 
            ff.loc[ant,'LL'].gains.values, 
            urr[ant],
            url[ant],
            ulr[ant],
            ull[ant],
            ant,
            ofile
        )
        fig.savefig (ofile, bbox_inches='tight', dpi=300)
        fig.clf()

        # plt.show ()
        # break

