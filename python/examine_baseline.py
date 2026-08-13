# coding: utf-8

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.ndimage import uniform_filter1d

####
"""
unapplied
applied 
rr, rl, lr, ll
of one baseline

data, gains
"""

ODIR   = "gwbexamine19"
PPTAG  = "SCAN_19"
DFPKL  = "scan19.bldata_df.pkl"
GGPKL  = "scan21.bldata_df.pkl"

PPRG   = f"{PPTAG}_r.gains"
PPLG   = f"{PPTAG}_l.gains"
##################################################
df       = pd.read_pickle(DFPKL)
ggf      = pd.read_pickle(GGPKL)

pprgains = pd.read_csv (PPRG,sep='\\s+').map(complex)
pplgains = pd.read_csv (PPLG,sep='\\s+').map(complex)
##################################################
# g for 130 or r
# h for 175 or l
"""
dmn_pq = gm_p * Mmn * conjugate(gn_q) 

ignore parallactic jones at this point
"""
##################################################
if False:
    if not os.path.exists ( ODIR ): os.mkdir ( ODIR )
    fig      = plt.figure('ebaseline')

    for (ant1, ant2), sdf in df.groupby(level=['ant1', 'ant2']):
        ##
        gp   = pprgains[ant1]
        gq   = pprgains[ant2]
        hp   = pplgains[ant1]
        hq   = pplgains[ant2]
        ##
        drr  = sdf.complex.loc[ant1, ant2, 'rr']
        drl  = sdf.complex.loc[ant1, ant2, 'rl']
        dlr  = sdf.complex.loc[ant1, ant2, 'lr']
        dll  = sdf.complex.loc[ant1, ant2, 'll']
        ##
        mrr  = drr / gp / np.conjugate ( gq )
        mrl  = drl / gp / np.conjugate ( hq )
        mlr  = dlr / hp / np.conjugate ( gq )
        mll  = dll / hp / np.conjugate ( hq )
        ##
        drr  = ggf.complex.loc[ant1, ant2, 'rr']
        drl  = ggf.complex.loc[ant1, ant2, 'rl']
        dlr  = ggf.complex.loc[ant1, ant2, 'lr']
        dll  = ggf.complex.loc[ant1, ant2, 'll']
        ##
        ##################################
        ofile= f"{ODIR}/{ant1}_{ant2}.png"
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

        fig.suptitle ( f"{ant1}&{ant2}\nGWB-applied=black applied=blue" )

        fig.savefig (ofile, bbox_inches='tight', dpi=300)
        fig.clf()
