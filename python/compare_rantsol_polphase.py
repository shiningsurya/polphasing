
from tqdm import tqdm

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#############################################
refant = "C00"
SRC   = "3C138"
PPTAG = "psols/pp_3323_0"
nchans= 2048
RSOLS = "rantsols/polphasing-main/rantsols/0"
QSOLS = f"qq_{SRC}_diag.csv"
QSOLF = f"qq_{SRC}_diag_df.pkl"
ODIR  = f"compare_rantsol_polphase_3323"
#############################################
## read rantsol
rslgains = pd.read_csv (RSOLS+"_175.dat",sep='\\s+', nrows=nchans).shift(1,axis='columns').map(lambda x : np.exp(1.0j*np.deg2rad(x))).drop(columns=['#'])
rsrgains = pd.read_csv (RSOLS+"_130.dat",sep='\\s+', nrows=nchans).shift(1,axis='columns').map(lambda x : np.exp(1.0j*np.deg2rad(x))).drop(columns=['#'])
## read polphase
pplgains = pd.read_csv (PPTAG+"_ll.gains", sep='\\s+').map(complex)
pprgains = pd.read_csv (PPTAG+"_rr.gains", sep='\\s+').map(complex)
## reference polphase
refll  = np.exp(1.0j*np.angle(pplgains[refant]))
refrr  = np.exp(1.0j*np.angle(pprgains[refant]))
##
if not os.path.exists ( ODIR ): os.mkdir ( ODIR )
#############################################
"""
two columns: L and R

focus on phase difference

rantsol and ref-polphase

polphase
"""

fig  = plt.figure ('polphase_rantsol', figsize=(5,5))

for ant in tqdm (pprgains.columns, desc='antennas'):
# for ant in ['C02']:

    ofile = os.path.join ( ODIR, f"{ant}.png" )
    axes = fig.subplots ( 2,2,sharex=True, sharey='col', gridspec_kw={'hspace':0., 'wspace':0.04}, )

    ## left

    axes[0,0].plot ( np.angle(rslgains[ant]), c='k', alpha=0.5 )
    axes[0,0].plot ( np.angle(refll*np.conjugate(pplgains[ant])), c='b', alpha=0.5 )

    axes[0,1].plot ( np.angle(rsrgains[ant]), c='k', alpha=0.5 )
    axes[0,1].plot ( np.angle(refrr*np.conjugate(pprgains[ant])), c='b', alpha=0.5 )

    axes[1,0].plot ( np.angle(np.conjugate(pplgains[ant])), c='b' )
    axes[1,1].plot ( np.angle(np.conjugate(pprgains[ant])), c='b' )

    axes[0,1].yaxis.tick_right()
    axes[0,1].yaxis.set_label_position('right')
    axes[1,1].yaxis.tick_right()
    axes[1,1].yaxis.set_label_position('right')

    axes[1,0].set_xlabel ('Channel')
    axes[1,1].set_xlabel ('Channel')

    axes[0,0].set_title ("L")
    axes[0,1].set_title ("R")

    axes[0,0].set_ylabel ('Phase / rad')
    axes[1,0].set_ylabel ('Phase / rad')

    axes[0,1].set_ylabel ('Phase / rad')
    axes[1,1].set_ylabel ('Phase / rad')

    fig.suptitle(f"{ant} refant={refant}\nrantsol=black polphase=blue")

    # plt.show ()
    fig.savefig (ofile, bbox_inches='tight', dpi=300)
    fig.clf()
