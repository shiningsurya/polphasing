
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


X = sys.argv[1]

df = pd.read_csv (X, sep='\\s+', names=['chan', 'ant1', 'ant2', 'corr', 'complex']).set_index(['corr','chan'])
df['complex'] = df['complex'].apply(complex)

drr  = df['complex'].loc['rr']
drl  = df['complex'].loc['rl']
dlr  = df['complex'].loc['lr']
dll  = df['complex'].loc['ll']

fig = plt.figure ('examine_one')

axes = fig.subplots ( 4,2,sharex=True, sharey='col', gridspec_kw={'hspace':0., 'wspace':0.04},  )

axes[0,0].semilogy ( np.abs(drr), c='k' )
axes[1,0].semilogy ( np.abs(drl), c='k' )
axes[2,0].semilogy ( np.abs(dlr), c='k' )
axes[3,0].semilogy ( np.abs(dll), c='k' )

axes[0,1].plot ( np.angle(drr), c='k' )
axes[1,1].plot ( np.angle(drl), c='k' )
axes[2,1].plot ( np.angle(dlr), c='k' )
axes[3,1].plot ( np.angle(dll), c='k' )

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

fig.suptitle (os.path.basename(X))

plt.show ()


