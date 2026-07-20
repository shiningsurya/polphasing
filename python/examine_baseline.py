# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.ndimage import uniform_filter1d

model  = pd.read_csv ("/home/shining/shit/gmrt_phase_test/polphasing/code/polphasing/models/3C138_2048_97656.25_550000000_model.txt", sep='\\s+')
model['rr']  = model['stokes_i']
model['ll']  = model['stokes_i']
model['lr']  = model['stokes_q'] - 1.0j*model['stokes_u']
model['rl']  = model['stokes_q'] + 1.0j*model['stokes_u']

DFPKL  = "xscan12_df.pkl"
PPRG   = "SCAN_12_r.gains"
PPLG   = "SCAN_12_l.gains"
RSRG   = "rantsols/12_130.dat"
RSLG   = "rantsols/12_175.dat"
##################################################
cf       = pd.read_pickle(DFPKL).set_index(['ant1','ant2'])

rslgains = pd.read_csv (RSLG,sep='\\s+', nrows=2048).shift(1,axis='columns').map(lambda x : np.exp(1.0j*np.deg2rad(x))).drop(columns=['#'])
rsrgains = pd.read_csv (RSRG,sep='\\s+', nrows=2048).shift(1,axis='columns').map(lambda x : np.exp(1.0j*np.deg2rad(x))).drop(columns=['#'])

pprgains = pd.read_csv (PPRG,sep='\\s+').map(complex)
pplgains = pd.read_csv (PPLG,sep='\\s+').map(complex)
##################################################
ant1,ant2= 'C00','C02'
# g for 130 or r
# h for 175 or l
def get_baseline ( ant1, ant2 ):
    gp,gq,hp,hq  = pprgains[ant1],pprgains[ant2],pplgains[ant1],pplgains[ant2]
    __row        = cf.loc[ant1, ant2]
    vpq_rr       = __row.rr
    vpq_rl       = __row.rl
    vpq_lr       = __row.lr
    vpq_ll       = __row.ll
    return vpq_rr, vpq_rl, vpq_lr, vpq_ll, gp, gq, hp, hq
####
rr, rl, lr, ll, gp, gq, hp, hq = get_baseline ( 'C00', 'C02' )
##################################################

import sys
sys.exit(0)
fig = plt.figure ('whypolphase-weird')
plt.plot ( np.angle(cf.rr.loc['C08','E03']) )
plt.plot ( np.angle(pprgains['C08']*cf.rr.loc['C08','E03']*np.conjugate(pprgains['E03']) ) )
plt.clf()
plt.plot ( np.angle(pprgains['C08']*cf.rr.loc['C08','E03']*np.conjugate(pprgains['E03']) ) )
plt.plot ( np.angle(rsrgains['C08']*cf.rr.loc['C08','E03']*np.conjugate(rsrgains['E03']) ) )
plt.clf()
plt.plot ( np.angle(rsrgains['C08']*cf.rr.loc['C08','E03']*np.conjugate(rsrgains['E03']) ) )
plt.clf()
plt.plot ( np.angle(rsrgains['C08']*cf.rl.loc['C08','E03']*np.conjugate(rslgains['E03']) ) )
plt.plot ( np.angle(pprgains['C08']*cf.rl.loc['C08','E03']*np.conjugate(pplgains['E03']) ) )
plt.clf()
plt.plot ( np.angle(pprgains['C08']*cf.rl.loc['C08','E03']*np.conjugate(pplgains['E03']) ) )
plt.clf()
plt.plot ( np.angle(rsrgains['C08']*cf.rl.loc['C08','E03']*np.conjugate(rslgains['E03']) ) )
plt.plot ( np.angle(rslgains['C08']*cf.lr.loc['C08','E03']*np.conjugate(rsrgains['E03']) ) )
plt.clf()
plt.plot ( np.angle(rslgains['C08']*cf.lr.loc['C08','E03']*np.conjugate(rsrgains['E03']) ) )
plt.clf()
plt.plot ( np.angle(rslgains['C08']*cf.lr.loc['C08','E03']*np.conjugate(rsrgains['E03']) ) )
plt.clf()
plt.clf()
plt.clf()
plt.plot ( np.angle(cf.rr.loc['C08','E03']) )
plt.plot ( np.angle(cf.rr.loc['C00','C01']) )
plt.plot ( np.angle(cf.rr.loc['C00','C02']) )
plt.clf()
plt.plot ( np.angle(cf.rl.loc['C00','C02']) )
plt.plot ( np.angle(cf.lr.loc['C00','C02']) )
plt.clf()
plt.plot ( np.angle(cf.lr.loc['C00','C02']) )
plt.plot ( np.angle(cf.rl.loc['C00','C02']) )
plt.clf()
plt.plot ( np.angle(cf.rl.loc['C00','C02']) )
plt.plot ( np.angle(rsrgains['C00']*cf.rl.loc['C00','C02']*np.conjugate(rslgains['C02']) ) )
plt.plot ( np.angle(pprgains['C00']*cf.rl.loc['C00','C02']*np.conjugate(pplgains['C02']) ) )
plt.clf()
plt.plot ( np.angle(pprgains['C00']*cf.rl.loc['C00','C02']*np.conjugate(pplgains['C02']) ) )
plt.plot ( np.angle(rsrgains['C00']*cf.rr.loc['C00','C02']*np.conjugate(rsrgains['C02']) ) )
plt.clf()
plt.plot ( np.angle(rsrgains['C00']*cf.rr.loc['C00','C02']*np.conjugate(rsrgains['C02']) ) )
plt.plot ( np.angle(rsrgains['C02']*cf.rr.loc['C00','C02']*np.conjugate(rsrgains['C00']) ) )
plt.clf()
plt.plot ( np.angle(rsrgains['C00']*cf.rr.loc['C00','C02']*np.conjugate(rsrgains['C02']) ) )
plt.clf()
plt.plot ( np.abs(rsrgains['C00']*cf.rr.loc['C00','C02']*np.conjugate(rsrgains['C02']) ) )
plt.plot ( np.abs(pprgains['C00']*cf.rr.loc['C00','C02']*np.conjugate(pprgains['C02']) ) )
plt.clf()
bl = cf.rr.loc['C00','C02']
vpq = cf.rr.loc['C00','C02']
gp,gq  = ppgains['C00'],ppgains['C02']
gp,gq  = pprgains['C00'],pprgains['C02']
plt.plot ( np.abs(vpq/gp/gq) )
plt.yscale('log')
plt.plot ( np.abs(vpq/gp/np.conjugate(gq)) )
plt.clf()
plt.plot ( np.abs(vpq/gp/np.conjugate(gq)) )
plt.yscale('log')
plt.clf()
plt.plot ( np.angle(vpq/gp/np.conjugate(gq)) )
model.head()
plt.clf()
plt.plot ( model['stokes_i'] )
plt.clf()
plt.plot ( np.angle(gp*model['stokes_i']*np.conjugate(gq)) )
plt.plot ( np.angle(vpq) )
plt.plot ( np.angle(np.conjugate(vpq)) )
plt.plot ( np.angle(gp*model['stokes_i']*np.conjugate(gq)) )
plt.clf()
plt.plot ( np.angle(hp*model['stokes_i']*np.conjugate(hq)) )
plt.plot ( np.angle(vpq_ll) )
plt.plot ( np.angle(np.conjugate(vpq_ll)) )
plt.clf()
plt.plot ( np.angle(gp*(model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(gp*(model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(hp*(model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(gq)) )
plt.plot ( np.angle(hq*(model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(gp)) )
plt.clf()
plt.plot ( np.angle(model['stokes_q'] - 1.0j*model['stokes_u']) )
plt.plot ( np.angle(model['stokes_q'] + 1.0j*model['stokes_u']) )
plt.clf()
plt.plot ( np.angle(hp) )
plt.plot ( np.angle(gq) )
plt.clf()
plt.plot ( np.angle(hp*(model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(gq)) )
plt.plot ( np.angle(vpq_lr) )
plt.plot ( np.angle(vpq_rl) )
plt.plot ( np.angle(hp*(model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(gq)) )
plt.clf()
plt.plot ( np.angle(hp*(model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(gq)) )
plt.plot ( np.angle(vpq_rl) )
plt.plot ( np.angle(hp*(model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(gq)) )
plt.plot ( np.angle(gp*(model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(gp*(model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(gp*(model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(np.conjugate(vpq_rl)) )
plt.clf()
plt.plot ( np.angle(np.conjugate(vpq_rl)) )
plt.plot ( np.angle(gp*(model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(gp*(model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(gp*(model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.clf()
plt.plot ( np.angle(gp*(model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(np.conjugate(vpq_rl)) )
plt.plot ( np.angle(np.conjugate(vpq_lr)) )
plt.clf()
plt.plot ( np.angle(np.conjugate(vpq_rl)) )
plt.plot ( np.angle(hq*(model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(gp)) )
plt.plot ( np.angle(hq*(model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(gp)) )
plt.plot ( np.angle(vpq_rl) )
plt.plot ( np.angle(hq*(model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(gp)) )
plt.clf()
plt.plot ( uniform_filter1d(np.angle(vpq_rl)) )
plt.plot ( uniform_filter1d(np.angle(vpq_rl),8) )
plt.plot ( uniform_filter1d(np.angle(vpq_rl),4) )
plt.plot ( uniform_filter1d(np.angle(vpq_rl),8) )
plt.plot ( np.angle(gp*(model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(gp*(model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.clf()
plt.plot ( uniform_filter1d(np.angle(np.conjugate(vpq_rl)),8) )
plt.plot ( uniform_filter1d(np.angle(np.conjugate(vpq_rl)),8) , c='k')
plt.plot ( np.angle(gp*(model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(gp*(model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(gp*(-model['stokes_q'] - 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.plot ( np.angle(gp*(-model['stokes_q'] + 1.0j*model['stokes_u'])*np.conjugate(hq)) )
plt.clf()
plt.plot ( uniform_filter1d(np.abs(np.conjugate(vpq_rl)),8) , c='k')
plt.yscale('log')
plt.plot ( uniform_filter1d(np.abs(np.conjugate(vpq_lr)),8) , c='k')
plt.plot ( uniform_filter1d(np.abs(np.conjugate(vpq_rr)),8) , c='k')
plt.plot ( uniform_filter1d(np.abs(np.conjugate(vpq_ll)),8) , c='k')
plt.clf()
plt.plot ( uniform_filter1d(np.abs(np.conjugate(vpq_ll)),8) , c='k')
plt.plot ( uniform_filter1d(np.abs(np.conjugate(vpq_rr)),8) , c='k')
plt.clf()
plt.plot ( uniform_filter1d(np.abs(np.conjugate(vpq_rr)),8) , c='k')
plt.plot ( uniform_filter1d(np.angle(np.conjugate(vpq_rr)),8) , c='k')
plt.clf()
plt.plot ( uniform_filter1d(np.angle(np.conjugate(vpq_rr)),8) )
plt.plot ( uniform_filter1d(np.angle((vpq_rr)),8) )
plt.plot ( np.angle(gp*(model['stokes_i'])*np.conjugate(gq)) )
plt.clf()
plt.plot ( np.abs(gp*(model['stokes_i'])*np.conjugate(gq)) )
plt.plot ( uniform_filter1d(np.abs((vpq_rr)),8) )
plt.clf()
plt.plot ( uniform_filter1d(np.abs((vpq_rr)),8) )
plt.yscale('log')
plt.plot ( np.abs(gp*(model['stokes_i'])*np.conjugate(gq)) )
plt.clf()
plwd
get_ipython().run_line_magic('pwd', '')
