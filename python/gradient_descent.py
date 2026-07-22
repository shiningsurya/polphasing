"""
polarimetric phasing 
but solving using gradient descent.

i am not happy with the preformance of LMsolver
"""
from collections import defaultdict

from tqdm import tqdm

import numpy as np
import pandas as pd

class PolPhasing:
    """
    Solves for complex gains for both hands (R and L) for every channel independently

    It solves using gradient descent method similar to rantsol.

    For Nant number of antennas, 
    there are 4 * Nant real input parameters and (Nant C 2) * 4 * 2 residuals.

    It is an over-determined system.

    objective and jacobian functions are provided.
    Jacobian matrix is sparse.

    input parameters are treated as 
    [R I R I R I R I .... R I R I] 
    |------| |-----| .... |------|
    [ Gain(hand=R,ant=i) Gain(hand=L,ant=i) ... Gain(hand=R,ant=k) Gain(hand=L,ant=k) ]
    real and imaginary parts 

    """
    def __init__ (self, dataframe, nbl, nchan):
        """
        dataframe contains
        ['chan', 'complex', 'model', 'ant1', 'ant2', 'band1', 'band2']
        """
        ## these cols are hardcoded for now
        self.df   = dataframe
        self.nbl  = nbl
        self.nchan= nchan
        ## also contains selfs
        self.ants = sorted(list(set(self.df['ant1']).union(set(self.df['ant2']))))
        self.nant = len(self.ants)
        self.ant2idx = {ant:inum for inum,ant in enumerate(self.ants)}
        ### compute b
        ## gains
        self.rgains    = {ant:np.zeros(self.nchan,dtype=np.complex64) for ant in self.ants}
        self.lgains    = {ant:np.zeros(self.nchan,dtype=np.complex64) for ant in self.ants}
        ## residuals
        self.cost  = np.zeros (self.nchan)
        ### 
        self.npar = self.nant * 4

    def gd_iterate ( self, usol, data, model, b1, b2,  alpha=0.45 ):
        """
        one iteration of gradient descent
        """
        g_nr  = np.zeros ( 2*self.nant, dtype=np.complex64 )
        g_dr  = np.zeros ( 2*self.nant, dtype=np.float32 )
        ###################################
        for ibl in range ( self.nbl ):
            ## for every polar baseline
            _b1  = b1[ibl]
            _b2  = b2[ibl]
            ##
            d    = data.iloc[ibl]
            m    = model.iloc[ibl]
            ##
            g1   = usol[2*_b1 + 0] + 1.0j*usol[2*_b1 + 1]
            g2   = usol[2*_b2 + 0] + 1.0j*usol[2*_b2 + 1]
            ##
            ### direct
            g_nr[_b1] += ( d * np.conjugate(m) * g2 )
            g_dr[_b1] += np.real( m*np.conjugate(m) * g2*np.conjugate(g2) )
            ### conj
            g_nr[_b2] += ( np.conjugate(d) * np.conjugate(m) * g1 )
            g_dr[_b2] += np.real( m*np.conjugate(m) * g1*np.conjugate(g1) )
        ###################################
        pg    = g_nr / g_dr 
        vsol  = np.array ( pg.view(np.float32) )
        ###################################
        ## do lerp
        return (1.0 - alpha)*usol + alpha*vsol
        # return usol + alpha*( vsol - usol )

    def solve_chan (self, ichan, max_iterations=100, alpha=0.45, delta=0.1):
        """
        solves for one channel
        """
        ## every ant
        ## ensure float32
        ## [ Gain(hand=R) Gains(hand=L) ...  ]
        jsol    = np.ones ( self.npar, dtype=np.float32 )

        chanf   = self.df.loc[ichan]
        data    = chanf['complex']
        model   = chanf['model']

        ant1,band1 = chanf['ant1'],chanf['band1']
        ant2,band2 = chanf['ant2'],chanf['band2']

        b1      = np.array ( [2*self.ant2idx[iant]+iband for iant,iband in zip(ant1,band1)] )
        b2      = np.array ( [2*self.ant2idx[iant]+iband for iant,iband in zip(ant2,band2)] )

        ##
        cost       = 0.0
        last_cost  = 0.0
        retcode    = 0
        for it in range ( max_iterations ):
            ## gradient descent iterate
            jsol   = self.gd_iterate ( jsol, data, model, b1, b2, alpha=alpha )
            ## forward pass
            pg1    = jsol[2*b1 + 0] + 1.0j*jsol[2*b1 + 1]
            pg2    = jsol[2*b2 + 0] + 1.0j*jsol[2*b2 + 1]
            pfm    = pg1 * model * np.conjugate(pg2)
            ## cost
            cost   = np.sum ( np.power ( np.abs ( data - pfm ),2.0 ) )

            ## termination condition
            if abs ( cost - last_cost ) < delta:
                retcode = 1
                break

            ## save last code
            last_cost = cost

        #
        self.cost[ichan] = cost
        found_sol        = jsol

        if retcode == 0:
            print ("did not coverge for chan={ichan:d}")

        ## save 
        for ant in self.ants:
            idx          = self.ant2idx [ ant ]
            self.rgains[ant][ichan] = found_sol[4*idx + 0] + 1.0j*found_sol[4*idx + 1]
            self.lgains[ant][ichan] = found_sol[4*idx + 2] + 1.0j*found_sol[4*idx + 3]

if __name__ == "__main__":

    if False:
        bldata  = pd.read_csv ("UNFLIP_SCAN12.bldata", sep='\\s+')
        bldata['complex'] = bldata['complex'].apply(complex)
        ## take nonself correlation products
        bl      = pd.DataFrame ( bldata [ (bldata['ant1'] != bldata['ant2']) ] )
        ## bands 
        bmap    = {'r':0, 'l':1}
        bl['band1']  = bl['correlation'].apply ( lambda x : bmap[x[0]] )
        bl['band2']  = bl['correlation'].apply ( lambda x : bmap[x[1]] )
        ## models
        model  = pd.read_csv ("/home/shining/shit/gmrt_phase_test/polphasing/code/polphasing/models/3C138_2048_97656.25_550000000_model.txt", sep='\\s+')
        model['rr']  = model['stokes_i']
        model['ll']  = model['stokes_i']
        model['lr']  = model['stokes_q'] - 1.0j*model['stokes_u']
        model['rl']  = model['stokes_q'] + 1.0j*model['stokes_u']

        bl['model']  = [model[icorr].iloc[ichan] for icorr,ichan in zip (bl['correlation'],bl['chan'])]

        bl.to_pickle ( "bldata_df.pkl" )
    else:
        bl      = pd.read_pickle("bldata_df.pkl").set_index(['chan'])


    nchan = 2048
    nbl   = 1012

    gder  = PolPhasing ( bl, nbl, nchan )
