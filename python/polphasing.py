"""
Polarimetric phasing using standard polarized quasar
"""

from collections import defaultdict

from tqdm import tqdm

import numpy as np
import pandas as pd

import scipy.optimize as so

class PolPhasing:
    """
    Solves for complex gains for both hands (R and L) for every channel independently

    It solves using scipy.optimize.least_squares(method='lm') which uses Levenburg Marquadt

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
    def __init__ (self, dataframe, refant='C00'):
        """
        dataframe: ant1, ant2, bl, {rr, rl, lr, ll} x {phased, unphased}

        only need col because we somehow need to find nchan

        dfarrrys
        """
        ## these cols are hardcoded for now
        cols      = ['rr_unphased', 'rl_unphased', 'lr_unphased', 'll_unphased']
        self.df   = dataframe
        self.nbl  = self.df.shape[0]
        self.nchan= self.df [ cols[0] ].iloc [ 0 ].size
        ## also contains selfs
        self.ants = sorted(list(set(self.df['ant1']).union(set(self.df['ant2']))))
        self.nant = len(self.ants)
        self.ant2idx = {ant:inum for inum,ant in enumerate(self.ants)}
        ### preget ants
        self.df['iant1'] = self.df['ant1'].apply ( lambda ant : self.ant2idx[ant] )
        self.df['iant2'] = self.df['ant2'].apply ( lambda ant : self.ant2idx[ant] )
        ## gains
        self.rgains    = {ant:np.zeros(self.nchan,dtype=np.complex64) for ant in self.ants}
        self.lgains    = {ant:np.zeros(self.nchan,dtype=np.complex64) for ant in self.ants}
        ## residuals
        self.cost  = np.zeros (self.nchan)
        ## precomputing
        ## reals and imags indices of solution
        ## GZ HZ 
        self.a1gr  = 4 * self.df['iant1'] 
        self.a1gi  = 4 * self.df['iant1'] + 1
        self.a1hr  = 4 * self.df['iant1'] + 2
        self.a1hi  = 4 * self.df['iant1'] + 3
        ###
        self.a2gr  = 4 * self.df['iant2']
        self.a2gi  = 4 * self.df['iant2'] + 1
        self.a2hr  = 4 * self.df['iant2'] + 2
        self.a2hi  = 4 * self.df['iant2'] + 3
        ## dfarrays for every frequency
        ## complex64 is to be ensured
        ## (chan, pol,baseline)
        self.dfarrays  = np.complex64 ( np.array ( self.df[ cols ].values.tolist() ) ).T
        ### XXX need to think of a better way to propagate nans
        ### least_squares does not like nans
        ### for now, replacing nans with zero
        self.dfarrays  = np.nan_to_num ( self.dfarrays )
        ### 
        self.npar = self.nant * 4
        self.nres = self.nbl * 4 * 2
        ## i have to ensure complex64
        self.res  = np.zeros ( self.nbl * 4, dtype=np.complex64 )
        ## the following indicies are for complex64
        ## which is then viewed as float32
        self.rirr = slice(0*self.nbl,1*self.nbl)
        self.rirl = slice(1*self.nbl,2*self.nbl)
        self.rilr = slice(2*self.nbl,3*self.nbl)
        self.rill = slice(3*self.nbl,4*self.nbl)
        ## rr are from 0*nbl:1*nbl in complex64 array
        ## re(rr) = 2*(0:nbl)
        ## im(rr) = 2*(0:nbl) + 1
        self.bl_rr_re = 2*np.r_[0*self.nbl:1*self.nbl]
        self.bl_rr_im = 2*np.r_[0*self.nbl:1*self.nbl] + 1
        ## ll are from 3*nbl:4*nbl in complex64 array
        self.bl_ll_re = 2*np.r_[3*self.nbl:4*self.nbl]
        self.bl_ll_im = 2*np.r_[3*self.nbl:4*self.nbl] + 1
        ## rl are from 1*nbl:2*nbl in complex64 array
        self.bl_rl_re = 2*np.r_[1*self.nbl:2*self.nbl]
        self.bl_rl_im = 2*np.r_[1*self.nbl:2*self.nbl] + 1
        ## lr are from 2*nbl:3*nbl
        self.bl_lr_re = 2*np.r_[2*self.nbl:3*self.nbl]
        self.bl_lr_im = 2*np.r_[2*self.nbl:3*self.nbl] + 1

    def solve_chan (self, ichan, model):
        """
        solves for one channel
        """
        ## every ant
        ## ensure float32
        ## [ Gain(hand=R) Gains(hand=L) ...  ]
        initial_sol      = np.ones ( self.npar, dtype=np.float32 )

        ## call minimize
        res = so.least_squares (
            self.objective,
            initial_sol,
            jac=self.jacobian,
            # jac='3-point',
            method='lm',
            x_scale='jac',
            verbose=1,
            args=(ichan, model),
        )

        #
        found_sol  = res.x
        # if not res.success:
        print(f" status={res.status:d} msg={res.message} cost={res.cost:.3e} nfev={res.nfev:d} njev={res.njev:d}")
        self.cost[ichan] = res.cost

        ## save 
        for ant in self.ants:
            idx          = self.ant2idx [ ant ]
            self.rgains[ant][ichan] = found_sol[4*idx + 0] + 1.0j*found_sol[4*idx + 1]
            self.lgains[ant][ichan] = found_sol[4*idx + 2] + 1.0j*found_sol[4*idx + 3]

    def objective (self, initial_sol, ichan, model ):
        """
        Function to be minimized

        residuals as R I R I 
        [rr .. rl .. lr .. ll]

        model should be 
            rr rl lr ll

        """
        ## ensure initial_sol is float32
        gp    = initial_sol [ self.a1gr ] + 1.0j*initial_sol [ self.a1gi ]
        hp    = initial_sol [ self.a1hr ] + 1.0j*initial_sol [ self.a1hi ]
        gqc   = initial_sol [ self.a2gr ] - 1.0j*initial_sol [ self.a2gi ]
        hqc   = initial_sol [ self.a2hr ] - 1.0j*initial_sol [ self.a2hi ]
        ## precompute for every frequency
        Xpq   = self.dfarrays[ichan]
        ## ensure self.dfarrays[ichan] is np.complex64
        self.res[self.rirr] = Xpq [0] - ( model[0] * gp * gqc )
        self.res[self.rirl] = Xpq [1] - ( model[1] * gp * hqc )
        self.res[self.rilr] = Xpq [2] - ( model[2] * hp * gqc )
        self.res[self.rill] = Xpq [3] - ( model[3] * hp * hqc )
        return self.res.view(np.float32)

    def jacobian (self, initial_sol, ichan, model):
        """
        returns jacobian
        ichan is here simply to match the function signatures

        model should be 
            rr rl lr ll

        rr = I + V
        ll = I - V
        rl = Q + jU
        lr = Q - jU
        there may be a sign flip with rl and lr
        XXX need to investigate further

        V is set to zero.

        jacobian is (#residuals,#parameters)
        residuals are [rr rl lr ll]
        real and imaginary parts
        """
        I, Q, U, V = model
        I     = model [ 0 ]
        ## somewhere sanity check 
        ## if model[0] == model[3] since V = 0
        Q     = np.real ( model[1] )
        U     = np.imag ( model[1] )
        ## ensure initial_sol is float32
        gpr   = initial_sol [ self.a1gr ] 
        gpi   = initial_sol [ self.a1gi ]
        hpr   = initial_sol [ self.a1hr ] 
        hpi   = initial_sol [ self.a1hi ]
        gqr   = initial_sol [ self.a2gr ] 
        gqi   = initial_sol [ self.a2gi ]
        hqr   = initial_sol [ self.a2hr ] 
        hqi   = initial_sol [ self.a2hi ]

        ## precompute for every frequency
        ## jacobian does not contain the data
        ## Xpq   = self.dfarrays[ichan]

        ret   = np.zeros ( ( self.nres, self.npar ), dtype=np.float32 )

        ###########################
        ## all should be negative
        ###########################

        ## re(rr) 
        ret[self.bl_rr_re, self.a1gr]  = -1.0 * (   I * gqr )
        ret[self.bl_rr_re, self.a1gi]  = -1.0 * (   I * gqi )
        ret[self.bl_rr_re, self.a2gr]  = -1.0 * (   I * gpr )
        ret[self.bl_rr_re, self.a2gi]  = -1.0 * (   I * gpi )

        ## im(rr)
        ret[self.bl_rr_im, self.a1gr]  = -1.0 * ( - I * gqi )
        ret[self.bl_rr_im, self.a1gi]  = -1.0 * (   I * gqr )
        ret[self.bl_rr_im, self.a2gr]  = -1.0 * (   I * gpi )
        ret[self.bl_rr_im, self.a2gi]  = -1.0 * ( - I * gpr )

        ## re(ll)
        ret[self.bl_ll_re, self.a1hr]  = -1.0 * (   I * hqr )
        ret[self.bl_ll_re, self.a1hi]  = -1.0 * (   I * hqi )
        ret[self.bl_ll_re, self.a2hr]  = -1.0 * (   I * hpr )
        ret[self.bl_ll_re, self.a2hi]  = -1.0 * (   I * hpi )

        ## im(ll)
        ret[self.bl_ll_im, self.a1hr]  = -1.0 * ( - I * hqi )
        ret[self.bl_ll_im, self.a1hi]  = -1.0 * (   I * hqr )
        ret[self.bl_ll_im, self.a2hr]  = -1.0 * (   I * hpi )
        ret[self.bl_ll_im, self.a2hi]  = -1.0 * ( - I * hpr )

        ## re(rl)
        ret[self.bl_rl_re, self.a1gr]  = -1.0 * ( Q*hqr + U*hqi )
        ret[self.bl_rl_re, self.a1gi]  = -1.0 * ( Q*hqi - U*hqr )
        ret[self.bl_rl_re, self.a2hr]  = -1.0 * ( Q*gpr - U*gpi )
        ret[self.bl_rl_re, self.a2hi]  = -1.0 * ( Q*gpi + U*gpr )

        ## im(rl)
        ret[self.bl_rl_im, self.a1gr]  = -1.0 * (-Q*hqi + U*hqr )
        ret[self.bl_rl_im, self.a1gi]  = -1.0 * ( Q*hqr + U*hqi )
        ret[self.bl_rl_im, self.a2hr]  = -1.0 * ( Q*gpi + U*gpr )
        ret[self.bl_rl_im, self.a2hi]  = -1.0 * (-Q*gpr + U*gpi )

        ## re(lr)
        ret[self.bl_lr_re, self.a1hr]  = -1.0 * ( Q*gqr - U*gqi )
        ret[self.bl_lr_re, self.a1hi]  = -1.0 * ( Q*gqi + U*gqr )
        ret[self.bl_lr_re, self.a2gr]  = -1.0 * ( Q*hpr + U*hpi )
        ret[self.bl_lr_re, self.a2gi]  = -1.0 * ( Q*hpi - U*hpr )

        ## im(lr)
        ret[self.bl_lr_im, self.a1hr]  = -1.0 * (-Q*gqi - U*gqr)
        ret[self.bl_lr_im, self.a1hi]  = -1.0 * ( Q*gqr - U*gqi)
        ret[self.bl_lr_im, self.a2gr]  = -1.0 * ( Q*hpi - U*hpr)
        ret[self.bl_lr_im, self.a2gi]  = -1.0 * (-Q*hpr - U*hpi)

        ##
        return ret

if __name__ == "__main__":
    pf = pd.read_pickle ("3C138_df.pkl")

    noself = pd.DataFrame ( pf [ pf['ant1'] != pf['ant2'] ] )

    phasing = PolPhasing ( noself )
    
    #model   = [ 1.0, 0.2 + 0.2j, 0.2 - 0.2j, 1.0 ]
    mf    = pd.read_csv ("3C138_band4_2K_model_df.csv")
    iref  = mf['I'].iloc[0] ## first chan
    for k in ['I','Q','U']: mf[k] = mf[k] / iref

    # for ichan in tqdm ( range(500,525), desc='channels', unit='chan' ):
    for ichan in tqdm ( range(2048), desc='channels', unit='chan' ):
        mi, mq, mu = mf['I'].iloc[ichan], mf['Q'].iloc[ichan], mf['U'].iloc[ichan]
        model      = [mi, mq+1.0j*mu, mq-1.0*mu, mi] 
        phasing.solve_chan ( ichan, model )




