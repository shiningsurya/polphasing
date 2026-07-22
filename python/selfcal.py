"""
We use solved complex gains and the data used to solve for the gains, to update the model.

This is like self-calibration.
We update the model.

i see that RR and LL fits well with the data.
But RL and LR do not.

i suspect it is because model for RL and LR (like Q and U) is not accurate.

so this selfcal step to make it 
"""

from collections import defaultdict

from tqdm import tqdm

import numpy as np
import pandas as pd

import scipy.optimize as so

class Selfcal:
    """
    Solves for model (IQU per channel) using solved complex gains and corresponding data.

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
    def __init__ (self, rrdf, rldf, lrdf, lldf, nchan=2048):
        """
        dataframe: ant1, ant2, bl, {rr, rl, lr, ll} x {phased, unphased}
        this is out dataframe version of output of :write_baselines:

        only need col because we somehow need to find nchan

        """
        ### data dataframes
        self.blrr = rrdf
        self.blrl = rldf
        self.bllr = lrdf
        self.blll = lldf
        ###  each is in exploded form
        ### take first chan, its shape is nbl
        self.nbl  = self.blrr.loc[0].shape[0]
        self.nchan= 2048
        ## residuals
        self.cost  = np.zeros (self.nchan, dtype=np.float32)
        ## complex64 is to be ensured
        ## (chan, pol,baseline)
        # self.dfarrays  = np.complex64 ( np.array ( self.df[ cols ].values.tolist() ) ).T
        ### XXX need to think of a better way to propagate nans
        ### least_squares does not like nans
        ### for now, replacing nans with zero
        ### 
        self.npar = 3
        ###
        self.li   = 0
        self.lq   = 1
        self.lu   = 2
        ###
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
        ## rl are from 1*nbl:2*nbl in complex64 array
        self.bl_rl_re = 2*np.r_[1*self.nbl:2*self.nbl]
        self.bl_rl_im = 2*np.r_[1*self.nbl:2*self.nbl] + 1
        ## lr are from 2*nbl:3*nbl
        self.bl_lr_re = 2*np.r_[2*self.nbl:3*self.nbl]
        self.bl_lr_im = 2*np.r_[2*self.nbl:3*self.nbl] + 1
        ## ll are from 3*nbl:4*nbl in complex64 array
        self.bl_ll_re = 2*np.r_[3*self.nbl:4*self.nbl]
        self.bl_ll_im = 2*np.r_[3*self.nbl:4*self.nbl] + 1
        ###
        self.fitted_i = np.zeros ( self.nchan, dtype=np.float32 )
        self.fitted_q = np.zeros ( self.nchan, dtype=np.float32 )
        self.fitted_u = np.zeros ( self.nchan, dtype=np.float32 )

    def solve_chan (self, ichan):
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
            args=(ichan,),
        )

        #
        found_sol  = res.x
        # if not res.success:
        print(f" status={res.status:d} msg={res.message} cost={res.cost:.3e} nfev={res.nfev:d} njev={res.njev:d}")
        self.cost[ichan] = res.cost

        self.fitted_i[ichan]  = found_sol[self.li]
        self.fitted_q[ichan]  = found_sol[self.lq]
        self.fitted_u[ichan]  = found_sol[self.lu]

    def objective (self, initial_sol, ichan ):
        """
        Function to be minimized

        residuals as R I R I 
        [rr .. rl .. lr .. ll]

        model should be 
            rr rl lr ll

        """
        ## ensure initial_sol is float32
        mi, mq, mu  = initial_sol
        model       = [mi, mq+1.0j*mu, mq-1.0*mu, mi] 
        ## precompute for every frequency
        Xrr         = self.blrr.loc [ ichan ]
        Xrl         = self.blrl.loc [ ichan ]
        Xll         = self.blll.loc [ ichan ]
        Xlr         = self.bllr.loc [ ichan ]
        ## ensure self.dfarrays[ichan] is np.complex64
        self.res[self.rirr] = Xrr['complex'] - ( model[0] * Xrr['g1'] * Xrr['g2c'] )
        self.res[self.rirl] = Xrl['complex'] - ( model[1] * Xrl['g1'] * Xrl['g2c'] )
        self.res[self.rilr] = Xlr['complex'] - ( model[2] * Xlr['g1'] * Xlr['g2c'] )
        self.res[self.rill] = Xll['complex'] - ( model[3] * Xll['g1'] * Xll['g2c'] )
        return self.res.view(np.float32)

    def jacobian (self, initial_sol, ichan):
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
        ## ensure initial_sol is float32
        mi, mq, mu  = initial_sol
        model       = [mi, mq+1.0j*mu, mq-1.0*mu, mi] 
        ## precompute for every frequency
        Xrr         = self.blrr.loc [ ichan ]
        Xrl         = self.blrl.loc [ ichan ]
        Xll         = self.blll.loc [ ichan ]
        Xlr         = self.bllr.loc [ ichan ]

        ret   = np.zeros ( ( self.nres, self.npar ), dtype=np.float32 )

        ###########################
        ## all should be negative
        ###########################

        tr    = np.real ( Xrr['g1'] )
        ti    = np.imag ( Xrr['g1'] )
        sr    = np.real ( Xrr['g2c'] )
        si    = -np.imag ( Xrr['g2c'] )
        ## g2c is pre-conjugated.
        ret[self.bl_rr_re, self.li]  = -1.0 * ( si*ti + sr*tr )
        ret[self.bl_rr_im, self.li]  = -1.0 * (-si*tr + sr*ti )

        tr    = np.real ( Xll['g1'] )
        ti    = np.imag ( Xll['g1'] )
        sr    = np.real ( Xll['g2c'] )
        si    = -np.imag ( Xll['g2c'] )
        ## g2c is pre-conjugated.
        ret[self.bl_ll_re, self.li]  = -1.0 * ( si*ti + sr*tr )
        ret[self.bl_ll_im, self.li]  = -1.0 * (-si*tr + sr*ti )

        tr    = np.real ( Xrl['g1'] )
        ti    = np.imag ( Xrl['g1'] )
        sr    = np.real ( Xrl['g2c'] )
        si    = -np.imag ( Xrl['g2c'] )
        ## g2c is pre-conjugated.
        ret[self.bl_rl_re, self.lq]  = -1.0 * ( si*ti + sr*tr )
        ret[self.bl_rl_im, self.lq]  = -1.0 * (-si*tr + sr*ti )
        ret[self.bl_rl_re, self.lu]  = -1.0 * ( si*tr - sr*ti )
        ret[self.bl_rl_im, self.lu]  = -1.0 * ( si*ti + sr*tr )

        tr    = np.real ( Xlr['g1'] )
        ti    = np.imag ( Xlr['g1'] )
        sr    = np.real ( Xlr['g2c'] )
        si    = -np.imag ( Xlr['g2c'] )
        ## g2c is pre-conjugated.
        ret[self.bl_lr_re, self.lq]  = -1.0 * ( si*ti + sr*tr )
        ret[self.bl_lr_im, self.lq]  = -1.0 * (-si*tr + sr*ti )
        ret[self.bl_lr_re, self.lu]  =  1.0 * ( si*tr - sr*ti )
        ret[self.bl_lr_im, self.lu]  =  1.0 * ( si*ti + sr*tr )
        ## the sign negation here because 
        ## lr is Q - jU
        ## there is already a minus sign there

        ##
        return ret

if __name__ == "__main__":
    """
    preparing the input files is a bit pain
    """
    nchan   = 2048
    ##
    rgains  = pd.read_csv ("UNFLIP_r.gains", sep='\\s+').map(complex)
    lgains  = pd.read_csv ("UNFLIP_l.gains", sep='\\s+').map(complex)

    bltag   = "UNFLIP_SCAN12"

    if False:
        bldata  = pd.read_csv (f"{bltag}.bldata", sep='\\s+')
        bldata['complex'] = bldata['complex'].apply(complex)
        ## take nonself correlation products
        blrr    = pd.DataFrame ( bldata [ (bldata['ant1'] != bldata['ant2']) & (bldata['correlation'] == 'rr') ] )
        blrl    = pd.DataFrame ( bldata [ (bldata['ant1'] != bldata['ant2']) & (bldata['correlation'] == 'rl') ] )
        bllr    = pd.DataFrame ( bldata [ (bldata['ant1'] != bldata['ant2']) & (bldata['correlation'] == 'lr') ] )
        blll    = pd.DataFrame ( bldata [ (bldata['ant1'] != bldata['ant2']) & (bldata['correlation'] == 'll') ] )
        ## add solved complex gains to the table
        def add_gains ( df, g1, g2 ):
            """
            df has 'chan', 'ant1' and 'ant2'
            ant2 gain is conjugated
            """
            df['g1']  = [ g1[iant].iloc[ichan] for ichan,iant in zip(df.chan, df.ant1) ]
            df['g2c'] = [ np.conjugate(g2[iant].iloc[ichan]) for ichan,iant in zip(df.chan, df.ant2) ]
            return df
        
        blrr  = add_gains ( blrr, rgains, rgains )
        blrl  = add_gains ( blrl, rgains, lgains )
        bllr  = add_gains ( bllr, lgains, rgains )
        blll  = add_gains ( blll, lgains, lgains )
        ##
        blrr.to_pickle ( f"{bltag}_rr_df.pkl" )
        blrl.to_pickle ( f"{bltag}_rl_df.pkl" )
        bllr.to_pickle ( f"{bltag}_lr_df.pkl" )
        blll.to_pickle ( f"{bltag}_ll_df.pkl" )
    else:
        blrr  = pd.read_pickle (f"{bltag}_rr_df.pkl").set_index(['chan','ant1', 'ant2'])
        blrl  = pd.read_pickle (f"{bltag}_rl_df.pkl").set_index(['chan','ant1', 'ant2'])
        bllr  = pd.read_pickle (f"{bltag}_lr_df.pkl").set_index(['chan','ant1', 'ant2'])
        blll  = pd.read_pickle (f"{bltag}_ll_df.pkl").set_index(['chan','ant1', 'ant2'])


    selfcal = Selfcal ( blrr, blrl, bllr, blll, nchan=nchan )

    for ichan in tqdm (range(nchan), desc='Channel', unit='chan'):
        selfcal.solve_chan ( ichan )

    
    


