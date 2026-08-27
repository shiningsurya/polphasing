"""
make dataframe for easier access
"""
import os
import numpy as np
import pandas as pd

def get_args():
    import argparse
    agp = argparse.ArgumentParser("make_solutions_df", description="Make complex gains")
    add = agp.add_argument
    add ('-p','--polphase', help='Tag of polphase', dest='polphase')
    add ('-r','--rantsol', help='Tag of rantsol', dest='rantsol')
    add ('-n','--nchan', help='Number of channels', dest='nchans', type=int, default=2048)
    return agp.parse_args()

def make_indexed_gaintable (rg, lg):
    """
    multindex dataframe with {ant, corr} as index

    it helps a lot down the line
    """
    ##
    trg   = rg.transpose()
    tlg   = lg.transpose()
    ###
    trg['gains']  = trg.values.tolist()
    tlg['gains']  = tlg.values.tolist()
    ###
    rr    = pd.DataFrame ( trg['gains'] )
    ll    = pd.DataFrame ( tlg['gains'] )
    ###
    gg    = pd.concat([rr, ll], keys=['r', 'l'], names=['corr', 'ant']).swaplevel().sort_index()
    gg['gains'] = gg['gains'].apply(np.array)
    return gg

if __name__ == "__main__":
    args = get_args()

    if args.polphase:
        og  = os.path.basename ( args.polphase )
        rg  = pd.read_csv ( args.polphase + "_r.gains", sep='\\s+' ).map(complex)
        lg  = pd.read_csv ( args.polphase + "_l.gains", sep='\\s+' ).map(complex)

        gg  = make_indexed_gaintable ( rg, lg )
        gg.to_pickle ( og + "_polphase_df.pkl" )

    if args.rantsol:
        og  = os.path.basename ( args.rantsol )
        rg  = pd.read_csv ( args.rantsol +".130.dat",sep='\\s+', nrows=args.nchans).shift(1,axis='columns').map(lambda x : np.exp(1.0j*np.deg2rad(x))).drop(columns=['#'])
        lg  = pd.read_csv ( args.rantsol +".175.dat",sep='\\s+', nrows=args.nchans).shift(1,axis='columns').map(lambda x : np.exp(1.0j*np.deg2rad(x))).drop(columns=['#'])

        gg  = make_indexed_gaintable ( rg, lg )
        gg.to_pickle ( og + "_rantsol_df.pkl" )

        # rg.to_pickle ( og + "_rantsol_r_df.pkl" )
        # lg.to_pickle ( og + "_rantsol_l_df.pkl" )

