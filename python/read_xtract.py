"""
xtract out 
either self or cross
"""
import os

from collections import defaultdict

import numpy as np
import pandas as pd

default_ant_list = "C00 C01 C02 C03 C04 C05 C06 C08 C09 C10 C11 C12 C13 C14 E02 E03 E04 E05 E06 S01 S02 S03 S04 S06 W01 W02 W03 W04 W05 W06 C07 S05"

def name_of ( p ):
    """
    we do not know this for sure yet
    """
    if p == "130": return 'r'
    elif p == "175": return 'l'
    else:
        raise ValueError (f"What pol? = {p}")

def read_out ( filename ):
    """
    returns pandas dataframe
    """
    kv   = dict()
    labs = defaultdict(list)
    with open (filename, "r") as _f:
        lines = [a.strip() for a in _f.readlines()]

    for l in lines:
        if l.startswith ("#"):
            toks = l[1:].split()
            ###
            if toks[0] in ['LTAT', 'FREQ00']:
                kv [ toks[0] ] = float ( toks[1] )
            elif toks[0] in ['NROWS', 'NCOLS', 'CHAN1', 'CHAN2', 'CHANINC']:
                kv [ toks[0] ] = int ( toks[1] )
            elif toks[0].startswith ("OBJECT"):
                kv [ 'SOURCE' ] = toks[1]
                kv [ 'RA' ] = toks[3]
                kv [ 'DEC' ] = toks[5]
                kv [ 'wavelength' ] = toks[7]
#OBJECT= SKY_ND RA=   5h34m0055.2s DEC=  23d02'000031" Lambda= 0.545077 DATE-OBS= Mon Mar 23 00:01:09 2026
# 0       1     2      3            4     5              6      7         8        9
            elif toks[0].startswith ("LABEL"):
                # print ( toks )
                i   = int ( toks[0][len('LABEL'):] )
                if i == 0:
                    ## this is channel
                    # labs['li'].append ( i )
                    continue
                tag = toks[1]
                ap  = toks[2]
                ###
                tant1,tant2  = tag.split(':')
                kant1,kant2  = tant1.split('-'), tant2.split('-')

                labs['li'].append ( i )
                labs['ant1'].append ( kant1[0] )
                labs['ant2'].append ( kant2[0] )
                labs['sb1'].append ( kant1[1] )
                labs['sb2'].append ( kant2[1] )
                labs['pol1'].append ( kant1[2] )
                labs['pol2'].append ( kant2[2] )
                labs['ap'].append ( ap )
            elif toks[0].startswith("End"):
                continue
            else:
                print ("unknown # line=", l)

    ## column labels
    colf = pd.DataFrame ( labs )
    ## sanity check
    if len(labs['li']) != kv['NCOLS']-1:
        raise RuntimeError (f" number of columns not matching... expected={kv['NCOLS']} got={len(labs['li'])}")
    # if not np.all ( colf['ant1'] == colf['ant2'] ):
        # raise RuntimeError (f" These are not self-correlations")
    if not np.all ( colf['sb1'] == colf['sb2'] ):
        raise RuntimeError (f" These are not same sideband")
    ## XXX we need same sideband because 
    ## otherwise it does not make sense
    ## and we do not consider it anymore
    ## 
    names = ['chan']
    blset = []
    for _,irow in colf.iterrows ():
        p1 = name_of ( irow.pol1 )
        p2 = name_of ( irow.pol2 )
        t  = f"{irow.ant1}_{irow.ant2}_{p1}{p2}_{irow.ap}"
        names.append ( t )
        blset.append ( (irow.ant1, irow.ant2) )

    ### unique (ant1,ant2) baselines
    blset = list(set(blset))

    df   = pd.read_csv ( filename, sep='\\s+', comment='#', names=names )
    #################
    ret  = {'ant1':[], 'ant2':[], 'rr':[], 'rl':[], 'lr':[], 'll':[]}
    for (ant1,ant2) in blset:
        ## column names 
        _rr  = df [f"{ant1}_{ant2}_rr_a"] * np.exp ( 1.0j * np.deg2rad ( df[f"{ant1}_{ant2}_rr_p"] ) )
        _rl  = df [f"{ant1}_{ant2}_rl_a"] * np.exp ( 1.0j * np.deg2rad ( df[f"{ant1}_{ant2}_rl_p"] ) )
        _lr  = df [f"{ant1}_{ant2}_lr_a"] * np.exp ( 1.0j * np.deg2rad ( df[f"{ant1}_{ant2}_lr_p"] ) )
        _ll  = df [f"{ant1}_{ant2}_ll_a"] * np.exp ( 1.0j * np.deg2rad ( df[f"{ant1}_{ant2}_ll_p"] ) )
        ## save
        ret['ant1'].append ( ant1 )
        ret['ant2'].append ( ant2 )
        ret['rr'].append   ( _rr )
        ret['rl'].append   ( _rl )
        ret['lr'].append   ( _lr )
        ret['ll'].append   ( _ll )
    #################
    cf   = pd.DataFrame ( ret ).sort_values(['ant1','ant2'])
    # return kv, colf, df, pd.DataFrame(rrf), pd.DataFrame(llf),pd.DataFrame(rlf)
    return kv, cf

def get_args ():
    import argparse
    agp = argparse.ArgumentParser ('read_out')
    add = agp.add_argument
    add ('ofile', help=':xtract: generated output file')
    return agp.parse_args()

if __name__ == "__main__":
    args = get_args ()

    odf  = os.path.basename(args.ofile) + "_df.pkl"

    kv, cf  = read_out ( args.ofile )

    cf.to_pickle(odf)
