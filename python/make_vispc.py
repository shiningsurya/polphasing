# coding: utf-8
"""

take bldata_df average over all baselines

write the result as dataframe
"""

import numpy as np
import pandas as pd

###############################
def get_args():
    import argparse
    agp = argparse.ArgumentParser("make_vispc", description="Makes frequency averaged PC from visi")
    add = agp.add_argument
    add ('bldata_df', help='Output of :make_bldata:')
    add ('-p','--polphase', help='polphase type :make_solutions_df:', dest='polphase')
    add ('-f','--fullpolphase', help='polphase type :make_solutions_df:', dest='fullpolphase')
    return agp.parse_args()

if __name__ == "__main__":
    args  = get_args()

    bn = args.bldata_df[:-len("df.pkl")]

    psr = pd.read_pickle( args.bldata_df )


    ############################
    if args.polphase:
        dd  = {'rr':[], 'rl':[], 'lr':[], 'll':[]}
        sols = pd.read_pickle(args.polphase)

        for row in psr.itertuples():
            ant1, ant2, corr = row.Index
            ## skip if self 
            if ant1 == ant2: continue
            ## extract corr
            corr1, corr2 = corr
            ## calibration
            ga = sols.gains.loc[ant1, corr1] 
            gb = sols.gains.loc[ant2, corr2]
            ###
            # normalize for fun
            # ga = np.exp(1.0j * np.angle(ga))
            # gb = np.exp(1.0j * np.angle(gb))
            ###

            caled = np.array ( row.complex  ) / ga / np.conjugate(gb)
            ## append
            dd[corr].append ( caled )

        dfs = dict()
        for k,v in dd.items():
            dfs[k] = np.mean(v, axis=0)

        dd  = pd.DataFrame(dfs)
        odf = bn + "vispc_caled_df.pkl"
        print ( f"Writing to {odf}")
        dd.to_pickle(odf)

    elif args.fullpolphase:
        raise NotImplementedError("uwuw")
    else:
        ### remove the self terms
        i1 = psr.index.get_level_values('ant1')
        i2 = psr.index.get_level_values('ant2')
        corr = psr.index.get_level_values('correlation')

        ## get flags
        irr = (i1 != i2) & (corr=='rr')
        irl = (i1 != i2) & (corr=='rl')
        ilr = (i1 != i2) & (corr=='lr')
        ill = (i1 != i2) & (corr=='ll')

        print ("averaging ... ", end='')

        ## do average
        rr  = np.array(list(psr.complex.loc[irr])).mean(0)
        rl  = np.array(list(psr.complex.loc[irl])).mean(0)
        lr  = np.array(list(psr.complex.loc[ilr])).mean(0)
        ll  = np.array(list(psr.complex.loc[ill])).mean(0)

        print ( "done" )

        ## save as dataframe
        dd  = pd.DataFrame ( {'rr':rr, 'rl':rl, 'lr':lr, 'll':ll} )

        odf = bn + "vispc_df.pkl"
        print ( f"Writing to {odf}")
        dd.to_pickle(odf)

