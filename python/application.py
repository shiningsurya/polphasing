"""
apply solutions
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm


from collections import defaultdict
###############################
def get_args():
    import argparse
    agp = argparse.ArgumentParser("application", description="Applies gain solutions to bldata")
    add = agp.add_argument
    add ('bldata_df', help='Output of :make_bldata_df:')
    add ('-p','--polphase', help='Tag of polphase', dest='polphase')
    add ('-f','--fullpolphase', help='Tag of full polphase', dest='fullpolphase')
    add ('-u','--use-full', help='Use full jones inversion, otherwise only use diagonal terms', action='store_true', dest='use_full')
    return agp.parse_args()

if __name__ == "__main__":
    args = get_args()

    DF   = args.bldata_df

    soltag  = ""
    if args.polphase: 
        soltag = args.polphase

        PPRG   = f"{args.polphase}_r.gains"
        PPLG   = f"{args.polphase}_l.gains"

        rgains = pd.read_csv (PPRG,sep='\\s+').map(complex)
        lgains = pd.read_csv (PPLG,sep='\\s+').map(complex)

    if args.fullpolphase: 
        soltag = args.fullpolphase

        rrgains = pd.read_csv (f"{args.fullpolphase}_rr.gains",sep='\\s+').map(complex)
        rlgains = pd.read_csv (f"{args.fullpolphase}_rl.gains",sep='\\s+').map(complex)
        lrgains = pd.read_csv (f"{args.fullpolphase}_lr.gains",sep='\\s+').map(complex)
        llgains = pd.read_csv (f"{args.fullpolphase}_ll.gains",sep='\\s+').map(complex)


    df   = pd.read_pickle(DF)

    gg   = defaultdict(list)

    if args.polphase:
        for row in df.itertuples():

            ## pull index
            ant1, ant2, (corr1, corr2) = row.Index

            gg['ant1'].append ( ant1 )
            gg['ant2'].append ( ant2 )
            gg['correlation'].append ( corr1+corr2 )

            gg['unapplied'].append ( row.complex )

            if corr1 == 'r' and corr2 == 'r':
                cal  = row.complex / rgains[ant1] / np.conjugate(rgains[ant2])
                gg['applied'].append ( cal )
            elif corr1 == 'r' and corr2 == 'l':
                cal  = row.complex / rgains[ant1] / np.conjugate(lgains[ant2])
                gg['applied'].append ( cal )
            elif corr1 == 'l' and corr2 == 'r':
                cal  = row.complex / lgains[ant1] / np.conjugate(rgains[ant2])
                gg['applied'].append ( cal )
            elif corr1 == 'l' and corr2 == 'l':
                cal  = row.complex / lgains[ant1] / np.conjugate(lgains[ant2])
                gg['applied'].append ( cal )
            else:
                raise RuntimeError("corr not understood")

    elif args.fullpolphase:

        gf   = df.groupby(level=['ant1','ant2'])

        conj = np.conjugate

        if args.use_full:
            print ("performing full inversion ... ")
            soltag += "_fullsol"

            for ((ant1,ant2), igf) in gf:

                ## gains
                gprr, gprl, gplr, gpll  = rrgains[ant1], rlgains[ant1], lrgains[ant1], llgains[ant1]
                gqrr, gqrl, gqlr, gqll  = rrgains[ant2], rlgains[ant2], lrgains[ant2], llgains[ant2]

                ## data
                dpqrr = np.array ( igf.complex.loc[ant1, ant2, 'rr'] )
                dpqrl = np.array ( igf.complex.loc[ant1, ant2, 'rl'] )
                dpqlr = np.array ( igf.complex.loc[ant1, ant2, 'lr'] )
                dpqll = np.array ( igf.complex.loc[ant1, ant2, 'll'] )

                """

                math from math_application

                """
                ## rr
                unapplied = dpqrr
                applied   = dpqll*gprl*conj(gqrl)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                - dpqlr*gprl*conj(gqll)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                - dpqrl*gpll*conj(gqrl)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                + dpqrr*gpll*conj(gqll)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl))
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "rr" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )

                unapplied  = dpqrl
                applied    = -dpqll*gprl*conj(gqrr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        + dpqlr*gprl*conj(gqlr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        + dpqrl*gpll*conj(gqrr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        - dpqrr*gpll*conj(gqlr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl))
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "rl" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )

                unapplied  = dpqlr
                applied    = -dpqll*gprr*conj(gqrl)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        + dpqlr*gprr*conj(gqll)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        + dpqrl*gplr*conj(gqrl)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        - dpqrr*gplr*conj(gqll)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl))
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "lr" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )

                unapplied  = dpqll
                applied    = dpqll*gprr*conj(gqrr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        - dpqlr*gprr*conj(gqlr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        - dpqrl*gplr*conj(gqrr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        + dpqrr*gplr*conj(gqlr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl))
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "ll" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )

        else:
            print ("performing only diagonal inversion with fullpolphase solutions ... ")
            soltag += "_diagsol"

            for ((ant1,ant2), igf) in gf:
                ## gains
                gprr, gpll  = rrgains[ant1], llgains[ant1]
                gqrr, gqll  = rrgains[ant2], llgains[ant2]

                ## data
                dpqrr = np.array ( igf.complex.loc[ant1, ant2, 'rr'] )
                dpqrl = np.array ( igf.complex.loc[ant1, ant2, 'rl'] )
                dpqlr = np.array ( igf.complex.loc[ant1, ant2, 'lr'] )
                dpqll = np.array ( igf.complex.loc[ant1, ant2, 'll'] )

                ## rr
                unapplied = dpqrr
                applied   = dpqrr / gprr / conj(gqrr)
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "rr" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )

                ## rl
                unapplied = dpqrl
                applied   = dpqrl / gprr / conj(gqll)
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "rl" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )

                ## lr
                unapplied = dpqlr
                applied   = dpqlr / gpll / conj(gqrr)
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "lr" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )

                ## ll
                unapplied = dpqll
                applied   = dpqll / gpll / conj(gqll)
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "ll" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )

    gg  = pd.DataFrame(gg).set_index(['ant1','ant2','correlation']).sort_index()

    odf  = DF[:-len("bldata_df.pkl")] + soltag + ".bldata_df.pkl"
    print(f"Writing to {odf}")

    gg.to_pickle(odf)
