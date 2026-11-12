"""
apply solutions
"""
import numpy as np
import numpy.linalg as nl
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
    add ('-g','--gwb', help='soldf, made with :make_solutions_df: with -g,--gwb', dest='gwb')
    add ('-f','--fullpolphase', help='Tag of full polphase', dest='fullpolphase')
    add ('-s','--soldf', help='Solutions df :make_solutions_df: or :normalizer:', dest='soldf')
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
        soltag += "_fullsol"

        rrgains = pd.read_csv (f"{args.fullpolphase}_rr.gains",sep='\\s+').map(complex)
        rlgains = pd.read_csv (f"{args.fullpolphase}_rl.gains",sep='\\s+').map(complex)
        lrgains = pd.read_csv (f"{args.fullpolphase}_lr.gains",sep='\\s+').map(complex)
        llgains = pd.read_csv (f"{args.fullpolphase}_ll.gains",sep='\\s+').map(complex)

    if args.soldf:
        args.use_full = True

        soltag = args.soldf[:-len("_df.pkl")]
        soltag += "_soldf"

        sf     = pd.read_pickle(args.soldf)

        rrgains = sf.gains.loc[:,'rr']
        rlgains = sf.gains.loc[:,'rl']
        lrgains = sf.gains.loc[:,'lr']
        llgains = sf.gains.loc[:,'ll']


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

    elif args.fullpolphase or args.soldf:

        gf   = df.groupby(level=['ant1','ant2'])

        conj = np.conjugate

        if args.use_full:
            print ("performing full inversion ... ")

            for ((ant1,ant2), igf) in gf:

                ## gains
                gprr, gprl, gplr, gpll  = rrgains[ant1], rlgains[ant1], lrgains[ant1], llgains[ant1]
                gqrr, gqrl, gqlr, gqll  = rrgains[ant2], rlgains[ant2], lrgains[ant2], llgains[ant2]

                ## data
                dpqrr = np.array ( igf.complex.loc[ant1, ant2, 'rr'] )
                dpqrl = np.array ( igf.complex.loc[ant1, ant2, 'rl'] )
                dpqlr = np.array ( igf.complex.loc[ant1, ant2, 'lr'] )
                dpqll = np.array ( igf.complex.loc[ant1, ant2, 'll'] )

                # VV    = np.array ( [[dpqrr, dpqrl], [dpqlr, dpqll]] )
                # JP    = np.array ( [[gprr, gprl],[gplr, gpll]] )
                # JQ    = np.array ( [[gqrr, gqrl],[gqlr, gqll]] )
                # AA    = np.zeros_like ( VV )

                # iJP   = np.zeros_like (JP)
                # iJQH  = np.zeros_like (JQ)
                # for i in range(JP.shape[2]):
                    # __iJp   = np.zeros_like (JP[...,2])
                    # __iJqH  = np.zeros_like (JP[...,2])
                    # try:
                        # __iJp   = nl.inv ( JP[...,i] )
                    # except:
                        # __iJp[0,0]  =  1.0
                        # __iJp[1,1]  =  1.0
                    # try:
                        # __iJqH  = nl.inv ( JQ[...,i] ).transpose().conjugate()
                    # except:
                        # __iJqH[0,0] = 1.0
                        # __iJqH[1,1] = 1.0
                    
                    # AA[...,i]   = nl.multi_dot ( (__iJp, VV[...,i], __iJqH) )

                # AA    = np.einsum('abk,bck,cdk -> adk', iJP, VV, iJQH)

                """

                math from math_application

                ignore the determinant terms

                | gpll  -gprl | | Drr  Drl | | cgqll -cgqlr |
                | -gplr  gprr | | Dlr  Dll | |-cgqrl  cgqrr |

                rr term
                | gpll*Drr - gprl*Dlr   gpll*Drl - gprl*Dll |

                """
                ## rr
                unapplied = dpqrr
                # applied   = AA[0,0]
                applied   = dpqll*gprl*conj(gqrl)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                - dpqlr*gprl*conj(gqll)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                - dpqrl*gpll*conj(gqrl)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                + dpqrr*gpll*conj(gqll)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl))
                model     = gprl*conj(gqrl) + gprr*conj(gqrr)
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "rr" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )
                gg['model'].append ( model )
                gg['error_fit'].append ( np.abs(unapplied - model)**2 )
                gg['error_cal'].append ( np.abs(applied - 1.0)**2 )

                unapplied  = dpqrl
                # applied    = AA[0,1]
                applied    = -dpqll*gprl*conj(gqrr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        + dpqlr*gprl*conj(gqlr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        + dpqrl*gpll*conj(gqrr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        - dpqrr*gpll*conj(gqlr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl))
                model      = gprl*conj(gqll) + gprr*conj(gqlr)
                residual  = applied 
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "rl" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )
                gg['model'].append ( model )
                gg['error_fit'].append ( np.abs(unapplied - model)**2 )
                gg['error_cal'].append ( np.abs(applied)**2 )

                unapplied  = dpqlr
                # applied    = AA[1,0]
                applied    = -dpqll*gprr*conj(gqrl)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        + dpqlr*gprr*conj(gqll)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        + dpqrl*gplr*conj(gqrl)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        - dpqrr*gplr*conj(gqll)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl))
                model      = gpll*conj(gqrl) + gplr*conj(gqrr)
                residual  = applied 
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "lr" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )
                gg['model'].append ( model )
                gg['error_fit'].append ( np.abs(unapplied - model)**2 )
                gg['error_cal'].append ( np.abs(applied)**2 )

                unapplied  = dpqll
                # applied    = AA[1,1]
                applied    = dpqll*gprr*conj(gqrr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        - dpqlr*gprr*conj(gqlr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        - dpqrl*gplr*conj(gqrr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl)) \
                        + dpqrr*gplr*conj(gqlr)/(gpll*gprr*conj(gqll)*conj(gqrr) - gpll*gprr*conj(gqlr)*conj(gqrl) - gplr*gprl*conj(gqll)*conj(gqrr) + gplr*gprl*conj(gqlr)*conj(gqrl))
                model      = gpll*conj(gqll) + gplr*conj(gqlr)
                residual  = applied  - 1.0
                gg['ant1'].append ( ant1 )
                gg['ant2'].append ( ant2 )
                gg['correlation'].append ( "ll" )
                gg['unapplied'].append ( unapplied )
                gg['applied'].append ( applied )
                gg['model'].append ( model )
                gg['error_fit'].append ( np.abs(unapplied - model)**2 )
                gg['error_cal'].append ( np.abs(applied - 1.0)**2 )

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

    if args.gwb:
        """
        gwb solutions are already inverted.
        """

        soltag = args.gwb[:-len("_df.pkl")]

        sf     = pd.read_pickle(args.gwb)

        rgains = sf.gains.loc[:,'r']
        lgains = sf.gains.loc[:,'l']
        for row in df.itertuples():

            ## pull index
            ant1, ant2, (corr1, corr2) = row.Index

            gg['ant1'].append ( ant1 )
            gg['ant2'].append ( ant2 )
            gg['correlation'].append ( corr1+corr2 )

            gg['unapplied'].append ( row.complex )

            if corr1 == 'r' and corr2 == 'r':
                cal  = row.complex * rgains[ant1] * np.conjugate(rgains[ant2])
                gg['applied'].append ( cal )
            elif corr1 == 'r' and corr2 == 'l':
                cal  = row.complex * rgains[ant1] * np.conjugate(lgains[ant2])
                gg['applied'].append ( cal )
            elif corr1 == 'l' and corr2 == 'r':
                cal  = row.complex * lgains[ant1] * np.conjugate(rgains[ant2])
                gg['applied'].append ( cal )
            elif corr1 == 'l' and corr2 == 'l':
                cal  = row.complex * lgains[ant1] * np.conjugate(lgains[ant2])
                gg['applied'].append ( cal )
            else:
                raise RuntimeError("corr not understood")

    gg  = pd.DataFrame(gg).set_index(['ant1','ant2','correlation']).sort_index()

    odf  = DF[:-len("bldata_df.pkl")] + soltag + ".bldata_df.pkl"
    print(f"Writing to {odf}")

    gg.to_pickle(odf)
