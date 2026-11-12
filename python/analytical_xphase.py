"""
analytical crosshand phase computation
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
    agp = argparse.ArgumentParser("analytical_xphase", description="Measures crosshandphase from formula")
    add = agp.add_argument
    add ('bldata_df', help='full gains applied on polarized quasar')
    add ('-m','--model', help='model for the polarized quasar', dest='model', required=True)
    add ('-p','--par', help='Parallactic angle gains _{pr,pl}.gains soldf', dest='pargains', required=True)
    return agp.parse_args()

if __name__ == "__main__":
    args = get_args()

    data  = pd.read_pickle( args.bldata_df )

    data  = data [ data.index.get_level_values('ant1') != data.index.get_level_values('ant2') ] 

    data  = data.loc[:,:,['rl','lr']]

    parf  = pd.read_pickle(args.pargains)
    model = pd.read_csv (args.model, sep='\\s+')

    gg   = defaultdict(list)

    MRL   = model['stokes_q'] + 1.0j*model['stokes_u']
    MLR   = model['stokes_q'] - 1.0j*model['stokes_u']

    for (ant1,ant2), row in data.applied.groupby(level=('ant1','ant2')):

        _pr   = parf.gains.loc[ant1, 'pr']
        _pl   = parf.gains.loc[ant1, 'pl']

        _cqr   = np.conjugate ( parf.gains.loc[ant2, 'pr'] )
        _cql   = np.conjugate ( parf.gains.loc[ant2, 'pl'] )

        # _pr   = np.conjugate ( parf.gains.loc[ant1, 'pr'] )
        # _pl   = np.conjugate ( parf.gains.loc[ant1, 'pl'] )

        # _cqr   =  parf.gains.loc[ant2, 'pr'] 
        # _cql   =  parf.gains.loc[ant2, 'pl'] 

        _mrl   = MRL * _pr * _cql
        _mlr   = MLR * _pl * _cqr

        _arl   = row.loc[ant1, ant2, 'rl'] * np.conjugate ( _mrl )
        _alr   = row.loc[ant1, ant2, 'lr'] * np.conjugate ( _mlr )

        _brl   = np.conjugate ( row.loc[ant1, ant2, 'rl'] ) * _mrl
        _blr   = np.conjugate ( row.loc[ant1, ant2, 'lr'] ) * _mlr

        gg['ant1'].append ( ant1 )
        gg['ant2'].append ( ant2 )
        gg['arl'].append ( np.array(_arl) )

        gg['mrl'].append ( np.array(_mrl) )
        gg['mlr'].append ( np.array(_mlr) )

        # gg['ant1'].append ( ant1 )
        # gg['ant2'].append ( ant2 )
        gg['alr'].append ( np.array( _alr ) )

        # gg['ant1'].append ( ant1 )
        # gg['ant2'].append ( ant2 )
        gg['brl'].append ( np.array( _brl ) )

        # gg['ant1'].append ( ant1 )
        # gg['ant2'].append ( ant2 )
        gg['blr'].append ( np.array( _blr ) )

    gg = pd.DataFrame (gg).set_index(['ant1','ant2'])
    gg.to_pickle("test_analytical_xphase.pkl")
