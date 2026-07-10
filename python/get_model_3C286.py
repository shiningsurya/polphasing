"""
produces pol model to be used for solving gains

we use CASA convention to fit and extract

"""

import numpy as np
from scipy.optimize import curve_fit
import pandas as pd

import matplotlib.pyplot as plt

CAL     = "3C286"

# the following three parameters 
# decide the frequency axis
# will be matched against input LTA
# file to decide which model file to pick
FEDGE   = 550_000_000 # 550 MHz
FBW     = 97_656.25   # 97.65625 kHz
NCHAN   = 2048        # nchannels
CFILE   = f"../models/{CAL}_{NCHAN}_{FBW}_{FEDGE}_casa.txt"
MFILE   = f"../models/{CAL}_{NCHAN}_{FBW}_{FEDGE}_model.txt"
PFILE   = f"../models/{CAL}_{NCHAN}_{FBW}_{FEDGE}_model.png"

REFFREQ = FEDGE
REFFREQ_= f"{FEDGE/1E6:.4}MHz"
REFFREQ_GHZ  = REFFREQ / 1E9

FAXIS_GHZ  = ( FEDGE + ( ( np.arange(NCHAN) + 0.5 ) * FBW )  ) / 1E9

# REFFREQ = 551562500.0 # Hz
# REFFREQ_= "551.5625MHz"
# REFFREQ_GHZ = REFFREQ / 1E9 # f0 in GHz
##################################################

def CASA_poly ( freqs, c0, c1, c2, c3 ):
    """ casa polynomial (freqs in GHz)  """
    ## fraction freq.
    _x  = ( freqs - REFFREQ_GHZ ) / REFFREQ_GHZ
    return c0 + (c1*_x) + (c2*_x*_x) + (c3*_x*_x*_x)

def PB_Flux (f, *an):
    """ Perley+Butler flux polynomial (GHz)"""
    return np.power (10, an[0] + an[1]*np.log10(f) + an[2]*(np.log10(f))**2)

def CASA_alphabeta (freqs, Sa, alpha, beta ):
    """ Sa flux at f0 (freqs in GHz) """
    _x  = freqs / REFFREQ_GHZ
    return Sa* _x ** (alpha + beta*np.log10(_x) )

## from some meerkat paper
def evpa_3c286 ( fghz ):
    """
    https://arxiv.org/abs/2603.09001
    """
    x = np.log10(fghz)
    # print ( fghz, x )
    padeg = np.polyval ( [3_790, 615, 57, 26], x )
    return np.rad2deg(np.arctan(np.tan(np.deg2rad(padeg))))

perley_butler_flux_an    = np.array ( [1.2481, -0.4507, -0.1798, 0.0357] )
##################################################
## 
## from Perley+Butler (2013 and 2017)
perley_butler = pd.DataFrame ( dict(
    freq_ghz = np.array ( [1.050, 1.450, 1.640, 1.950] ),
    lp_percent = np.array ( [8.6, 9.5, 9.9, 10.1] ),
    pa_deg   = np.array ( [33., 33., 33., 33.] )
) )

perley_butler['lp_frac'] = perley_butler['lp_percent'] / 100.
perley_butler['flux_jy'] = perley_butler['freq_ghz'].apply ( lambda x : PB_Flux ( x, *perley_butler_flux_an ) )
##############################################################
## combined list
freq_ghz  = list(perley_butler['freq_ghz'])
flux_jy   = list(perley_butler['flux_jy'])
lp_frac   = list(perley_butler['lp_frac'])
# pa_deg    = list(perley_butler['pa_deg']) 

ALPHABETA,_    = curve_fit ( CASA_alphabeta, freq_ghz, flux_jy )
POL_INDICES,_  = curve_fit ( CASA_poly, freq_ghz, lp_frac )

POL_ANGLES,_   = curve_fit ( CASA_poly, FAXIS_GHZ, np.deg2rad(evpa_3c286(FAXIS_GHZ)) )

with open (CFILE, 'w') as f:
    f.write ( f"CAL\t{CAL}\n" )
    f.write ( f"REF_FREQ\t{REFFREQ_}\n" )
    f.write ( f"POL_INDICES\t{POL_INDICES}\n" )
    f.write ( f"POL_ANGLES\t{POL_ANGLES}\n" )
    f.write ( f"ALPHABETA\t{ALPHABETA}\n" )

print ('CAL', CAL, sep='\t')
print ('REF_FREQ', REFFREQ_, sep='\t')
print ('POL_INDICES', POL_INDICES, sep='\t')
print ('POL_ANGLES', POL_ANGLES, sep='\t')
print ('ALPHABETA', ALPHABETA, sep='\t')

########################################################
cf = dict()

_freq_mhz  = FAXIS_GHZ * 1E3

cf['freqs_mhz'] = _freq_mhz
cf['ichan']     = np.arange(NCHAN)

mflux   = CASA_alphabeta ( _freq_mhz*1E-3, ALPHABETA[0], *ALPHABETA[1:] )
mlp     = CASA_poly ( _freq_mhz*1E-3, *POL_INDICES )
mpa     = np.deg2rad ( evpa_3c286 (_freq_mhz*1E-3) )

cf['I'] = mflux
cf['V'] = mflux * 0.
cf['Q'] = mflux * mlp * np.cos (2.0*mpa)
cf['U'] = mflux * mlp * np.sin (2.0*mpa)

# cf      = pd.DataFrame(cf)
# cf.to_csv (MFILE, index=False)
with open(MFILE, 'w') as f:
    ## header
    st  = "{0:<9} {1:<9} {2:<9} {3:<9}\n".format ( "freqs", "stokes_i", "stokes_q", "stokes_u" )
    # print ( st, end='' )
    f.write ( st )
    for ichan in range (NCHAN):
        st  = "{0: ^9.6f} {1: ^9.6f} {2: ^9.6f} {3: ^9.6f}\n".format (
            _freq_mhz[ichan],
            cf['I'][ichan],
            cf['Q'][ichan],
            cf['U'][ichan],
        )
        # print ( st,end='' )
        f.write ( st )
########################################################
_pb  = dict(marker='D', c='k')
_cx  = dict(marker='o', c='blue')

fig  = plt.figure ('model')

axflux, axpfrac, axpa = fig.subplots ( 3, 1, sharex=True )


#####

axflux.plot ( _freq_mhz, mflux )

# axflux.scatter ( perley_butler.freq_ghz*1E3, perley_butler.flux_jy, **_pb )

axpfrac.plot ( _freq_mhz, mlp  )
# axpfrac.scatter ( perley_butler.freq_ghz*1E3, perley_butler.lp_frac, **_pb)

axpa.plot ( _freq_mhz, np.rad2deg ( CASA_poly ( _freq_mhz*1E-3, *POL_ANGLES ) ) )
# axpa.plot ( _freq_mhz, evpa_3c286 (_freq_mhz*1E-3) )
# axpa.scatter ( perley_butler.freq_ghz*1E3, perley_butler.pa_deg, **_pb)

axpa.set_ylabel ('PA / deg')
axflux.set_ylabel ('Flux / Jy')
axpfrac.set_ylabel ('Pfraction')
axpa.set_xlabel ('Freq / MHz')

f0  = REFFREQ_GHZ

axpa.axvline ( f0 * 1E3, ls='--', c='r' )
axpfrac.axvline ( f0 * 1E3, ls='--', c='r' )
axflux.axvline ( f0 * 1E3, ls='--', c='r' )

# plt.show ()
fig.savefig (PFILE, dpi=300, bbox_inches='tight')



