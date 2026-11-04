"""
math for levenburg marquardt
writing J to give to lmder
solving for everything at once

iterative
but iterating one gain type at a time
"""
from collections import defaultdict

import sys
import pandas as pd
import sympy as sp
from sympy import (cse, simplify, symbols, Matrix, I, re, im, expand, conjugate, Derivative, latex, cos, sin, exp, conjugate)
from sympy.printing.pycode import PythonCodePrinter
# from sympy.printing.ccode import CCodePrinter

"""
jacobian is size complex(number of polar baselines, number of gain terms)
or twice of each dimension in real.

The polar baseline order is same as LTA file, since we will be reading in the same order.
The gain terms ordering is {rr, rl, lr, ll}. 

For each gain
u denotes real part of gain
v denotes imag part of gain
p denotes first antenna in the baseline
q denotes second antenna in the baseline

Each term appearing in residual is of the form, 
{u,v}{p,q}{rr,rl,lr,ll}

So, uqrr is the real (rr gain of the second antenna)
"""

## re/im of gains
uprr, vprr     = symbols('uprr vprr', real=True)
uprl, vprl     = symbols('uprl vprl', real=True)
uplr, vplr     = symbols('uplr vplr', real=True)
upll, vpll     = symbols('upll vpll', real=True)

uqrr, vqrr     = symbols('uqrr vqrr', real=True)
uqrl, vqrl     = symbols('uqrl vqrl', real=True)
uqlr, vqlr     = symbols('uqlr vqlr', real=True)
uqll, vqll     = symbols('uqll vqll', real=True)

gainterms = {
    'p_real_rr' : uprr,
    'p_imag_rr' : vprr,

    'p_real_rl' : uprl,
    'p_imag_rl' : vprl,

    'p_real_lr' : uplr,
    'p_imag_lr' : vplr,

    'p_real_ll' : upll,
    'p_imag_ll' : vpll,

    'q_real_rr' : uqrr,
    'q_imag_rr' : vqrr,

    'q_real_rl' : uqrl,
    'q_imag_rl' : vqrl,

    'q_real_lr' : uqlr,
    'q_imag_lr' : vqlr,

    'q_real_ll' : uqll,
    'q_imag_ll' : vqll,
}

rrgainterms = {
    'p_real_rr' : uprr,
    'p_imag_rr' : vprr,

    'q_real_rr' : uqrr,
    'q_imag_rr' : vqrr,
}

llgainterms = {
    'p_real_ll' : upll,
    'p_imag_ll' : vpll,

    'q_real_ll' : uqll,
    'q_imag_ll' : vqll,
}
rlgainterms = {
    'p_real_rl' : uprl,
    'p_imag_rl' : vprl,

    'q_real_rl' : uqrl,
    'q_imag_rl' : vqrl,
}
lrgainterms = {
    'p_real_lr' : uplr,
    'p_imag_lr' : vplr,

    'q_real_lr' : uqlr,
    'q_imag_lr' : vqlr,
}

rllrgainterms = {
    'p_real_rl' : uprl,
    'p_imag_rl' : vprl,

    'p_real_lr' : uplr,
    'p_imag_lr' : vplr,

    'q_real_rl' : uqrl,
    'q_imag_rl' : vqrl,

    'q_real_lr' : uqlr,
    'q_imag_lr' : vqlr,
}

## parang of p,q
chip, chiq     = symbols('chi_p chi_q', real=True)

zp  = cos(chip) + I*sin(chip)
zq  = cos(chiq) + I*sin(chiq)
czp = cos(chip) - I*sin(chip)
czq = cos(chiq) - I*sin(chiq)

## gains
gprr   = uprr + I*vprr
gprl   = uprl + I*vprl
gplr   = uplr + I*vplr
gpll   = upll + I*vpll

gqrr   = uqrr + I*vqrr
gqrl   = uqrl + I*vqrl
gqlr   = uqlr + I*vqlr
gqll   = uqll + I*vqll
## conjugates of gq
cgqrr   = uqrr - I*vqrr
cgqrl   = uqrl - I*vqrl
cgqlr   = uqlr - I*vqlr
cgqll   = uqll - I*vqll

"""
-----MODELFORWARDPASS----
RR
gprl*zp*conj(gqrl)*conj(zq) + gprr*zq*conj(gqrr)*conj(zp)
RL
gprl*zp*conj(gqll)*conj(zq) + gprr*zq*conj(gqlr)*conj(zp)
LR
gpll*zp*conj(gqrl)*conj(zq) + gplr*zq*conj(gqrr)*conj(zp)
LL
gpll*zp*conj(gqll)*conj(zq) + gplr*zq*conj(gqlr)*conj(zp)
-----MODELFORWARDPASS----

rr: RR, RL, LR
rl: RR, RL, LR
lr:     RL, LR, LL
ll:     RL, LR, LL

i am only considering the model terms here

the negative sign here because of 
data - model

ibl is polarbaseline
jbl is baseline floor (ibl / 4)
"""

res_rr  = -gprl*zp*cgqrl*czq - gprr*zq*cgqrr*czp

res_rl  = -gprl*zp*cgqll*czq - gprr*zq*cgqlr*czp

res_lr  = -gpll*zp*cgqrl*czq - gplr*zq*cgqrr*czp

res_ll  = -gpll*zp*cgqll*czq - gplr*zq*cgqlr*czp

cresterms = {
    'res_rr' : res_rr,
    'res_rl' : res_rl,
    'res_lr' : res_lr,
    'res_ll' : res_ll,
}

resterms = dict()
for t,e in cresterms.items():
    _real, _imag = re(e), im(e)
    resterms[f"{t}_real"] = _real
    resterms[f"{t}_imag"] = _imag

def action( t, d ):
    """
    expand, take derivative, simplify
    """
    return simplify ( expand(t).diff(d) )

## derivaties

ibl, ldfjac, ip, iq = symbols('jbl ldfjac ip iq', real=True)
"""
residual start index is jbl

| jbl, jbl+1 | jbl+2, jbl+3 | jbl+4, jbl+5 |

every jacobian is shape ( 2 * 3 * baselines, 2 * nantenna )

three correlations. each correlation is complex. 
every antenna has one complex gain. 
"""

def populate_jac ( restdict, gaindict ):
    jindices, jexps = list(), list()
    for _res, _resexp in restdict.items():
        ## for every residual term
        ## _res is res_{real,imag}_{rr,rl,lr,ll}
        """
        jac.shape = (2*3*nbaselines, 2*nant)
        """
        # print (f"----------------------------------")
        # print (f"           {_res}                 ")
        # print (f"----------------------------------")
        rmap = {k:ibl+i for i,k in enumerate ( restdict.keys() )}
        # print ("serialized +1", rmap, sep='\n')

        ridx = rmap [ _res ]

        ## for every gainterm
        ## _gain is {p,q}_{real,imag}_{rr,rl,lr,ll}
        """
        p antenna index -> ip
        q antenna index -> iq
              | re | im | 
        2ip + |  0 |  1 |
        2iq + |  0 |  1 |

        but in iterative, we need (rr), (ll) and (rl,lr)
        """
        gmap = {
            'p_real_rr' : 2*ip,
            'p_imag_rr' : 2*ip + 1,

            'p_real_rl' : 2*ip,
            'p_imag_rl' : 2*ip + 1,

            'p_real_lr' : 2*ip,
            'p_imag_lr' : 2*ip + 1,

            'p_real_ll' : 2*ip,
            'p_imag_ll' : 2*ip + 1,

            'q_real_rr' : 2*iq,
            'q_imag_rr' : 2*iq + 1,

            'q_real_rl' : 2*iq,
            'q_imag_rl' : 2*iq + 1,

            'q_real_lr' : 2*iq,
            'q_imag_lr' : 2*iq + 1,

            'q_real_ll' : 2*iq,
            'q_imag_ll' : 2*iq + 1,
        }
        for _gain, _gainexp in gaindict.items():
            gidx  = gmap [ _gain ]
            
            ##
            _dexp = action ( _resexp, _gainexp )
            ##
            # print ( _res, _gain, _dexp )
            # print ( _res, _gain, _dexp )
            if _dexp.is_zero : continue

            """
            jac shape according to MINPACK is (npar, ndata)
            or
            (gainterms, residual terms)
            fastest index is polarbaseline index
            """

            fidx = expand(ridx + ldfjac*gidx)
            # print ( fidx, _dexp, sep='   ' )
            jindices.append ( fidx )
            jexps.append ( _dexp )
    ## do CSE and find 
    jidx  = cse(jindices, symbols=sp.numbered_symbols(prefix='idx'))
    jacx  = cse(jexps, symbols=sp.numbered_symbols(prefix='jac'))

    # print (f"-------expressions-------")
    print ()
    for a,b in jidx[0]:
        print (f"const int {a} ( {b} ) ;")
    print ()
    for a,b in jacx[0]:
        print (f"const real_type {a} ( {b} ) ;")
    # print (f"-----------fjac----------")
    print ()
    mat  = 'fjac'
    for ji, ja in zip ( jidx[1], jacx[1] ):
        print (f"{mat}[ {ji} ] = {ja} ;")

print (f"##########################################")
print (f"######  res_rr_real res_rr_imag    #######")
print (f"##########################################")
populate_jac ( {k:resterms[k] for k in ['res_rr_real', 'res_rr_imag', 'res_rl_real', 'res_rl_imag', 'res_lr_real', 'res_lr_imag'] }, rrgainterms )

print (f"##########################################")
print (f"###   res_rl_real, res_rl_imag         ###")
print (f"##########################################")
populate_jac ( {k:resterms[k] for k in ['res_rr_real', 'res_rr_imag', 'res_rl_real', 'res_rl_imag', 'res_lr_real', 'res_lr_imag'] }, rlgainterms )

print (f"##########################################")
print (f"###   res_lr_real, res_lr_imag         ###")
print (f"##########################################")
populate_jac ( {k:resterms[k] for k in ['res_rl_real', 'res_rl_imag', 'res_lr_real', 'res_lr_imag', 'res_ll_real', 'res_ll_imag'] }, lrgainterms )


print (f"##########################################")
print (f"######  res_ll_real res_ll_imag    #######")
print (f"##########################################")
populate_jac ( {k:resterms[k] for k in ['res_rl_real', 'res_rl_imag', 'res_lr_real', 'res_lr_imag', 'res_ll_real', 'res_ll_imag'] }, llgainterms )

# print (f"##########################################")
# print (f"###  res_rl,lr_real res_rl,lr_imag     ###")
# print (f"##########################################")
# populate_jac ( {k:resterms[k] for k in ['res_rl_real', 'res_rl_imag', 'res_lr_real', 'res_lr_imag'] }, rllrgainterms )

