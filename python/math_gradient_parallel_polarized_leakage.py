"""
math to compute gradient
with only parallel complex gains

p,q are baselines
m,n are polarization hands
"""
from collections import defaultdict
import numpy as np

import sys
import pandas as pd
from sympy import (symbols, Matrix, I, re, im, expand, conjugate, Derivative, latex, exp)
from sympy.printing.pycode import PythonCodePrinter

## gains
gp_rr, gp_ll = symbols('g_p^{rr} g_p^{ll}')
gq_rr, gq_ll = symbols('g_q^{rr} g_q^{ll}')
## xphase
phi          = symbols('phi', real=True)
## model
## this will be parang compensated
m_rr, m_rl, m_lr, m_ll     = symbols('m^{rr} m^{rl} m^{lr} m^{ll}')
## data
d_rr, d_rl, d_lr, d_ll     = symbols('d_{pq}^{rr} d_{pq}^{rl} d_{pq}^{lr} d_{pq}^{ll}')
## leakage    
pleakrl, pleaklr  = symbols('p_d^{rl} p_d^{lr}')
qleakrl, qleaklr  = symbols('q_d^{rl} q_d^{lr}')

### testsub
### populate random complex numbers into everything except phi
rx      = np.random.randn(40)
testsub = dict()
testsub[gp_rr] = rx[0] + 1.0j*rx[1]
testsub[gp_ll] = rx[6] + 1.0j*rx[7]
testsub[gq_rr] = rx[8] + 1.0j*rx[9]
testsub[gq_ll] = rx[14] + 1.0j*rx[15]
testsub[m_rr]  = rx[16] + 1.0j*rx[17]
testsub[m_rl]  = rx[18] + 1.0j*rx[19]
testsub[m_lr]  = rx[20] + 1.0j*rx[21]
testsub[m_ll]  = rx[22] + 1.0j*rx[23]
testsub[d_rr]  = rx[24] + 1.0j*rx[25]
testsub[d_rl]  = rx[26] + 1.0j*rx[27]
testsub[d_lr]  = rx[28] + 1.0j*rx[29]
testsub[d_ll]  = rx[30] + 1.0j*rx[31]
testsub[pleakrl]  = rx[32] + 1.0j*rx[33]
testsub[pleaklr]  = rx[34] + 1.0j*rx[35]
testsub[qleakrl]  = rx[36] + 1.0j*rx[37]
testsub[qleaklr]  = rx[38] + 1.0j*rx[39]

## matrices
gp  = Matrix([[gp_rr, 0],[0, gp_ll]])
gq  = Matrix([[gq_rr, 0],[0, gq_ll]])
xphase  = Matrix([[exp(I*phi), 0], [0, 1]])
mm  = Matrix ([[m_rr, m_rl], [m_lr, m_ll]])
pleak  = Matrix([[1, pleakrl],[pleaklr, 1]])
qleak  = Matrix([[1, qleakrl],[qleaklr, 1]])

## forward

herm  = lambda x : x.transpose().conjugate()

jones = gp * pleak * xphase

fff = gp  * pleak * xphase * mm * herm(xphase) * herm(qleak) * herm(gq)

frr  = expand ( fff[0,0] )
frl  = expand ( fff[0,1] )
flr  = expand ( fff[1,0] )
fll  = expand ( fff[1,1] )

## residuals
err  = d_rr - frr
erl  = d_rl - frl
elr  = d_lr - flr
ell  = d_ll - fll

## norm
srr  = expand ( conjugate(err) * err )
srl  = expand ( conjugate(erl) * erl )
slr  = expand ( conjugate(elr) * elr )
sll  = expand ( conjugate(ell) * ell )

### sum of all 
sss   = srr + srl + slr + sll
# sss  = expand ( conjugate(eee) * eee )
# arger = lambda t,a,b,c,d : sss.diff(t).collect( (Derivative(conjugate(t),t), conjugate(a), conjugate(b), conjugate(c), conjugate(d)), evaluate=False)

## i want to print cxxcode
## using python because the expressing is simple
## using std::conj;
print_settings  = PythonCodePrinter (settings={'user_functions':{'conjugate':'conj'}})
# printer     = lambda p : print_settings.doprint(p).replace("_","").replace("^","")
def printer (p):
    pp = print_settings.doprint(p)
    pp = pp.replace("_","")
    pp = pp.replace("^","")
    pp = pp.replace('math.exp(1j*phi)','zp')
    pp = pp.replace('math.exp(2*1j*phi)','zp2')
    pp = pp.replace('math.exp(-1j*phi)','czp')
    pp = pp.replace('math.exp(-2*1j*phi)','czp2')
    pp = pp.replace('1j','I')
    return pp

print ("---FORWARD---")
print ("RR", printer(frr), sep='\n' )
print ("RL", printer(frl), sep='\n' )
print ("LR", printer(flr), sep='\n' )
print ("LL", printer(fll), sep='\n' )
print ("---FORWARD---")

## replace sub/super scripts because my variable names do not have them
def arger (t, terms, popdkey=True): 
    """
    Derivative(conjugate(z),z) term is ignored
    Wirtinger derivative

    m_rr, m
    """
    cerms = [conjugate(t) for t in terms]

    ###
    dkey  = Derivative(conjugate(t),t)
    kv = sss.diff(t).collect( 
        (
            dkey,
            *cerms
        ), 
    evaluate=False)
    if popdkey: kv.pop ( dkey )
    ###
    tv    = {conjugate(k):conjugate(v) for k,v in kv.items()}
    lv    = "+".join([latex(k*v) for k,v in tv.items()])
    dv    = {printer(k):printer(v) for k,v in tv.items()}
    return lv,dv,tv

## collect terms
terms = [d_rr, d_rl, d_lr, d_ll, gp_rr, gp_ll]
p_rr  = arger ( gp_rr, terms )
p_ll  = arger ( gp_ll, terms )

terms = [conjugate(d_rr), conjugate(d_rl), conjugate(d_lr), conjugate(d_ll), gq_rr, gq_ll]
q_rr  = arger ( gq_rr, terms )
q_ll  = arger ( gq_ll, terms )

terms = [d_rr, d_rl, d_lr, d_ll, conjugate(d_rr), conjugate(d_rl), conjugate(d_lr), conjugate(d_ll), exp(I*phi), exp(-I*phi) ]
dphi  = arger ( phi, terms, popdkey=False )

### latex printing
## switching it off when generating code
print ("GPRR", p_rr[0], sep='\n')
print ("GPLL", p_ll[0], sep='\n')
print ("GQRR", q_rr[0], sep='\n')
print ("GQLL", q_ll[0], sep='\n')
print ("PHI",  dphi[0], sep='\n')

### coefficient printing
locs = defaultdict(list)
catamap = {
    'conj(dpqrr)':'dqprr',
    'conj(dpqrl)':'dqprl',
    'conj(dpqlr)':'dqplr',
    'conj(dpqll)':'dqpll',
}
def codeprintaction ( k, v, tag ):
    """
    k : dlr
    v : -gqrr*conj(gplr)

    tag : Drr

    Drr_coeff_{k} = v

    """
    ptag = tag
    pk   = k
    if pk in catamap.keys():
        pk = catamap[pk]

    stmt = f"const complex_type {ptag}_coeff_{pk} = {v} ; "

    locs[pk].append (stmt)

    # print ( stmt )

    # print (f"const complex_type {swapper(tag)}_coeff_{swapper(k)} = {swapper(v)} ; ")

# print ("----------   RR   ------------")
for k,v in p_rr[1].items(): codeprintaction ( k, v,"Dgprr" )
for k,v in q_rr[1].items(): codeprintaction ( k, v,"Dgqrr" )
# print ("----------   LL   ------------")
for k,v in p_ll[1].items(): codeprintaction ( k, v,"Dgpll" )
for k,v in q_ll[1].items(): codeprintaction ( k, v,"Dgqll" )
# print ("----------   PHI  ------------")
for k,v in dphi[1].items(): codeprintaction ( k, v,"Dphi" )

##########################

for k,v in locs.items():
    print ("-------------------------")
    # print (k)
    for iv in v: print (iv)

"""
derivative of sss wrt phi is real. 
we can test this by subs testsub (which puts random complex values for everything)
looking at the imag part.
It is numerically zero, the value is about ~1E-15 or so.
"""
