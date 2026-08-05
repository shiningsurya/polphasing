"""
math to compute gradient
with only parallel complex gains

p,q are baselines
m,n are polarization hands
"""
from collections import defaultdict

import sys
import pandas as pd
from sympy import (symbols, Matrix, I, re, im, expand, conjugate, Derivative, latex)
from sympy.printing.pycode import PythonCodePrinter

gp_rr, gp_ll = symbols('g_p^{rr} g_p^{ll}')
gq_rr, gq_ll = symbols('g_q^{rr} g_q^{ll}')
m_rr, m_rl, m_lr, m_ll     = symbols('m^{rr} m^{rl} m^{lr} m^{ll}')
d_rr, d_rl, d_lr, d_ll     = symbols('d_{pq}^{rr} d_{pq}^{rl} d_{pq}^{lr} d_{pq}^{ll}')

# einstein summation

frr  = gp_rr * m_rr * conjugate(gq_rr)
frl  = gp_rr * m_rl * conjugate(gq_ll)
flr  = gp_ll * m_lr * conjugate(gq_rr)
fll  = gp_ll * m_ll * conjugate(gq_ll)

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
printer     = lambda p : print_settings.doprint(p).replace("_","").replace("^","")
## replace sub/super scripts because my variable names do not have them
def arger (t, terms): 
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
    kv.pop ( dkey )
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

### latex printing
## switching it off when generating code
print ("GPRR", p_rr[0], sep='\n')
print ("GPLL", p_ll[0], sep='\n')
print ("GQRR", q_rr[0], sep='\n')
print ("GQLL", q_ll[0], sep='\n')

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

##########################

for k,v in locs.items():
    print ("-------------------------")
    # print (k)
    for iv in v: print (iv)

