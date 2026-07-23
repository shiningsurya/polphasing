"""
full leakages math
for gradient descent update equations
"""

import sys
import pandas as pd
from sympy import (symbols, Matrix, I, re, im, expand, conjugate, Derivative, latex)

gp_rr, gp_rl, gp_lr, gp_ll = symbols('g_p^{rr} g_p^{rl} g_p^{lr} g_p^{ll}')
gq_rr, gq_rl, gq_lr, gq_ll = symbols('g_q^{rr} g_q^{rl} g_q^{lr} g_q^{ll}')
m_rr, m_rl, m_lr, m_ll     = symbols('m^{rr} m^{rl} m^{lr} m^{ll}')
d_rr, d_rl, d_lr, d_ll     = symbols('d^{rr} d^{rl} d^{lr} d^{ll}')

# einstein summation
#mn         ma     ab                bn
frr  =  (gp_rr * m_rr * conjugate(gq_rr)) +\
        (gp_rr * m_rl * conjugate(gq_lr)) +\
        (gp_rl * m_lr * conjugate(gq_rr)) +\
        (gp_rl * m_ll * conjugate(gq_lr))

#mn         ma     ab                bn
frl  =  (gp_rr * m_rr * conjugate(gq_rl)) +\
        (gp_rr * m_rl * conjugate(gq_ll)) +\
        (gp_rl * m_lr * conjugate(gq_rl)) +\
        (gp_rl * m_ll * conjugate(gq_ll))

#mn         ma     ab                bn
flr  =  (gp_lr * m_rr * conjugate(gq_rr)) +\
        (gp_lr * m_rl * conjugate(gq_lr)) +\
        (gp_ll * m_lr * conjugate(gq_rr)) +\
        (gp_ll * m_ll * conjugate(gq_lr))

#mn         ma     ab                bn
fll  =  (gp_lr * m_rr * conjugate(gq_rl)) +\
        (gp_lr * m_rl * conjugate(gq_ll)) +\
        (gp_ll * m_lr * conjugate(gq_rl)) +\
        (gp_ll * m_ll * conjugate(gq_ll))

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
# arger = lambda t,a,b,c,d : sss.diff(t).collect( (Derivative(conjugate(t),t), conjugate(a), conjugate(b), conjugate(c), conjugate(d)), evaluate=False)
def arger (t,a,b,c,d): 
    """
    c,d are the data terms
    a,b are gain terms that we need to estimate

    Derivative(conjugate(z),z) term is ignored
    Wirtinger derivative
    """
    ca = conjugate (a)
    cb = conjugate (b)
    cc = conjugate (c)
    cd = conjugate (d)
    ###
    kv = sss.diff(t).collect( (Derivative(conjugate(t),t), ca, cb, cc, cd), evaluate=False)
    ###
    aa = conjugate ( kv [ ca ] )
    bb = conjugate ( kv [ cb ] )
    cc = conjugate ( kv [ cc ] )
    dd = conjugate ( kv [ cd ] )
    ###
    term  = a*aa + b*bb + c*cc + d*dd
    return term

## collect terms
trr   = arger ( gp_rr, gp_rr, gp_rl, d_rr, d_rl ) 
trl   = arger ( gp_rl, gp_rr, gp_rl, d_rr, d_rl ) 
tlr   = arger ( gp_lr, gp_ll, gp_lr, d_ll, d_lr ) 
tll   = arger ( gp_ll, gp_ll, gp_lr, d_ll, d_lr ) 


print ("RR", latex(trr), sep='\n')
print ("RL", latex(trl), sep='\n')
print ("LR", latex(tlr), sep='\n')
print ("LL", latex(tll), sep='\n')

