"""
full leakages math
to compute model from solved Jones matrices and data
"""

import sys
import pandas as pd
from sympy import (symbols, Matrix, I, re, im, expand, conjugate, Derivative, latex)
from sympy.printing.pycode import PythonCodePrinter

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

## i want to print cxxcode
## using python because the expressing is simple
## using std::conj;
print_settings  = PythonCodePrinter (settings={'user_functions':{'conjugate':'conj'}})
printer     = lambda p : print_settings.doprint(p).replace("_","").replace("^","")
## replace sub/super scripts because my variable names do not have them
def arger (t): 
    """
    Derivative(conjugate(z),z) term is ignored
    Wirtinger derivative

    m_rr, m
    """
    terms = [m_rr, m_rl, m_lr, m_ll, d_rr, d_rl, d_lr, d_ll]
    cerms = [conjugate(t) for t in terms]

    ###
    kv = sss.diff(t).collect( 
        (
            Derivative(conjugate(t),t), 
            *cerms
        ), 
    evaluate=False)
    dkey  = Derivative(conjugate(t),t)
    kv.pop ( dkey )
    ###
    tv    = {conjugate(k):conjugate(v) for k,v in kv.items()}
    lv    = "+".join([latex(k*v) for k,v in tv.items()])
    dv    = {printer(k):printer(v) for k,v in tv.items()}
    return lv,dv

trr, vrr   = arger ( m_rr ) 
trl, vrl   = arger ( m_rl ) 
tlr, vlr   = arger ( m_lr ) 
tll, vll   = arger ( m_ll ) 


print ("============== LATEX PRINTING ==============")

print ("RR", trr, sep='\n')
print ("RL", trl, sep='\n')
print ("LR", tlr, sep='\n')
print ("LL", tll, sep='\n')

print ("============== CODE PRINTING ==============")

def codeprintaction ( k, v, tag ):
    """
    k : dlr
    v : -gqrr*conj(gplr)

    tag : Drr

    Drr_coeff_{k} = v
    """
    print (f"const complex_type {tag}_coeff_{k} = {v} ; ")


print ("----------   RR   ------------")
for k,v in vrr.items():
    codeprintaction ( k, v, "Drr" )
print ("----------   RL   ------------")
for k,v in vrl.items():
    codeprintaction ( k, v, "Drl" )
print ("----------   LR   ------------")
for k,v in vlr.items():
    codeprintaction ( k, v, "Dlr" )
print ("----------   LL   ------------")
for k,v in vll.items():
    codeprintaction ( k, v, "Dll" )
