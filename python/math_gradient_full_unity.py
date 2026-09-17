"""
full leakages math
to compute gradient

With unpolarized source 

assume unity I and zero pol-stokes

p,q are baselines
m,n are polarization hands
"""

from collections import defaultdict

import sys
import pandas as pd
from sympy import (symbols, Matrix, I, re, im, expand, conjugate, Derivative, latex, cos, sin)
from sympy.printing.pycode import PythonCodePrinter
# from sympy.printing.ccode import CCodePrinter

class CustomPythonPrinter(PythonCodePrinter):
    def _print_Pow(self, expr):
        base, exp = expr.as_base_exp()
        # Intercept when the exponent is exactly 2
        if exp == 2:
            base_str = self._print(base)
            return f"({base_str} * {base_str})"
        elif exp == 3:
            base_str = self._print(base)
            return f"({base_str} * {base_str} * {base_str})"

        # Fall back to default behavior for all other powers
        return super()._print_Pow(expr)

gp_rr, gp_rl, gp_lr, gp_ll = symbols('g_p^{rr} g_p^{rl} g_p^{lr} g_p^{ll}')
gq_rr, gq_rl, gq_lr, gq_ll = symbols('g_q^{rr} g_q^{rl} g_q^{lr} g_q^{ll}')
# m_rr, m_rl, m_lr, m_ll     = symbols('m^{rr} m^{rl} m^{lr} m^{ll}')
# iunpol                     = symbols('I_{up}', real=True)
d_rr, d_rl, d_lr, d_ll     = symbols('d_{pq}^{rr} d_{pq}^{rl} d_{pq}^{lr} d_{pq}^{ll}')

## parang
## taken from math_parang.py
# t1, t2                     = symbols('theta_1 theta_2', real=True)
# c1,c2         = cos(t1), cos(t2)
# s1,s2         = sin(t1), sin(t2)

# c1, c2    = symbols('c_1 c_2', real=True)
# s1, s2    = symbols('s_1 s_2', real=True)

# z1   = c1 + I*s1
# z2   = c2 + I*s2

z1, z2 = symbols('z_1 z_2')

z1c  = conjugate(z1)
z2c  = conjugate(z2)

# z1c  = c1 - I*s1
# z2c  = c2 - I*s2

r1   = Matrix([[z1c, 0],[0, z1]])
r2   = Matrix([[z2c, 0],[0, z2]])

r2h  = r2.transpose().conjugate()

gp   = Matrix([[gp_rr, gp_rl], [gp_lr, gp_ll]])
gq   = Matrix([[gq_rr, gq_rl], [gq_lr, gq_ll]])
# mm   = Matrix([[m_rr, m_rl], [m_lr, m_ll]])
# mm   = Matrix([[iunpol, 0], [0, iunpol]])
mm   = r1 * Matrix([[1, 0], [0, 1]]) * r2h
gqh  = gq.transpose().conjugate()

fff  = gp * mm * gqh

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
# arger = lambda t,a,b,c,d : sss.diff(t).collect( (Derivative(conjugate(t),t), conjugate(a), conjugate(b), conjugate(c), conjugate(d)), evaluate=False)

## i want to print cxxcode
## using python because the expressing is simple
## using std::conj;
print_settings  = CustomPythonPrinter (settings={'user_functions':{'conjugate':'conj','sin':'sin', 'cos':'cos'}})
# print_settings  = CCodePrinter (settings={'user_functions':{'conjugate':'conj'}})
printer     = lambda p : print_settings.doprint(p).replace("_","").replace("^","")
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
    if popdkey:
        kv.pop ( dkey )
    ###
    tv    = {conjugate(k):conjugate(v) for k,v in kv.items()}
    lv    = "+".join([latex(k*v) for k,v in tv.items()])
    dv    = {printer(k):printer(v) for k,v in tv.items()}
    return lv,dv,tv

## collect terms
terms = [d_rr, d_rl, d_lr, d_ll, gp_rr, gp_rl, gp_lr, gp_ll]
p_rr  = arger ( gp_rr, terms )
p_rl  = arger ( gp_rl, terms )
p_lr  = arger ( gp_lr, terms )
p_ll  = arger ( gp_ll, terms )

terms = [conjugate(d_rr), conjugate(d_rl), conjugate(d_lr), conjugate(d_ll), gq_rr, gq_rl, gq_lr, gq_ll]
q_rr  = arger ( gq_rr, terms )
q_rl  = arger ( gq_rl, terms )
q_lr  = arger ( gq_lr, terms )
q_ll  = arger ( gq_ll, terms )

### latex printing
## switching it off when generating code
print ("GPRR", p_rr[0], sep='\n')
print ("GPRL", p_rl[0], sep='\n')
print ("GPLR", p_lr[0], sep='\n')
print ("GPLL", p_ll[0], sep='\n')
print ("GQRR", q_rr[0], sep='\n')
print ("GQRL", q_rl[0], sep='\n')
print ("GQLR", q_lr[0], sep='\n')
print ("GQLL", q_ll[0], sep='\n')

### model forward pass
print ("-----MODELFORWARDPASS----")

print ( "RR",  printer ( frr ), sep='\n' )
print ( "RL",  printer ( frl ), sep='\n' )
print ( "LR",  printer ( flr ), sep='\n' )
print ( "LL",  printer ( fll ), sep='\n' )

print ("-----MODELFORWARDPASS----")

### coefficient printing
locs = defaultdict(list)
catamap = {
    'conj(dpqrr)':'dqprr',
    'conj(dpqrl)':'dqprl',
    'conj(dpqlr)':'dqplr',
    'conj(dpqll)':'dqpll',
}

groups = {
    'g1': ['dpqrr', 'dqprr', 'gprr', 'gqrr'], 
    'g2': ['dpqrl', 'dqprl', 'gprl', 'gqrl'], 
    'g3': ['dpqlr', 'dqplr', 'gplr', 'gqlr'], 
    'g4': ['dpqll', 'dqpll', 'gpll', 'gqll'], 
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
    # if pk.startswith("conj(") and pk.endswith(")"):
        # swapper = lambda t : t.translate ( str.maketrans({'p':'q','q':'p'}) )
        # pk      = swapper(pk[len('conj('):-len(')')])

    stmt = f"const complex_type {ptag}_coeff_{pk} = {v} ; "

    for k,v in groups.items():
        if pk in v:
            locs[k].append ( stmt )
            continue


    # locs[pk].append (stmt)
    # locs[ptag].append (stmt)

    # print ( stmt )

    # print (f"const complex_type {swapper(tag)}_coeff_{swapper(k)} = {swapper(v)} ; ")

# print ("----------   RR   ------------")
for k,v in p_rr[1].items(): codeprintaction ( k, v,"Dgprr" )
for k,v in q_rr[1].items(): codeprintaction ( k, v,"Dgqrr" )
# print ("----------   RL   ------------")
for k,v in p_rl[1].items(): codeprintaction ( k, v,"Dgprl" )
for k,v in q_rl[1].items(): codeprintaction ( k, v,"Dgqrl" )
# print ("----------   LR   ------------")
for k,v in p_lr[1].items(): codeprintaction ( k, v,"Dgplr" )
for k,v in q_lr[1].items(): codeprintaction ( k, v,"Dgqlr" )
# print ("----------   LL   ------------")
for k,v in p_ll[1].items(): codeprintaction ( k, v,"Dgpll" )
for k,v in q_ll[1].items(): codeprintaction ( k, v,"Dgqll" )

##########################

for k,v in locs.items():
    print ("-------------------------")
    # print (k)
    for iv in v: print (iv)
