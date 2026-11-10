"""
parallel leakages math
to compute gradient

With polarized source 

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
    def _print_conjugate(self, expr):
        ist = self._print(expr.args[0])
        return f"c{ist}"

        # Fall back to default behavior for all other powers
        return super()._print_Pow(expr)

gp_rr, gp_ll = symbols('g_p^{rr} g_p^{ll}')
gq_rr, gq_ll = symbols('g_q^{rr} g_q^{ll}')

lp_rl, lp_lr = symbols('l_p^{rl} l_p^{lr}')
lq_rl, lq_lr = symbols('l_q^{rl} l_q^{lr}')

gp_rl  = gp_rr * lp_rl
gp_lr  = gp_ll * lp_lr

gq_rl  = gq_rr * lq_rl
gq_lr  = gq_ll * lq_lr

m_rr, m_rl, m_lr, m_ll     = symbols('m^{rr} m^{rl} m^{lr} m^{ll}')
d_rr, d_rl, d_lr, d_ll     = symbols('d_{pq}^{rr} d_{pq}^{rl} d_{pq}^{lr} d_{pq}^{ll}')

gp   = Matrix([[gp_rr, gp_rl], [gp_lr, gp_ll]])
gq   = Matrix([[gq_rr, gq_rl], [gq_lr, gq_ll]])
gqh  = gq.transpose().conjugate()

mm   = Matrix([[m_rr, m_rl], [m_lr, m_ll]])

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
terms = [d_rr, d_rl, d_lr, d_ll, gp_rr, gp_ll]
p_rr  = arger ( gp_rr, terms )
p_ll  = arger ( gp_ll, terms )

terms = [conjugate(d_rr), conjugate(d_rl), conjugate(d_lr), conjugate(d_ll), gq_rr, gq_ll]
q_rr  = arger ( gq_rr, terms )
q_ll  = arger ( gq_ll, terms )

### latex printing
## switching it off when generating code
# print ("GPRR", p_rr[0], sep='\n')
# print ("GPRL", p_rl[0], sep='\n')
# print ("GPLR", p_lr[0], sep='\n')
# print ("GPLL", p_ll[0], sep='\n')
# print ("GQRR", q_rr[0], sep='\n')
# print ("GQRL", q_rl[0], sep='\n')
# print ("GQLR", q_lr[0], sep='\n')
# print ("GQLL", q_ll[0], sep='\n')

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
    'cdpqrr':'dqprr',
    'cdpqrl':'dqprl',
    'cdpqlr':'dqplr',
    'cdpqll':'dqpll',
}

groups = {
    'g1': ['dpqrr', 'dqprr', 'gprr', 'gqrr'], 
    'g2': ['dpqrl', 'dqprl', 'gprl', 'gqrl'], 
    'g3': ['dpqlr', 'dqplr', 'gplr', 'gqlr'], 
    'g4': ['dpqll', 'dqpll', 'gpll', 'gqll'], 
}
def codeprintaction ( k, v, c, tag ):
    """

    k: key
    v: val
    c: collect terms

    """
    ec   = v.collect(c, evaluate=False)

    tv   = ""
    for _ek, _ev in ec.items():

        __evar  = f"{tag}_{printer(k)}_{printer(_ek)}"
        stmt = f"const complex_type {__evar} = {printer(_ev)} ; "

        tv += f"{__evar}*{printer(_ek)} "

        print ( stmt )

    print (tv)

print ("######################################")
print ("need to edit norm(g[p,q][rr,ll]) terms")
print ("######################################")

print ("\n---------- GPRR   ------------\n")
# for k,v in p_rr[2].items(): codeprintaction ( k, v, "Dgprr" )

codeprintaction (d_rr, p_rr[2][d_rr], (gq_rr, gq_ll), "Dgprr")
codeprintaction (d_rl, p_rr[2][d_rl], (gq_rr, gq_ll), "Dgprr")
codeprintaction (gp_rr, p_rr[2][gp_rr], (gq_rr*conjugate(gq_rr), gq_ll*conjugate(gq_ll)), "Dgprr")

print ("\n---------- GQRR   ------------\n")

codeprintaction (conjugate(d_rr), q_rr[2][conjugate(d_rr)], (gp_rr, gp_ll), "Dgqrr")
codeprintaction (conjugate(d_lr), q_rr[2][conjugate(d_lr)], (gp_rr, gp_ll), "Dgqrr")
codeprintaction (gq_rr, q_rr[2][gq_rr], (gp_rr*conjugate(gp_rr), gp_ll*conjugate(gp_ll)), "Dgqrr")

print ("\n---------- GPLL   ------------\n")

codeprintaction (d_lr, p_ll[2][d_lr], (gq_rr, gq_ll), "Dgpll")
codeprintaction (d_ll, p_ll[2][d_ll], (gq_rr, gq_ll), "Dgpll")
codeprintaction (gp_ll, p_ll[2][gp_ll], (gq_rr*conjugate(gq_rr), gq_ll*conjugate(gq_ll)), "Dgpll")

print ("\n---------- GQLL   ------------\n")

codeprintaction (conjugate(d_rl), q_ll[2][conjugate(d_rl)], (gp_rr, gp_ll), "Dgqll")
codeprintaction (conjugate(d_ll), q_ll[2][conjugate(d_ll)], (gp_rr, gp_ll), "Dgqll")
codeprintaction (gq_ll, q_ll[2][gq_ll], (gp_rr*conjugate(gp_rr), gp_ll*conjugate(gp_ll)), "Dgqll")

##########################

# for k,v in locs.items():
    # print ("-------------------------")
    # print (k)
    # for iv in v: print (iv)


