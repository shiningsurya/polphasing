"""
math to do parallactic angle correction
and leakage compensation

we require this as we solve for leakages separately and for parallel gains separately

p,q are baselines
m,n are polarization hands
"""

import sys
import pandas as pd
from sympy import (symbols, Matrix, I, re, im, expand, conjugate, Derivative, latex, cos, sin)
from sympy.printing.pycode import PythonCodePrinter

m_rr, m_rl, m_lr, m_ll     = symbols('m^{rr} m^{rl} m^{lr} m^{ll}')
###########################
# t1, t2                     = symbols('theta_1 theta_2', real=True)
# c1,c2         = cos(t1), cos(t2)
# s1,s2         = sin(t1), sin(t2)

# z1   = c1 + I*s1
# z2   = c2 + I*s2

# z1c  = c1 - I*s1
# z2c  = c2 - I*s2

z1, z2   = symbols('z_p z_q')

z1c  = z1.conjugate()
z2c  = z2.conjugate()

r1   = Matrix([[z1c, 0],[0, z1]])
r2   = Matrix([[z2c, 0],[0, z2]])

r2h  = r2.transpose().conjugate()

##########################
## leakage terms
dpr, dpl, dqr, dql = symbols('l_p^r l_p^l l_q^r l_q^l')
leak_p  = Matrix([[1, dpr],[dpl, 1]])
leak_q  = Matrix([[1, dqr],[dql, 1]])

leak_q_h = leak_q.transpose().conjugate()

##########################
mm   = Matrix([[m_rr, m_rl],[m_lr, m_ll]])

## leakage_p * par_p * model * par_q_H * leakage_q_H
out  = leak_p * r1 * mm * r2h * leak_q_h

##########################

orr  = expand(out[0,0])
orl  = expand(out[0,1])
olr  = expand(out[1,0])
oll  = expand(out[1,1])

## using python because the expressing is simple
## using std::conj;
print_settings  = PythonCodePrinter (settings={'user_functions':{'conjugate':'conj'}})
printer     = lambda p : print_settings.doprint(p).replace("_","").replace("^","")
## replace sub/super scripts because my variable names do not have them

print ( "RR", orr, sep='\n' )
print ( "RL", orl, sep='\n' )
print ( "LR", olr, sep='\n' )
print ( "LL", oll, sep='\n' )

print("-----------------")
print ("RR", printer(orr), sep='\n' )
print ("RL", printer(orl), sep='\n' )
print ("LR", printer(olr), sep='\n' )
print ("LL", printer(oll), sep='\n' )
