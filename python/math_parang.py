"""
math to do parallactic angle correction

p,q are baselines
m,n are polarization hands
"""

import sys
import pandas as pd
from sympy import (symbols, Matrix, I, re, im, expand, conjugate, Derivative, latex, cos, sin)
from sympy.printing.pycode import PythonCodePrinter

m_rr, m_rl, m_lr, m_ll     = symbols('m^{rr} m^{rl} m^{lr} m^{ll}')
t1, t2                     = symbols('theta_1 theta_2', real=True)
c1,c2         = cos(t1), cos(t2)
s1,s2         = sin(t1), sin(t2)

###########################
z1   = c1 + I*s1
z2   = c2 + I*s2

z1c  = c1 - I*s1
z2c  = c2 - I*s2

r1   = Matrix([[z1c, 0],[0, z1]])
r2   = Matrix([[z2c, 0],[0, z2]])

r2h  = r2.transpose().conjugate()

mm   = Matrix([[m_rr, m_rl],[m_lr, m_ll]])

out  = r1 * mm * r2h


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
