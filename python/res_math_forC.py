"""
same math as jac_math
The math is the same, but the implementation differs

in C++, we are running a loop over polar baselines. 
so, we want to compute jacobian given any two gains

so in one row there will only be 4 complex terms
and it should be consistent

compare with remath in objectives.C

"""
import sys
import pandas as pd
from sympy import (symbols, Matrix, I, re, im, expand)

### data
dr, di  = symbols('d^r d^i', real=True)

### model
mr, mi  = symbols('m^r m^i', real=True)

### gain of first polarbaseline
pr, pi  = symbols('p^r p^i', real=True)

### gain of second polarbaseline
qr, qi  = symbols('q^r q^i', real=True)
"""
no need to have p1 or q2. 
it is implied p is first and q is second
"""
## make them complex
dc  = dr + I*di
mc  = mr + I*mi
pc  = pr + I*pi
qc  = qr + I*qi
## helps to define conjugate of qc
qcc = qr - I*qi
##
"""

imodel = pc * mc * qcc
res    = data - imodel

directly write re and im parts
"""
imodel  = pc * mc * qcc
res     = dc - imodel

rres = re ( expand(res) )
ires = im ( expand(res) )


print (" real(res) ", rres, sep='\t')
print (" imag(res) ", ires, sep='\t')

"""
real(res) 	d^r + m^i*p^i*q^r - m^i*p^r*q^i - m^r*p^i*q^i - m^r*p^r*q^r
imag(res) 	d^i - m^i*p^i*q^i - m^i*p^r*q^r - m^r*p^i*q^r + m^r*p^r*q^i

"""

