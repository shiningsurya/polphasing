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

### model
mr, mi  = symbols('m^r m^i', real=True)

### gain of first polarbaseline
tr, ti  = symbols('t^r t^i', real=True)

### gain of second polarbaseline
sr, si  = symbols('s^r s^i', real=True)
"""
no need to have t1 or s2. 
it is implied t is first and s is second
"""
## make them complex
mc  = mr + I*mi
tc  = tr + I*ti
sc  = sr + I*si
## helps to define conjugate of sc
scc = sr - I*si
##
"""
error = data - model

after taking derivative wrt gain terms,
data vanishes,
there is a minus that is implied and will be written while implementing
"""
res = tc * mc * scc

rres = re ( expand(res) )
ires = im ( expand(res) )


print (" real (res) / tr ", rres.diff(tr), sep='\t')
print (" imag (res) / tr ", ires.diff(tr), sep='\t')

print (" real (res) / ti ", rres.diff(ti), sep='\t')
print (" imag (res) / ti ", ires.diff(ti), sep='\t')

print (" real (res) / sr ", rres.diff(sr), sep='\t')
print (" imag (res) / sr ", ires.diff(sr), sep='\t')

print (" real (res) / si ", rres.diff(si), sep='\t')
print (" imag (res) / si ", ires.diff(si), sep='\t')

print (" real (res) / mr ", rres.diff(mr), sep='\t')
print (" imag (res) / mr ", ires.diff(mr), sep='\t')

print (" real (res) / mi ", rres.diff(mi), sep='\t')
print (" imag (res) / mi ", ires.diff(mi), sep='\t')

"""
 real (res) / tr 	m^i*s^i + m^r*s^r
 imag (res) / tr 	m^i*s^r - m^r*s^i
 real (res) / ti 	-m^i*s^r + m^r*s^i
 imag (res) / ti 	m^i*s^i + m^r*s^r
 real (res) / sr 	-m^i*t^i + m^r*t^r
 imag (res) / sr 	m^i*t^r + m^r*t^i
 real (res) / si 	m^i*t^r + m^r*t^i
 imag (res) / si 	m^i*t^i - m^r*t^r

 real (res) / mr 	s^i*t^i + s^r*t^r
 imag (res) / mr 	-s^i*t^r + s^r*t^i
 real (res) / mi 	s^i*t^r - s^r*t^i
 imag (res) / mi 	s^i*t^i + s^r*t^r

"""

