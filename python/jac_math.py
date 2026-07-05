# coding: utf-8
import sys
import pandas as pd
from sympy import (symbols, Matrix, I, re, im, expand)
# import sympy as sp

## gains for both hands for two elements
# gp,gq,hp,hq = symbols('g_p g_q h_p h_q')
gpr,gpi,gqr,gqi = symbols('g_p^r g_p^i g_q^r g_q^i',real=True)
hpr,hpi,hqr,hqi = symbols('h_p^r h_p^i h_q^r h_q^i',real=True)
## complex gains
gp    = gpr + I*gpi
gq    = gqr + I*gqi
hp    = hpr + I*hpi
hq    = hqr + I*hqi
## gains matrix
gainp = Matrix([[gp,0],[0,hp]])
gainq = Matrix([[gq,0],[0,hq]])
## stokes def
SI,V,Q,U = symbols('I_S V Q U',real=True)
## coherence matrix
rho = Matrix([[SI+V,Q+I*U],[Q-I*U,SI-V]])
## i am going to set V=0
rho = rho.subs({V:0})
## if unpolarized source
# rho = rho.subs({V:0,Q:0,U:0})
## Radio interferometer measurement equation
rime = gainp @ rho @ gainq.conjugate().transpose()
####
print ( rime )
## elements of the observed visibility matrix
err  = expand ( rime[0,0] )
erl  = expand ( rime[0,1] )
elr  = expand ( rime[1,0] )
ell  = expand ( rime[1,1] )

dt   = [gpr, gpi, hpr, hpi, gqr, gqi, hqr, hqi]
dn   = ["gpr", "gpi", "hpr", "hpi","gqr", "gqi", "hqr", "hqi"]

from collections import defaultdict
df   = defaultdict(list)

def action (ee, tag):
    """

    """
    df['term'].append ( tag )
    for d,n in zip(dt,dn):
        df[n].append ( str( ee.diff(d).simplify() ) )

action ( re(err), "re(rr)")
action ( im(err), "im(rr)")
action ( re(erl), "re(rl)")
action ( im(erl), "im(rl)")
action ( re(elr), "re(lr)")
action ( im(elr), "im(lr)")
action ( re(ell), "re(ll)")
action ( im(ell), "im(ll)")

df = pd.DataFrame (df).set_index('term')

df.to_pickle ("jac_df.pkl")
df.to_html ("jac.html")

print ( df )

sys.exit (0)

rime
rime.diff(g_p)
rime.diff(gp)
re(rime)
re(rime[0,0])
re(rime[0,0]).diff(re(gp))
rime[0,0].diff(gp)
re(rime[0,0].diff(gp))
rime
rime[0,0]
rime[0,0].diff(gp)
rime[0,0].diff(gq)
rime[0,0].diff(gq).simplify()
rime[0,0]
rime[0,0].subs({gp:gpr+I*gpi,gq:gqr+I*gqi})
gpr
rime[0,0].subs({gp:gpr+I*gpi,gq:gqr+I*gqi})
r00 = rime[0,0].subs({gp:gpr+I*gpi,gq:gqr+I*gqi})
r00.simplify()
r00
r00 = r00.simplify()
r00
r00 = expand(r00)
r00
r00r, r00i = r00.as_real_imag()
r00r
r00i
r00r.diff(gpr)
r00r.diff(gpi)
r00i.diff(gpr)
r00i.diff(gpi)
r00 = expand(rime[0,0].subs({gp:gpr+I*gpi,gq:gqr+I*gqi}))
r00
rime
r00 = expand(rime[0,0].subs({gp:gpr+I*gpi,gq:gqr+I*gqi,hp:hpr+I*hpi,hq:hqr+I*hqi}))
r00
r01 = expand(rime[0,1].subs({gp:gpr+I*gpi,gq:gqr+I*gqi,hp:hpr+I*hpi,hq:hqr+I*hqi}))
r01
r01.diff(gpr)
re(r01).diff(gpr)
re(r01).diff(gpi)
re(r01).diff(hpr)
re(r01).diff(hpi)
re(r01).diff(hqr)
re(r01).diff(hqi)
im(r01)
r00 = expand(rime[0,0].subs({gp:gpr+I*gpi,gq:gqr+I*gqi,hp:hpr+I*hpi,hq:hqr+I*hqi}))
r00
rime
r00
re(r00)
re(r00).diff(gpr)
im(r00).diff(gpr)
re(r00).diff(gpr).simplify()
im(r00).diff(gpr).simplify()
re(r00).diff(gpi).simplify()
im(r00).diff(gpi).simplify()
re(r00).diff(gqr).simplify()
re(r00).diff(gqi).simplify()
im(r00).diff(gqr).simplify()
re(r00).diff(gqi).simplify()
im(r00).diff(gqi).simplify()
re(r00).diff(hpr).simplify()
r00
r00 = expand(rime[0,0].subs({gp:gpr+I*gpi,gq:gqr+I*gqi,hp:hpr+I*hpi,hq:hqr+I*hqi}))
re(r00)
re(r00).diff(gpr)
re(r00).diff(gpr).simplify()
re(r00).diff(gpi).simplify()
re(r00).diff(hpr).simplify()
re(r00).diff(hpi).simplify()
re(r00).diff(gqr).simplify()
re(r00).diff(gqi).simplify()
re(r00).diff(hpr).simplify()
re(r00).diff(hpi).simplify()
re(r00).diff(hqr).simplify()
re(r00).diff(hqi).simplify()
get_ipython().run_line_magic('pwd', '')
