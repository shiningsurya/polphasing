"""

inverse applying

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

        # Fall back to default behavior for all other powers
        return super()._print_Pow(expr)

gp_rr, gp_rl, gp_lr, gp_ll = symbols('g_p^{rr} g_p^{rl} g_p^{lr} g_p^{ll}')
gq_rr, gq_rl, gq_lr, gq_ll = symbols('g_q^{rr} g_q^{rl} g_q^{lr} g_q^{ll}')
m_rr, m_rl, m_lr, m_ll     = symbols('m^{rr} m^{rl} m^{lr} m^{ll}')
iunpol                     = symbols('I_{up}', real=True)
d_rr, d_rl, d_lr, d_ll     = symbols('d_{pq}^{rr} d_{pq}^{rl} d_{pq}^{lr} d_{pq}^{ll}')

Jp   = Matrix([[gp_rr, gp_rl],[gp_lr, gp_ll]])
Jq   = Matrix([[gq_rr, gq_rl],[gq_lr, gq_ll]])

D    = Matrix([[d_rr, d_rl],[d_lr, d_ll]])

JqH  = Jq.transpose().conjugate()

BB   = Jp.inv() * D * JqH.inv()

brr  = expand (BB[0,0])
brl  = expand (BB[0,1])
blr  = expand (BB[1,0])
bll  = expand (BB[1,1])

class CustomPythonPrinter(PythonCodePrinter):
    def _print_Pow(self, expr):
        base, exp = expr.as_base_exp()
        # Intercept when the exponent is exactly 2
        if exp == 2:
            base_str = self._print(base)
            return f"({base_str} * {base_str})"

        # Fall back to default behavior for all other powers
        return super()._print_Pow(expr)

print_settings  = CustomPythonPrinter (settings={'user_functions':{'conjugate':'conj','sin':'sin', 'cos':'cos'}})
# print_settings  = CCodePrinter (settings={'user_functions':{'conjugate':'conj'}})
printer     = lambda p : print_settings.doprint(p).replace("_","").replace("^","")

print ("-- RR", printer(brr), "", sep='\n')
print ("-- RL", printer(brl), "", sep='\n')
print ("-- LR", printer(blr), "", sep='\n')
print ("-- LL", printer(bll), "", sep='\n')

