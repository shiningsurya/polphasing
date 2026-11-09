"""
https://journals.uwyo.edu/index.php/ela/article/view/677

explicit polar decomposition

A = [[a,b],[c,d]] 

a,b,c,d are complex
"""
import sympy as sp

a,b,c,d = sp.symbols('a b c d')
A       = sp.Matrix([[a,b],[c,d]])

da      = sp.symbols('D_A', real=True)

theta   = sp.symbols('theta', real=True)

dagger  = lambda d : d.conjugate().transpose()

et      = sp.exp ( sp.I * theta )

# et      = A.det() / da

adjA_D  = dagger(A).adjugate()

zt      = A + et*adjA_D

zs      = (2*da + a*a.conjugate() + d*d.conjugate() + b*b.conjugate() + c.conjugate()*c) ** -0.5
