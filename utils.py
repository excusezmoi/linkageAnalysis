import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tabulate import tabulate

def getParameters():
    params = input("Enter the parameters, in the sequence of r2, r1, r4, r5, r6, r7, w2, a2(separate by space): ").split()
    keys = ["rv2", "rv1", "rv4", "rv5", "rv6", "rv7", "wv2", "av2"]
    return dict(zip(keys, map(float, params)))

def findValueT5T6(rv1,rv2,rv4,rv5,rv6,rv7):
    #find the desired initial values for t5 and t6
    
    x, y = sym.symbols('x, y')
    
    phi = sym.atan(rv2/rv1)
    xright = rv4*sym.sin(phi)
    yright = rv4*sym.cos(phi) - rv1
    xleft = -rv7
    _yleft = 0
    b = ((xright - xleft)**2 + yright**2)**0.5
    eq1 = sym.Eq((x+b)**2 + y**2 - rv6**2,0)
    eq2 = sym.Eq(x**2 + y**2 - rv5**2,0)
    result = sym.solve([eq1,eq2],(x,y))

    for sol in result:
        # print(sol[0], sol[1])
        if sol[0].is_real and sol[1].is_real:
            x, y = sol[0], abs(sol[1])
            break

    tp1 = sym.atan(y/x)

    #There are two circumstances. If the try one want work, the other will
    try:
        if x*y < 0:
            tp1 = sym.N(tp1+sym.pi)
            tp2 = sym.atan(y/(x+b))
        if y/(x+b) < 0:
            tp2 = sym.N(tp2+sym.pi)
        theta = sym.atan(yright/(xright-xleft))
        tp1 += theta
        tp2 += theta
        tp2 += sym.N(sym.pi)
        # print("try")

    except:
        t3min56 = sym.N(sym.pi - sym.atan(rv1 / rv7) - sym.acos((rv1**2 + rv7**2 + rv4**2 - (rv5 + rv6)**2) / (2 * sym.sqrt(rv1**2 + rv7**2) * rv4)))
        xright = rv4 * sym.cos(t3min56)
        yright = rv4 * sym.sin(t3min56) - rv1
        tp1 = sym.atan(yright/(xright + rv7))
        tp1 = sym.N(sym.pi + tp1)
        tp2 = tp1 + 0.1
        # print("except")

    return tp1, tp2