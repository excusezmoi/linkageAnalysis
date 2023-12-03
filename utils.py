import sympy as sym
import numpy as np

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
    # print(result)

    for sol in result:
        # print(sol[0], sol[1])
        if sol[0].is_real and sol[1].is_real:
            x, y = sol[0], abs(sol[1])
            break

    t5Ini = sym.atan(y/x)

    #There are two circumstances. If the try one won't work, the other will
    try:
        if x*y < 0:
            t5Ini = sym.N(t5Ini+sym.pi)
            t6Ini = sym.atan(y/(x+b))
        if y/(x+b) < 0:
            t6Ini = sym.N(t6Ini+sym.pi)
        theta = sym.atan(yright/(xright-xleft))
        t5Ini += theta
        t6Ini += theta
        t6Ini += sym.N(sym.pi)
        # print("try")

    except:
        t3min56 = sym.N(sym.pi - sym.atan(rv1 / rv7) - sym.acos((rv1**2 + rv7**2 + rv4**2 - (rv5 + rv6)**2) / (2 * sym.sqrt(rv1**2 + rv7**2) * rv4)))
        xright = rv4 * sym.cos(t3min56)
        yright = rv4 * sym.sin(t3min56) - rv1
        t5Ini = sym.atan(yright/(xright + rv7))
        t5Ini = sym.N(sym.pi + t5Ini)
        t6Ini = t5Ini - 0.1
        # print("except")

    return t5Ini, t6Ini

def t3Limits(rv1, rv2, rv4, rv5, rv6, rv7):
    t3max123 = sym.N(sym.pi / 2 + sym.asin(rv2 / rv1))
    t3min123 = sym.N(sym.pi / 2 - sym.asin(rv2 / rv1))

    if rv5 == rv6:
        t3max56 = sym.N(sym.pi - sym.atan(rv1 / rv7))
    else:
        t3max56 = sym.N(sym.pi - sym.atan(rv1 / rv7) - sym.acos((rv1**2 + rv7**2 + rv4**2 - (rv5 - rv6)**2) / (2 * sym.sqrt(rv1**2 + rv7**2) * rv4)))

    t3min56 = sym.N(sym.pi - sym.atan(rv1 / rv7) - sym.acos((rv1**2 + rv7**2 + rv4**2 - (rv5 + rv6)**2) / (2 * sym.sqrt(rv1**2 + rv7**2) * rv4)))

    maxt3 = min(t3max123, t3max56) if t3max56.is_real else t3max123
    mint3 = max(t3min123, t3min56) if t3min56.is_real else t3min123

    return maxt3, mint3

def invalidAngle(tv3, maxt3, mint3):
    return round(tv3, 4) > round(maxt3,4) or round(tv3,4) < round(mint3,4)

def loopOneEquations(rv1, rv2, r3, tv2, t3, delr3, delt3):
    # r3, t3, t5, t6, delr3, delt3, delt5, delt6 = sym.symbols('r3, t3, t5, t6, delr3, delt3, delt5, delt6')
    f1 = rv2 * sym.cos(tv2) - r3 * sym.cos(t3)
    f2 = rv1 + rv2 * sym.sin(tv2) - r3 * sym.sin(t3)
    eq1 = sym.Eq(sym.diff(f1, r3) * delr3 + sym.diff(f1, t3) * delt3, -f1)
    eq2 = sym.Eq(sym.diff(f2, r3) * delr3 + sym.diff(f2, t3) * delt3, -f2)
    return f1, f2, eq1, eq2

def loopTwoEquations(rv4, rv5, rv6, rv7, rv8, tv3, t5, t6, delt5, delt6):
    # r3, t3, t5, t6, delr3, delt3, delt5, delt6 = sym.symbols('r3, t3, t5, t6, delr3, delt3, delt5, delt6')
    f3 = rv4 * sym.cos(tv3) + rv5 * sym.cos(t5) + rv6 * sym.cos(t6) + rv7
    f4 = rv4 * sym.sin(tv3) + rv5 * sym.sin(t5) + rv6 * sym.sin(t6) - rv8
    eq3 = sym.Eq(sym.diff(f3, t5) * delt5 + sym.diff(f3, t6) * delt6, -f3)
    eq4 = sym.Eq(sym.diff(f4, t5) * delt5 + sym.diff(f4, t6) * delt6, -f4)
    return f3, f4, eq3, eq4

def zeroToTwoPi(angle):
    return angle % (2 * sym.pi)

# def NMsolveEquations(fA, fB, eqA, eqB, X, Y, vX, vY, delX, delY, normalize=[False, False]):
def NMsolveEquations(fA, fB, eqA, eqB, X, Y, vX, vY, delX, delY, normalize=[False, False]):
    result = sym.solve([eqA, eqB], (delX, delY))
    delvX, delvY = 0, 0
    while abs(sym.N(fA.subs({X: vX, Y: vY}))) > 0.0001 or abs(sym.N(fB.subs({X: vX, Y: vY}))) > 0.0001:
        delvX = sym.N(result[delX].subs({X: vX, Y: vY}))
        delvY = sym.N(result[delY].subs({X: vX, Y: vY}))
        vX += delvX
        vY += delvY
        vX, vY = [zeroToTwoPi(var) if norm else var for var, norm in zip([vX, vY], normalize)]
    return vX, vY

def calculatePositions(rv1, rv3, rv4, rv5, tv3, tv5):
    xr3 = rv3 * sym.cos(tv3)
    yr3 = -rv1 + rv3 * sym.sin(tv3)
    xr5 = rv4 * sym.cos(tv3)  # crankshaft x-position
    yr5 = rv4 * sym.sin(tv3) - rv1  # crankshaft y-position
    xr6 = rv4 * sym.cos(tv3) + rv5 * sym.cos(tv5)  # connecting rod x-position
    yr6 = rv4 * sym.sin(tv3) + rv5 * sym.sin(tv5) - rv1  # connecting rod y-position
    xr3, yr3, xr5, yr5, xr6, yr6 = [*map(lambda x: np.float64(sym.N(x)), [xr3, yr3, xr5, yr5, xr6, yr6])]
    return xr3, yr3, xr5, yr5, xr6, yr6