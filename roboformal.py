#This program uses symbolic calculation, while the iterative numerical methods to find the solution of non-linear equations are used.
import sympy as sym
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

#the symbols needed later, theta is represented as t, delxx is the calibrating value of the unknowns (xx) we want to estimate
#v stands for velocity, w stands for angular velocity (rad/s), a is the angular acceleration, except for ar3 being the acceleration of vector r3
r3, t3, t5, t6, delr3, delt3, delt5, delt6, v3, w3, w5, w6, ar3, at3, a5, a6, t, x, y= sym.symbols('r3, t3, t5, t6, delr3, delt3, delt5, delt6, v3, w3, w5, w6, ar3, at3, a5, a6, t, x, y')

# It is important to distinguish between the symbols and values in a symbolic calculation. The line above are the symbols, which are the unknows we aim to solve or we solved in the way;
# the line below is the ones being numeric.
rv2, rv1, rv4, rv5, rv6, rv7, wv2, av2 = [float(i) for i in input("Enter the parameters, in the sequence of r2, r1, r4, r5, r6, r7, w2, a2: ").split()] #Input of the given parameters
which = int(input("enter 1 if you want the angle position relation plot; enter 2 if you want the position of point D & E; enter 3 if you want velocity relation plot; else if you want acceleration plot: "))
rv8 = rv1
# a = np.full( (10,1),0.00001)
t2plot6 = np.full((36,1),0.00001)
t3plot = np.full((36,1),0.00001)
r3plot = np.full((36,1),0.00001)
t5plot = np.full((36,1),0.00001)
t6plot = np.full((36,1),0.00001)
dxplot = np.full((36,1),0.00001)
dyplot = np.full((36,1),0.00001)
explot = np.full((36,1),0.00001)
eyplot = np.full((36,1),0.00001)
evxplot = np.full((36,1),0.00001)
evyplot = np.full((36,1),0.00001)
eaxplot = np.full((36,1),0.00001)
eayplot = np.full((36,1),0.00001)
table = list()

#find the desired initial values for t5 and t6
phi = sym.atan(rv2/rv1)
xright = rv4*sym.sin(phi)
yright = rv4*sym.cos(phi) - rv1
xleft = -rv7
yleft = 0
b = ((xright - xleft)**2 + yright**2)**0.5
eq1 = sym.Eq((x+b)**2 + y**2 - rv6**2,0)
eq2 = sym.Eq(x**2 + y**2 - rv5**2,0)
result = sym.solve([eq1,eq2],(x,y))

for i in result:
    if i[0].is_real and i[1].is_real:
        x = i[0]
        y = abs(i[1])
tp1 = sym.atan(y/x)

#There are two circumstances.  If the first method won’t work, the second will
try:
    if x*y <0:
        tp1 = sym.N(tp1+sym.pi)
        tp2 = sym.atan(y/(x+b))
    if y/(x+b) <0:
        tp2 = sym.N(tp2+sym.pi)
    theta = sym.atan(yright/(xright-xleft))
    tp1 += theta
    tp2 += theta
    tp2 += sym.N(sym.pi)

except:
    t3min56 = sym.N(sym.pi - sym.atan(rv1 / rv7) - sym.acos((rv1**2 + rv7**2 + rv4**2 - (rv5 + rv6)**2) / (2 * sym.sqrt(rv1**2 + rv7**2) * rv4)))
    xright = rv4 * sym.cos(t3min56)
    yright = rv4 * sym.sin(t3min56) - rv1
    tp1 = sym.atan(yright/(xright + rv7))
    tp1 = sym.N(sym.pi + tp1)
    tp2 = tp1 + 0.1
#The initial guesses of the positions in Loop one are all ones.
rv3 = 1.0
tv3 = 1.0
tv5 = tp1
tv6 = tp2

#Whether the input link r2 can revolute 360 degree is depend on maximum and minimum angle of r3
#The maximum and minimum angles have two restricts: one by Loop1 and one by Loop2
#t3max123 is the maximum angle allowed by Loop 1, same for t3min123
t3max123 = sym.N(sym.pi / 2 + sym.asin(rv2 / rv1))
t3min123 = sym.N(sym.pi / 2 - sym.asin(rv2 / rv1))

#t3max56 is the maximum angle allowed by Loop 2, same for t3min56
if rv5 == rv6:
    t3max56 = sym.N(sym.pi - sym.atan(rv1 / rv7))
else:
    t3max56 = sym.N(sym.pi - sym.atan(rv1 / rv7) - sym.acos((rv1**2 + rv7**2 + rv4**2 - (rv5 - rv6)**2) / (2 * sym.sqrt(rv1**2 + rv7**2) * rv4)))
t3min56 = sym.N(sym.pi - sym.atan(rv1 / rv7) - sym.acos((rv1**2 + rv7**2 + rv4**2 - (rv5 + rv6)**2) / (2 * sym.sqrt(rv1**2 + rv7**2) * rv4)))

# this yields the possible moving range of r3
#Because of the numbers in sympy are all complex, the complex ones should be excluded
if t3max56.is_real:
    maxt3 = min(t3max123, t3max56)
else:
    maxt3 = t3max123

if t3min56.is_real:
    mint3 = max(t3min123, t3min56)
else:
    mint3 = t3min123
# print("maxt3:", maxt3, "mint3:", mint3)

#Note that if and only if maxt3 == t3max123 and mint3 ==t3min123, r3 can revolute 360 degree!

#take 35 values of input angles with an equal interval of 10 degree for the for loop
for i in range(1,37):
    tv2 = 10*i/180*sym.pi
    print("tv2 =", 10*i)
#     t2plot[i-1] = np.float64(10*i)

    eqt = sym.Eq(0.5 * av2 * t**2 + wv2 * t - sym.N(10/180*sym.pi),0)
    resultt = sym.solve(eqt,t)
    tv = max(resultt)
    wv2 = wv2 + av2 * tv
    # print("f")
    #Loop1
    #x-direction, we want to find solutions for f1 == 0 
    f1 = rv2 * sym.cos(tv2) - r3 * sym.cos(t3) # eq1
    #y-direction
    f2 = rv1 + rv2 * sym.sin(tv2) - r3 * sym.sin(t3) #eq2
    eq1 = sym.Eq(sym.diff(f1, r3) * delr3 + sym.diff(f1, t3) * delt3, -f1)
    eq2 = sym.Eq(sym.diff(f2, r3) * delr3 + sym.diff(f2, t3) * delt3, -f2)
    result1 = sym.solve([eq1, eq2], (delr3, delt3)) # the result contain the symbolic representation of delr3 and delt3
    delr3v = 0
    delt3v = 0
    # print("s")
    #everytime we adjust the value of r3 and t3, we plug them into eq1 and eq2 to find if f1 and f2 are both lesser than or equal to 0.0001
    while abs(sym.N(f1.subs({t3: tv3, r3: rv3}))) > 0.0001 or abs(sym.N(f2.subs({t3: tv3, r3: rv3}))) > 0.0001:
        delr3v = sym.N(result1[delr3].subs({t3:tv3, r3: rv3}))
        delt3v = sym.N(result1[delt3].subs({t3:tv3, r3: rv3}))    
        rv3 += delr3v
        tv3 += delt3v
    
    #if the resulted values of r3 or t3 are not in the limits, it cannot occur.
    if round(tv3, 4) > round(maxt3,4) or round(tv3,4) < round(mint3,4): #the round() is used to avoid the wrong result from the calculation above
        print(f"This input angle {10*i} degree cannot occur!")
        table.append([10*i,"N","N","N","N","N","N","N","N","N"])
        continue
    else:
        t2plot6[i-1] = np.float64(10*i)
        t3plot[i-1] = np.float64(tv3/sym.N(sym.pi)*180)
        r3plot[i-1] = np.float64(rv3*100)
        table.append([10*i])
#         print("rv3", rv3,"tv3:", tv3)
    # print("t")
    #Loop 2, basicallly the same as above
    f3 = rv4 * sym.cos(tv3) + rv5 * sym.cos(t5) + rv6 * sym.cos(t6) + rv7
    f4 = rv4 * sym.sin(tv3) + rv5 * sym.sin(t5) + rv6 * sym.sin(t6) - rv8
    eq3 = sym.Eq(sym.diff(f3, t5) * delt5 + sym.diff(f3, t6) * delt6, -f3)
    eq4 = sym.Eq(sym.diff(f4, t5) * delt5 + sym.diff(f4, t6) * delt6, -f4)
    result2 = sym.solve([eq3, eq4], (delt5, delt6))
    delt5v = 0
    delt6v = 0
    while abs(sym.N(f3.subs({t5: tv5, t6: tv6}))) > 0.0001 or abs(sym.N(f4.subs({t5: tv5, t6: tv6}))) > 0.0001:
        delt5v = sym.N(result2[delt5].subs({t5:tv5, t6: tv6}))
        delt6v = sym.N(result2[delt6].subs({t5:tv5, t6: tv6}))    
        tv5 += delt5v
        tv6 += delt6v
        
        #let the angles yielded always be in 0 to 2*pi
        while tv5 < 0:
            tv5 += sym.N(2 * sym.pi)
        while tv6 < 0:
            tv6 += sym.N(2 * sym.pi)
        while tv5 > 2 * sym.pi:
            tv5 -= sym.N(2 * sym.pi)
        while tv6 > 2 * sym.pi:
            tv6 -= sym.N(2 * sym.pi)
    t5plot[i-1] = np.float64(tv5/sym.N(sym.pi)*180)
    t6plot[i-1] = np.float64(tv6/sym.N(sym.pi)*180)
    # print("tv5 = ", tv5, "tv6 = ", tv6)

#     print("tv5:",tv5,"tv6:", tv6)
#     print(f"position D:({sym.N(rv4*sym.cos(tv3))},{sym.N(rv4*sym.sin(tv3)-rv1)})", f"position E:({sym.N(rv4*sym.cos(tv3)+rv5*sym.cos(tv5))},{sym.N(rv4*sym.sin(tv3)+rv5*sym.sin(tv5)-rv1)})")
    dx = round(sym.N(rv4*sym.cos(tv3)),4)
    dy = round(sym.N(rv4*sym.sin(tv3)-rv1), 4)
    ex = round(sym.N(rv4*sym.cos(tv3)+rv5*sym.cos(tv5)), 4)
    ey = round(sym.N(rv4*sym.sin(tv3)+rv5*sym.sin(tv5)-rv1), 4)
    dxplot[i-1] = np.float64(dx)
    dyplot[i-1] = np.float64(dy)
    explot[i-1] = np.float64(ex)
    eyplot[i-1] = np.float64(ey)
#     print("ey = ", ey)
    table[i-1].extend((dx, dy, ex, ey))
    # print("five")
#velocity
    eq1v = sym.Eq(sym.N(-rv2 * sym.sin(tv2) * wv2) - v3 * sym.cos(tv3) + rv3 * sym.sin(tv3) * w3, 0)
    eq2v = sym.Eq(sym.N(rv2 * sym.cos(tv2) * wv2) - v3 * sym.sin(tv3) - rv3 * sym.cos(tv3) * w3, 0)
    resultv1 = sym.solve([eq1v,eq2v],(v3,w3))
    wv3 = resultv1[w3]
    vv3 = resultv1[v3]
    eq3v = sym.Eq(-rv4 * sym.sin(tv3) * wv3 - rv5 * sym.sin(tv5) * w5 - sym.N(rv6 * sym.sin(tv6)) * w6, 0)
    eq4v = sym.Eq(rv4 * sym.cos(tv3) * wv3 + rv5 * sym.cos(tv5) * w5 + rv6 * sym.N(sym.cos(tv6)) * w6, 0)
    resultv2 = sym.solve([eq3v,eq4v],(w5,w6))
    wv5 = resultv2[w5]
    wv6 = resultv2[w6]
#     print(f"Ve:({-rv4 * sym.sin(tv3) * wv3 - rv5 * sym.sin(tv5) * wv5},{rv4 * sym.cos(tv3) * wv3 + rv5 * sym.cos(tv5) * wv5})")
    table[i-1].extend((round(-rv4 * sym.sin(tv3) * wv3 - rv5 * sym.sin(tv5) * wv5,4), round(rv4 * sym.cos(tv3) * wv3 + rv5 * sym.cos(tv5) * wv5,4)))
    if -rv4 * sym.sin(tv3) * wv3 - rv5 * sym.sin(tv5) * wv5 == 0:
        evxadd = 0.00001
        evxplot[i-1] = np.float64(evxadd)
    else:
        evxplot[i-1] = np.float64(round(-rv4 * sym.sin(tv3) * wv3 - rv5 * sym.sin(tv5) * wv5,4))
    if rv4 * sym.cos(tv3) * wv3 + rv5 * sym.cos(tv5) * wv5 == 0:
        evyadd = 0.00001
        evyplot[i-1] = np.float64(evyadd)
    else:
        evyplot[i-1] = np.float64(round(rv4 * sym.cos(tv3) * wv3 + rv5 * sym.cos(tv5) * wv5,4))
    # print("six")
#acceleration
    eq1a = sym.Eq(-rv2 * sym.N(sym.cos(tv2)) * wv2 ** 2 - rv2 * sym.N(sym.sin(tv2)) * av2 - ar3 * sym.cos(tv3) + 2 * (vv3 * sym.N(sym.sin(tv3)) * wv3) + rv3 * sym.N(sym.cos(tv3)) * wv3 ** 2 + rv3 * sym.sin(tv3) * at3, 0)
    eq2a = sym.Eq(-rv2 * sym.N(sym.sin(tv2)) * wv2 ** 2 + rv2 * sym.N(sym.cos(tv2)) * av2 - ar3 * sym.sin(tv3) - 2 * (vv3 * sym.N(sym.cos(tv3)) * wv3) + rv3 * sym.N(sym.sin(tv3)) * wv3 ** 2 - rv3 * sym.cos(tv3) * at3, 0)
    resulta1 = sym.solve([eq1a,eq2a],(ar3,at3))
    atv3 = resulta1[at3]
    arv3 = resulta1[ar3]
    eq3a = sym.Eq(-rv4 * sym.cos(tv3) * wv3 ** 2 - rv4 * sym.sin(tv3) * atv3 - rv5 * sym.cos(tv5) * wv5 ** 2 - rv5 * sym.sin(tv5) * a5 - rv6 * sym.cos(tv6) * wv6 ** 2 - rv6 * sym.sin(tv6) * a6, 0)
    eq4a = sym.Eq(-rv4 * sym.sin(tv3) * wv3 ** 2 + rv4 * sym.cos(tv3) * atv3 - rv5 * sym.sin(tv5) * wv5 ** 2 - rv5 * sym.cos(tv5) * a5 - rv6 * sym.sin(tv6) * wv6 ** 2 + rv6 * sym.cos(tv6) * a6, 0)
    resultav2 = sym.solve([eq3a,eq4a],(a5,a6))
    av5 = resultav2[a5]
    av6 = resultav2[a6]
#     print(f"Ae:({-rv4 * sym.cos(tv3) * wv3 ** 2 - rv4 * sym.sin(tv3) * atv3 - rv5 * sym.cos(tv5) * wv5 ** 2 - rv5 * sym.sin(tv5) * av5},{-rv4 * sym.sin(tv3) * wv3 ** 2 + rv4 * sym.cos(tv3) * atv3 - rv5 * sym.sin(tv5) * wv5 ** 2 - rv5 * sym.cos(tv5) * av5})")
    table[i-1].extend((round(-rv4 * sym.cos(tv3) * wv3 ** 2 - rv4 * sym.sin(tv3) * atv3 - rv5 * sym.cos(tv5) * wv5 ** 2 - rv5 * sym.sin(tv5) * av5,4),round( -rv4 * sym.sin(tv3) * wv3 ** 2 + rv4 * sym.cos(tv3) * atv3 - rv5 * sym.sin(tv5) * wv5 ** 2 - rv5 * sym.cos(tv5) * av5, 4)))
    table[i-1].append(round(wv6, 4))
    eaxplot[i-1] = np.float64(round(-rv4 * sym.cos(tv3) * wv3 ** 2 - rv4 * sym.sin(tv3) * atv3 - rv5 * sym.cos(tv5) * wv5 ** 2 - rv5 * sym.sin(tv5) * av5,4))
    eayplot[i-1] = np.float64(round( -rv4 * sym.sin(tv3) * wv3 ** 2 + rv4 * sym.cos(tv3) * atv3 - rv5 * sym.sin(tv5) * wv5 ** 2 - rv5 * sym.cos(tv5) * av5, 4))
    # print("seven")
t2plot6 = t2plot6[t2plot6 != 0.00001]
t3plot = t3plot[t3plot != 0.00001]
r3plot = r3plot[r3plot != 0.00001]
t5plot = t5plot[t5plot != 0.00001]
t6plot = t6plot[t6plot != 0.00001]
dxplot = dxplot[dxplot != 0.00001]
dyplot = dyplot[dyplot != 0.00001]
explot = explot[explot != 0.00001]
eyplot = eyplot[eyplot != 0.00001]
evxplot = evxplot[evxplot != 0.00001]
evyplot = evyplot[evyplot != 0.00001]
eaxplot = eaxplot[eaxplot != 0.00001]
eayplot = eayplot[eayplot != 0.00001]
if which == 1:
    # print(which)
    plt.plot(t2plot6, r3plot, label = "r3")
    plt.plot(t2plot6, t3plot, label = "theta 3")
    plt.plot(t2plot6, t5plot, label = "theta 5")
    plt.plot(t2plot6, t6plot, label = "theta 6")
    plt.xlabel("theta2(degree)")
    plt.ylabel("theta(degree)\n100 times r")
    plt.legend()
    plt.show()
else:
    fig = plt.figure()
    ax = fig.add_subplot(projection = "3d")
    if which == 2:
        ax.plot(t2plot6, dxplot, dyplot, label='dp')
        ax.plot(t2plot6, explot, eyplot, label='ep')
    elif which == 3:
        ax.plot(t2plot6, evxplot, evyplot, label='ev')
    else:
        ax.plot(t2plot6, eaxplot, eayplot, label='ea')
    ax.plot
    ax.set_xlabel("theta2(degree)")
    ax.set_ylabel("x")
    ax.set_zlabel("y")
    ax.legend()
    plt.show()
    
print(tabulate(table, headers = ["theta2","Dpx","Dpy", "Epx", "Epy", "Evx", "Evy", "Eax", "Eay", "w6"],floatfmt=".4f"))

#1 2 4 1.5 1.5 4.1 1 1