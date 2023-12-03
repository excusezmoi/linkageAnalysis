#This program uses symbolic calculation, while the iterative numerical methods to find the solution of non-linear equations are used.
import sympy as sym
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utils import *

def link3VelocityEquations(rv2, rv3, tv2, tv3, wv2, w3, v3):
    eq1v = sym.Eq(-rv2 * sym.sin(tv2) * wv2 - v3 * sym.cos(tv3) + rv3 * sym.sin(tv3) * w3, 0)
    eq2v = sym.Eq(rv2 * sym.cos(tv2) * wv2 - v3 * sym.sin(tv3) - rv3 * sym.cos(tv3) * w3, 0)
    return eq1v, eq2v

def link56AngularVelocityEquations(rv4, rv5, rv6, tv3, tv5, tv6, wv3, w5, w6):
    eq3v = sym.Eq(-rv4 * sym.sin(tv3) * wv3 - rv5 * sym.sin(tv5) * w5 - sym.N(rv6 * sym.sin(tv6)) * w6, 0)
    eq4v = sym.Eq(rv4 * sym.cos(tv3) * wv3 + rv5 * sym.cos(tv5) * w5 + rv6 * sym.N(sym.cos(tv6)) * w6, 0)
    return eq3v, eq4v

def link3AccelerationEquations(rv2, rv3, tv2, tv3, vv3, wv2, wv3, av2, ar3, at3):
    eq1a = sym.Eq(-rv2 * sym.N(sym.cos(tv2)) * wv2 ** 2 - rv2 * sym.N(sym.sin(tv2)) * av2 - ar3 * sym.cos(tv3) + 2 * (vv3 * sym.N(sym.sin(tv3)) * wv3) + rv3 * sym.N(sym.cos(tv3)) * wv3 ** 2 + rv3 * sym.sin(tv3) * at3, 0)
    eq2a = sym.Eq(-rv2 * sym.N(sym.sin(tv2)) * wv2 ** 2 + rv2 * sym.N(sym.cos(tv2)) * av2 - ar3 * sym.sin(tv3) - 2 * (vv3 * sym.N(sym.cos(tv3)) * wv3) + rv3 * sym.N(sym.sin(tv3)) * wv3 ** 2 - rv3 * sym.cos(tv3) * at3, 0)
    return eq1a, eq2a

def link56AngularAccelerationEquations(rv4, rv5, rv6, tv3, tv5, tv6, wv3, wv5, wv6, a5, a6, atv3):
    eq3a = sym.Eq(-rv4 * sym.cos(tv3) * wv3 ** 2 - rv4 * sym.sin(tv3) * atv3 - rv5 * sym.cos(tv5) * wv5 ** 2 - rv5 * sym.sin(tv5) * a5 - rv6 * sym.cos(tv6) * wv6 ** 2 - rv6 * sym.sin(tv6) * a6, 0)
    eq4a = sym.Eq(-rv4 * sym.sin(tv3) * wv3 ** 2 + rv4 * sym.cos(tv3) * atv3 - rv5 * sym.sin(tv5) * wv5 ** 2 - rv5 * sym.cos(tv5) * a5 - rv6 * sym.sin(tv6) * wv6 ** 2 + rv6 * sym.cos(tv6) * a6, 0)
    return eq3a, eq4a

def symSolveEquations(equations, variables):
    result = sym.solve(equations, variables)
    variables = [result[variable] for variable in variables]
    return variables

#take 35 values of input angles with an equal interval of 10 degree for the for loop
def inputIncrementAngles(num = 37) -> np.ndarray:
    return np.linspace(0, 360, num)[1:]/180*sym.pi

def npFullColumn(num = 36, value = 0.00001) -> np.ndarray:
    return np.full((num,1), value)

class PlotData:
    def __init__(self):
        (self.t2plot6, self.t3plot, self.r3plot, self.t5plot, self.t6plot, self.dxplot, 
        self.dyplot, self.explot, self.eyplot, self.evxplot, self.evyplot, self.eaxplot, 
        self.eayplot) = [npFullColumn() for _ in range(13)]

    def addData(self, step, rv1, rv3, rv4, rv5, tv3, tv5, tv6, wv3, wv5, av5, atv3):
        f = np.float64
        r = lambda x: round(x, 4)

        self.t2plot6[step] = f(10*(step+1))
        self.t3plot[step] = f(tv3/sym.N(sym.pi)*180)
        self.r3plot[step] = f(rv3*100)
        self.t5plot[step] = f(tv5/sym.N(sym.pi)*180)
        self.t6plot[step] = f(tv6/sym.N(sym.pi)*180)

        dx, dy = rv4 * sym.cos(tv3), rv4 * sym.sin(tv3) - rv1
        ex, ey = rv4 * sym.cos(tv3) + rv5 * sym.cos(tv5), rv4 * sym.sin(tv3) + rv5 * sym.sin(tv5) - rv1

        # for val, array in zip([dx, dy, ex, ey], [dxplot, dyplot, explot, eyplot]):
        #     array[step] = f(r(sym.N(val)))
        
        self.dxplot[step], self.dyplot[step], self.explot[step], self.eyplot[step] = [*map(lambda val: f(r(sym.N(val))), [dx, dy, ex, ey])]
        
        evx = -rv4 * sym.sin(tv3) * wv3 - rv5 * sym.sin(tv5) * wv5
        evy = rv4 * sym.cos(tv3) * wv3 + rv5 * sym.cos(tv5) * wv5

        self.evxplot[step] = f(r(evx) if evx != 0 else 0.00001)
        self.evyplot[step] = f(r(evy) if evy != 0 else 0.00001)

        eax = -rv4 * sym.cos(tv3) * wv3 ** 2 - rv4 * sym.sin(tv3) * atv3 - rv5 * sym.cos(tv5) * wv5 ** 2 - rv5 * sym.sin(tv5) * av5
        eay = -rv4 * sym.sin(tv3) * wv3 ** 2 + rv4 * sym.cos(tv3) * atv3 - rv5 * sym.sin(tv5) * wv5 ** 2 - rv5 * sym.cos(tv5) * av5

        self.eaxplot[step], self.eayplot[step] = f(r(eax)), f(r(eay))
        
        return dx, dy, ex, ey, evx, evy, eax, eay
    
    def filterUnchangedValue(self, array: np.ndarray, num: float) -> np.ndarray:
        return array[array != num]

    def cleanInvalidValues(self):
        plot_attributes = ['t2plot6', 't3plot', 'r3plot', 't5plot', 't6plot', 
                        'dxplot', 'dyplot', 'explot', 'eyplot', 'evxplot', 
                        'evyplot', 'eaxplot', 'eayplot']

        for attr in plot_attributes:
            setattr(self, attr, self.filterUnchangedValue(getattr(self, attr), 0.00001))
    
    def display(self, option):
        if option == 1:
            # print(option)
            plt.plot(self.t2plot6, self.r3plot, label = "r3")
            plt.plot(self.t2plot6, self.t3plot, label = "theta 3")
            plt.plot(self.t2plot6, self.t5plot, label = "theta 5")
            plt.plot(self.t2plot6, self.t6plot, label = "theta 6")
            plt.xlabel("theta2(degree)")
            plt.ylabel("theta(degree)\n100 times r")
            plt.legend()
            plt.show()
        else:
            fig = plt.figure()
            ax = fig.add_subplot(projection = "3d")
            if option == 2:
                ax.plot(self.t2plot6, self.dxplot, self.dyplot, label='dp')
                ax.plot(self.t2plot6, self.explot, self.eyplot, label='ep')
            elif option == 3:
                ax.plot(self.t2plot6, self.evxplot, self.evyplot, label='ev')
            else:
                ax.plot(self.t2plot6, self.eaxplot, self.eayplot, label='ea')
            ax.plot
            ax.set_xlabel("theta2(degree)")
            ax.set_ylabel("x")
            ax.set_zlabel("y")
            ax.legend()
            plt.show()

class Table:
    def __init__(self):
        self.table = pd.DataFrame(columns = ["theta2(deg)","Dpx","Dpy", "Epx", "Epy", "Evx", "Evy", "Eax", "Eay", "w6"])

    def addRow(self, step, valid, *values):
        if valid:
            self.table.loc[step] = [10*(step + 1), *values]
        else:
            self.table.loc[step] = [10*(step + 1), *["N"]*9]

    def display(self):
        print(self.table)

def kinematic(rv1,rv2,rv4,rv5,rv6,rv7,wv2,av2):

    #symbolic variables to be calculated
    r3, t3, t5, t6, delr3, delt3, delt5, delt6 = sym.symbols('r3, t3, t5, t6, delr3, delt3, delt5, delt6')
    v3, w3, w5, w6, ar3, at3, a5, a6 = sym.symbols('v3, w3, w5, w6, ar3, at3, a5, a6')

    #prepare to record data
    table = Table()
    plotData = PlotData()

    #Initial values
    rv3 = 1.0
    tv3 = 1.0
    tv5, tv6 = findValueT5T6(rv1,rv2,rv4,rv5,rv6,rv7)
    rv8 = rv1

    #find input angle limitations
    maxt3, mint3 = t3Limits(rv1, rv2, rv4, rv5, rv6, rv7)

    #find input angles
    inputAngles = inputIncrementAngles(num = 37)

    for step, tv2 in enumerate(inputAngles):
        print("tv2 =", 360/(37-1)*(step + 1))
        rv3, tv3 = NMsolveEquations(*loopOneEquations(rv1, rv2, r3, tv2, t3, delr3, delt3), r3, t3, rv3, tv3, delr3, delt3, normalize=[False, True])
        
        if invalidAngle(tv3, maxt3, mint3):
            print(f"This input angle {sym.N(tv2/sym.pi*180)} degree cannot occur!")
            table.addRow(step, False)
            continue

        #numerical methods
        tv5, tv6 = NMsolveEquations(*loopTwoEquations(rv4, rv5, rv6, rv7, rv8, tv3, t5, t6, delt5, delt6), t5, t6, tv5, tv6, delt5, delt6, normalize=[True, True])

        #symbolic methods
        vv3, wv3 = symSolveEquations(link3VelocityEquations(rv2, rv3, tv2, tv3, wv2, w3, v3), (v3, w3))
        wv5, wv6 = symSolveEquations(link56AngularVelocityEquations(rv4, rv5, rv6, tv3, tv5, tv6, wv3, w5, w6), (w5, w6))
        arv3, atv3 = symSolveEquations(link3AccelerationEquations(rv2, rv3, tv2, tv3, vv3, wv2, wv3, av2, ar3, at3), (ar3, at3))
        av5, av6 = symSolveEquations(link56AngularAccelerationEquations(rv4, rv5, rv6, tv3, tv5, tv6, wv3, wv5, wv6, a5, a6, atv3), (a5, a6))
    #     print(f"Ve:({-rv4 * sym.sin(tv3) * wv3 - rv5 * sym.sin(tv5) * wv5},{rv4 * sym.cos(tv3) * wv3 + rv5 * sym.cos(tv5) * wv5})")
    #     print(f"Ae:({-rv4 * sym.cos(tv3) * wv3 ** 2 - rv4 * sym.sin(tv3) * atv3 - rv5 * sym.cos(tv5) * wv5 ** 2 - rv5 * sym.sin(tv5) * av5},{-rv4 * sym.sin(tv3) * wv3 ** 2 + rv4 * sym.cos(tv3) * atv3 - rv5 * sym.sin(tv5) * wv5 ** 2 - rv5 * sym.cos(tv5) * av5})")
        dx, dy, ex, ey, evx, evy, eax, eay = plotData.addData(step, rv1, rv3, rv4, rv5, tv3, tv5, tv6, wv3, wv5, av5, atv3)
        table.addRow(step, True, dx, dy, ex, ey, evx, evy, eax, eay, wv6)
    plotData.cleanInvalidValues()
    return plotData, table

#1 2 4 1.5 1.5 4.1 1 1
##1 2 4 3.5 3.5 4.1 1 1

if __name__ == "__main__":
    params = getParameters()
    plotData, table = kinematic(**params)
    # plotData.display(int(input("enter 1 if you want the angle position relation plot; enter 2 if you want the position of point D & E; enter 3 if you want velocity relation plot; else if you want acceleration plot: ")))
    plotData.display(1)
    plotData.display(2)
    plotData.display(3)
    plotData.display(4)
    table.display()