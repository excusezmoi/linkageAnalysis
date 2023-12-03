import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from utils import *

def getFrameDuration():
    dt = float(input("Enter desired duration for each frame(in sec): "))
    return dt

def getIterate():
    iterate = bool(int(input("Do you want the animation to repeat?(1/0): ")))
    return iterate

def frameDurationToMs(dt):
    inter = dt * 1000 #display time in millisecond
    return inter

#Find the input angle to be calculated
def generateInputAngles(wv2, av2, dt):
    inputAngles = []
    nwv2 = wv2
    tv2 = 0
    while tv2 < sym.N(2*sym.pi):
        inputAngles.append(tv2)
        tv2 += nwv2*dt + 0.5 * av2 * dt**2
        nwv2 += av2*dt
    return inputAngles

def calculateAllData(rv1,rv2,rv4,rv5,rv6,rv7,wv2,av2,dt):

    r3, t3, t5, t6, delr3, delt3, delt5, delt6 = sym.symbols('r3, t3, t5, t6, delr3, delt3, delt5, delt6')

    #Initial values
    rv3 = 1.0
    tv3 = 1.0
    tv5, tv6 = findValueT5T6(rv1,rv2,rv4,rv5,rv6,rv7)
    rv8 = rv1

    maxt3, mint3 = t3Limits(rv1, rv2, rv4, rv5, rv6, rv7)
    inputAngles = generateInputAngles(wv2, av2, dt)

    xr3, yr3, xr5, yr5, xr6, yr6 = [np.zeros(len(inputAngles)) for _ in range(6)]

    for step, tv2 in enumerate(inputAngles):
        rv3, tv3 = solveEquations(*loopOneEquations(rv1, rv2, r3, tv2, t3, delr3, delt3), r3, t3, rv3, tv3, delr3, delt3, normalize=[False, True])
        if invalidAngle(tv3, maxt3, mint3):
            print(f"This input angle {sym.N(tv2/sym.pi*180)} degree cannot occur!")
            continue
        tv5, tv6 = solveEquations(*loopTwoEquations(rv4, rv5, rv6, rv7, rv8, tv3, t5, t6, delt5, delt6), t5, t6, tv5, tv6, delt5, delt6, normalize=[True, True])
        xr3[step], yr3[step], xr5[step], yr5[step], xr6[step], yr6[step] = calculatePositions(rv1, rv3, rv4, rv5, tv3, tv5)

    return xr3, yr3, xr5, yr5, xr6, yr6

def linkageAnimation(params,dt,iterate):

    rv1 = params['rv1']
    rv7 = params['rv7']

    xr2, yr2 = 0, 0
    xr4, yr4 = 0, -rv1
    xrf, yrf = -rv7, 0

    xr3, yr3, xr5, yr5, xr6, yr6 = list(map(lambda arr: arr[arr != 0], calculateAllData(**params, dt = dt)))

    def findLimitXY(rv1, rv7, xr5=xr5, yr5=yr5, xr6=xr6, yr6=yr6):

        xMax = float(max(*xr5, 0)) + 1
        xMin = float(min(*xr6, -rv7)) - 1

        yMax = float(max(*yr5, *yr6)) + 1
        yMin = -rv1 - 1

        return {"xlim":(xMin, xMax), "ylim":(yMin, yMax)}


    # set up the figure and subplot
    fig = plt.figure()
    ax = fig.add_subplot(
        111, aspect="equal", autoscale_on=False, **findLimitXY(rv1, rv7)
    )
    # add grid lines, title and take out the axis tick labels
    ax.grid(alpha=0.5)
    ax.set_title("Strange mechanism")
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    (line, ) = ax.plot(
        [], [], "o-", lw=5, color="#2b8cbe"
    )  # color from: http://colorbrewer2.org/

    def init():
        line.set_data([], [])
        return line,

    # animation function
    def animate(step):
        x_points = [xr3[step], xr2, xr4, xr5[step], xr6[step], xrf]
        y_points = [yr3[step], yr2, yr4, yr5[step], yr6[step], yrf]

        line.set_data(x_points, y_points)
        return line,

    # call the animation
    _ani = animation.FuncAnimation(
        fig, animate, init_func=init, frames=len(xr3), interval=frameDurationToMs(dt), blit=True, repeat=iterate
    )

    # show the animation
    plt.show()
    ##1 2 4 3.5 3.5 4.1 1 1
    ## 1 2 4 3.5 3.5 4.1 3.1415926 0

if __name__ == "__main__":

    params = getParameters()
    linkageAnimation(params, dt = getFrameDuration(),iterate = getIterate())