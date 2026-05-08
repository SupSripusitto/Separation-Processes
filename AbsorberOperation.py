import numpy as np
import matplotlib.pyplot as plt

''' User input part'''

# Output display (1 for needed, 0 for not needed)
graph = 1
report = 1

# Solute mole ratio (or approximately mole fraction) in absorbent feed
X0 = 0

# Gas feed solute mole ratio
YN1 = 0.0067

# Gas outlet spec. (put negative value for recovery of solute in absorbent, no need to put here if you know the exact L/V ratio)
Y1 = -0.67

# Molar gas (without solute) flow rate
V = 1

# Molar liquid (without solute) flow rate (put negative value for multiples of the minimum value, if negative, you need to put Y1.)
L = -1.5

# Equilibrium constant/Partitioning coefficient (allow only for constant value)
K = 67

# Number of stages
N = 3

''' The end of user inputs, the next part is only the algorithms.'''

# Find operating line slope
slope = None

if Y1 < 0:
    Y1 = (1+Y1)*YN1

if L > 0:
    slope = L/V
elif Y1 > 0:
    slope = -L*(YN1-Y1)/(YN1/(K*(YN1+1)-YN1)-X0)

# Concentration profile from guessed Y1 function
def Conc_Profile(y1):
    XYonline = np.zeros((2*N+1,2))
    XYonline[0,:] = [X0,y1]

    for i in range(1,2*N+1,2):
        yn = XYonline[i-1,1]
        XYonline[i,:] = [yn/K,yn]

        xn = XYonline[i,0]
        XYonline[i+1,:] = [xn,y1+slope*(xn-X0)]
    return XYonline

# Quasi-Newton Solver to get feed composition the same as YN1
Tol = 1e-7
err = 1

if X0 > Y1/K:
    Tol = 1
    report = 0
    graph = 0
    print("!!! Too high absorbent feed composition")

while err > Tol:
    XY = Conc_Profile(Y1)
    YN1c = XY[-1,1]
    der = 1e7*(Conc_Profile(Y1+1e-7)[-1,1]-YN1c)
    Y1 -= (YN1c-YN1)/der
    err = np.abs(YN1-YN1c)/YN1

XN = (YN1-Y1)/slope+X0

if report:
    print("===== Calculation Report =====")
    print("Gas outlet mole ratio:       ",Y1)
    print("Liquid outlet mole ratio:    ",XN)
    print("Solute recovery in absorbent:",(YN1-Y1)/YN1)
    print("===== End of the report =====")

if graph:
    plt.plot(XY[:,0],XY[:,1])
    plt.plot([0,YN1*1.05/K],[0,YN1*1.05])
    plt.plot([X0,XN],[Y1,YN1])
    plt.legend(["Staircase","Equilibrium line","Operating line"])
    plt.show()