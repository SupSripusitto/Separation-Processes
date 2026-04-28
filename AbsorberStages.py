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

# Gas outlet spec. (put negative value for recovery of solute in absorbent)
Y1 = -0.67

# Molar gas (without solute) flow rate
V = 1

# Molar liquid (without solute) flow rate (put negative value for multiples of the minimum value)
L = -1.5

# Equilibrium constant/Partitioning coefficient (allow only for constant value)
K = 67

# Maximum number of stages (to prevent infinite loop calculation, no need to change the value from 150 stages, if it isn't necessary)
Nm = 150

''' The end of user inputs, the next part is only the algorithms.'''

# Find operating line slope
slope = None

if L > 0:
    slope = L/V
elif Y1 > 0:
    slope = -K*(YN1-Y1)/YN1*L
else:
    slope = K*Y1*L
    Y1 = (1+Y1)*YN1

# Find the ending criterion
XN = X0 + (YN1-Y1)/slope

XYonline = np.zeros((2*Nm,2))

XYonline[0,:] = [X0,Y1]

i = 0

while XYonline[i,0] < XN:
    Yn = XYonline[i,1]
    i += 1
    XYonline[i,:] = [Yn/K,Yn]
    
    Xn = XYonline[i,0]
    i += 1
    XYonline[i,:] = [Xn,Y1 + slope*(Xn-X0)]

    if i == 2*Nm-2:
        print("Too extreme condition, the number of stages is over the limit.")
        break

if report:
    print("===== Calculation Report =====")
    print("Number of stages:                 ", i//2)
    print("Outlet liquid mole raio:          ", XYonline[i,0])
    print("Maximum feed mole ratio capable:  ", XYonline[i,1])
    print("Liquid to Feed ratio:             ", slope)
    print("===== End of the report =====")
if graph:
    plt.plot(XYonline[:i+1,0],XYonline[:i+1,1])
    plt.plot([0,1.05*Xn],[0,1.05*K*Xn])
    plt.plot([X0,XN],[Y1,YN1])
    plt.legend(["Staircase","Equilibrium line","Operating line"])
    plt.show()