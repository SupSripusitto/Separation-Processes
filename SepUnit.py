import numpy as np
import matplotlib.pyplot as plt

def AbsorberOperation(X0 : float, YN1 : float, Y1 : float, V : float, L : float, K : float, N : int, report : bool = True, graph : bool = True):
    '''
    Calculate the operating conditions from spec. & design
    Arg:
        X0 : Inlet liquid mole ratio
        YN1 : Inlet gas mole ratio
        Y1 : Outlet gas spec. (put negative value for recovery fraction)
        V : Molar gas (w/o solute) flow rate
        L : Molar liquid (w/o solute) flow rate (put neg. for times of minimum liquid flow rate)
        K : Equilibrium constant / Partitioning coefficient
        N : Number of stages
        report (optional) : show report or not
        graph (optional) : show graph or not
    '''
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
    elif Y1 > YN1:
        Tol = 1
        report = 0
        graph = 0
        print("!!! Outlet mole ratio must less than the inlet.")

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

def AbsorberDesign(X0 : float, YN1 : float, Y1 : float, V : float, L : float, K : float, Nm : int = 150, report : bool = True, graph : bool = True):
    '''
    Calculate number of stages from spec.
    Arg :
        X0 : Inlet liquid mole ratio
        YN1 : Inlet gas mole ratio
        Y1 : Outlet gas spec. (put negative value for recovery fraction)
        V : Molar gas (w/o solute) flow rate
        L : Molar liquid (w/o solute) flow rate (put neg. for times of minimum liquid flow rate)
        K : Equilibrium constant / Partitioning coefficient
        Nm (optional) : Maximum number of stages (150 by default)
        report (optional) : show report or not
        graph (optional) : show graph or not
    '''
    slope = None

    if Y1 < 0:
        Y1 = (1+Y1)*YN1

    if L > 0:
        slope = L/V
    elif Y1 > 0:
        slope = -L*(YN1-Y1)/(YN1/(K*(YN1+1)-YN1)-X0)

    # Find the ending criterion
    XN = X0 + (YN1-Y1)/slope

    XYonline = np.zeros((2*Nm,2))

    XYonline[0,:] = [X0,Y1]

    i = 0

    if X0 > Y1/K:
        XYonline[0,:] = np.ones((1,2))
        print("!!! Too high absorbent feed composition")
        report = 0
        graph = 0
    elif Y1 > YN1:
        XYonline[0,:] = np.ones((1,2))
        report = 0
        graph = 0
        print("!!! Outlet mole ratio must less than the inlet.")

    while XYonline[i,0] < XN:
        Yn = XYonline[i,1]
        i += 1
        XYonline[i,:] = [Yn/K,Yn]
        
        Xn = XYonline[i,0]
        i += 1
        XYonline[i,:] = [Xn,Y1 + slope*(Xn-X0)]

        if i == 2*Nm-2:
            print("!!! The number of stages is over the limit. (Check your spec. or the parameter you put in)")
            report = 0
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