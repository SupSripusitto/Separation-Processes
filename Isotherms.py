import numpy as np
import matplotlib.pyplot as plt

'''Use input section'''

# Partial pressure & adsorbed amount data
p = np.array([5e-4,1e-3,2e-3,5e-3,1e-2,2e-2])
q = np.array([6.7e-5,11.2e-5,18e-5,33e-5,51e-5,78e-5])

# Isotherm name (has built-in "Langmuir", "Freundlich" or put "Custom" for the custom one)
isot = "Freundlich"

# Number of parameters (2 for Langmuir & Freundlich)
n_param = 2

# (Optional) Custom isotherm
def Custom(p,param):
    return param[0]*p

'''The end of the user input, the next part is only the algorithm'''

def Langmuir(p,param):
    return param[0,0]*param[1,0]*p/(1+param[0,0]*p)

def Freundlich(p,param):
    return param[0,0]*p**(1/param[1,0])

# Mapping string to function by dictionary
Isot = {"Freundlich":Freundlich, "Langmuir":Langmuir, "Custom":Custom}

param = np.ones((n_param,1))
Tol = 1e-11
grad = np.ones((n_param,1))
mse = 0

# Loss function defining
def MSE(param):
    q_est = np.zeros((1,len(q)))
    for i in range(len(q)):
        q_est[0,i] = Isot[isot](p[i],param)

    return np.sum((q-q_est)**2)/len(q)

lr = 1000

# Gradient descent
while np.abs(np.max(grad)) > Tol:
    # Find gradient
    for i in range(n_param):
        d = np.zeros((n_param,1))
        d[i] += 1e-10
        grad[i,0] = 1e10*(MSE(param + d)-MSE(param))
    
    # Line search
    while MSE(param - lr*grad) > MSE(param):
        lr /= 1.5
    
    param -= lr*grad

pspan = np.linspace(0,q[-1]*1.2,100)

print("RMSE:      ",np.sqrt(MSE(param)))
print("Parameters:",param.T)

pspan = np.linspace(0,p[-1]*1.2)

plt.scatter(p,q)
plt.plot(pspan,[Isot[isot](i,param) for i in pspan])
plt.show()