import numpy as np
from jaxtyping import Array,Float
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math as mt

def exponencial_decay(t:Float[Array,"dim"],u:Float[Array,"dim"]) -> Float[Array,"2"]:
    return -.5*u
#la solucion es a *e^-t/2
e = np.e
if __name__=='__main__':
    sol=solve_ivp(fun = exponencial_decay,t_span=(0,10),y0=(2,4,6,8),t_eval=np.linspace(0,10),dense_output=True,)
    print(sol)
    x = np.linspace(0,10)
    print(sol.y)
    coeficientes = sol.y[:,0]
    print(coeficientes)
    plt.figure()
    for i in coeficientes:
        Y = i*e**(-x/2)
        plt.plot(x,Y,label=f'u(0)={i}')
        
    plt.xlabel('Time t')
    plt.ylabel('u(t)')
    plt.show()
    