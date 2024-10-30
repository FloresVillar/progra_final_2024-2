import matplotlib.pyplot as plt
import numpy as np
import math as mt
from mpl_toolkits.mplot3d import Axes3D

def grafica_simple():
    lista = np.array([-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
    cuadrado = lista**2
    plt.plot(lista,cuadrado,'o')
    plt.show()

def f1(x):
    return (mt.cos(x))**2

def f2(x):
    return mt.sin(x)*mt.cos(x)

def f3(x):
    try :
        mt.sin(x)/x 
    except ZeroDivisionError:
        return float('nan')
    return mt.sin(x)/x 

def grafica_de_f12():
    N = 100
    min = -2*mt.pi
    max = -min 
    delta = (max- min )/(N - 1)
    X = [0]*N
    Y1 =[0]*N
    Y2 =[0]*N
    for i in range(N):
        x = min + i*delta
        X[i] = x 
        Y1[i] = f1(x)
        Y2[i] = f2(x)
    plt.plot(X,Y1,"r-",'o',label="f1")
    plt.plot(X,Y2,"g-",'o',label="f2")
    plt.xlabel("$x$")
    plt.ylabel("f1 = cos(x)^2  f2= sin(x)*cos(x)")
    plt.title("f1 y f2")
    plt.legend()
    plt.show() 

def grafica_de_f3():
    N=301
    min = -15
    max = -min
    delta  = (max-min)/(N -1)
    X = [0]*N
    Y = [0]*N
    for i in range(N):
        x = min + i*delta
        X[i] = x
        Y[i] = f3(x)
    plt.plot(X,Y)
    plt.show()

def funcion_normal(x,u,sg):
    e = mt.e
    pi = mt.pi
    return (e**(-((x-u)**2)/(2*sg**2)))/(sg*(2*pi)**.5)

def grafica_de_funcion_normal():
    N = 1000
    min = -5
    max= -min
    X = np.linspace(min,max,N)
    u = 0
    sg = 1
    Y =[0]*N 
    for i,x in enumerate(X):
        Y[i] = funcion_normal(x,u,sg)
    plt.plot(X,Y,linestyle = '--',color="y")
    plt.title('Distribucion de Gauss')
    plt.xlabel("$x$")
    plt.ylabel('$N({},{})$'.format(u,sg))    
    plt.show()

def grafica_de_fxy():
    fig = plt.figure ()
    ax = fig.add_subplot(111,projection = "3d")
    x = np.linspace(-5,5,50)
    y = np.linspace(-5,5,50)
    X,Y = np.meshgrid(x,y)
    Z = X**2 + Y**2
    ax.plot_surface(X,Y,Z,cmap = "viridis")
    plt.show()

if __name__=='__main__':
    grafica_simple()
    #grafica_de_f12()
    #grafica_de_f3()
    #grafica_de_funcion_normal()
    #grafica_de_fxy()
#https://jorgedelossantos.github.io/apuntes-python/Matplotlib.html