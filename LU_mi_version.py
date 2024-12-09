import numpy as np
import math as mt

def lu_iniciales(l,u,n):
    for i in range(n):
        for j in range(n):
            if(i==j):
                l[i,j]=1
    for i in range(n):
        for j in range(n):
            if(i>j):
                u[i,j] = 0

def doolitle(a,l,u):
    n = len(a)
    for i in range(n):
        for j in range(i):
            suma = 0
            for k in range(j):
                suma = suma + l[i,k]*u[k,j] 
            l[i,j] = (a[i,j]-suma)/u[j,j] 
        for j in range(i,n):
            suma = 0
            for k in range(i):
                suma = suma +l[i,k]*u[k,j] 
            u[i,j] = (a[i,j] -suma)/l[i,i]

def sustitucion_progresiva(l,b):
    n = len(l)
    y = np.zeros((n,1),float)
    for i in range(n):
        suma = 0
        for j in range(0,i):
            suma = suma + l[i,j]*y[j,0]
        y[i,0] = (b[i]-suma)/l[i,i]
    return y

def sustitucion_regresiva(u,y):
    n = len(u)
    x = np.zeros((n,1),float)
    for i in range(n-1,-1,-1):
        suma = 0 
        for j in range(i+1,n):
            suma = suma + u[i,j]*x[j,0]
        x[i] = (y[i,0]-suma)/u[i,i]
    return x

if __name__ =='__main__':
    A = np.array([[12.0,1.0,4.0,4.0],[6.0,10.0,15.0,18.0],[4.0,5.0,8.0,7.0],[4.0,5.0,7.0,1.0]])
    b = np.array([[0.0],[20.0],[9.0],[50.0]])
    n = len(A)
    l=np.zeros((n,n),float)
    u=np.zeros((n,n),float)
    lu_iniciales(l,u,n)
    print(l)
    print(u)
    doolitle(A,l,u)
    print(l)
    print(u)
    #Ly=b progresiva 
    y = sustitucion_progresiva(l,b)
    print(y)
    print(sustitucion_regresiva(u,y))
#A=LU
#Ax=b
#LUx=b
# y
#Ly=b  L es inferior-sus progresiva
#Ux=y  U es superior-sus regresiva
#incluir restricciones luego 