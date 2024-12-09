import numpy as np 
import math as mt 
from LU_mi_version import sustitucion_progresiva
from LU_mi_version import sustitucion_regresiva
def isRegular(a):
    n = len(a)
    for i in range(n):
        a_i = a[:i,:i]
        if np.linalg.det(a_i)==0:
            return False
    return True 

def isSymetric(a):
    n = len(a)
    for i in range(n):
        for j in  range(n):
            if(a[i,j]!=a[j,i]):
                return False
    return True

def ld_iniciales():
    for i in range(n):
        for j in range(n):
            if(j> i ):
                l[i,j]=0
            if(i==j):
                l[i,j]=1
    for i in range(n):
        for j in range(n):
            if(i !=j ):
                d[i,j]=0.0
            else:
                d[i,j]=1.0          

def ldlt(a,l,d):
    for k in range(n):
        suma = 0
        for p in range(k):
            suma = suma +(l[k,p]**2)*d[p,p]
        d[k,k] = a[k,k] - suma
        if(d[k,k] ==0):
            break
        for i in range(k+1,n):
            suma = 0
            for p in range(k):
                suma = suma + l[i,p]*l[k,p]*d[p,p]
            l[i,k] = (a[i,k] - suma)/d[k,k]

if __name__=='__main__':
    a = np.array([[12.0,1.0,4.0,4.0],[6.0,10.0,15.0,18.0],[4.0,5.0,8.0,7.0],[4.0,5.0,7.0,1.0]])
    b = np.array([[0.0],[20.0],[9.0],[50.0]])
    n = len(a)
    if(not isSymetric(a)):
        b = np.transpose(a)@b 
        a = np.transpose(a)@a
    l = np.zeros((n,n),float)
    d = np.zeros((n,n),float)
    ld_iniciales()
    ldlt(a,l,d)
    print(l)
    print(d)
    z = sustitucion_progresiva(l,b)    #lz = b  l es inferior-progresiva
    y = np.zeros((n,1),float)
    for i in range(n):                 #Dy = z   D es diagonal
        y[i,0]=z[i,0]/d[i,i]
    x = sustitucion_regresiva(np.transpose(l),y)#l^t x = y es superior-regresiva
    print(x)
