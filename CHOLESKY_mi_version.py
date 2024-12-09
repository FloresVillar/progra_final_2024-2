import numpy as np     
import math as mat 
from LDL_mi_version import isSymetric 
from LU_mi_version import sustitucion_regresiva
from LU_mi_version import sustitucion_progresiva

def cholesky(a,g):
    n = len(a)
    for j in range(n):
        suma = 0
        for k in range(j):
            suma = suma + g[j,k]**2
        g[j,j] = mat.pow(a[j,j]-suma,0.5)
        for i in range(j+1,n):
            suma = 0
            for k in range(j):
                suma = suma + g[i,k]*g[j,k]
            g[i,j] = (a[i,j]-suma)/g[j,j]

def l_inicial(g):
    n = len(g)
    for j in range(n):
        for i in range(n):
            if(i==j):
                g[i,j]=1.0

if __name__=='__main__':
    a = np.array([[12.0,1.0,4.0,4.0],[6.0,10.0,15.0,18.0],[4.0,5.0,8.0,7.0],[4.0,5.0,7.0,1.0]])
    b = np.array([[0.0],[20.0],[9.0],[50.0]])
    n = len(a)
    if(not isSymetric(a)):
        b = np.transpose(a)@b
        a = np.transpose(a)@a
    g = np.zeros((n,n),float)
    l_inicial(g)
    cholesky(a,g)
    z = sustitucion_progresiva(g,b) #gz= b' es inferior- progresiva
    x = np.zeros((n,1),float)       #g^t x = z   superior - regresiva    
    x = sustitucion_regresiva(np.transpose(g),z)
    print(x)
"""
    si n = 12    NB= 4      k = 3  j : 1 → 3       
    j=1
    a11      a12     a13     a14        a15     a16     a17     a18          a19      a110       a111       a112
    a21      a22     a23     a24        a25     a26     a27     a28          a29      a210       a211       a212
A11 a31      a32     a33     a34        a35     a36     a37     a38          a39      a310       a311       a312
    a41      a42     a43     a44        a45     a46     a47     a48          a49      a410       a411       a412
    
    a51      a52     a53     a54        a55     a56     a57     a58          a59      a510       a511       a512
A21 a61      a62     a63     a64        a65     a66     a67     a68          a69      a610       a611       a612
    a71      a72     a73     a74        a75     a76     a77     a78          a79      a710       a711       a712
    a81      a82     a83     a84        a85     a86     a87     a88          a89      a810       a811       a812
                                                    A22
    a91      a92     a93     a94        a95     a96     a97     a98          a99      a910       a911       a912
A31 a101     a102    a103    a104       a105    a106    a107    a108         a109    a1010       a1011      a1012    
    a111     a112    a113    a114       a115    a116    a117    a118         a119    a1110       a1112      a1112
    a121     a122    a123    a124       a125    a126    a127    a128         a129    a1210       a1211      a1212
                                                    A32                                             A33
  
    L11 = Cholesky(A11)   luego actualizar  para A21 : L21 = A21 L11^-t 
                                            para A31 : L31 = A31 L11^-t
                        luego actualizar A22 y A33
                                            A22 = A22 - L21 L21^t
                                            A33 = A33 - L31 L21^t


                                                j=2
    a11      a12     a13     a14        a15     a16     a17     a18          a19      a110       a111       a112
    a21      a22     a23     a24        a25     a26     a27     a28          a29      a210       a211       a212
A11 a31      a32     a33     a34        a35     a36     a37     a38          a39      a310       a311       a312
    a41      a42     a43     a44        a45     a46     a47     a48          a49      a410       a411       a412
    
    a51      a52     a53     a54        a55     a56     a57     a58          a59      a510       a511       a512
A21 a61      a62     a63     a64        a65     a66     a67     a68          a69      a610       a611       a612
    a71      a72     a73     a74        a75     a76     a77     a78          a79      a710       a711       a712
    a81      a82     a83     a84        a85     a86     a87     a88          a89      a810       a811       a812
                                                    A22
    a91      a92     a93     a94        a95     a96     a97     a98          a99      a910       a911       a912
A31 a101     a102    a103    a104       a105    a106    a107    a108         a109    a1010       a1011      a1012    
    a111     a112    a113    a114       a115    a116    a117    a118         a119    a1110       a1112      a1112
    a121     a122    a123    a124       a125    a126    a127    a128         a129    a1210       a1211      a1212
                                                    A32                                             A33
  
    L22 = Cholesky(A22)   luego actualizar  para A32 : L32 = A32 L22^-t 
                        luego actualizar A33
                                            A33 = A33 - L32 L32^t


                                                                                        j=3
    a11      a12     a13     a14        a15     a16     a17     a18          a19      a110       a111       a112
    a21      a22     a23     a24        a25     a26     a27     a28          a29      a210       a211       a212
A11 a31      a32     a33     a34        a35     a36     a37     a38          a39      a310       a311       a312
    a41      a42     a43     a44        a45     a46     a47     a48          a49      a410       a411       a412
    
    a51      a52     a53     a54        a55     a56     a57     a58          a59      a510       a511       a512
A21 a61      a62     a63     a64        a65     a66     a67     a68          a69      a610       a611       a612
    a71      a72     a73     a74        a75     a76     a77     a78          a79      a710       a711       a712
    a81      a82     a83     a84        a85     a86     a87     a88          a89      a810       a811       a812
                                                    A22
    a91      a92     a93     a94        a95     a96     a97     a98          a99      a910       a911       a912
A31 a101     a102    a103    a104       a105    a106    a107    a108         a109    a1010       a1011      a1012    
    a111     a112    a113    a114       a115    a116    a117    a118         a119    a1110       a1112      a1112
    a121     a122    a123    a124       a125    a126    a127    a128         a129    a1210       a1211      a1212
                                                    A32                                             A33
  
    L33 = Cholesky(A33)    





    a11 a12 a13 a14     a15 a16 a17 a18      a19 a110 a111 a112
    a21 a22 a23 a24     a25 a26 a27 a28      a29 a210 a211 a212
A11 a31 a32 a33 a34     a35 a36 a37 a38      a39 a310 a311 a312
    a41 a42 a43 a44     a45 a46 a47 a48      a49 a410 a411 a412
    
    a51 a52 a53 a54     a55 a56 a57 a58      a59 a510 a511 a512
    a61 a62 a63 a64     a65 a66 a67 a68      a69 a610 a611 6a12
A21 a71 a72 a73 a74     a75 a76 a77 a78      a79 a710 a711 a712
    a81 a82 a83 a84     a85 a86 a87 a88      a89 a810 a811 a812
                            A22
    a91  a92  a93  a94   a95  a96  a97  a98   a99  a910  a911  a912
    a101 a102 a103 a104  a105 a106 a107 a108  a109 a1010 a1011 a1012
A31 a111 a112 a113 a114  a115 a116 a117 a118  a119 a1110 a1111 a1112
    a121 a122 a123 a124  a125 a126 a127 a128  a129 a1210 a1211 a1212
                            A32                         A33
"""                     