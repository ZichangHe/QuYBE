import sys
import numpy as np
import numpy.linalg as npla
from math import *
    
I = np.eye(2)
CX = np.eye(4)
CX[2,2] = CX[3,3] = 0
CX[2,3] = CX[3,2] = 1

def Rx(t):
    m = np.zeros((2,2),dtype=complex)
    m[0,0] = m[1,1] = np.cos(t/2)
    m[1,0] = m[0,1] = -1j*np.sin(t/2)
    return m

def Rz(t):
    m = np.zeros((2,2),dtype=complex)
    m[0,0] = np.exp(-1j*t/2)
    m[1,1] = np.exp(1j*t/2)
    return m

def G(t1,t2,opt):
    if opt == 1:
        G = np.kron(Rx(pi/2),Rx(pi/2))@CX@np.kron(Rx(t1),Rz(t2))@CX@np.kron(Rx(-pi/2),Rx(-pi/2))
    elif opt == 2:
        G = CX@np.kron(Rx(t1),Rz(t2))@CX
    elif opt == 3:
        G = np.kron(Rz(pi/2),Rz(pi/2))@CX@np.kron(Rx(t1),Rz(t2))@CX@np.kron(Rz(-pi/2),Rz(-pi/2))
    return G

def YBE_V2A(t_init):
    '''
    YBE transform for Heisenberg XY/YZ/XZ model (no external field)
    V-shape to A-shape
    input: six rotations
    output: transfomed rotations
    '''

    def find_rotation(s,c):
        '''
        based on the signs of sin(r) and cos(r), find the right rotation 'r'.
        '''
        r = 0.0
        if s > 0:
            if c > 0:
                r = np.arctan(s/c)
            elif c < 0:
                r = np.pi + np.arctan(s/c)
            else:
                r = 0.5*np.pi
        elif s < 0:
            if c > 0:
                r = np.arctan(s/c)
            elif c < 0:
                r = np.pi + np.arctan(s/c)
            else:
                r = 1.5*np.pi
        else:
            if c > 0:
                r = 0.0
            elif c < 0:
                r = np.pi
            else:
                # cos(t10) = 0
                print('xxx')
                r = np.random.rand()
        return r

    t_trans = np.zeros(len(t_init))
    # L.H.S. of Eqs (1-16) in analysis.tex
    Eq = np.zeros(16)
    Eq[0]  =  np.sin(0.5*t_init[2])*np.cos(0.5*(t_init[0]-t_init[4]))*np.cos(0.5*(t_init[1]-t_init[5]))*np.sin(0.5*t_init[3])
    Eq[1]  =  np.cos(0.5*t_init[2])*np.cos(0.5*(t_init[0]-t_init[4]))*np.cos(0.5*(t_init[1]+t_init[5]))*np.sin(0.5*t_init[3])
    Eq[2]  = -np.sin(0.5*t_init[2])*np.cos(0.5*(t_init[0]+t_init[4]))*np.sin(0.5*(t_init[1]-t_init[5]))*np.cos(0.5*t_init[3])
    Eq[3]  =  np.cos(0.5*t_init[2])*np.cos(0.5*(t_init[0]+t_init[4]))*np.sin(0.5*(t_init[1]+t_init[5]))*np.cos(0.5*t_init[3])
    Eq[4]  =  np.sin(0.5*t_init[2])*np.cos(0.5*(t_init[0]+t_init[4]))*np.cos(0.5*(t_init[1]-t_init[5]))*np.cos(0.5*t_init[3])
    Eq[5]  =  np.cos(0.5*t_init[2])*np.cos(0.5*(t_init[0]+t_init[4]))*np.cos(0.5*(t_init[1]+t_init[5]))*np.cos(0.5*t_init[3])
    Eq[6]  = -np.sin(0.5*t_init[2])*np.cos(0.5*(t_init[0]-t_init[4]))*np.sin(0.5*(t_init[1]-t_init[5]))*np.sin(0.5*t_init[3])
    Eq[7]  =  np.cos(0.5*t_init[2])*np.cos(0.5*(t_init[0]-t_init[4]))*np.sin(0.5*(t_init[1]+t_init[5]))*np.sin(0.5*t_init[3])
    Eq[8]  =  np.sin(0.5*t_init[2])*np.sin(0.5*(t_init[0]+t_init[4]))*np.cos(0.5*(t_init[1]-t_init[5]))*np.cos(0.5*t_init[3])
    Eq[9]  =  np.cos(0.5*t_init[2])*np.sin(0.5*(t_init[0]+t_init[4]))*np.cos(0.5*(t_init[1]+t_init[5]))*np.cos(0.5*t_init[3])
    Eq[10] =  np.sin(0.5*t_init[2])*np.sin(0.5*(t_init[0]-t_init[4]))*np.sin(0.5*(t_init[1]-t_init[5]))*np.sin(0.5*t_init[3])
    Eq[11] = -np.cos(0.5*t_init[2])*np.sin(0.5*(t_init[0]-t_init[4]))*np.sin(0.5*(t_init[1]+t_init[5]))*np.sin(0.5*t_init[3])
    Eq[12] = -np.sin(0.5*t_init[2])*np.sin(0.5*(t_init[0]-t_init[4]))*np.cos(0.5*(t_init[1]-t_init[5]))*np.sin(0.5*t_init[3])
    Eq[13] = -np.cos(0.5*t_init[2])*np.sin(0.5*(t_init[0]-t_init[4]))*np.cos(0.5*(t_init[1]+t_init[5]))*np.sin(0.5*t_init[3])
    Eq[14] = -np.sin(0.5*t_init[2])*np.sin(0.5*(t_init[0]+t_init[4]))*np.sin(0.5*(t_init[1]-t_init[5]))*np.cos(0.5*t_init[3])
    Eq[15] =  np.cos(0.5*t_init[2])*np.sin(0.5*(t_init[0]+t_init[4]))*np.sin(0.5*(t_init[1]+t_init[5]))*np.cos(0.5*t_init[3])

    # an auxiliary matrix with elements being the L.H.S. of the 16 Eqs. see Eq (45) in analysis.tex
    # the matrix is only used for demonstration of sign assignments
    Mat = np.zeros((4,4))
    Mat[0,0] = Eq[0] 
    Mat[0,1] = Eq[1] 
    Mat[1,2] = Eq[2] 
    Mat[1,3] = Eq[3] 
    Mat[1,0] = Eq[4] 
    Mat[1,1] = Eq[5] 
    Mat[0,2] = Eq[6] 
    Mat[0,3] = Eq[7] 
    Mat[3,0] = Eq[8] 
    Mat[3,1] = Eq[9] 
    Mat[2,2] = Eq[10]
    Mat[2,3] = Eq[11]
    Mat[2,0] = Eq[12]
    Mat[2,1] = Eq[13]
    Mat[3,2] = Eq[14]
    Mat[3,3] = Eq[15]

    # a,b are constant arrays correspond to Eqs. (17-32) in analysis.tex
    a = np.zeros(8)
    b = np.zeros(8)

    m = 0
    for i in range(0,16,2):
        a[m] = Eq[i]*Eq[i] + Eq[i+1]*Eq[i+1]
        m += 1

    m = 0
    for i in range(0,4):
        b[m] = Eq[i]*Eq[i] + Eq[i+4]*Eq[i+4]
        m += 1

    for i in range(8,12):
        b[m] = Eq[i]*Eq[i] + Eq[i+4]*Eq[i+4]
        m += 1

    # follow Eqs. (33) and (35) in analysis.tex
    t_trans[2] = 2*np.arcsin(np.sqrt(np.sum(a[4:8])))
    t_trans[3] = 2*np.arcsin(np.sqrt(np.sum(a[1:8:2])))

    # col and row vectors correspond to the L.H.S. vectors of Eq. (45) in analysis.tex
    col = np.zeros(4)
    row = np.zeros(4)

    # follow Eqs. (37-44) to fill col and row vectors
    row[0] = np.sign(Mat[0,0])*np.sqrt(b[0]+b[4])
    row[1] = np.sign(Mat[0,1])*np.sqrt(b[1]+b[5])
    row[2] = np.sign(Mat[0,2])*np.sqrt(b[2]+b[6])
    row[3] = np.sign(Mat[0,3])*np.sqrt(b[3]+b[7])
    col[0] = np.sign(Mat[0,0])*np.sqrt(a[0]+a[3])
    col[1] = np.sign(Mat[1,0])*np.sqrt(a[1]+a[2])
    col[2] = np.sign(Mat[2,0])*np.sqrt(a[5]+a[6])
    col[3] = np.sign(Mat[3,0])*np.sqrt(a[4]+a[7])

    Mtmp = np.outer(col,row)
    if np.allclose(Mat,-1*Mtmp):
        row = -1*row

    # find t_trans[0] +/- t_trans[5] and t_trans[1] +/- t_trans[6]
    # take Eq. (61) in analysis.tex as reference, t7 <-> t_trans[0], t8 <-> t_trans[1], ..., t12 <-> t_trans[5]
    t0p4 = find_rotation(row[0],row[1])
    t0m4 = find_rotation(row[2],row[3])
    t1p5 = find_rotation(col[0],col[1])
    t1m5 = find_rotation(col[2],col[3])

    t_trans[0] = (t0p4 + t0m4)
    t_trans[4] = (t0p4 - t0m4)
    t_trans[1] = (t1p5 + t1m5)
    t_trans[5] = (t1p5 - t1m5)

    return t_trans

def YBE_A2V(t_init):
    '''
    YBE transform for Heisenberg XY/YZ/XZ model (no external field) 
    A-shape to V-shape
    input: six rotations
    output: transfomed rotations
    '''

    def find_rotation(s,c):
        '''
        based on the signs of sin(r) and cos(r), find the right rotation 'r'.
        '''
        r = 0.0
        if s > 0:
            if c > 0:
                r = np.arctan(s/c)
            elif c < 0:
                r = np.pi + np.arctan(s/c)
            else:
                r = 0.5*np.pi
        elif s < 0:
            if c > 0:
                r = np.arctan(s/c)
            elif c < 0:
                r = np.pi + np.arctan(s/c)
            else:
                r = 1.5*np.pi
        else:
            if c > 0:
                r = 0.0
            elif c < 0:
                r = np.pi
            else:
                # cos(t10) = 0
                r = np.random.rand()
        return r

    t_trans = np.zeros(len(t_init))
    # R.H.S. of Eqs (1-16) in analysis.tex
    Eq = np.zeros(16)
    Eq[0]  =  np.cos(0.5*t_init[2])*np.sin(0.5*(t_init[0]+t_init[4]))*np.sin(0.5*(t_init[1]+t_init[5]))*np.cos(0.5*t_init[3])
    Eq[1]  =  np.cos(0.5*t_init[2])*np.cos(0.5*(t_init[0]+t_init[4]))*np.sin(0.5*(t_init[1]+t_init[5]))*np.cos(0.5*t_init[3])
    Eq[2]  =  np.cos(0.5*t_init[2])*np.sin(0.5*(t_init[0]-t_init[4]))*np.cos(0.5*(t_init[1]+t_init[5]))*np.sin(0.5*t_init[3])
    Eq[3]  =  np.cos(0.5*t_init[2])*np.cos(0.5*(t_init[0]-t_init[4]))*np.cos(0.5*(t_init[1]+t_init[5]))*np.sin(0.5*t_init[3])
    Eq[4]  =  np.cos(0.5*t_init[2])*np.sin(0.5*(t_init[0]+t_init[4]))*np.cos(0.5*(t_init[1]+t_init[5]))*np.cos(0.5*t_init[3])
    Eq[5]  =  np.cos(0.5*t_init[2])*np.cos(0.5*(t_init[0]+t_init[4]))*np.cos(0.5*(t_init[1]+t_init[5]))*np.cos(0.5*t_init[3])
    Eq[6]  =  np.cos(0.5*t_init[2])*np.sin(0.5*(t_init[0]-t_init[4]))*np.sin(0.5*(t_init[1]+t_init[5]))*np.sin(0.5*t_init[3])
    Eq[7]  =  np.cos(0.5*t_init[2])*np.cos(0.5*(t_init[0]-t_init[4]))*np.sin(0.5*(t_init[1]+t_init[5]))*np.sin(0.5*t_init[3])
    Eq[8]  =  np.sin(0.5*t_init[2])*np.sin(0.5*(t_init[0]+t_init[4]))*np.cos(0.5*(t_init[1]-t_init[5]))*np.cos(0.5*t_init[3])
    Eq[9]  =  np.sin(0.5*t_init[2])*np.cos(0.5*(t_init[0]+t_init[4]))*np.cos(0.5*(t_init[1]-t_init[5]))*np.cos(0.5*t_init[3])
    Eq[10] =  np.sin(0.5*t_init[2])*np.sin(0.5*(t_init[0]-t_init[4]))*np.sin(0.5*(t_init[1]-t_init[5]))*np.sin(0.5*t_init[3])
    Eq[11] =  np.sin(0.5*t_init[2])*np.cos(0.5*(t_init[0]-t_init[4]))*np.sin(0.5*(t_init[1]-t_init[5]))*np.sin(0.5*t_init[3])
    Eq[12] =  np.sin(0.5*t_init[2])*np.sin(0.5*(t_init[0]+t_init[4]))*np.sin(0.5*(t_init[1]-t_init[5]))*np.cos(0.5*t_init[3])
    Eq[13] =  np.sin(0.5*t_init[2])*np.cos(0.5*(t_init[0]+t_init[4]))*np.sin(0.5*(t_init[1]-t_init[5]))*np.cos(0.5*t_init[3])
    Eq[14] =  np.sin(0.5*t_init[2])*np.sin(0.5*(t_init[0]-t_init[4]))*np.cos(0.5*(t_init[1]-t_init[5]))*np.sin(0.5*t_init[3])
    Eq[15] =  np.sin(0.5*t_init[2])*np.cos(0.5*(t_init[0]-t_init[4]))*np.cos(0.5*(t_init[1]-t_init[5]))*np.sin(0.5*t_init[3])

    # an auxiliary matrix with elements being the L.H.S. of the 16 Eqs. see Eq (74) in analysis.tex
    # the matrix is only used for demonstration of sign assignments
    Mat = np.zeros((4,4))
    Mat[0,0] =  Eq[0] 
    Mat[0,1] =  Eq[1] 
    Mat[1,2] = -Eq[2] 
    Mat[1,3] =  Eq[3] 
    Mat[1,0] =  Eq[4] 
    Mat[1,1] =  Eq[5] 
    Mat[0,2] = -Eq[6] 
    Mat[0,3] =  Eq[7] 
    Mat[3,0] =  Eq[8] 
    Mat[3,1] =  Eq[9] 
    Mat[2,2] =  Eq[10]
    Mat[2,3] = -Eq[11]
    Mat[2,0] = -Eq[12]
    Mat[2,1] = -Eq[13]
    Mat[3,2] = -Eq[14]
    Mat[3,3] =  Eq[15]

    # sys.exit(1)
    # e,f are constant arrays correspond to Eqs. (46-61) in analysis.tex
    e = np.zeros(8)
    f = np.zeros(8)

    e[0] = Eq[0]*Eq[0] + Eq[6]*Eq[6]
    e[1] = Eq[1]*Eq[1] + Eq[7]*Eq[7]
    e[2] = Eq[2]*Eq[2] + Eq[4]*Eq[4]
    e[3] = Eq[3]*Eq[3] + Eq[5]*Eq[5]
    e[4] = Eq[8]*Eq[8] + Eq[14]*Eq[14]
    e[5] = Eq[9]*Eq[9] + Eq[15]*Eq[15]
    e[6] = Eq[10]*Eq[10] + Eq[12]*Eq[12]
    e[7] = Eq[11]*Eq[11] + Eq[13]*Eq[13]
    
    m = 0
    for i in range(0,4):
        f[m] = Eq[i]*Eq[i] + Eq[i+12]*Eq[i+12]
        m += 1

    for i in range(4,8):
        f[m] = Eq[i]*Eq[i] + Eq[i+4]*Eq[i+4]
        m += 1

    # follow Eqs. (62) and (64) in analysis.tex
    t_trans[2] = 2*np.arcsin(np.sqrt(e[0]+e[2]+e[4]+e[6]))
    t_trans[3] = 2*np.arcsin(np.sqrt(e[0]+e[1]+e[6]+e[7]))

    # col and row vectors correspond to the L.H.S. vectors of Eq. (74) in analysis.tex
    col = np.zeros(4)
    row = np.zeros(4)

    # follow Eqs. (66-73) to fill col and row vectors
    row[0] = np.sign(Mat[0,0])*np.sqrt(f[0]+f[4])
    row[1] = np.sign(Mat[0,1])*np.sqrt(f[1]+f[5])
    row[2] = np.sign(Mat[0,2])*np.sqrt(f[2]+f[6])
    row[3] = np.sign(Mat[0,3])*np.sqrt(f[3]+f[7])
    col[0] = np.sign(Mat[0,0])*np.sqrt(e[0]+e[1])
    col[1] = np.sign(Mat[1,0])*np.sqrt(e[2]+e[3])
    col[2] = np.sign(Mat[2,0])*np.sqrt(e[6]+e[7])
    col[3] = np.sign(Mat[3,0])*np.sqrt(e[4]+e[5])
    
    Mtmp = np.outer(col,row)
    if np.allclose(Mat,-1*Mtmp):
        row = -1*row

    # find t_trans[0] +/- t_trans[5] and t_trans[1] +/- t_trans[6]
    # take Eq. (61) in analysis.tex as reference, t7 <-> t_trans[0], t8 <-> t_trans[1], ..., t12 <-> t_trans[5]
    t1p5 = find_rotation(row[3],row[1])
    t1m5 = find_rotation(row[2],row[0])
    t0p4 = find_rotation(col[3],col[1])
    t0m4 = find_rotation(col[2],col[0])

    # print(t1p5,t1m5,t0p4,t0m4)
    t_trans[0] = (t0p4 + t0m4)
    t_trans[4] = (t0p4 - t0m4)
    t_trans[1] = (t1p5 + t1m5)
    t_trans[5] = (t1p5 - t1m5)

    return t_trans

'''
verify 
--|theta_init[0]|-------------------|theta_init[4]|-- 
--|theta_init[1]|--|theta_init[2]|--|theta_init[5]|-- = 
-------------------|theta_init[3]|-------------------   

--------------------|theta_trans[2]|--------------------
--|theta_trans[0]|--|theta_trans[3]|--|theta_trans[4]|--
--|theta_trans[1]|--------------------|theta_trans[5]|--
'''

