import sys
from QuYBE.core_utils import *

if __name__ == "__main__":
    sys.argv.append('XY')
    if len(sys.argv) < 2:
        print("A pilot code to generate and verify YBE relation for Heisenberg XY/XZ/YZ")
        print("usage: python3 YBE_updated.py [ XY | XZ | YZ ]")
        sys.exit(1)
    
    ### input rotations
    theta_init = np.zeros(6)
    for i in range(6):
        theta_init[i] = np.random.rand()
    
    print('Inital rotations: ',theta_init)
    
    ## verify for H = (1) Hx + Hy, (2) Hx + Hz, or (3) Hy + Hz
    
    Hmodel = sys.argv[1]
    if Hmodel == 'XY':
        opt = 1
    elif Hmodel == 'XZ':
        opt = 2
    elif Hmodel == 'YZ':
        opt = 3
    else:
        print("The code only supports obtaining YBE relations for Heisenbergy XY/XZ/YZ models")
        sys.exit(1)
    
    A = np.kron(G(theta_init[0],theta_init[1],opt),I)
    B = np.kron(I,G(theta_init[2],theta_init[3],opt))
    C = np.kron(G(theta_init[4],theta_init[5],opt),I)
    Left = A@B@C
    
    theta_trans = YBE_V2A(theta_init)
    print('V2A YBE transformed rotations: ',theta_trans)
    
    A = np.kron(I,G(theta_trans[0],theta_trans[1],opt))
    B = np.kron(G(theta_trans[2],theta_trans[3],opt),I)
    C = np.kron(I,G(theta_trans[4],theta_trans[5],opt))
    Right = A@B@C
    
    print('Verifying V2A YBE for ',Hmodel,' model, Left => Right?',np.allclose(Left,Right))
    
    theta_trans2 = YBE_A2V(theta_trans)
    print('A2V YBE transformed rotations: ',theta_trans2)
    
    A = np.kron(G(theta_trans2[0],theta_trans2[1],opt),I)
    B = np.kron(I,G(theta_trans2[2],theta_trans2[3],opt))
    C = np.kron(G(theta_trans2[4],theta_trans2[5],opt),I)
    Left = A@B@C
    
    print('Verifying A2V YBE for ',Hmodel,' model, Right => Left? ',np.allclose(Left,Right))