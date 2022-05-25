import sys
import numpy as np
from math import *
from YBE_updated import *

class block:
    def __init__(self, theta=None, opt=1, row='12'):
        if theta is None:
            theta = np.random.rand(2)
            self.theta = theta
        else:
            assert len(theta) == 2
            self.theta = theta
        if row == '12':
            self.m = np.kron(G(self.theta[0],self.theta[1],opt),I)
        elif row == '23':
            self.m = np.kron(I,G(self.theta[0],self.theta[1],opt))
        self.row=row
 
    def merge(self, block2):
        assert self.row == block2.row
        self.theta += block2.theta
        if self.row == '12':
            self.m = np.kron(G(self.theta[0],self.theta[1],opt),I)
        elif self.row == '23':
            self.m = np.kron(I,G(self.theta[0],self.theta[1],opt))
    
sys.argv.append('XY')
if len(sys.argv) < 2:
    print("A pilot code to generate and verify YBE relation for Heisenberg XY/XZ/YZ")
    print("usage: python3 YBE_updated.py [ XY | XZ | YZ ]")
    sys.exit(1)
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
    
### Generate original chain
p = 4
row12 = []
row23 = []
for _ in range(p):
    row12.append(block(row='12'))
    row23.append(block(row='23'))
ori = np.kron(np.kron(I,I),I)
for i in range(p):
    ori = ori@row12[i].m@row23[i].m
 
### Apply YBE transformation
new_row12 = row12[:]
new_row23 = row23[:]
YBE_count = 0
while(len(new_row12)>1):
    # end up with a single A shape circuit
    if len(new_row12) >= len(new_row23):
        # when equal size, by default start from V shape
        # do the V2A transform
        if len(new_row12) >= 2 and len(new_row23) >= 1: 
            theta_init=np.concatenate((np.concatenate((new_row12[0].theta,new_row23[0].theta),axis=None), \
                                      new_row12[1].theta), axis=None)
            theta_trans = YBE_V2A(theta_init)
            YBE_count+=1
            del new_row12[:2] #more efficient than pop
            new_row12 = [block(theta=theta_trans[2:4],row='12')] + new_row12
            del new_row23[:1]
            new_row23 = [block(theta=theta_trans[0:2],row='23'), block(theta=theta_trans[4:6],row='23')] + new_row23

            # merge two blocks in row23
            if len(new_row23) >= 3:
                new_row23[1].merge(new_row23[2])
                del new_row23[2]
    else:
        # do the A2V transofrm
        if len(new_row12) >= 1 and len(new_row23) >= 2:
            theta_init=np.concatenate((np.concatenate((new_row23[0].theta,new_row12[0].theta),axis=None), \
                                      new_row23[1].theta), axis=None)
            theta_trans = YBE_A2V(theta_init)
            YBE_count+=1
            # new_row12.pop(0)
            del new_row12[:1] 
            new_row12 = [block(theta=theta_trans[0:2],row='12'), block(theta=theta_trans[4:6],row='12')] + new_row12
            del new_row23[:2] 
            new_row23 = [block(theta=theta_trans[2:4],row='23')] + new_row23
            
            # merge two blocks in row12
            if len(new_row12) >= 3:
                new_row12[1].merge(new_row12[2])
                del new_row12[2]
                
assert (len(new_row12)==1 and len(new_row23)==2)              
new = new_row23[0].m@new_row12[0].m@new_row23[1].m
print(f'Verify for 3-qubit {p}-depth chain, {YBE_count}# YBE, {np.allclose(ori,new)}')

