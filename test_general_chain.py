import sys
import numpy as np
from math import *
import copy
from YBE_updated import *

class block:
    def __init__(self, theta=None, opt=1, row=0, col=None, N=6):
        # N is the # of qubits
        assert (row < N-1 and row >= 0)
        if theta is None:
            theta = np.random.rand(2)
            self.theta = theta
        else:
            assert len(theta) == 2
            self.theta = theta
        self.col = col
        self.row = row
        self.ind = {'row':self.row, 'col':self.col}
        self.N=N
        temp1, temp2 = 1, 1
        for _ in range(row):
            temp1 = np.kron(I,temp1)
        for _ in range(row+1,N-1):
            temp2 = np.kron(I,temp2)
        self.m = np.kron(np.kron(temp1,G(self.theta[0],self.theta[1],opt)),temp2)
        
 
    def merge(self, block2):
        assert self.row == block2.row
        self.theta += block2.theta
        if self.row == 1:
            self.m = np.kron(G(self.theta[0],self.theta[1],opt),I)
        elif self.row == 2:
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
N = 10 #current N is even
p = N//2 # not merge yet
YBE_count = 0
rows = [[] for _ in range(N-1)] #[[]] * N is wrong
for j in range(2*p):
    for i in range(N-1):
        if j % 2 == 0:
            if i % 2 ==0:
                rows[i].append(block(row=i,col=j,N=N))
            else:
                rows[i].append([])
        else:
            if i % 2 == 0:
                rows[i].append([])
            else:
                rows[i].append(block(row=i,col=j,N=N))
print(f'Generate N={N}, depth={p} XY circuits')
ori = 1
ori2 = 1
for _ in range(N):
    ori=np.kron(ori,I)
    ori2=np.kron(ori2,I)
for j in range(2*p):
    for i in range(N-1):
        if rows[i][j] != []:
            ori=ori@rows[i][j].m 
### Test other origins         
# for j in range(2*p):        
#     for i in [0,2,4,1,3]: #permutation is allowed
#         if rows[i][j] != []:
#             ori2=ori2@rows[i][j].m 
# print(f'Ori == Ori2: {np.allclose(ori,ori2)}')      

#### Step 1
def move_first_row(new_rows,k=0): 
    global YBE_count
    for s in range((N-k)//2-1): 
        ind_ini = 2*s+2 
        ind = ind_ini #rightmost ind in V shape
        for row in range(k):
            new_rows[row][ind+1:ind+1] = [[],[]] 
        new_rows[k][ind+1:ind+1] = [[],[]] 
        new_rows[k+1][ind:ind] = [[],[]] 
        new_rows[k+2][ind-1:ind-1] = [[],[]] 
        for row in range(1,N-1-2-k): 
            new_rows[k+2+row][ind-1+row:ind-1+row] = [[],[]]
        p_temp = max(len(x) for x in new_rows)
        print(f'#{k+1} row: Start moving the #{s+1} block with {p_temp} layers') 
        ####
        for i in range(k,N-ind): 
            theta_init=np.concatenate((np.concatenate((new_rows[i][ind-2].theta,new_rows[i+1][ind-1].theta),axis=None), \
                                                  new_rows[i][ind].theta), axis=None)
            theta_trans = YBE_V2A(theta_init)
            YBE_count += 1
            if i == k:
                new_rows[i][ind-2] = []
                new_rows[i][ind] = []
                new_rows[i][ind] = block(theta=theta_trans[2:4],row=i,N=N) 
                new_rows[i+1][ind-1] = []
                new_rows[i+1][ind-1] = block(theta=theta_trans[0:2],row=i+1,N=N) 
                new_rows[i+1][ind+1] = block(theta=theta_trans[4:6],row=i+1,N=N)
                ind+=3
            else:
                new_rows[i][ind-2] = []
                new_rows[i][ind] = []
                new_rows[i][ind-2] = block(theta=theta_trans[2:4],row=i,N=N) 
                new_rows[i+1][ind-1] = []
                new_rows[i+1][ind-3] = block(theta=theta_trans[0:2],row=i+1,N=N) 
                new_rows[i+1][ind-1] = block(theta=theta_trans[4:6],row=i+1,N=N)
                ind+=1
            ####
            new = 1
            for _ in range(N):
                new=np.kron(new,I)
            for col in range(p_temp): 
                for row in range(N-1):
                    if new_rows[row][col] != []:
                        new=new@new_rows[row][col].m 
            print(f'#{i-k+1} YBE, Ori == New: {np.allclose(ori,new)}') 
        print(f'#{k+1} row: Finish moving the #{s+1} blocK')
        
        # move back some blocks
        for row in range(N-1): 
            for col in range(row-k+ind_ini+4, p_temp): 
                new_rows[row][col-2], new_rows[row][col] = new_rows[row][col], new_rows[row][col-2]

        for col in range(p_temp-1,-1,-1):
            temp=[row[col] for row in new_rows]
            if not any(temp):
                for row in range(N-1):
                    del new_rows[row][col]
        p_temp = max(len(x) for x in new_rows)
            
        new = 1
        for _ in range(N):
            new=np.kron(new,I)
        for col in range(p_temp): 
            for row in range(N-1):
                if new_rows[row][col] != []:
                    new=new@new_rows[row][col].m 
        print(f'Move back with layer {p_temp}, Ori == New: {np.allclose(ori,new)}\n')  
        
    return new_rows

def move_second_row(new_rows,k=1): 
    global YBE_count
    for s in range((N-k)//2-1): 
        ind_ini = 2*s+3
        ind = ind_ini #rightmost ind in V shape
        for row in range(k):
            new_rows[row][ind:ind] = [[],[]] 
        new_rows[k][ind+1:ind+1] = [[],[]] 
        new_rows[k+1][ind:ind] = [[],[]] 
        new_rows[k+2][ind-1:ind-1] = [[],[]] 
        for row in range(1,N-1-2-k): 
            new_rows[k+2+row][ind-1+row:ind-1+row] = [[],[]]
        p_temp = max(len(x) for x in new_rows)
        print(f'#{k+1} row: Start moving the #{s+1} block with {p_temp} layers') 
        ####
        for i in range(k,N-ind): 
            theta_init=np.concatenate((np.concatenate((new_rows[i][ind-2].theta,new_rows[i+1][ind-1].theta),axis=None), \
                                                  new_rows[i][ind].theta), axis=None)
            theta_trans = YBE_V2A(theta_init)
            YBE_count += 1
            if i == k:
                new_rows[i][ind-2] = []
                new_rows[i][ind] = []
                new_rows[i][ind] = block(theta=theta_trans[2:4],row=i,N=N) 
                new_rows[i+1][ind-1] = []
                new_rows[i+1][ind-1] = block(theta=theta_trans[0:2],row=i+1,N=N) 
                new_rows[i+1][ind+1] = block(theta=theta_trans[4:6],row=i+1,N=N)
                ind+=3
            else:
                new_rows[i][ind-2] = []
                new_rows[i][ind] = []
                new_rows[i][ind-2] = block(theta=theta_trans[2:4],row=i,N=N) 
                new_rows[i+1][ind-1] = []
                new_rows[i+1][ind-3] = block(theta=theta_trans[0:2],row=i+1,N=N) 
                new_rows[i+1][ind-1] = block(theta=theta_trans[4:6],row=i+1,N=N)
                ind+=1
            ####
            new = 1
            for _ in range(N):
                new=np.kron(new,I)
            for col in range(p_temp): 
                for row in range(N-1):
                    if new_rows[row][col] != []:
                        new=new@new_rows[row][col].m 
            print(f'#{i-k+1} YBE, Ori == New: {np.allclose(ori,new)}') 
        print(f'#{k+1} row: Finish moving the #{s+1} block')
        
        ### move back some blocks
        for row in range(N-1): 
            for col in range(row-k+ind_ini+4, p_temp): 
                new_rows[row][col-2], new_rows[row][col] = new_rows[row][col], new_rows[row][col-2]        
    
        for col in range(p_temp-1,-1,-1):
            temp=[row[col] for row in new_rows]
            if not any(temp):
                for row in range(N-1):
                    del new_rows[row][col]
        p_temp = max(len(x) for x in new_rows)
            
        new = 1
        for _ in range(N):
            new=np.kron(new,I)
        for col in range(p_temp): 
            for row in range(N-1):
                if new_rows[row][col] != []:
                    new=new@new_rows[row][col].m 
        print(f'Move back with layer {p_temp}, Ori == New: {np.allclose(ori,new)}\n')  
        
    return new_rows

def move_last_diag(new_rows,k=1): #fill N-1-k row
    global YBE_count
    for s in range(ceil(k/2)):
        ind_ini_row = N-2-k+2*s 
        if k==1:
            ind_ini_col = 2*k-1 #1
        else:
            ind_ini_col = k+1
        ind_row = ind_ini_row #top ind in A shape
        ind_col = ind_ini_col #top ind in A shape
        
        for row in range(ind_row+1):
            new_rows[ind_row-row][ind_col+row:ind_col+row] = [[],[]] 
        for row in range(ind_row+1,N-1): 
            new_rows[row][ind_col-1:ind_col-1] = [[],[]]
            
        p_temp = max(len(x) for x in new_rows)
        print(f'#{k+1} reverse-diagnal: Start moving the #{s+1} block with {p_temp} layers') 
        ####
        for i in range(ind_ini_row+1,N-2-k,-1): 
            # i is the bottom row of A
            
            if i < ind_ini_row+1: 
                ## update blank blocks
                ind_row-=1
                ind_col-=1
                for row in range(ind_row+1):
                    new_rows[ind_row-row][ind_col+row:ind_col+row] = [[],[]] 
                for row in range(ind_row+1,N-1): 
                    new_rows[row][ind_col-1:ind_col-1] = [[],[]]
            
                
            p_temp = max(len(x) for x in new_rows)    
            ###
            theta_init=np.concatenate((np.concatenate((new_rows[i][ind_col+2-1].theta,new_rows[i-1][ind_col+2].theta),axis=None), \
                                                  new_rows[i][ind_col+2+1].theta), axis=None)
            theta_trans = YBE_V2A(theta_init)
            YBE_count += 1
            
            new_rows[i][ind_col+2+1] = []
            new_rows[i][ind_col+2-1] = block(theta=theta_trans[2:4],row=i,N=N)
            new_rows[i-1][ind_col+2-2] = []
            new_rows[i-1][ind_col+2-2] = block(theta=theta_trans[0:2],row=i-1,N=N)
            new_rows[i-1][ind_col+2] = block(theta=theta_trans[4:6],row=i-1,N=N) 
            
            if i < ind_ini_row+1:
                ## move back the blank blocks in the last i
                last_ind_row = ind_row + 1 
                last_ind_col = ind_col + 1 + 2
                for row in range(last_ind_row+1):
                    for col in range(last_ind_col+row+2,p_temp):
                        new_rows[last_ind_row-row][col], new_rows[last_ind_row-row][col-2] = \
                        new_rows[last_ind_row-row][col-2], new_rows[last_ind_row-row][col]
                for row in range(last_ind_row+1,N-1): 
                    for col in range(last_ind_col+1,p_temp): #ind_col+(row-ind_ini_row-1)
                        new_rows[row][col-2], new_rows[row][col]= \
                        new_rows[row][col], new_rows[row][col-2]
                if i == N-2-k+1:
                    for row in range(ind_row+1):
                        for col in range(ind_col,p_temp):
                            new_rows[row][col-2], new_rows[row][col]= \
                            new_rows[row][col], new_rows[row][col-2]
                    for row in range(ind_row+1, N-1):
                        for col in range(ind_col+1,p_temp):
                            new_rows[row][col-2], new_rows[row][col]= \
                            new_rows[row][col], new_rows[row][col-2]
            if s>0 and i==ind_ini_row+1:
                for row in range(ind_ini_row+2,N-1):
                    for col in range(ind_col+2+2,ind_col+2+2+(row-ind_ini_row-2)+1):
                        new_rows[row][col-2], new_rows[row][col]= \
                        new_rows[row][col], new_rows[row][col-2]
            ####
            new = 1
            for _ in range(N):
                new=np.kron(new,I)
            for col in range(p_temp): 
                for row in range(N-1):
                    if new_rows[row][col] != []:
                        new=new@new_rows[row][col].m 
            print(f'#{ind_ini_row+1-i+1} YBE, Ori == New: {np.allclose(ori,new)}') 
        print(f'#{k+1} reverse-diagnal: Finish moving the #{s+1} blocK')
        
        # move back some blocks
        if s == 0 and k==1:
            for row in range(ind_ini_row+1):
                for col in range(ind_col+row,p_temp):
                    new_rows[ind_ini_row-row][col], new_rows[ind_ini_row-row][col-1] = \
                    new_rows[ind_ini_row-row][col-1], new_rows[ind_ini_row-row][col]
            for row in range(ind_ini_row+1,N-1): 
                for col in range(ind_col,p_temp):
                    new_rows[row][col-1], new_rows[row][col]= \
                    new_rows[row][col], new_rows[row][col-1]
        elif s == 0 and k>=2:
            for row in range(ind_ini_row+1):
                for col in range(ind_col+row,p_temp):
                    new_rows[ind_row-row][col], new_rows[ind_ini_row-row][col-2] = \
                    new_rows[ind_row-row][col-2], new_rows[ind_ini_row-row][col]
            for row in range(ind_ini_row+1,N-1): 
                for col in range(ind_col+1,p_temp): #ind_col+(row-ind_ini_row-1)
                    new_rows[row][col-2], new_rows[row][col]= \
                    new_rows[row][col], new_rows[row][col-2]
            for row in range(ind_ini_row+2,N-1):
                for col in range(ind_col+2,ind_col+2+(row-ind_ini_row-2)+1):
                    new_rows[row][col-2], new_rows[row][col]= \
                    new_rows[row][col], new_rows[row][col-2]
                
        for col in range(p_temp-1,-1,-1):
            temp=[row[col] for row in new_rows]
            if not any(temp):
                for row in range(N-1):
                    del new_rows[row][col]
            
        p_temp = max(len(x) for x in new_rows)
            
        new = 1
        for _ in range(N):
            new=np.kron(new,I)
        for col in range(p_temp): 
            for row in range(N-1):
                if new_rows[row][col] != []:
                    new=new@new_rows[row][col].m 
        print(f'Move back with layer {p_temp}, Ori == New: {np.allclose(ori,new)}\n')  
        
    return new_rows


new_rows = copy.deepcopy(rows) 
print(f'Step 1: Start constructing a trangular-like circuit\n' )
k=0
while( N-k > 3): 
    # If N-k > 4, then move back from k=3
    # If N-k > 3, then move back from k=1
    if k % 2 == 0:
        new_rows=move_first_row(new_rows,k=k)
    else:
        new_rows=move_second_row(new_rows,k=k)
    k+=1
print(f'Finish constructing a trangular-like circuit with #{YBE_count} YBE\n' )
# new_rows1 = move_first_row(new_rows,k=0)
# new_rows2 = move_second_row(new_rows1,k=1)
# new_rows3 = move_first_row(new_rows2,k=2)
# new_rows4 = move_second_row(new_rows3,k=3)
# m = [sublist[:2*p-2] for sublist in new_rows[1:N-1]] #slicing lst of lst

#### Step 2
print(f'Step 2: Start moving back a trangular-like circuit\n' )
# new_rows1 = move_last_diag(new_rows,k=1)
# new_rows2 = move_last_diag(new_rows1,k=2)
# new_rows3 = move_last_diag(new_rows2,k=3)
# new_rows4 = move_last_diag(new_rows3,k=4)
# new_rows5 = move_last_diag(new_rows4,k=5)
# new_rows6 = move_last_diag(new_rows5,k=6)
k=1
while( N-k >= 2): 
    new_rows = move_last_diag(new_rows,k=k)
    k+=1
print(f'Finish transformation with #{YBE_count} YBE\n' )

