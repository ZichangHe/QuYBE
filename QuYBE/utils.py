import copy
from .core_utils import *

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
        # self.col = col
        self.row = row
        # self.ind = {'row':self.row, 'col':self.col}
        self.N=N
        self.opt=opt
        temp1, temp2 = 1, 1
        for _ in range(row):
            temp1 = np.kron(I,temp1)
        for _ in range(row+1,N-1):
            temp2 = np.kron(I,temp2)
        self.m = np.kron(np.kron(temp1,G(self.theta[0],self.theta[1],opt)),temp2)
        
    def merge(self, block2):
        assert self.row == block2.row
        self.theta += block2.theta
        temp1, temp2 = 1, 1
        for _ in range(self.row):
            temp1 = np.kron(I,temp1)
        for _ in range(self.row+1,self.N-1):
            temp2 = np.kron(I,temp2)
        self.m = np.kron(np.kron(temp1,G(self.theta[0],self.theta[1],self.opt)),temp2)

def generate_circuit(N, p, theta=None):
    rows = [[] for _ in range(N-1)] 
    for j in range(p):
        for i in range(N-1):
            theta = theta
            # if theta is None, generate random angles by default
            if j % 2 == 0:
                if i % 2 ==0:
                    rows[i].append(block(theta=theta,row=i,col=j,N=N))
                else:
                    rows[i].append([])
            else:
                if i % 2 == 0:
                    rows[i].append([])
                else:
                    rows[i].append(block(theta=theta,row=i,col=j,N=N))
    print(f'Generate N={N}, depth={p}, XY circuits')
    return rows

def sv_simulation(rows):
    output = 1
    N = len(rows)+1
    for _ in range(N):
        output = np.kron(output,I)
    for col in range(len(rows[0])): 
        for row in range(N-1):
            if rows[row][col] != []:
                if N % 2 == 0:
                    output = rows[row][col].m @ output
                else: ### To match the results from qiskit simulator
                    output = output @ rows[row][col].m 
    return output

def J2theta(Jx,Jy,delta_t):
    ''' From dynamics parameters to rows' parameters'''
    # Jz = 0
    hbar = 0.658212 
    gamma = (Jx+Jy)*delta_t/hbar
    delta = (Jx-Jy)*delta_t/hbar
    theta = np.array([-(gamma+delta),-(gamma-delta)])
    return theta

def merge_to_minimal(rows,ori=None,verbose=False):
    YBE_count = 0
    N = len(rows)+1
    p = len(rows[0])
    new_rows = [sublist[:N] for sublist in rows] #slicing lst of lst
    k = N
    reflect_count = 0
    while k < p:
        print(f'======Start merging #{reflect_count+1} layer======\n')
        to_merge = [sublist[k] for sublist in rows]
        if reflect_count % 2 == 0:
            new_rows, new_YBE_count = transform_V2A(new_rows,verbose)
            YBE_count += new_YBE_count
        else:
            new_rows, new_YBE_count = transform_A2V(new_rows,verbose)
            YBE_count += new_YBE_count
        for i in range(N-1):
            if to_merge[i]!=[]:
                new_rows[i][-1].merge(to_merge[i])
        k+=1
        reflect_count+=1
    
    ## generate state vector from compressed circuits 
    new = sv_simulation(new_rows)
      
    ## generate state vector from original circuits           
    if ori is None: 
        ori = sv_simulation(rows)
        
    print(f'Finish merging, Ori == New: {np.allclose(ori,new)}\n') #{k+1} layer
    return new_rows, YBE_count              
  
              
def move_first_row(new_rows,k=0,mode='V2A',verbose=False,debug=True): 
    '''
    The first row means that the first block is NOT empty
    '''
    YBE_count = 0
    N = len(new_rows) + 1
    if mode == 'V2A': #when it is used in transform_V2A
        if N % 2 == 0:
            s_size = (N-k)//2-1
        else:
            s_size = (N-k-1)//2
    elif mode == 'A2V': #when it is used in transform_A2V
        if N % 2 == 0:    
            s_size = (N-1-k)//2 
        else:
            s_size = (N-k)//2-1
    for s in range(s_size): 
        ind_ini = 2*s+2 
        ind = ind_ini #rightmost ind in V shape
        for row in range(k):
            new_rows[row][ind+1:ind+1] = [[],[]] 
        new_rows[k][ind+1:ind+1] = [[],[]] 
        new_rows[k+1][ind:ind] = [[],[]] 
        if k+2 <= N-2:
            new_rows[k+2][ind-1:ind-1] = [[],[]] 
            for row in range(1,N-1-2-k): 
                new_rows[k+2+row][ind-1+row:ind-1+row] = [[],[]]
        p_temp = max(len(x) for x in new_rows)
        if verbose is True:
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
            if debug is True:
                new = sv_simulation(new_rows)
                print(f'#{i-k+1} YBE, Ori == New: {np.allclose(ori,new)}') 
        if verbose is True:
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
        if debug is True:    
            new = sv_simulation(new_rows)
            print(f'Move back with layer {p_temp}, Ori == New: {np.allclose(ori,new)}\n')  
        
    return new_rows, YBE_count


def move_second_row(new_rows,k=1,mode='V2A',verbose=False,debug=False): 
    '''
    The second row means that the first block is empty
    '''
    YBE_count = 0
    N = len(new_rows) + 1
    if mode=='V2A': #when it is used in transform_V2A
        if N % 2 == 0:
            s_size = (N-k)//2-1
        else:
            s_size = (N-1-k)//2
    elif mode == 'A2V': #when it is used in transform_A2V
        if N % 2==0:
            s_size = (N-1-k)//2
        else:
            s_size = (N-k)//2-1 
            
    for s in range(s_size): 
        ind_ini = 2*s+3
        ind = ind_ini #rightmost ind in V shape
        for row in range(k):
            new_rows[row][ind:ind] = [[],[]] 
        new_rows[k][ind+1:ind+1] = [[],[]] 
        new_rows[k+1][ind:ind] = [[],[]] 
        if k+2 <= N-2:
            new_rows[k+2][ind-1:ind-1] = [[],[]] 
            for row in range(1,N-1-2-k): 
                new_rows[k+2+row][ind-1+row:ind-1+row] = [[],[]]
        p_temp = max(len(x) for x in new_rows)
        if verbose is True:
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
            if debug is True:
                new = sv_simulation(new_rows)
                print(f'#{i-k+1} YBE, Ori == New: {np.allclose(ori,new)}') 
        if verbose is True:
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
        if debug is True:    
            new = sv_simulation(new_rows)
            print(f'Move back with layer {p_temp}, Ori == New: {np.allclose(ori,new)}\n')  
        
    return new_rows, YBE_count

def move_diag_to23(new_rows,k=1,verbose=False,debug=False): #fill N-1-k row
    '''
    --------------------|theta_trans[2]|--------------------
    --|theta_trans[0]|--|theta_trans[3]|--|theta_trans[4]|--
    --|theta_trans[1]|--------------------|theta_trans[5]|--
    23: The first block of the first row is empty
    '''
    YBE_count = 0
    N = len(new_rows) + 1
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
        if verbose is True:
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
            theta_trans = YBE_A2V(theta_init)
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
            if debug is True:
                new = sv_simulation(new_rows)
                print(f'#{ind_ini_row+1-i+1} YBE, Ori == New: {np.allclose(ori,new)}') 
        if verbose is True:
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
        if debug is True:    
            new = sv_simulation(new_rows)
            print(f'Move back with layer {p_temp}, Ori == New: {np.allclose(ori,new)}\n')  
        
    return new_rows, YBE_count


def move_diag_to12(new_rows,k=2,verbose=False,debug=False): #fill N-1-k row
    '''
    --|theta_init[0]|-------------------|theta_init[4]|-- 
    --|theta_init[1]|--|theta_init[2]|--|theta_init[5]|-- = 
    -------------------|theta_init[3]|-------------------   
    12: The first block of the first row is NOT empty
    '''
    YBE_count = 0
    N = len(new_rows) + 1
    for s in range(int(k/2)):
        ind_ini_row = N-k-1+2*s   
        ind_ini_col = k+1
        
        ind_row = ind_ini_row #top ind in A shape
        ind_col = ind_ini_col #top ind in A shape
        
        for row in range(ind_row+1):
            new_rows[ind_row-row][ind_col+row:ind_col+row] = [[],[]] 
        for row in range(ind_row+1,N-1): 
            new_rows[row][ind_col-1:ind_col-1] = [[],[]]
            
        p_temp = max(len(x) for x in new_rows)
        if verbose == True:
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
            theta_trans = YBE_A2V(theta_init)
            YBE_count += 1
            new_rows[i][ind_col+2+1] = []
            new_rows[i][ind_col+2-1] = block(theta=theta_trans[2:4],row=i,N=N)
            new_rows[i-1][ind_col+2-2] = []
            new_rows[i-1][ind_col+2-2] = block(theta=theta_trans[0:2],row=i-1,N=N)
            new_rows[i-1][ind_col+2] = block(theta=theta_trans[4:6],row=i-1,N=N) 
            
            ###
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

            if i == ind_ini_row+1:
                for row in range(ind_ini_row+2,N-1):
                    for col in range(ind_col+2+2,ind_col+2+2+(row-ind_ini_row-2)+1):
                        new_rows[row][col-2], new_rows[row][col]= \
                        new_rows[row][col], new_rows[row][col-2]
            ###
            if debug is True:
                new = 1
                for _ in range(N):
                    new=np.kron(new,I)
                for col in range(p_temp): 
                    for row in range(N-1):
                        if new_rows[row][col] != []:
                            new=new@new_rows[row][col].m 
                print(f'#{ind_ini_row+1-i+1} YBE, Ori == New: {np.allclose(ori,new)}') 
        if verbose is True:
            print(f'#{k+1} reverse-diagnal: Finish moving the #{s+1} blocK')
        
        for col in range(p_temp-1,-1,-1):
            temp=[row[col] for row in new_rows]
            if not any(temp):
                for row in range(N-1):
                    del new_rows[row][col]
            
        p_temp = max(len(x) for x in new_rows)
        if debug is True:  
            new = sv_simulation(new_rows)
            print(f'Move back with layer {p_temp}, Ori == New: {np.allclose(ori,new)}\n')  
        
    return new_rows, YBE_count


def transform_V2A(rows,verbose=False,debug=False):
    N = len(rows) + 1
    YBE_count = 0
    new_rows = copy.deepcopy(rows) 
    print(f'====== Step 1: Start constructing a trangular-like circuit ======\n' )
    k = 0
    if N % 2 == 0:
        end = 3
    else:
        end = 2
    while( N-k > end): 
        if k % 2 == 0:
            new_rows, new_YBE_count=move_first_row(new_rows,k=k,verbose=verbose,debug=debug)
            YBE_count += new_YBE_count
        else:
            new_rows, new_YBE_count=move_second_row(new_rows,k=k,verbose=verbose,debug=debug)
            YBE_count += new_YBE_count
        k+=1
    #print(f'Finish constructing a trangular-like circuit with #{YBE_count} YBE\n' )

    #### Step 2
    print(f'====== Step 2: Start moving back a trangular-like circuit ======\n' )
    k = 1
    while( N-k >= 2): 
        if N % 2 == 0:
            new_rows, new_YBE_count = move_diag_to23(new_rows,k=k,verbose=verbose,debug=debug)
            YBE_count += new_YBE_count
        else:
            new_rows, new_YBE_count = move_diag_to12(new_rows,k=k,verbose=verbose,debug=debug)
            YBE_count += new_YBE_count
        k+=1
    #print(f'Finish transformation with #{YBE_count} YBE\n' )
    return new_rows, YBE_count

def transform_A2V(rows,verbose=False,debug=False):
    N = len(rows) + 1
    YBE_count = 0
    new_rows = copy.deepcopy(rows) 
    #### Step 1: 
    print(f'====== Step 1: Start constructing a trangular-like circuit ======\n' )
    k = 0
    if N % 2 == 0:
        end = 3
    else:
        end = 2
    while( N-k >= end): 
        if k % 2 == 0:
            new_rows, new_YBE_count = move_second_row(new_rows,k=k,mode='A2V',verbose=verbose,debug=debug) 
            YBE_count += new_YBE_count
        else:
            new_rows, new_YBE_count = move_first_row(new_rows,k=k,mode='A2V',verbose=verbose,debug=debug)
            YBE_count += new_YBE_count
        k+=1
    #print(f'Finish constructing a trangular-like circuit with #{YBE_count} YBE\n' )

    #### Step 2: need a different move back
    print(f'====== Step 2: Start moving back a trangular-like circuit ======\n' )
    k = 1
    while( N-k >= 2): 
        if N%2 == 0:
            new_rows, new_YBE_count = move_diag_to12(new_rows,k=k,verbose=verbose,debug=debug)
            YBE_count += new_YBE_count
        else:
            new_rows, new_YBE_count = move_diag_to23(new_rows,k=k,verbose=verbose,debug=debug)
            YBE_count += new_YBE_count
        k+=1
    #print(f'Finish transformation with #{YBE_count} YBE\n' )
    return new_rows, YBE_count

