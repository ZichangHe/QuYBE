import qiskit 
from qiskit import QuantumRegister, QuantumCircuit, execute, Aer
import numpy as np

def state_str2num(basis_state_as_str):
    return int(basis_state_as_str, 2)

def state_reverse(basis_state_as_num, nqubits):
    basis_state_as_str = state_num2str(basis_state_as_num, nqubits)
    new_str = basis_state_as_str[::-1]
    return state_str2num(new_str)

def get_adjusted_state(state):
    """
    Convert qubit ordering invert for state vector
    https://github.com/rsln-s/QAOA_tutorial/blob/main/Hands-on.ipynb
    """
    nqubits = np.log2(state.shape[0])
    if nqubits % 1:
        raise ValueError("Input vector is not a valid statevector for qubits.")
    nqubits = int(nqubits)

    adjusted_state = np.zeros(2**nqubits, dtype=complex)
    for basis_state in range(2**nqubits):
        adjusted_state[state_reverse(basis_state, nqubits)] = state[basis_state]
    return adjusted_state

def input_state_circ(n,circ):
    for i in range (0,n):
        if (i % 2) == 0:
         circ.x(i)
        else:
         circ.id(i)
    
    return circ

def trotter_circuit(n,Jx,Jy,Jz,delta_t,circ):
    hbar = 0.658212 
    if (Jz==0.0):
        gamma = (Jx+Jy)*delta_t/hbar
        delta = (Jx-Jy)*delta_t/hbar

    # Odd Layer
    for i in range (0,n-1,2):
      circ.rx(np.pi/2,i)
      circ.rx(np.pi/2,i+1)
      circ.cx(i,i+1)
      circ.rx(-(gamma+delta),i)
      circ.rz(-(gamma-delta),i+1)
      circ.cx(i,i+1)
      circ.rx(-np.pi/2,i)
      circ.rx(-np.pi/2,i+1)
      
    # Even Layer
    for i in range (1,n-1,2):
      circ.rx(np.pi/2,i)
      circ.rx(np.pi/2,i+1)
      circ.cx(i,i+1)
      circ.rx(-(gamma+delta),i)
      circ.rz(-(gamma-delta),i+1)
      circ.cx(i,i+1)
      circ.rx(-np.pi/2,i)
      circ.rx(-np.pi/2,i+1)
          
    return circ

def trotter_circuit_half(n,Jx,Jy,Jz,delta_t,circ):
    hbar = 0.658212 
    if (Jz==0.0):
        gamma = (Jx+Jy)*delta_t/hbar
        delta = (Jx-Jy)*delta_t/hbar

    # Odd Layer
    for i in range (0,n-1,2):
      circ.rx(np.pi/2,i)
      circ.rx(np.pi/2,i+1)
      circ.cx(i,i+1)
      circ.rx(-(gamma+delta),i)
      circ.rz(-(gamma-delta),i+1)
      circ.cx(i,i+1)
      circ.rx(-np.pi/2,i)
      circ.rx(-np.pi/2,i+1)
      
    return circ

def row_to_trotter_circuit(rows,circ):
    if rows[0][0] == []:
        opt = rows[0][1].opt
    else:
        opt = rows[0][0].opt
    if opt == 1:
        n = len(rows) + 1
        for col in range(len(rows[0])):
            for i in range(0, n-1):
                if rows[i][col] != []:
                    angle = rows[i][col].theta
                    circ.rx(np.pi/2,i)
                    circ.rx(np.pi/2,i+1)
                    circ.cx(i,i+1)
                    circ.rx(angle[0],i)
                    circ.rz(angle[1],i+1)
                    circ.cx(i,i+1)
                    circ.rx(-np.pi/2,i)
                    circ.rx(-np.pi/2,i+1)
    if opt == 2:
        n = len(rows) + 1
        for col in range(len(rows[0])):
            for i in range (0,n-1):
                if rows[i][col] != []:
                    angle = rows[i][col].theta
                    circ.cx(i,i+1)
                    circ.rx(angle[0],i)
                    circ.rz(angle[1],i+1)
                    circ.cx(i,i+1)
         
    if opt == 3:
        n = len(rows) + 1
        for col in range(len(rows[0])):
            for i in range (0,n-1):
                if rows[i][col] != []:
                    angle = rows[i][col].theta
                    circ.rz(np.pi/2,i)
                    circ.rz(np.pi/2,i+1)
                    circ.cx(i,i+1)
                    circ.rx(angle[0],i)
                    circ.rz(angle[1],i+1)
                    circ.cx(i,i+1)
                    circ.rz(-np.pi/2,i)
                    circ.rz(-np.pi/2,i+1)
              
    return circ
