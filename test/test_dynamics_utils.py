from QuYBE.qiskit_utils import input_state_circ,trotter_circuit,trotter_circuit_half,row_to_trotter_circuit
from QuYBE.utils import merge_to_minimal,generate_circuit,sv_simulation, J2theta
from qiskit import QuantumRegister, QuantumCircuit, execute, Aer
import numpy as np

if __name__ == "__main__":
    N = 5 # N can be both even and odd
    step = 10
    p = step*2 # p>=N
    
    ### get Initial state
    backend = Aer.get_backend('statevector_simulator')
    circ = QuantumCircuit(N)
    circ = input_state_circ(N,circ)
    result = execute(circ, backend).result()
    ini_state = result.get_statevector()
    
    ### Generate XY dynamics circuits 
    Jx = -0.8
    Jy = -0.2
    delta_t=0.05
    circ = QuantumCircuit(N)
    circ = input_state_circ(N,circ)
    backend = Aer.get_backend('statevector_simulator')
    # backend = Aer.get_backend('unitary_simulator')
    for _ in range(int(p//2)):
        trotter_circuit(N,Jx,Jy,0,delta_t,circ)
    if p % 2 != 0:
        trotter_circuit_half(N,Jx,Jy,0,delta_t,circ)
    state = execute(circ, backend).result().get_statevector()
    #circ_unitary = execute(circ, backend).result().get_unitary(circ)
    
    
    ### Generate XY circuit blocks/rows
    theta = J2theta(Jx,Jy,delta_t)
    rows = generate_circuit(N, p, theta=theta)
    ### Ground truth simulation result
    ori_unitary = sv_simulation(rows)
    ori_state = ori_unitary @ ini_state
    ### Verify the uncompressed circuit and rows
    print(np.allclose(state, ori_state))
    
    
    ### YBE compression
    merged, YBE_count = merge_to_minimal(rows,ori_unitary,verbose=False)
    print(f'Finish all transformation with #{YBE_count} YBE')
    circ = QuantumCircuit(N)
    circ = input_state_circ(N,circ)
    ### Verify the compressed circuit and rows
    circ = row_to_trotter_circuit(merged,circ)
    YBE_state = execute(circ, backend).result().get_statevector()
    print(np.allclose(state, YBE_state))
    
    
