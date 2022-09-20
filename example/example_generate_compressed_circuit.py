#Sahil Gulania & Zichang He#
#ANL & USCB 2022#
# Generating quantum circuit for quantum tim evolution of
# XY, YZ, ZX Hamiltonian both with compression and without 
# Compression
import numpy as np
import scipy as sp
import sys
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, IBMQ
from QuYBE.qiskit_utils import *
from QuYBE.utils import *

N_spin = 5
Jx = -0.8
Jy = -0.2
Jz = 0
del_t = 0.05 #np.pi/100
total_step = 10

for i in range (0,total_step):
    qr = QuantumRegister(N_spin, 'q')
    cr = ClassicalRegister(N_spin, 'c')
    circ = QuantumCircuit(qr, cr)
    # Initial state
    circ.barrier(qr)
    input_state_circ(N_spin,circ)
    circ.barrier(qr)
    # Initial state

    # Time evolution using trotter approximation
    for k in range(0,i):
        trotter_circuit(N_spin,Jx,Jy,Jz,del_t,circ)

    circ.barrier(qr)

    # Measurement 
    circ.measure(qr, cr)

    # Write QASM code to a file
    with open("qtd_"+str(N_spin)+"_qubits_trotter_"+str(i)+".qasm", "a") as o:
        o.write(circ.qasm())

    # YBE Compressed QASM Code 
    step_YBE = i
    p = step_YBE*2 # p>=N
    ### Generate XY circuit blocks/rows
    rows = generate_circuit(N_spin, p, Jx, Jy, Jz, del_t, debug=False)
    
    merged, YBE_count = merge_to_minimal(rows,debug=False,verbose=False)
    #print(f'Finish all transformation with #{YBE_count} YBE')
    YBE_circ = QuantumCircuit(N_spin)
    YBE_circ = input_state_circ(N_spin,YBE_circ)
    YBE_circ = row_to_trotter_circuit(merged,YBE_circ)
    YBE_circ.measure_all(True)
    with open("qtd_"+str(N_spin)+"_qubits_trotter_YBE_"+str(i)+".qasm", "a") as o_ybe:
        o_ybe.write(YBE_circ.qasm())


