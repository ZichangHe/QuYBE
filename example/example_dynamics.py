from QuYBE.qiskit_utils import *
from QuYBE.utils import *
import pickle

N = 5
Jx = -0.8
Jy = 0.1
Jz = 0
delta_t = 0.050 
step = 10
p = step*2 # p>=N
# backend = Aer.get_backend('statevector_simulator')
backend = Aer.get_backend('unitary_simulator')
save = False
debug = False

### Generate XY circuit blocks/rows
rows = generate_circuit(N, p, Jx, Jy, Jz, delta_t, debug=debug)
if debug == True:
    ori_unitary = sv_simulation(rows)

### Compression   
import time
start = time.time()
merged, YBE_count = merge_to_minimal(rows, debug=debug, verbose=False)
end = time.time()
print(f'Finish all transformation with #{YBE_count} YBE')
print(f'The time of execution of above program is :{(end-start)}')
circ = QuantumCircuit(N)
# circ = input_state_circ(N,circ)

### Verify the compressed circuit and rows
circ = row_to_trotter_circuit(merged,circ) 
if debug == True:
    # YBE_state = execute(circ, backend).result().get_statevector()
    YBE_unitary = execute(circ, backend).result().get_unitary(circ)
    print(np.allclose(ori_unitary, YBE_unitary))

### Print the qasm circuit
# print(circ.qasm())

if save is True:
    f = open(
            f"N{N}_Jx{Jx}_Jy{Jy}_delta_t{delta_t}_step{step}_circuit.pckl",
            "wb",
        )
    pickle.dump(
        [circ],
        f,
    )
    f.close()
    
    