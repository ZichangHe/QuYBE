# QuYBE

Exactly simulate 1D-Heisenberg dynamics with a compressed circuit. This is the implementation of the algorithm in *[Quantum time dynamics employing the Yang-Baxter equation for circuit compression](https://journals.aps.org/pra/pdf/10.1103/PhysRevA.106.012412)*

## Installation

QuYBE is written in Python 3, with only dependency on numpy and qiskit (optional). You can enter the fold and install with 

```
pip install .
```

## Example

Run the compression XY-model in analytical operators by 
    
```
python example/example_YBE.py
```

Run the compression of XY Heisenberg dynamics model in qiskit circuits by 

```
python example/example_dynamics.py
```

## Cite
If you find this repository useful, please cite us at ...


