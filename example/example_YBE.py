from QuYBE.utils import generate_circuit, sv_simulation, transform_V2A, transform_A2V, merge_to_minimal

if __name__ == "__main__":
    N = 5 # N can be both even and odd
    p = 5 # p>=N
    ### Generate XY circuit blocks
    rows = generate_circuit(N, p)
    ### Ground truth simulation result
    ori = sv_simulation(rows)
    
    if p == N:
        verbose = True
        YBE_count = 0
        new23, new_YBE_count = transform_V2A(rows,verbose)
        YBE_count += new_YBE_count
        new12, new_YBE_count = transform_A2V(new23,verbose)
        YBE_count += new_YBE_count
        print(f'Finish all transformation with #{YBE_count} YBE')
    else:
        verbose = True
        merged, YBE_count = merge_to_minimal(rows,ori,verbose)
        print(f'Finish all transformation with #{YBE_count} YBE')
        
        