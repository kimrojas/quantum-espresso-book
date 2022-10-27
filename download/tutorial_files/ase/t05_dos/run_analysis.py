with open('./run_files/run_calc.dat', 'w', buffering=1) as f:
    print("", file=f)
    print("LOWDIN CHARGE", "- "*10, file=f)
    
    with open('./run_files/projwfc.out', 'r') as fout:
        lines = fout.read().splitlines()
    
    copy = False
    for line in lines:
        if 'Lowdin Charges' in line:
            copy = True
            continue
        
        if 'PROJWFC' in line:
            copy = False
            continue

        if copy:
            print(line, file=f) 
