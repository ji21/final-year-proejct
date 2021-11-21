row_list = [["center", "linkers", "porphyrins", "hl gap", "homo", "lumo"]]

linker_carbons = [2,4,6]
#add 6 later
repeating_units = [1,2,3,4,5]
centers = ["Zn+2", "Fe+2", "Ni+2", "Co+2", "Cu+2"]

for r in repeating_units:
    for c in linker_carbons:
        for center in centers:
            fname = f'{center}_{c}_{r+1}_out' 
