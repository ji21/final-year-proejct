import os.path

linker_carbons = [2,4,6]
#repeating_units = [1,2,3,4,5,6,7,8]
repeating_units = [1,2,3,4,5,6]
centers = ["", "Zn+2", "Fe+2", "Ni+2", "Cu+2", "Co+2"]
missing_files = []
total = 0
for center in centers:
    for r in repeating_units:
        for c in linker_carbons:
            total += 1
            fname = f'{center}_{c}_{r+1}_out/xtbopt.xyz'
            #print(fname)
            if not os.path.isfile(fname):
                missing_files.append((center,c,r+1))

try:
    afsdf
except Exception as e:
    print(e)

missing_files.sort()
print("optimimsations that failed silently: ", len(missing_files))
print("list of directories that wasn't optimised: ")
for file in missing_files:
    print(f'metal center: {file[0]}, linker_carbons: {file[1]}, number_of_porphyrins: {file[2]}')

print("expected optimisations: ", total)
print("optimisations that actually worked: ", total-len(missing_files))
