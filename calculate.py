import pandas as pd
import numpy as np
import os
import stk
import stko
import csv
import stk
import stko
import os
class Polymer():
    def __init__(self, center, linker_carbons, repeating_units, charge):
        self.dir = f'{center}_{linker_carbons}_{repeating_units+1}_out'
        if center != "":
            center = "[" + center + "]"
        self.center = center
        self.linker_carbons = linker_carbons
        self.repeating_units = repeating_units
        self.charge = charge
        self.polymer = None
        self.optimised = False
        self.result = None
        self.repeating_block = None
        self.mid_bb = None
        self.bbs = []
        print("init ok, will output to: ", self.dir)
        
    def get_linker_smiles(self):
        linker_carbons = self.linker_carbons-2
        linker = "C#C"
        for i in range(0, linker_carbons, 2):
            linker += "C#C"
        linker = "Br" + linker + "I"
        return linker
    
    def get_porphyrin_smiles(self, functional_group_1="", functional_group_2=""):
        if self.center == "":
            return f'C1=CC2=C{functional_group_1}C3=CC=C(N3)C=C4C=CC(=N4)C{functional_group_2}=C5C=CC(N5)=CC1(=N2)'
        return f'C1=CC2=C{functional_group_1}C3=CC=C([N]3)C=C4C=CC(=N4)C{functional_group_2}=C5C=CC(=N5)C=C1[N]2'

    def set_building_blocks_without_center(self):
        Br = "(Br)"
        I = "(I)"
        self.linker = stk.BuildingBlock(self.get_linker_smiles(), functional_groups=[stk.BromoFactory(), stk.IodoFactory()])
        self.first_bb = stk.BuildingBlock(smiles=self.get_porphyrin_smiles(Br), functional_groups=[stk.BromoFactory()])
        self.last_bb = stk.BuildingBlock(smiles=self.get_porphyrin_smiles(functional_group_2=I), functional_groups=[stk.IodoFactory()])
        self.mid_bb = stk.BuildingBlock(smiles=self.get_porphyrin_smiles(Br, I), functional_groups=[stk.BromoFactory(), stk.IodoFactory()])
        self.bbs = [self.first_bb]
        for i in range(1, self.repeating_units):
            self.bbs.append(self.linker)
            self.bbs.append(self.mid_bb)
        self.bbs.append(self.linker)
        self.bbs.append(self.last_bb)

    def set_building_blocks_with_center(self):
        Br = "(Br)"
        I = "(I)"
        if "Zn" in self.center:
            atom =stk.SingleAtom(stk.Zn(0, charge=2))
        elif "Fe" in self.center:
            atom =stk.SingleAtom(stk.Fe(0, charge=2))
        elif "Mg" in self.center:
            atom =stk.SingleAtom(stk.Mg(0, charge=2))
        elif "Mn" in self.center:
            atom =stk.SingleAtom(stk.Mn(0, charge=2))
        
        self.linker = stk.BuildingBlock(self.get_linker_smiles(), functional_groups=[stk.BromoFactory(), stk.IodoFactory()])
        self.first_bb = self.build_porphyrin_with_center(Br)
                #add building block with center
        self.first_bb = stk.BuildingBlock.init_from_molecule(self.first_bb, functional_groups = [stk.BromoFactory()])
        self.last_bb = self.build_porphyrin_with_center(I)
        self.last_bb = stk.BuildingBlock.init_from_molecule(self.last_bb, functional_groups = [stk.IodoFactory()])
        self.mid_bb = self.build_porphyrin_with_center(Br, I)
        self.mid_bb = stk.BuildingBlock.init_from_molecule(self.mid_bb, functional_groups = [stk.BromoFactory(), stk.IodoFactory()])
        self.bbs = [self.first_bb]
        for i in range(1, self.repeating_units):
            self.bbs.append(self.linker)
            self.bbs.append(self.mid_bb)
        self.bbs.append(self.linker)
        self.bbs.append(self.last_bb)
            
            
    def build(self):
        #builds molecule using stk
        if self.center == "":
            self.set_building_blocks_without_center()
        else:
            self.set_building_blocks_with_center()
        pattern = "A"
        if self.repeating_units > 1:
            pattern += (self.repeating_units-1)*"BC"
            pattern += "BA"
        else:
            pattern += "BC"
        print("patter :", pattern)
        if self.center == "":
            self.polymer = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=tuple(self.bbs),
                    repeating_unit=pattern,
                    num_repeating_units=1,
                    optimizer=stk.Collapser(scale_steps=False)
                )      
            )
        else:
            self.polymer = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=tuple(self.bbs),
                    repeating_unit=pattern,
                    num_repeating_units=1,
                    optimizer=stk.Collapser(scale_steps=False)
                )
            )

        print("finished building polymer")
        return

    def calculate(self):
        if self.center == "":
            fname = './free_base/' + self.dir + '/xtbopt.xyz'
            out = './free_base/' + self.dir + "energy"
        else:
            fname = self.dir + '/xtbopt.xyz'
            out = self.dir + "energy"
        p = self.polymer.with_structure_from_file(path=fname)
        print(fname)

        xtb = stko.XTBEnergy(
                xtb_path='/home/jeff/miniconda3/bin/xtb',
                unlimited_memory=True,
                output_dir=out,
        )
        r = xtb.get_results(p)
        try:
            self.gap = r.get_homo_lumo_gap()[0]
            hl_orbitals = r.get_homo_lumo_orbitals()
            self.homo = hl_orbitals[0]['HOMO'][2]
            self.lumo = hl_orbitals[0]['LUMO'][2]
            self.energy = r.get_total_energy()[0]
            self.fermi = r.get_fermi_level()[0]
            print(r.get_full_dipole_moments())
            print(hl_orbitals)
            print(r.get_homo_lumo_gap())
        except Exception as e:
            print(e)
            return False
        return True


    def build_porphyrin_with_center(self, functional_group_1="", functional_group_2=""):
        if "Zn" in self.center:
            atom =stk.SingleAtom(stk.Zn(0, charge=2))
        elif "Fe" in self.center:
            atom =stk.SingleAtom(stk.Fe(0, charge=2))
        elif "Co" in self.center:
            atom =stk.SingleAtom(stk.Co(0, charge=2))
        elif "Cu" in self.center:
            atom =stk.SingleAtom(stk.Cu(0, charge=2))
        elif "Ni" in self.center:
            atom =stk.SingleAtom(stk.Ni(0, charge=2))
        elif "Mn" in self.center:
            atom =stk.SingleAtom(stk.Mn(0, charge=2))



        metal = stk.BuildingBlock(
            smiles=self.center,
            functional_groups=(
                atom
                for i in range(4)
            ),
            position_matrix=[[0, 0, 0]],
        )
        
        porphyrin = stk.BuildingBlock(
            smiles=(
                self.get_porphyrin_smiles(functional_group_1, functional_group_2)
            ),
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7]~[#6]',
                    bonders=(1,),
                    deleters=(),
                ),
            ],
        )
        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.Porphyrin(
                metals=metal,
                ligands=porphyrin,
            ),
        )
        return complex
        
    def check_polymer(self):
        if self.polymer == None:
            raise ValueError("Polymer has not been built yet") 
        return
    
    def check_optimised(self):
        self.check_polymer()
        if self.optimised != True:
            raise ValueError("Polymer has not been optimised yet")
                
    
    def optimise(self):
        #calculates energy using stko
        if self.center == "":
            opt = stko.OptimizerSequence(
                stko.ETKDG(),
                stko.XTB(xtb_path='/home/jeff/miniconda3/bin/xtb',
                    unlimited_memory=True,
                    num_cores=1,
                    output_dir=self.dir,
                )
            )
        else:
            opt = stko.OptimizerSequence(
                stko.XTB(xtb_path='/home/jeff/miniconda3/bin/xtb',
                    unlimited_memory=True,
                    num_cores=1,
                    output_dir=self.dir,
                )
            )
        print("optimising")
        try:
            self.polymer = opt.optimize(mol=self.polymer)
        except Exception as e:
            print(e)
        print("finished optimising")
        return

row_list = [["center", "linkers", "porphyrins", "gap", "homo", "lumo", "energy", "fermi"]]
linker_carbons = [2,4,6]
#add 6 later
repeating_units = [1,2,3,4,5,6,7,8]
centers = ["", "Zn+2", "Fe+2", "Ni+2", "Co+2", "Cu+2"]

#polymer = Polymer(center="Co+2", linker_carbons=2, repeating_units=5, charge=0)
#polymer.build()
#polymer.calculate()


#print(fadfa)


for r in repeating_units:
    for c in linker_carbons:
        for center in centers:
            if center == "":
                fname = f'./free_base/{center}_{c}_{r+1}_out/xtbopt.xyz'
            else:
                fname = f'{center}_{c}_{r+1}_out/xtbopt.xyz'
            if os.path.isfile(fname):
                print(fname)
                polymer = Polymer(center=center, linker_carbons=c, repeating_units=r, charge=0)
                polymer.build()
                if polymer.calculate():
                    features = [center, c, r+1, polymer.gap, polymer.homo, polymer.lumo, polymer.energy, polymer.fermi]
                    row_list.append(features)
        print(row_list)
             

with open('data.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(row_list)
print("SUCCESS")

#polymer = Polymer(center="Zn+2", linker_carbons=6, repeating_units=4, charge=0)
#polymer.build()
#polymer.calculate()
