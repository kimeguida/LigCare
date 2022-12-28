

# Copyright (c) 2022 LIT
# For academic use only.
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.




import sys
import os
import numpy as np
import networkx as nx
import argparse
from time import strftime, localtime
import copy






def lig_mol2_info(mol2_):
    with open(mol2_, 'r') as f:
        mol2 = f.read()
        
    # get atom block
    mol2_atom_block = mol2.split("@<TRIPOS>ATOM\n")[1].split("@<TRIPOS>")[0]
    mol2_atom_block = mol2_atom_block.split("\n")
    del mol2_atom_block[-1]

    atom_types = {}
    coordinates = {}
    for l in mol2_atom_block:
        mol2lines = l.split()
        try:
            indice = int(mol2lines[0]) # mol2 indice starts at 1
            type_atm = mol2lines[5]
            x = float(mol2lines[2])
            y = float(mol2lines[3])
            z = float(mol2lines[4])
            atom_types[indice] = type_atm
            coordinates[indice] = [x, y, z]
        except ValueError:
            continue
    mol2_bond_block = mol2.split("@<TRIPOS>BOND\n")[1].split("@<TRIPOS>")[0]
    mol2_bond_block = mol2_bond_block.split("\n")
    del mol2_bond_block[-1]
    connectivity = {}
    bond_types = {}
    for l in mol2_bond_block:
        mol2line = l.split()
        try:
            idx_1 = int(mol2line[1])
            idx_2 = int(mol2line[2])
            bond_type = mol2line[3]
            bond_types[(idx_1, idx_2)] = bond_type
        except ValueError:
            continue


        if idx_1 not in connectivity:
            connectivity[idx_1] = [[atom_types[idx_2].split('.')[0].upper()],
                                 [idx_2]]
        else:
            connectivity[idx_1][0].append(atom_types[idx_2].split('.')[0].upper())
            connectivity[idx_1][1].append(idx_2)
 
        if idx_2 not in connectivity:
            connectivity[idx_2] = [[atom_types[idx_1].split('.')[0].upper()],
                                 [idx_1]]
        else:
            connectivity[idx_2][0].append(atom_types[idx_1].split('.')[0].upper())
            connectivity[idx_2][1].append(idx_1)
    
    return coordinates, connectivity, atom_types, bond_types




def add_new_ph4(coordinates_, index_, new_ph4_, coordinates_plus_):
    new_index = len(coordinates_) + len(coordinates_plus_) + 1
    coords = copy.deepcopy(coordinates_[index_][:3])
    #coordinates_plus_[new_index] = coords
    coordinates_plus_[new_index] = coords
    coordinates_plus_[new_index].append(new_ph4_)




def lig_ph4(coordinates_, connectivity_, atom_types_, duplicate_=False):
    
    ACCEPTED_ATOMS = ['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Du']
    HETEROATOMS = ['O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I']
    HYDROGENS = ['H']

    coordinates_h = []
    #print('Total atoms:', len(coordinates_))
    coordinates_plus = {}

    for i in coordinates_:

        if atom_types_[i].split('.')[0] not in ACCEPTED_ATOMS:
            continue
        
        # Hydrophobic rules                   
        if atom_types_[i] in ['C.2', 'C.3', 'S.2', 'S.3']:
            if all([atom not in HETEROATOMS for atom in connectivity_[i][0]]):
                coordinates_[i].append('CA')
            else:
                coordinates_[i].append(None)
           
        elif atom_types_[i] == 'F':
            # discussion with Didier: not hydrophobic, maybe weak H-bond acceptor
            coordinates_[i].append(None)
        
        elif atom_types_[i] in ['Cl', 'Br', 'I']:
            coordinates_[i].append('CA')

        


        # Aromatic rules
        elif atom_types_[i] == 'C.ar':
            coordinates_[i].append('CZ')
            # if duplicate, added after ring processing

        elif atom_types_[i] == 'C.1':
            if all([atom not in HETEROATOMS for atom in connectivity_[i][0]]):
                coordinates_[i].append('CZ')  # before was CA
            # Discussion with Didier, no duplicate into CA
            else:
                coordinates_[i].append(None)




        
        # H-bond rules
        elif atom_types_[i] == 'N.ar':
            if any([atom in HYDROGENS for atom in connectivity_[i][0]]):
                coordinates_[i].append('N')
                # can interact with protein OH (-->OG), O- (-->NZ)
                if duplicate_:
                    add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                    add_new_ph4(coordinates_, i, 'NZ', coordinates_plus)

       
            else:
                coordinates_[i].append('O')
                # can interact with protein OH (-->OG), +NH (-->OD1)
                if duplicate_:
                    add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                    add_new_ph4(coordinates_, i, 'OD1', coordinates_plus)

        
        elif atom_types_[i] == 'N.am':
            if any([atom in HYDROGENS for atom in connectivity_[i][0]]):
                coordinates_[i].append('N')
                # can interact with protein OH (-->OG), O- (-->NZ)
                if duplicate_:
                    add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                    add_new_ph4(coordinates_, i, 'NZ', coordinates_plus)
            else:
                coordinates_[i].append(None)
        

        elif atom_types_[i] == 'N.pl3':
            if any([atom in HYDROGENS for atom in connectivity_[i][0]]):
                coordinates_[i].append('N')
                # can interact with protein OH (-->OG), O- (-->NZ)
                if duplicate_:
                    add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                    add_new_ph4(coordinates_, i, 'NZ', coordinates_plus)
            else:
                coordinates_[i].append(None)
       
                
        elif atom_types_[i] == 'N.2':
            if any([atom in HYDROGENS for atom in connectivity_[i][0]]):
                coordinates_[i].append('N') # e.g. imine
                # can interact with protein OH (-->OG), O- (-->NZ)
                if duplicate_:
                    add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                    add_new_ph4(coordinates_, i, 'NZ', coordinates_plus)
            else:
                coordinates_[i].append('O')
                # can interact with protein OH (-->OG), +NH (-->OD1)
                if duplicate_:
                    add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                    add_new_ph4(coordinates_, i, 'OD1', coordinates_plus)
            

        elif atom_types_[i] == 'N.3':
            if any([atom in HYDROGENS for atom in connectivity_[i][0]]):
                coordinates_[i].append('N')
                # can interact with protein OH (-->OG), O- (-->NZ)
                if duplicate_:
                    add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                    add_new_ph4(coordinates_, i, 'NZ', coordinates_plus)
            else:
                # N sp3 likely to be protonated, cannot be acceptor
                coordinates_[i].append(None)
           
                
        elif atom_types_[i] == 'O.3':
            if any([atom in HYDROGENS for atom in connectivity_[i][0]]):
                coordinates_[i].append('OG')
                # can interact with protein O (-->N),
                # O- (-->NZ), NH (-->O), +NH (-->OD1)
                if duplicate_:
                    add_new_ph4(coordinates_, i, 'N', coordinates_plus)
                    add_new_ph4(coordinates_, i, 'O', coordinates_plus)
                    add_new_ph4(coordinates_, i, 'NZ', coordinates_plus)
                    add_new_ph4(coordinates_, i, 'OD1', coordinates_plus)
            else:
                coordinates_[i].append('O')
                # can interact with protein OH (-->OG), +NH (-->OD1)
                if duplicate_:
                    add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                    add_new_ph4(coordinates_, i, 'OD1', coordinates_plus)

        
        elif atom_types_[i] == 'O.2':
            coordinates_[i].append('O')
            # can interact with protein OH (-->OG), +NH (-->OD1)
            if duplicate_:
                add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                add_new_ph4(coordinates_, i, 'OD1', coordinates_plus)
            
    
        elif atom_types_[i] == 'N.1':
            coordinates_[i].append('O')
            # can interact with protein OH (-->OG), +NH (-->OD1)
            if duplicate_:
                add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                add_new_ph4(coordinates_, i, 'OD1', coordinates_plus)
            
        
        # ion rules
        elif atom_types_[i] == 'O.co2':
            coordinates_[i].append('OD1')
            # can interact with protein OH (-->OG), NH (-->O)
            if duplicate_:
                add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                add_new_ph4(coordinates_, i, 'O', coordinates_plus)


        elif atom_types_[i] == 'N.4':
            coordinates_[i].append('NZ')
            # can interact with protein OH (-->OG). if hydrogen, C=O (-->N)
            if duplicate_:
                add_new_ph4(coordinates_, i, 'OG', coordinates_plus)
                if any([atom in HYDROGENS for atom in connectivity_[i][0]]):
                    add_new_ph4(coordinates_, i, 'N', coordinates_plus)

        
        # These in ligands ???? --> untreated
        elif atom_types_[i] == 'Mg':
            coordinates_[i].append('NZ')
                
        elif atom_types_[i] == 'Mn':
            coordinates_[i].append('NZ')
                
        elif atom_types_[i] == 'Zn':
            coordinates_[i].append('NZ')
                
        elif atom_types_[i] == 'Ca':
            coordinates_[i].append('NZ')
                
        elif atom_types_[i] == 'Fe':
            coordinates_[i].append('NZ')
                
        elif atom_types_[i] == 'Cu':
            coordinates_[i].append('NZ')
        


        # rules for H atoms
        elif atom_types_[i] == 'H':
            coordinates_[i].append(None)


        elif atom_types_[i] == 'Du':
            coordinates_[i].append(None)
        
        # rules for dummy/other atoms
        else:
            print('Untreated atom type -->', atom_types_[i], coordinates_[i])
            coordinates_[i].append(None)


        # missing
        # P.3
        # C.cat, etc.
        # Co.oh
        # H.spc
        # S.o
        # S.o2
        # O.spc
        # H.t3p
        # O.t3p
        # Cr.oh
        # Cr.th
        # Se
        # Mo
        # Si
        # li, Al, K, Na


    if duplicate_:
        print("Coordinates were duplicated with new ph4:")
        for idx, coords in coordinates_plus.items():
            print(idx, coords)
    #else:
    #    print("Coordinates were not duplicated")
    #print(coordinates_plus)

    return coordinates_plus






# bond list
# bond type

def find_aromatic_rings(bonds_, atom_types_):
    # aromacity by mol2 atom typing
    ar_cycles = []
    aromatic_bonds = []
    for key, value in bonds_.items():
        if value == 'ar':
            aromatic_bonds.append(key)
    if len(aromatic_bonds) > 2:
        graph = nx.Graph()
        graph.add_edges_from(aromatic_bonds)
        ar_cycles += nx.cycle_basis(graph)
    
    #print("ar_cycles", ar_cycles)
    sorted_ar_cycles = []
    for ac in ar_cycles:
        sorted_ar_cycles.append(sorted(ac))
    #print("ar_cycles sorted", sorted_ar_cycles)




    ### this implementation might be incorrect in a few cases, eg. P in cycles ###

    # when mol2 typing is not informative: no ar typing
    # aromacity rules
    # 1) cyclic
    # 2) conjugated
    # 3) Huckel rule, 4n+2 pi electrons
    # 4) planar
    
    bonds_atoms = []
    for atom_1, atom_2 in bonds_:
        bonds_atoms.append(atom_1)
        bonds_atoms.append(atom_2)
    
    all_cycles = find_rings(list(bonds_.keys()))
    sorted_all_cycles = []
    for c in all_cycles:
        sorted_all_cycles.append(sorted(c))

    #print("cycles all sorted", sorted_all_cycles)
    cycles = [c for c in sorted_all_cycles if c not in sorted_ar_cycles]
    #print("cycles", cycles)

    for c in cycles:

        pi_electrons = []

        for atom_id in c:
            hybridization = atom_types_[atom_id].split('.')[-1]
            #print("hybridization", hybridization)
            if hybridization == '2' or hybridization == '1' or hybridization == 'ar':
                pi_electrons.append(1)

            # if not sp2, sp1, does it have lone pairs ?
            # or is it a tricycle with carbocation ?
            else:
                lone_electrons = find_lone_electrons(bonds_atoms, atom_id,
                                            atom_types_[atom_id].split('.')[0])
                pi_electrons.append(lone_electrons)
        #print("pi_electrons", pi_electrons)

        if len(c) == 3 and pi_electrons.count(1) == 2 and pi_electrons.count(-1) == 1 :
            # case of three carbon ring with one carbocation: aromatic
            ar_cycles.append(c)

        elif 0 in pi_electrons:
            # not aromatic
            continue
        
        elif -1 in pi_electrons: # cation in cycle
            continue

        elif pi_electrons.count(2) > 1: # more than two lone pairs
            continue

        elif len(c) > 3 and (sum(pi_electrons)-2) % 4 == 0:
            if all_multi_bond_inside_cycle(c, bonds_, sorted_all_cycles, sorted_ar_cycles):
                #print("!!!!!!!!!!!!!! multi_bond_inside_cycle")
                ar_cycles.append(c)


    #print("ar_cycles", ar_cycles)
    return ar_cycles





def find_lone_electrons(bonds_atoms_, atom_, symbol_):
    ground_state_bonds = {'N': 3,
                          'C': 4, 
                          'O': 2, 
                          'S': 2,
                          } 
    heteroatoms_with_lone = ['N', 'O', 'S']

    # hybridization are only 3, pl3 etc. sp2/sp1/ar are excluded
    # in a cycle, an atom with double bonds will be sp2
    
    if symbol_ in ground_state_bonds:
        missing_bonds = (ground_state_bonds[symbol_] - bonds_atoms_.count(atom_))
        
        if missing_bonds < 0:
            lone_electrons = -1 # missing 
        
        elif missing_bonds > 0:
            lone_electrons = 2
        
        else:
            lone_electrons = 0
            if symbol_ in heteroatoms_with_lone:
                lone_electrons = 2
    else:
        lone_electrons = -1 # continue
    return lone_electrons





def all_multi_bond_inside_cycle(cycle_, bonds_, all_cycles_, ar_cycles_):
    conjugation = True
    for atom in cycle_:
        for atom_pair in bonds_:
            if atom == atom_pair[0] and atom_pair[1] not in cycle_ and\
                                    bonds_[atom_pair] in ['2', '3', 'ar']: # not '1', 'nc', 'un', 'du'
                if atom_pair[1] in [a for c in ar_cycles_ for a in c]:
                    pass
                elif is_fused_to_aromatic_ring(atom, bonds_, cycle_, all_cycles_): # belong to another ring that is aromatic
                    pass
                else:
                    conjugation = False
                    break
            elif atom == atom_pair[1] and atom_pair[0] not in cycle_ and\
                                    bonds_[atom_pair] in ['2', '3', 'ar']: # not '1', 'nc', 'un', 'du'
                if atom_pair[0] in [a for c in ar_cycles_ for a in c]:
                    pass
                elif is_fused_to_aromatic_ring(atom, bonds_, cycle_, all_cycles_): # belong to another ring that is aromatic
                    pass
                else:
                    conjugation = False
                    break
            
    return conjugation

            



def is_fused_to_aromatic_ring(atom_, bonds_, cycle_, all_cycles_):
    fused_to_aromatic = False
    for c in all_cycles_:
        bond_types_in_cycle = []
        if c != cycle_ and atom_ in c:
            for i in range(len(c)-1):
                if (c[i], c[i+1]) in bonds_:
                    bond_types_in_cycle.append(bonds_[(c[i], c[i+1])])
                elif (c[i+1], c[i]) in bonds_:
                    bond_types_in_cycle.append(bonds_[(c[i+1], c[i])])
            if len(c)%2 == 0 and (bond_types_in_cycle.count('2') + bond_types_in_cycle.count('3') +
                    bond_types_in_cycle.count('pl3') + bond_types_in_cycle.count('ar')) == len(c)/2:
                fused_to_aromatic = True
            elif len(c)%2 == 0 and (bond_types_in_cycle.count('2') + bond_types_in_cycle.count('3') +
                    bond_types_in_cycle.count('pl3') + bond_types_in_cycle.count('ar')) >= (len(c)-1)/2:
                fused_to_aromatic = True

    return fused_to_aromatic





def find_rings(bonds_list_):
    if len(bonds_list_) > 2:
        graph = nx.Graph()
        graph.add_edges_from(bonds_list_)
        cycles = nx.cycle_basis(graph)
    else:
        cycles = []
    return cycles




def add_new_ph4_ring(coordinates_, ring_coordinates_, ring_index_, new_ph4_,
                        coordinates_plus_):
    new_index = len(coordinates_) + len(coordinates_plus_) + 1
    coords = copy.deepcopy(ring_coordinates_[ring_index_][:3])
    coordinates_plus_[new_index] = coords
    coordinates_plus_[new_index].append(new_ph4_)




def process_lig_ph4(coordinates_ph4_, aromatic_rings_, duplicate_=False):
    new_coordinates = {}
    processed = []
    for r in aromatic_rings_:
        #print(r)
        #print(coordinates_ph4_)
        ring_ph4 = [coordinates_ph4_[i][3] for i in r]
        new_ph4 = 'CZ'
        new_x = round(np.mean([coordinates_ph4_[i][0] for i in r]), 4)
        new_y = round(np.mean([coordinates_ph4_[i][1] for i in r]), 4)
        new_z = round(np.mean([coordinates_ph4_[i][2] for i in r]), 4)
        new_coordinates[tuple(r)] = [new_x, new_y, new_z, new_ph4]
        processed += r


    if duplicate_:
        coordinates_plus = {}
        for r in new_coordinates:
            # also CA
            add_new_ph4_ring(coordinates_ph4_, new_coordinates, r, 'CA', coordinates_plus)
        for idx in coordinates_plus:
            new_coordinates[idx] = coordinates_plus[idx]

        #print("Coordinates in ring were duplicated with new ph4:")
        #for idx, coords in coordinates_plus.items():
        #    print(idx, coords)


    #print(coordinates_ph4_)
    for idx in coordinates_ph4_:
        if idx in processed and coordinates_ph4_[idx][3] == 'CZ':
            continue
        if coordinates_ph4_[idx][3] is None:
            continue
        new_coordinates[idx] = coordinates_ph4_[idx]

    return new_coordinates




def write_mol2(ofile_, coordinates_ph4_, macromol_="SMALL"):

    ATOM_TYPE = {"OG":"O.3",
                  "N":"N.am",
                  "O":"O.2",
                  "NZ":"N.4",
                  "CZ":"C.ar",
                  "CA":"C.3",
                  "DU":"H",
                  "OD1":"O.co2",
                  "NOP":"H",
                  "H":"H"}

    RESIDUE = {"OG":"SER",
                "N":"ALA",
                "O":"ALA",
                "NZ":"LYS",
                "CZ":"PHE",
                "CA":"GLY",
                "DU":"CUB",
                "OD1":"ASP",
                "NOP":"CUB",
                "H":"CUB"}


    if ofile_[-5:] != '.mol2':
        ofile_ += '.mol2'
    name = os.path.basename(ofile_)
    name = os.path.splitext(name)[0]
    

    of_string = ""
    of_string += "# Modified by LigCare\n"
    of_string += "# Modification time: {}\n".format(
                            strftime("%a %d %b %Y %H:%M:%S", localtime()))
    of_string += "# Name: {}.mol2\n\n".format(name)

    of_string += "@<TRIPOS>MOLECULE\n"
    of_string += "{}\n".format(name)
    of_string += "{:>5}{:>6}{:>6}{:>6}{:>6}\n".format(
                                            len(coordinates_ph4_), 0, 0, 0, 0)
    of_string += "{}\n".format(macromol_)
    of_string += "NO_CHARGES\n"
    of_string += "@<TRIPOS>ATOM"




    for i, point in enumerate(coordinates_ph4_.values()):
        #print(point)
        x, y, z, ph4 = point
        #if ph4 == "NOP":
        #    ph4 = "DU"
        if ph4 is None:
            continue
        of_string += ("\n{:>7} {:<8} {:>9.4f} {:>9.4f} {:>9.4f} "
                      "{:<5} {:>5} {:<8} {:>9}".format(i+1,
                                                 ph4,
                                                 x,
                                                 y,
                                                 z,
                                                 ATOM_TYPE[ph4],
                                                 i+1,
                                                 RESIDUE[ph4]+str(i+1),
                                                 0.0000
                                                ))
    of_string += "\n@<TRIPOS>BOND"
    
    with open(ofile_, 'w') as of:
        of.write(of_string)
    print("written mol2 to {}".format(ofile_))
    return ofile_





def xtract_ligands(file_):
    ligands = []
    with open(file_, 'r') as f:
        ligands = f.read().split('\n')
    ligands.remove('')
    return ligands



def main(input_, output_, duplicate_):

    coordinates, connectivity, atom_types, bond_types = lig_mol2_info(input_)
    coordinates_plus = lig_ph4(coordinates, connectivity, atom_types, duplicate_)
    #print("coordinates_plus", coordinates_plus)
    coordinates = {**coordinates, **coordinates_plus}
    #print(coordinates)
    #for c, v in coordinates.items():
    #    print(c, v)

    aromatic_rings = find_aromatic_rings(bond_types, atom_types)
    #print(aromatic_rings)
    coordinates = process_lig_ph4(coordinates, aromatic_rings, duplicate_)
    #for c, v in coordinates.items():
    #    print(c, v)

    write_mol2(output_, coordinates, macromol_="SMALL")





if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument('--single', required=False, action='store_true',
                help='True or False, for single input, used with -i, -o options')

    parser.add_argument('-i', '--input', type=str, required=False,
                help='ligand mol2')

    parser.add_argument('-o', '--output', type=str, required=False,
                help='output file for ligand ph4 mol2')


    parser.add_argument('--database', required=False, action='store_true',
                help=('True or False, for list of ligands input,'
                        ' used with -l, -dir, -odir options'))

    parser.add_argument('-l', '--list', type=str, required=False,
                help='list of ligand mol2')

    parser.add_argument('-idir', '--inputdir', type=str, required=False,
                help='output directory for ligand ph4 mol2')

    parser.add_argument('-odir', '--outputdir', type=str, required=False,
                help='output directory for ligand ph4 mol2')


    parser.add_argument('--duplicate', required=False, action='store_true',                        
                help='duplicate ph4 by compatibility rules in mol2')


    args = parser.parse_args()



    duplicate = False
    if args.duplicate:
        duplicate = True



    if args.single:
        if args.input is None:
            print('required -i option: input ligand')
            sys.exit(1)

        if args.output is None:
            ofile = os.path.splitext(os.path.basename(args.input))[0] + '_ligph4.mol2'
        else:
            ofile = args.output

        main(args.input, ofile, duplicate)

    

    elif args.database:
        if args.list is None:
            print('required -l option: list of input ligands')
            sys.exit(1)

        if args.inputdir is None or args.outputdir is None:
            print('required -l option: list of input ligands')
            print('required -idir, -odir options: input and output directories')
            sys.exit(1)

     
        ligands = xtract_ligands(args.list)
        skipped = 0
        for i, lig in enumerate(ligands):
            print('processing.....................', lig, f'{(i+1)}/{len(ligands)}')
            ofile = os.path.splitext(os.path.basename(lig))[0]
            ofile = ofile.replace('_ligand', '') + '_ligph4.mol2'
            try:
                main(f'{args.inputdir}/{lig}', f'{args.outputdir}/{ofile}', duplicate)
            except IOError:
                print('incorrect path:', f'{args.inputdir}/{lig}', f'{args.outputdir}/{ofile}')
                skipped += 1
            except:
                print('skiping ligand................. ', f'{args.inputdir}/{lig}')
                with open('ligph4_skipped.txt', 'a') as sf:
                    sf.write(f'{lig}\n')
                skipped += 1

        if skipped > 0:
            print('############################################')
            print('skipped entries are written to: ligph4_skipped.txt')


    else:
        print('required --single or --database options')
        sys.exit(1)






