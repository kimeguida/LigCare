

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



import numpy as np
import argparse
import os
from time import strftime, localtime




def read_cav_mol2(mol2_):
    with open(mol2_, 'r') as f:
        mol2 = f.read()
    atom_area = mol2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0]
    atoms_tmp = atom_area.split('\n')
    del atoms_tmp[-1]

    coords = []
    ph4_features = []
    for atm in atoms_tmp:
        if atm == '':
            continue
        cols = atm.split()
        x = (round(float(cols[2]), 4))
        y = (round(float(cols[3]), 4))
        z = (round(float(cols[4]), 4))
        coords.append(np.array([x, y, z]))
        ph4_features.append(cols[1])

    return np.array(coords), ph4_features




def read_prot_mol2(mol2_):
    with open(mol2_, 'r') as f:
        mol2 = f.read()
    atom_area = mol2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0]
    atoms_tmp = atom_area.split('\n')
    del atoms_tmp[-1]

    coords = []
    for atm in atoms_tmp:
        if atm == '':
            continue
        cols = atm.split()
        atm = cols[1]
        if atm.startswith('H'):
            continue
        x = (round(float(cols[2]), 4))
        y = (round(float(cols[3]), 4))
        z = (round(float(cols[4]), 4))
        coords.append(np.array([x, y, z]))
    return np.array(coords)





def read_prot_mol2_features(mol2_):
    
    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'HID',
    'HIE', 'HIP', 'ILE', 'LEU', 'LYS', 'LYN', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
    'TYR', 'VAL', 'HOH', 'MG', 'FE', 'MN', 'ZN', 'CO', 'CA', 'NA']

    with open(mol2_, 'r') as f:
        mol2 = f.read()
    atom_area = mol2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0]
    atoms_tmp = atom_area.split('\n')
    del atoms_tmp[-1]

    coords = []
    res_coords_dict = {}
    features = []
    for atm in atoms_tmp:
        if atm == '':
            continue
        cols = atm.split()
        atm_name = cols[1]
        res = cols[7][:3]
        res_num = cols[7][3:]

        if res not in amino_acids:
            res = res[:2]
            res_num = cols[7][2:]
        if res not in amino_acids:
            continue

        try:
            res_num = int(res_num)
        except ValueError:
            continue

        x = (round(float(cols[2]), 4))
        y = (round(float(cols[3]), 4))
        z = (round(float(cols[4]), 4))
        
        res_coords_dict[(res, res_num, atm_name)] = np.array([x, y, z])

        if not atm_name.startswith('H'):
            coords.append(np.array([x, y, z]))
            features.append((res, res_num, atm_name))


    
    return np.array(coords), features, res_coords_dict





def read_mol2_features(mol2_):

    with open(mol2_, 'r') as f:
        mol2 = f.read()
    atom_area = mol2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0]
    atoms_tmp = atom_area.split('\n')
    del atoms_tmp[-1]

    coords = []
    features = []
    for atm in atoms_tmp:
        if atm == '':
            continue
        cols = atm.split()
        atm_symbol = cols[5].split('.')[0]

        x = (round(float(cols[2]), 4))
        y = (round(float(cols[3]), 4))
        z = (round(float(cols[4]), 4))

        coords.append(np.array([x, y, z]))
        features.append(atm_symbol)


    return np.array(coords), features




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
    name = os.path.splitext(ofile_)[0]

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




    for i, point in enumerate(coordinates_ph4_):
        if len(point) == 3:
            x, y, z = point
            ph4 = 'DU'
        elif len(point) == 4:
            x, y, z, ph4 = point
        else:
            print('ERROR: cannot read coordinates', point)
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
    #print("written mol2: {}".format(ofile_))
    return ofile_





def write_pruned_mol2(indexes_, mol2_, ofile_):
    coords, ph4_features = read_cav_mol2(mol2_)
    pruned_coordinates_ph4 = []
    for i in indexes_:
        pruned_coordinates_ph4.append([coords[i][0], coords[i][1], coords[i][2],
                            ph4_features[i]])

    write_mol2(ofile_, pruned_coordinates_ph4)






def read_ints_mol2(mol2_):
    """
    IChem interactions mol2
    """
    with open(mol2_, 'r') as f:
        mol2 = f.read()
    atom_area = mol2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0]
    atoms_tmp = atom_area.split('\n')
    del atoms_tmp[-1]

    coords = []
    ints_features = []
    for atm in atoms_tmp:
        if atm == '':
            continue
        cols = atm.split()
        res = cols[7]
        if res[2] != 'L':
            continue
        x = (round(float(cols[2]), 4))
        y = (round(float(cols[3]), 4))
        z = (round(float(cols[4]), 4))
        coords.append(np.array([x, y, z]))
        ints_features.append(cols[1])

    return np.array(coords), ints_features

