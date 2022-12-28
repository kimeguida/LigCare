



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



# beta version
# some cases need to be implemented



import os
import numpy as np
import argparse
from time import strftime, localtime
from sklearn.neighbors import NearestNeighbors



parser = argparse.ArgumentParser()
parser.add_argument('-p', '--protein', type=str, required=True,
            help='protein mol2')
parser.add_argument('-c', '--cavity', type=str, required=True,
            help='cavity mol2')
parser.add_argument('-o', '--output', type=str, required=True,
            help='output mol2')
args = parser.parse_args()


METALS = ['ZN', 'MG', 'CO', 'CA', 'FE']

def xtract_protein_coordinates(mol2_):
    with open(mol2_, 'r') as f:
        mol2 = f.read()
    mol2_atom_block = mol2.split("@<TRIPOS>ATOM\n")[1].split("@<TRIPOS>")[0]
    mol2_atom_block = mol2_atom_block.split("\n")
    del mol2_atom_block[-1]
    coordinates = []
    coordinates_h = {}
    atom_infos = []
    for line in mol2_atom_block:
        cols = line.split()
        idx = int(cols[0])
        atom = cols[1]
        if atom == 'HNCO' or atom == 'H' or atom == 'HN(CA)CO' or \
        atom == 'HNCACB' or atom == 'CBCA(CO)NH' or atom == 'HNCA' or \
        atom == 'CA(CO)NH' or atom == 'HNHA' or atom == 'HA(CO)NH' or \
        atom == 'HNcaCO' or atom == 'HNcoCA' or atom == 'HNcoCACB' or \
        atom == 'HNcoHA':
            atom = 'H'
        x = float(cols[2])
        y = float(cols[3])
        z = float(cols[4])
        res_info = cols[7]
        #print(res_info)
        if res_info[:2] in METALS:
            res_name = res_info[:2]
            res_num = int(res_info[2:])
        else:
            res_name = res_info[:3]
            if '_' not in res_info[3:]:
                res_num = int(res_info[3:])
            else:
                res_num = int(res_info[3:].replace('_', ''))*-1
        mol2_res_num = int(cols[6])
        if atom[0] == 'H':
            coordinates_h[(atom, mol2_res_num, res_name, res_num)] = np.array([x, y, z])
        else:
            coordinates.append(np.array([x, y, z]))
            atom_infos.append((atom, mol2_res_num, res_name, res_num))
            coordinates_h[(atom, mol2_res_num, res_name, res_num)] = np.array([x, y, z])
    return coordinates, atom_infos, coordinates_h


def xtract_cavity_coordinates(mol2_):
    with open(mol2_, 'r') as f:
        mol2 = f.read()
    mol2_atom_block = mol2.split("@<TRIPOS>ATOM\n")[1].split("@<TRIPOS>")[0]
    mol2_atom_block = mol2_atom_block.split("\n")
    del mol2_atom_block[-1]
    coordinates = []
    for line in mol2_atom_block:
        cols = line.split()
        x = float(cols[2])
        y = float(cols[3])
        z = float(cols[4])
        coordinates.append(np.array([x, y, z]))
    return coordinates



def get_binding_site_coordinates(coords_prot_, coords_cav_):
    #print(len(coords_cav_))
    #print(len(coords_prot_))
    neigh = NearestNeighbors(n_neighbors=len(coords_prot_), radius=3.5, algorithm='ball_tree').fit(coords_prot_)
    distances, indices = neigh.radius_neighbors(coords_cav_)
    #print(len(indices))

    bs_atoms = []
    for i in range(len(coords_cav_)):
        for j in indices[i]:
            bs_atoms.append(j)
    bs_atoms = set(bs_atoms)
    #print(len(bs_atoms))
    return bs_atoms



def towards_binding_site(root_, exit_, bs_center_):
    vec_1 = np.array(exit_)-np.array(root_)
    vec_2 = np.array(bs_center_)-np.array(root_)
    cos_angle = np.dot(vec_1, vec_2)/(np.linalg.norm(vec_1)*np.linalg.norm(vec_2))
    if 0 <= cos_angle <= 1:
        in_bs = True
    else:
        in_bs = False
    return in_bs



def atom_position(root_, exit_, dist_):
    vec_exit = np.array(exit_)-np.array(root_)
    k = dist_/np.linalg.norm(vec_exit)
    x = k * vec_exit[0] + root_[0]
    y = k * vec_exit[1] + root_[1]
    z = k * vec_exit[2] + root_[2]
    #return np.array([x, y, z])
    return (x, y, z)




def ring_atom_position(ring_, dist_):
    center_x = np.mean([p[0] for p in ring_])
    center_y = np.mean([p[1] for p in ring_])
    center_z = np.mean([p[2] for p in ring_])
    center_ring = np.array([center_x, center_y, center_z])
    #print(center_ring)
    vec_plan_1 = np.array(ring_[0])-np.array(center_ring)
    vec_plan_2 = np.array(ring_[1])-np.array(center_ring)
    vec_exit_1 = np.cross(vec_plan_1, vec_plan_2)
    vec_exit_2 = -1 * vec_exit_1
    #print('vec_exit_1', vec_exit_1)
    #print('vec_exit_2', vec_exit_2)
    projection_from_center = []
    for vec in [vec_exit_1, vec_exit_2]:
        k = dist_/np.linalg.norm(vec)
        x = k * vec[0] + center_ring[0]
        y = k * vec[1] + center_ring[1]
        z = k * vec[2] + center_ring[2]
        projection_from_center.append((x, y, z))
    return center_ring, projection_from_center



def prot_ph4(bs_atoms_, coords_, atom_infos_, coords_h_):
    center_x = np.mean([coords_[idx][0] for idx in bs_atoms_])
    center_y = np.mean([coords_[idx][1] for idx in bs_atoms_])
    center_z = np.mean([coords_[idx][2] for idx in bs_atoms_])
    center = np.array([center_x, center_y, center_z])
    coordinates_ph4 = []
    for idx in bs_atoms:
        atom, mol2_res_num, res_name, res_num = atom_infos_[idx]
        if atom == 'O' and res_name != 'HOH':
            root_atom = coords_h_[('C', mol2_res_num, res_name, res_num)]
            exit_atom = coords_[idx]
            if towards_binding_site(root_atom, exit_atom, center):
                coord = atom_position(root_atom, exit_atom, 3.5)
                ph4 = 'N'
                coordinates_ph4.append((coord, ph4))
        elif atom == 'N':
            root_atom = coords_[idx]
            exit_atom = coords_h_[('H', mol2_res_num, res_name, res_num)]
            if towards_binding_site(root_atom, exit_atom, center):
                coord = atom_position(root_atom, exit_atom, 3.5)
                ph4 = 'O'
                coordinates_ph4.append((coord, ph4))
        elif atom == 'CA':
            continue
        elif atom == 'C':
            continue
        elif atom == 'CB':
            continue
        elif res_name == 'GLY':
            continue
        elif res_name == 'ALA':
            continue

        elif res_name == 'VAL':
            if atom == 'CG1' or atom == 'CG2':
                exit_atom = coords_[idx]
                root_atom = coords_h_[('CB', mol2_res_num, res_name, res_num)]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'CA'
                    coordinates_ph4.append((coord, ph4))

        elif res_name == 'LEU':
            if atom == 'CG':
                exit_atom = coords_[idx]
                root_atom = coords_h_[('CB', mol2_res_num, res_name, res_num)]
            elif atom == 'CD1' or atom == 'CD2':
                exit_atom = coords_[idx]
                root_atom = coords_h_[('CG', mol2_res_num, res_name, res_num)]
            if towards_binding_site(root_atom, exit_atom, center):
                coord = atom_position(root_atom, exit_atom, 3.5)
                ph4 = 'CA'
                coordinates_ph4.append((coord, ph4))

        elif res_name == 'ILE':
            if atom == 'CG1' or atom == 'CG2':
                exit_atom = coords_[idx]
                root_atom = coords_h_[('CB', mol2_res_num, res_name, res_num)]
            elif atom == 'CD1':
                root_atom = coords_h_[('CG1', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
            if towards_binding_site(root_atom, exit_atom, center):
                coord = atom_position(root_atom, exit_atom, 3.5)
                ph4 = 'CA'
                coordinates_ph4.append((coord, ph4))

        elif res_name == 'PRO':
            if atom == 'CG':
                root_atom = coords_h_[('CB', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'CA'
                    coordinates_ph4.append((coord, ph4))
                root_atom = coords_h_[('CD', mol2_res_num, res_name, res_num)]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    coordinates_ph4.append((coord, ph4))

        elif res_name == 'MET':
             if atom == 'SD':
                root_atom = coords_h_[('CG', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'CA'
                    coordinates_ph4.append((coord, ph4))
                root_atom = coords_h_[('CE', mol2_res_num, res_name, res_num)]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    coordinates_ph4.append((coord, ph4))


        elif res_name == 'CYS':
            if atom == 'SG':
                root_atom = coords_h_[('CB', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'CA'
                    coordinates_ph4.append((coord, ph4))

        elif res_name == 'SER':
            if atom == 'OG':
                root_atom = coords_[idx]
                exit_atom = coords_h_[('HG', mol2_res_num, res_name, res_num)]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'OG'
                    coordinates_ph4.append((coord, ph4))

        elif res_name == 'THR':
            if atom == 'OG1':
                root_atom = coords_[idx]
                exit_atom = coords_h_[('HG1', mol2_res_num, res_name, res_num)]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'OG'
                    coordinates_ph4.append((coord, ph4))
            elif atom == 'CG2':
                root_atom = coords_h_[('CB', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'CA'
                    coordinates_ph4.append((coord, ph4))

        elif res_name == 'ASN':
            if atom == 'OD1':
                root_atom = coords_h_[('CG', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'N'
                    coordinates_ph4.append((coord, ph4))
            elif atom == 'ND2':
                root_atom = coords_[idx]
                exit_atom = coords_h_[('HD21', mol2_res_num, res_name, res_num)]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'O'
                    coordinates_ph4.append((coord, ph4))
                exit_atom = coords_h_[('HD22', mol2_res_num, res_name, res_num)]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    coordinates_ph4.append((coord, ph4))

        elif res_name == 'GLN':
            if atom == 'OE1':
                root_atom = coords_h_[('CD', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'N'
                    coordinates_ph4.append((coord, ph4))
            elif atom == 'NE2':
                root_atom = coords_[idx]
                exit_atom = coords_h_[('HE21', mol2_res_num, res_name, res_num)]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'O'
                    coordinates_ph4.append((coord, ph4))
                exit_atom = coords_h_[('HE22', mol2_res_num, res_name, res_num)]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'O'
                    coordinates_ph4.append((coord, ph4))
            elif atom == 'CG':
                root_atom = coords_h_[('CB', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    ph4 = 'CA'
                    coordinates_ph4.append((coord, ph4))

        if res_name == 'ASP':
            if atom == 'OD1' or atom == 'OD2':
                root_atom = coords_h_[('CG', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                ph4 = 'NZ'
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    coordinates_ph4.append((coord, ph4))
                exit_atom = center
                coord = atom_position(root_atom, exit_atom, 3.5)
                coordinates_ph4.append((coord, ph4))

        if res_name == 'GLU':
            if atom == 'OE1' or atom == 'OE2':
                root_atom = coords_h_[('CD', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                ph4 = 'NZ'
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    coordinates_ph4.append((coord, ph4))
                exit_atom = center
                coord = atom_position(root_atom, exit_atom, 4)
                coordinates_ph4.append((coord, ph4))

            elif atom == 'CG':
                root_atom = coords_h_[('CB', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 4)
                    ph4 = 'CA'
                    coordinates_ph4.append((coord, ph4))

        elif res_name == 'PHE':
            if atom == 'CD1' or atom == 'CD2' or \
                atom == 'CE1' or atom == 'CE2' or atom == 'CZ':
                ring = [coords_h_[('CG', mol2_res_num, res_name, res_num)],
                        coords_h_[('CD1', mol2_res_num, res_name, res_num)],
                        coords_h_[('CD2', mol2_res_num, res_name, res_num)],
                        coords_h_[('CE1', mol2_res_num, res_name, res_num)],
                        coords_h_[('CE2', mol2_res_num, res_name, res_num)],
                        coords_h_[('CZ', mol2_res_num, res_name, res_num)],
                        ]
                center_ring, projection_from_center = ring_atom_position(ring, 3.5)
                ph4 = 'CZ'
                if projection_from_center != []:
                    for coord in projection_from_center:
                        if towards_binding_site(center_ring, coord, center):
                            coordinates_ph4.append((coord, ph4))
                            
                root_atom = center_ring
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 4)
                    coordinates_ph4.append((coord, ph4))
                    
                    
        elif res_name == 'TYR':
            if atom == 'CD1' or atom == 'CD2' or \
                atom == 'CE1' or atom == 'CE2' or atom == 'CZ':
                ring = [coords_h_[('CG', mol2_res_num, res_name, res_num)],
                        coords_h_[('CD1', mol2_res_num, res_name, res_num)],
                        coords_h_[('CD2', mol2_res_num, res_name, res_num)],
                        coords_h_[('CE1', mol2_res_num, res_name, res_num)],
                        coords_h_[('CE2', mol2_res_num, res_name, res_num)],
                        coords_h_[('CZ', mol2_res_num, res_name, res_num)],
                        ]
                center_ring, projection_from_center = ring_atom_position(ring, 3.5)
                ph4 = 'CZ'
                if projection_from_center != []:
                    for coord in projection_from_center:
                        if towards_binding_site(center_ring, coord, center):
                            coordinates_ph4.append((coord, ph4))
                            
                root_atom = center_ring
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 4)
                    coordinates_ph4.append((coord, ph4))

            if atom == 'OH':
                root_atom = coords_[idx]
                exit_atom = coords_h_[('HH', mol2_res_num, res_name, res_num)]
                ph4 = 'OG'
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    coordinates_ph4.append((coord, ph4))
                
                    
        elif res_name == 'TRP':
            if atom == 'CD1' or atom == 'CD2' or \
                atom == 'CE2' or atom == 'CE3' or \
                atom == 'NE1' or atom == 'CZ2' or\
                atom == 'CZ3' or atom == 'CH2':
                ring = [coords_h_[('CD2', mol2_res_num, res_name, res_num)],                     
                        coords_h_[('CE2', mol2_res_num, res_name, res_num)],
                        coords_h_[('CE3', mol2_res_num, res_name, res_num)],
                        coords_h_[('CZ2', mol2_res_num, res_name, res_num)],
                        coords_h_[('CZ3', mol2_res_num, res_name, res_num)],
                        coords_h_[('CH2', mol2_res_num, res_name, res_num)],
                        ]       
                center_ring, projection_from_center = ring_atom_position(ring, 3.5)
                ph4 = 'CZ'
                
                if projection_from_center != []:
                    for coord in projection_from_center:
                        if towards_binding_site(center_ring, coord, center):
                            coordinates_ph4.append((coord, ph4))
                    
                root_atom = center_ring
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 4)
                    coordinates_ph4.append((coord, ph4))
                    
                    
                ring = [coords_h_[('CG', mol2_res_num, res_name, res_num)],                     
                        coords_h_[('CD1', mol2_res_num, res_name, res_num)],
                        coords_h_[('CD2', mol2_res_num, res_name, res_num)],
                        coords_h_[('NE1', mol2_res_num, res_name, res_num)],
                        coords_h_[('CE2', mol2_res_num, res_name, res_num)],
                        ]       
                center_ring, projection_from_center = ring_atom_position(ring, 3.5)             
                if projection_from_center != []:
                    for coord in projection_from_center:
                        if towards_binding_site(center_ring, coord, center):
                            coordinates_ph4.append((coord, ph4))
                    
                root_atom = center_ring
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 4)
                    coordinates_ph4.append((coord, ph4))
                    
                

        elif res_name == 'HIS':
            if atom == 'ND1' or atom == 'CD2' or atom == 'CE1' or atom == 'NE2':
                ring = [coords_h_[('CG', mol2_res_num, res_name, res_num)],
                        coords_h_[('ND1', mol2_res_num, res_name, res_num)],
                        coords_h_[('CD2', mol2_res_num, res_name, res_num)],
                        coords_h_[('CE1', mol2_res_num, res_name, res_num)],
                        coords_h_[('NE2', mol2_res_num, res_name, res_num)],
                        ]
                center_ring, projection_from_center = ring_atom_position(ring, 3.5)
                ph4 = 'CZ'
                if projection_from_center != []:
                    for coord in projection_from_center:
                        if towards_binding_site(center_ring, coord, center):
                            coordinates_ph4.append((coord, ph4))
                
            if atom == 'CD2' or atom == 'CE1':
                root_atom = center_ring
                exit_atom = coords_[idx]
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 4)
                    coordinates_ph4.append((coord, ph4))
                    
            if atom == 'NE2' and ('HE2', mol2_res_num, res_name, res_num) in coords_h_:
                root_atom = coords_[idx]
                exit_atom = coords_h_[('HE2', mol2_res_num, res_name, res_num)]
                ph4 = 'O'
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    coordinates_ph4.append((coord, ph4))
            
            if atom == 'ND1' and ('HD1', mol2_res_num, res_name, res_num) in coords_h_:
                root_atom = coords_[idx]
                exit_atom = coords_h_[('HD1', mol2_res_num, res_name, res_num)]
                ph4 = 'O'
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    coordinates_ph4.append((coord, ph4))
                
            
        elif res_name == 'LYS':     
            if atom == 'CG':
                exit_atom = coords_[idx]
                root_atom = coords_h_[('CB', mol2_res_num, res_name, res_num)]
                ph4 = 'CA'
            elif atom == 'CD':
                root_atom = coords_h_[('CG', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                ph4 = 'CA'
            elif atom == 'NZ':
                root_atom = coords_h_[('CE', mol2_res_num, res_name, res_num)]
                exit_atom = coords_[idx]
                ph4 = 'OD1'

            if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 3.5)
                    coordinates_ph4.append((coord, ph4))
                

        elif res_name == 'ARG':
            if atom == 'CG':
                exit_atom = coords_[idx]
                root_atom = coords_h_[('CB', mol2_res_num, res_name, res_num)]
            elif atom == 'NE' or atom == 'NH1' or atom == 'NH2':
                ring = [coords_h_[('NE', mol2_res_num, res_name, res_num)],
                        coords_h_[('NH1', mol2_res_num, res_name, res_num)],
                        coords_h_[('NH2', mol2_res_num, res_name, res_num)],
                        ]
                center_ring, projection_from_center = ring_atom_position(ring, 3.5)
                ph4 = 'CZ'
                if projection_from_center != []:
                    for coord in projection_from_center:
                        if towards_binding_site(center_ring, coord, center):
                            coordinates_ph4.append((coord, ph4))
                
                root_atom = center_ring
                exit_atom = coords_[idx]
                ph4 = 'OD1'
                if towards_binding_site(root_atom, exit_atom, center):
                    coord = atom_position(root_atom, exit_atom, 4)
                    coordinates_ph4.append((coord, ph4))
    
    #print(coordinates_ph4)
    unique_coordinates_ph4 = list(set(coordinates_ph4))
        
    return unique_coordinates_ph4




def clean_clashes(coords_prot_, coords_ph4_):
    coords = []
    for c, ph4 in coords_ph4_:
        coords.append(c)
    neigh = NearestNeighbors(n_neighbors=len(coords_prot_), radius=1.5, algorithm='ball_tree').fit(coords_prot_)
    distances, indices = neigh.radius_neighbors(coords)
    #print(len(indices))

    clean_coords_ph4 = []
    for i in range(len(coords)):
        if not indices[i].size == 0:
            clean_coords_ph4.append(coords_ph4_[i])
        
    return clean_coords_ph4




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




    for i, value in enumerate(coordinates_ph4_):
        coord, ph4 = value
        x, y, z = coord
        #if ph4 == "NOP":
        #    ph4 = "DU"
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



coords_protein, atom_infos_protein, coords_h_protein = xtract_protein_coordinates(args.protein)
coords_cavity = xtract_cavity_coordinates(args.cavity)


bs_atoms = get_binding_site_coordinates(coords_protein, coords_cavity)
#for idx in bs_atoms:
#    print(idx, atom_infos_protein[idx])

coordinates_ph4 = prot_ph4(bs_atoms, coords_protein, atom_infos_protein, coords_h_protein)
coordinates_ph4_clean = clean_clashes(coords_cavity, coordinates_ph4)
print('# points:', len(coordinates_ph4_clean))
write_mol2(args.output, coordinates_ph4_clean)
