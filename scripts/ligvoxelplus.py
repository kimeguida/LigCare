
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
# test version of lig_to_ph4 module
# some cases need to be implemented


#from numba import njit, prange
import argparse
import numpy as np
from sklearn.neighbors import NearestNeighbors
import os
from time import strftime, localtime
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d
from matplotlib.lines import Line2D



parser = argparse.ArgumentParser()
parser.add_argument('-i', "--iligand", type=str, required=True) 
parser.add_argument('-r', "--resolution", type=float, required=False,
                    default=1.5)
parser.add_argument('-o', "--oligand", type=str, required=False,
                    default='ligvoxelplus.mol2')
args = parser.parse_args()


# lig_voxel transforms ligands atoms in voxels pseudo atoms
# lig_voxel_plus augments data by adding pseudo atoms in the 6 face-adjacent voxels



def lig_mol2_info(ifile):
    with open(ifile, 'r') as f:
        mol2 = f.read()
        
    # get atom block
    mol2_atom_block = mol2.split("@<TRIPOS>ATOM\n")[1].split("@<TRIPOS>")[0]
    mol2_atom_block = mol2_atom_block.split("\n")
    del mol2_atom_block[-1]

    atom_types = {}
    coordinates = {}
    coords = [[], [], [], []]
    points = []
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
            coords[0].append(x)
            coords[1].append(y)
            coords[2].append(z)
            coords[3].append(indice)
            points.append([x, y, z])

        except ValueError:
            continue
        
    mol2_bond_block = mol2.split("@<TRIPOS>BOND\n")[1].split("@<TRIPOS>")[0]
    mol2_bond_block = mol2_bond_block.split("\n")
    del mol2_bond_block[-1]
    bonds = {}
    for l in mol2_bond_block:
        mol2lines = l.split()
        try:
            indice1 = int(mol2lines[1])
            indice2 = int(mol2lines[2])
        except ValueError:
            continue
        if indice1 not in bonds:
            bonds[indice1] = [[atom_types[indice2].split('.')[0].upper()],
                                 [indice2]]
        else:
            bonds[indice1][0].append(atom_types[indice2].split('.')[0].upper())
            bonds[indice1][1].append(indice2)
 
        if indice2 not in bonds:
            bonds[indice2] = [[atom_types[indice1].split('.')[0].upper()],
                                 [indice1]]
        else:
            bonds[indice2][0].append(atom_types[indice1].split('.')[0].upper())
            bonds[indice2][1].append(indice1)
    return coordinates, bonds, atom_types, coords, points



def lig_voxel(coord_, resol_=1.5):
    x_barycenter = np.mean(coord_[0])
    x_min = min(coord_[0])
    x_max = max(coord_[0])
    
    y_barycenter = np.mean(coord_[1])
    y_min = min(coord_[1])
    y_max = max(coord_[1])
    
    z_barycenter = np.mean(coord_[2])
    z_min = min(coord_[2])
    z_max = max(coord_[2])
    
    
    ############## compute x_cells
    x_cells = [x_barycenter]
    c = True
    i = 1
    while c :
        x_cells.append(x_barycenter + resol_*i)
        i += 1
        if x_cells[-1] > x_max:
            c = False
            
    c = True
    i = 1        
    while c :
        x_cells.append(x_barycenter - resol_*i)
        i += 1
        if x_cells[-1] < x_min:
            c = False
            
    x_cells = sorted(x_cells)
    #print("x_cells", x_cells, x_min, x_max)
    
    
    ############## compute y_cells
    y_cells = [y_barycenter]
    c = True
    i = 1
    while c :
        y_cells.append(y_barycenter + resol_*i)
        i += 1
        if y_cells[-1] > y_max:
            c = False
            
    c = True
    i = 1        
    while c :
        y_cells.append(y_barycenter - resol_*i)
        i += 1
        if y_cells[-1] < y_min:
            c = False
            
    y_cells = sorted(y_cells)
    #print("y_cells", y_cells, y_min, y_max)
    
    
    ############## compute z_cells
    z_cells = [z_barycenter]
    c = True
    i = 1
    while c :
        z_cells.append(z_barycenter + resol_*i)
        i += 1
        if z_cells[-1] > z_max:
            c = False
            
    c = True
    i = 1        
    while c :
        z_cells.append(z_barycenter - resol_*i)
        i += 1
        if z_cells[-1] < z_min:
            c = False
            
    z_cells = sorted(z_cells)
    #print("z_cells", z_cells, z_min, z_max)
    #print(coord_[0], x_min, x_max)
    
    ###################### get voxels
    voxels = []
    for x, y, z in zip(coord_[0], coord_[1], coord_[2]):
        #print(x, y, z)
        for i in range(len(x_cells)-1):
            for j in range(len(y_cells)-1):
                for k in range(len(z_cells)-1):
                    if x_cells[i] < x <= x_cells[i+1] and\
                        y_cells[j] < y <= y_cells[j+1] and\
                        z_cells[k] < z <= z_cells[k+1]:
                        point_x = (x_cells[i] + x_cells[i+1])/2
                        point_y = (y_cells[j] + y_cells[j+1])/2
                        point_z = (z_cells[k] + z_cells[k+1])/2
                        if [point_x, point_y, point_z] not in voxels:
                            voxels.append([point_x, point_y, point_z])
    #print(len(voxels), len(coord_[0]))
    return voxels
        
        
#@njit(parallel=True)        
def lig_voxel_plus(coord_, resol_=1.5):
    x_barycenter = np.mean(coord_[0])
    x_min = min(coord_[0])
    x_max = max(coord_[0])
    
    y_barycenter = np.mean(coord_[1])
    y_min = min(coord_[1])
    y_max = max(coord_[1])
    
    z_barycenter = np.mean(coord_[2])
    z_min = min(coord_[2])
    z_max = max(coord_[2])
    
    
    ############## compute x_cells
    x_cells = [x_barycenter]
    c = True
    i = 1
    while c :
        x_cells.append(x_barycenter + resol_*i)
        i += 1
        if x_cells[-1] > x_max:
            c = False
            
    c = True
    i = 1        
    while c :
        x_cells.append(x_barycenter - resol_*i)
        i += 1
        if x_cells[-1] < x_min:
            c = False
            
    x_cells = sorted(x_cells)
    #print("x_cells", x_cells, x_min, x_max)
    
    
    ############## compute y_cells
    y_cells = [y_barycenter]
    c = True
    i = 1
    while c :
        y_cells.append(y_barycenter + resol_*i)
        i += 1
        if y_cells[-1] > y_max:
            c = False
            
    c = True
    i = 1        
    while c :
        y_cells.append(y_barycenter - resol_*i)
        i += 1
        if y_cells[-1] < y_min:
            c = False
            
    y_cells = sorted(y_cells)
    #print("y_cells", y_cells, y_min, y_max)
    
    
    ############## compute z_cells
    z_cells = [z_barycenter]
    c = True
    i = 1
    while c :
        z_cells.append(z_barycenter + resol_*i)
        i += 1
        if z_cells[-1] > z_max:
            c = False
            
    c = True
    i = 1        
    while c :
        z_cells.append(z_barycenter - resol_*i)
        i += 1
        if z_cells[-1] < z_min:
            c = False
            
    z_cells = sorted(z_cells)
    #print("z_cells", z_cells, z_min, z_max)
    #print(coord_[0], x_min, x_max)
    
    ###################### get voxels
    voxels = []
    
    for x, y, z in zip(coord_[0], coord_[1], coord_[2]):
        #print(x, y, z)
        for i in range(len(x_cells)-1):
            for j in range(len(y_cells)-1):
                for k in range(len(z_cells)-1):
                    if x_cells[i] < x <= x_cells[i+1] and\
                        y_cells[j] < y <= y_cells[j+1] and\
                        z_cells[k] < z <= z_cells[k+1]:
                        point_x = (x_cells[i] + x_cells[i+1])/2
                        point_y = (y_cells[j] + y_cells[j+1])/2
                        point_z = (z_cells[k] + z_cells[k+1])/2

                        if [point_x, point_y, point_z] not in voxels:
                            voxels.append([point_x, point_y, point_z])

                        if [point_x+resol_, point_y, point_z] not in voxels:
                            voxels.append([point_x+resol_, point_y, point_z])
                        if [point_x-resol_, point_y, point_z] not in voxels:
                            voxels.append([point_x-resol_, point_y, point_z])
                        
                        if [point_x, point_y+resol_, point_z] not in voxels:
                            voxels.append([point_x, point_y+resol_, point_z])
                        if [point_x, point_y-resol_, point_z] not in voxels:
                            voxels.append([point_x, point_y-resol_, point_z])

                        if [point_x, point_y, point_z+resol_] not in voxels:
                            voxels.append([point_x, point_y, point_z+resol_])
                        if [point_x, point_y, point_z-resol_] not in voxels:
                            voxels.append([point_x, point_y, point_z-resol_])
                        
    #print(len(voxels), len(coord_[0]))
    return voxels





def lig_ph4(coordinates_, bonds_, atom_types_):
    
    HETEROATOMS = ['O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I']
    HYDROGENS = ['H']
    
    for i in coordinates_:
        
        # Hydrophobic rules
        if atom_types_[i] == 'C.1':
            if all([atom not in HETEROATOMS for atom in bonds_[i][0]]):
                coordinates_[i].append('CA')
            else:
                coordinates_[i].append('DU') # NOP
                
        elif atom_types_[i] == 'C.2':
            if all([atom not in HETEROATOMS for atom in bonds_[i][0]]):
                coordinates_[i].append('CA')
            else:
                coordinates_[i].append('DU') # NOP
                
        elif atom_types_[i] == 'C.3':
            if all([atom not in HETEROATOMS for atom in bonds_[i][0]]):
                #print(i, "HYD")
                coordinates_[i].append('CA')
            else:
                #print(i, 'NOP')
                coordinates_[i].append('DU') # NOP
                
        elif atom_types_[i] == 'S':
            if all([atom not in HETEROATOMS for atom in bonds_[i][0]]):
                coordinates_[i].append('CA')
            else:
                coordinates_[i].append('DU') # NOP
        
        elif atom_types_[i] == 'S.2':
            if all([atom not in HETEROATOMS for atom in bonds_[i][0]]):
                coordinates_[i].append('CA')
            else:
                coordinates_[i].append('DU') # NOP
                
        elif atom_types_[i] == 'S.3':
            if all([atom not in HETEROATOMS for atom in bonds_[i][0]]):
                coordinates_[i].append('CA')
            else:
                coordinates_[i].append('DU') # NOP
        
        elif atom_types_[i] == 'F':
            coordinates_[i].append('CA')
        
        elif atom_types_[i] == 'Cl':
            coordinates_[i].append('CA')
        
        elif atom_types_[i] == 'Br':
            coordinates_[i].append('CA')
            
        elif atom_types_[i] == 'I':
            coordinates_[i].append('CA')
        
        
        # Aromatic rules
        elif atom_types_[i] == 'C.ar':
            coordinates_[i].append('CZ')
        
        # H-bond rules
        elif atom_types_[i] == 'N.ar':
            if any([atom in HYDROGENS for atom in bonds_[i][0]]):
                coordinates_[i].append('N')
            else:
                coordinates_[i].append('O')
        
        elif atom_types_[i] == 'N.am':
            if any([atom in HYDROGENS for atom in bonds_[i][0]]):
                coordinates_[i].append('N')
            else:
                coordinates_[i].append('O')
        
        elif atom_types_[i] == 'N.pl3':
            if any([atom in HYDROGENS for atom in bonds_[i][0]]):
                coordinates_[i].append('N')
            else:
                coordinates_[i].append('O')
                
        elif atom_types_[i] == 'N.2':
            if any([atom in HYDROGENS for atom in bonds_[i][0]]):
                coordinates_[i].append('N')
            else:
                coordinates_[i].append('O')
            
        elif atom_types_[i] == 'N.3':
            if any([atom in HYDROGENS for atom in bonds_[i][0]]):
                coordinates_[i].append('N')
            else:
                coordinates_[i].append('O')
                
        elif atom_types_[i] == 'O.3':
            if any([atom in HYDROGENS for atom in bonds_[i][0]]):
                coordinates_[i].append('OG')
            else:
                coordinates_[i].append('O')
        
        elif atom_types_[i] == 'O.2':
                coordinates_[i].append('O')
        
        elif atom_types_[i] == 'N.1':
                coordinates_[i].append('O')
        
        # ion rules
        elif atom_types_[i] == 'O.co2':
                coordinates_[i].append('OD1')
            
        elif atom_types_[i] == 'N.4':
                coordinates_[i].append('NZ')
        
        # These in ligands ????
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
            coordinates_[i].append('H')
        
        else:
            coordinates_[i].append('DU')


        # missing
        # P.3
        # C.cat, etc.
        # H.spc


def lig_ph4_hydrogens(coordinates_, bonds_, atom_types_):
    for i in coordinates_:
        if coordinates_[i][3] == 'H':
            ind_heavy_atm = bonds_[i][1]
            if len(ind_heavy_atm) != 1:
                print("problem in hydrogens")
            ind_heavy_atm = ind_heavy_atm[0]
            ph4_heavy_atm = coordinates_[ind_heavy_atm][3]

            # hydrogens takes the properties of their heavy atom
            # hygrogens as an extention of their heavy atom
            if ph4_heavy_atm == 'DU':
                coordinates_[i][3] = 'DU'
            if ph4_heavy_atm == 'CA':
                coordinates_[i][3] = 'CA'
            if ph4_heavy_atm == 'CZ':
                coordinates_[i][3] = 'CZ'
            if ph4_heavy_atm == 'N':
                coordinates_[i][3] = 'N'
            if ph4_heavy_atm == 'NZ':
                coordinates_[i][3] = 'NZ'
            if ph4_heavy_atm == 'OG':
                coordinates_[i][3] = 'OG'
            
            # unlikely
            # H-bon acceptors or anion won't bond to H
            if ph4_heavy_atm == 'O':
                coordinates_[i][3] = 'DU' # NOP
            if ph4_heavy_atm == 'OD1':
                coordinates_[i][3] = 'DU' # NOP
            if ph4_heavy_atm == 'NOP':
                coordinates_[i][3] = 'DU' # NOP




def lig_voxel_ph4(lig_coords_ph4_, ligvoxelplus_):
    lig_coords = [c[:3] for c in lig_coords_ph4_.values()]
    ligvoxelplus_ph4 = []
    knn = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(lig_coords)
    distances, indices = knn.kneighbors(ligvoxelplus_)

    for i in range(len(indices)):
        dist = distances[i][0]
        if dist <= 1.5:
            ph4_ind = indices[i][0]+1 # indices start with 1

            ligvoxelplus_ph4.append([ligvoxelplus_[i][0],
                                        ligvoxelplus_[i][1],
                                        ligvoxelplus_[i][2],
                                        lig_coords_ph4_[ph4_ind][3]])
        
    #print(indices, distances, ligvoxelplus_ph4)

    return ligvoxelplus_ph4, distances



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


    clean_coordinates_ph4 = []
    for point in (coordinates_ph4_):
        x, y, z, ph4 = point
        if ph4 == "NOP" or ph4 == "DU" or ph4 == "H":
            continue
        clean_coordinates_ph4.append(point)


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
                                            len(clean_coordinates_ph4), 0, 0, 0, 0)
    of_string += "{}\n".format(macromol_)
    of_string += "NO_CHARGES\n"
    of_string += "@<TRIPOS>ATOM"




    for i, point in enumerate(clean_coordinates_ph4):
        x, y, z, ph4 = point
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




if __name__ == '__main__':

    coordinates, bonds, atom_types, coords, points = lig_mol2_info(args.iligand)

    #print("COORD", coordinates)
    #print("BONDS", bonds)
    #print("ATOM TYPES", atom_types)
    #print("COORD", coords)

    ligvoxelplus = lig_voxel(coords, args.resolution)
    #ligvoxelplus_1 = lig_voxel_plus(coords, 1)

    #print(ligvoxelplus)
    
    #print(len(ligvoxelplus))

    #print(len(coords))


    unique_ligvoxelplus = []
    dup = 0
    for c in ligvoxelplus:
        if c in unique_ligvoxelplus:
            dup += 1
        unique_ligvoxelplus.append(c)

    #print(dup)

    lig_ph4(coordinates, bonds, atom_types)
    #lig_ph4_hydrogens(coordinates, bonds, atom_types)
    #print(coordinates)

    #ligvoxelplus_ph4, distances = lig_voxel_ph4(coordinates, ligvoxelplus)
    ligvoxelplus_ph4, distances = lig_voxel_ph4(coordinates, points)

    write_mol2(args.oligand, ligvoxelplus_ph4, macromol_="SMALL")











    """


    # Vizualisation
    plt.rc('axes', labelsize=12)
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)
    mpl.rcParams["figure.dpi"] = 200

    ax = plt.axes(projection='3d')

    
    for x, y, z in [[c[0], c[1], c[2]] for c in ligvoxelplus_ph4 if c[3] == 'CA']:
        ax.scatter(x, y, z, c='g')
    for x, y, z in [[c[0], c[1], c[2]] for c in ligvoxelplus_ph4 if c[3] == 'CZ']:
        ax.scatter(x, y, z, c='b')
    for x, y, z in [[c[0], c[1], c[2]] for c in ligvoxelplus_ph4 if c[3] == 'O']:
        ax.scatter(x, y, z, c='r')
    for x, y, z in [[c[0], c[1], c[2]] for c in ligvoxelplus_ph4 if c[3] == 'OD1']:
        ax.scatter(x, y, z, c='y')
    for x, y, z in [[c[0], c[1], c[2]] for c in ligvoxelplus_ph4 if c[3] == 'OG']:
        ax.scatter(x, y, z, c='orange')
    for x, y, z in [[c[0], c[1], c[2]] for c in ligvoxelplus_ph4 if c[3] == 'N']:
        ax.scatter(x, y, z, c='purple')
    for x, y, z in [[c[0], c[1], c[2]] for c in ligvoxelplus_ph4 if c[3] == 'NZ']:
        ax.scatter(x, y, z, c='pink')
    for x, y, z in [[c[0], c[1], c[2]] for c in ligvoxelplus_ph4 if c[3] == 'DU']:
        ax.scatter(x, y, z, c='gray')
    for x, y, z in [[c[0], c[1], c[2]] for c in ligvoxelplus_ph4 if c[3] == 'NOP']:
        ax.scatter(x, y, z, c='k')


    ph4_names = ['CA', 'CZ', 'O', 'OD1', 'OG', 'N', 'NZ', 'DU', 'NOP']
    colors = ['g', 'b', 'r', 'y', 'orange', 'purple', 'magenta', 'gray', 'k']
    legend_labels = [Line2D([], [], marker='.', linestyle='', color=c, label=ph4)
                        for c, ph4 in zip(colors, ph4_names)]
    ax.legend(handles=legend_labels)
    
 



    cavity = []
    colors_cav = [plt.cm.get_cmap('Greens')(0.5),
                  plt.cm.get_cmap('Blues')(0.5),
                  plt.cm.get_cmap('Reds')(0.5),
                  plt.cm.get_cmap('YlOrBr')(0.3),
                  plt.cm.get_cmap('Oranges')(0.5),
                  plt.cm.get_cmap('Purples')(0.5),
                  plt.cm.get_cmap('RdPu')(0.4),
                  plt.cm.get_cmap('Greys')(0.3),
                  plt.cm.get_cmap('Greys')(0.7)
                ]

    with open('2rh1_1_cavity6.mol2', 'r') as f:
        mol2_atom_block = f.read().split("@<TRIPOS>ATOM\n")[1].split("@<TRIPOS>")[0]
        mol2_atom_block = mol2_atom_block.split("\n")
    del mol2_atom_block[-1]

    for l in mol2_atom_block:
        cols = l.split()
        cav_ph4 = cols[1]
        x = float(cols[2])
        y = float(cols[3])
        z = float(cols[4])
        cavity.append([x, y, z, cav_ph4])

    for x, y, z in [[c[0], c[1], c[2]] for c in cavity if c[3] == 'CA']:
        ax.scatter(x, y, z, c=colors_cav[0])
    for x, y, z in [[c[0], c[1], c[2]] for c in cavity if c[3] == 'CZ']:
        ax.scatter(x, y, z, c=colors_cav[1])
    for x, y, z in [[c[0], c[1], c[2]] for c in cavity if c[3] == 'O']:
        ax.scatter(x, y, z, c=colors_cav[2])
    for x, y, z in [[c[0], c[1], c[2]] for c in cavity if c[3] == 'OD1']:
        ax.scatter(x, y, z, c=colors_cav[3])
    for x, y, z in [[c[0], c[1], c[2]] for c in cavity if c[3] == 'OG']:
        ax.scatter(x, y, z, c=colors_cav[4])
    for x, y, z in [[c[0], c[1], c[2]] for c in cavity if c[3] == 'N']:
        ax.scatter(x, y, z, c=colors_cav[5])
    for x, y, z in [[c[0], c[1], c[2]] for c in cavity if c[3] == 'NZ']:
        ax.scatter(x, y, z, c=colors_cav[6])
    for x, y, z in [[c[0], c[1], c[2]] for c in cavity if c[3] == 'DU']:
        ax.scatter(x, y, z, c=colors_cav[7])
    for x, y, z in [[c[0], c[1], c[2]] for c in cavity if c[3] == 'NOP']:
        ax.scatter(x, y, z, c=colors_cav[8])

    ph4_names = ['CA', 'CZ', 'O', 'OD1', 'OG', 'N', 'NZ', 'DU', 'NOP']
    legend_labels += [Line2D([], [], marker='.', linestyle='', color=c, label='{}_cav'.format(ph4))
                        for c, ph4 in zip(colors_cav, ph4_names)]
    ax.legend(handles=legend_labels)

    plt.show()

    """