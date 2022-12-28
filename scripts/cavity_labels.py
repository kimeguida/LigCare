

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




import os
import numpy as np
import argparse
import copy
from sklearn.neighbors import NearestNeighbors
import iomol2, iodesc






def assign_ligph4_ints(ligph4_coords_, ints_coords_):
    # lig ph4 was placed as 
    neigh = NearestNeighbors(n_neighbors=1, radius=0.1,
                               algorithm='ball_tree').fit(ints_coords_)
    distances, indices = neigh.radius_neighbors(ligph4_coords_)
    return indices



def assign_cav_ligph4(cav_coords_, ligph4_coords_, d_=1.51):
    neigh = NearestNeighbors(n_neighbors=len(ligph4_coords_), radius=d_,
                               algorithm='ball_tree').fit(ligph4_coords_)
    distances, indices = neigh.radius_neighbors(cav_coords_)
    return indices




def label_cav(cav_ligph4_indices_, ligph4_ints_indices_, cav_features_,
                                ligph4_features_, ints_features_):

    apolar_int = {'CA'}
    polar_int = {'CZ', 'O', 'N', 'NZ', 'OD1', 'Zn'}

    labels = []
    for cav_idx in range(len(cav_ligph4_indices_)):
        #### labels:
        # 0 point index
        # 1 near any ligand ph4
        # 2 near ligand ph4 of same type
        # 3 near interacting ligand atom
        # 4 near polar-interacting ligand atom
        # ligand interaction type: can be many

        label = np.zeros(5)
        label[0] = cav_idx

        if not cav_ligph4_indices_[cav_idx].size == 0:
            label[1] = 1

            if cav_features_[cav_idx] in [ligph4_features_[ligph4_idx] 
                                for ligph4_idx in cav_ligph4_indices_[cav_idx]]:
                label[2] = 1
            
            for ligph4_idx in cav_ligph4_indices_[cav_idx]:
                if not ligph4_ints_indices_[ligph4_idx].size == 0:
                    label[3] = 1
                    polar_int_around = polar_int.intersection([ints_features_[i]
                                        for i in ligph4_ints_indices_[ligph4_idx]])
                    if len(polar_int_around) > 0:
                        label[4] = 1
                        break
        labels.append(label)


    return np.array(labels)
                    







if __name__ == '__main__':
    


    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cavity', type=str, required=True,
                help='volsite cavity mol2')
    parser.add_argument('-lph', '--ligph4', type=str, required=True,
                help='ligand ph4 mol2')
    parser.add_argument('-int', '--interaction', type=str, required=True,
                help='ints mol2')
    parser.add_argument('-o', '--ofile', type=str, required=True,
                help='output numpy file')
    args = parser.parse_args()



    ints_coords, ints_features, = iomol2.read_ints_mol2(args.interaction)
    #print(ints_coords, ints_features)

    cav_coords, cav_features = iomol2.read_cav_mol2(args.cavity)
    #print()

    ligph4_coords, ligph4_features = iomol2.read_cav_mol2(args.ligph4)
    #print(ligph4_coords, ligph4_features)


    ligph4_ints_indices = assign_ligph4_ints(ligph4_coords, ints_coords)
    #print(ligph4_ints_indices)


    cav_ligph4_indices = assign_cav_ligph4(cav_coords, ligph4_coords)
    #print(cav_ligph4_indices)


    cav_labels = label_cav(cav_ligph4_indices, ligph4_ints_indices, cav_features,
                                    ligph4_features, ints_features)

    #for l in cav_labels:
    #    print(l)

    iodesc.write_labels(cav_labels, args.ofile)

    #new_cav_labels = iodesc.read_labels(args.ofile)
    #print(new_cav_labels)