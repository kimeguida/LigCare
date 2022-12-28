

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




import argparse
import numpy as np
from sklearn.neighbors import NearestNeighbors

# from utils import bsa, iomol2
import bsa, iomol2, iodesc







def calc_query_ph4_fingerprint(indexes_, ph4_features_):
    positions = {'CA': 0,
                 'CZ': 1,
                 'N': 2,
                 'NZ': 3,
                 'O': 4,
                 'OD1': 5,
                 'OG': 6,
                 'DU': 7,
                 }
    fingerprint = {}
    for point_idx in indexes_:
        point_fingerprint = np.zeros(len(positions))
        ph4_feature = ph4_features_[point_idx]
        if ph4_feature == 'CA':
            point_fingerprint[positions['CA']] = 1

        elif ph4_feature == 'CZ':
            point_fingerprint[positions['CZ']] = 1
            #point_fingerprint[positions['CA']] = 1

        elif ph4_feature == 'O':
            point_fingerprint[positions['O']] = 1

        elif ph4_feature == 'OD1':
            point_fingerprint[positions['OD1']] = 1
            #point_fingerprint[positions['O']] = 1

        elif ph4_feature == 'N':
            point_fingerprint[positions['N']] = 1

        elif ph4_feature == 'NZ':
            point_fingerprint[positions['NZ']] = 1
            #point_fingerprint[positions['N']] = 1

        elif ph4_feature == 'OG':
            #point_fingerprint[positions['O']] = 1
            #point_fingerprint[positions['N']] = 1
            point_fingerprint[positions['OG']] = 1

        elif ph4_feature == 'DU':
            point_fingerprint[positions['DU']] = 1   
        
        fingerprint[point_idx] = point_fingerprint

    return fingerprint






def init_bsa_bins(max_bsa_, step_):
    bsa_bounds = np.arange(0, max_bsa_+step_, step_)
    return bsa_bounds




def calc_query_bsa_fingerprint(indexes_, bsa_, bsa_bounds_):
    l = len(bsa_bounds_) - 1
    fingerprint = {}
    for point_idx in indexes_:
        bsa = bsa_[point_idx]
        point_fingerprint = np.zeros(l)
        for i in range(len(bsa_bounds_)-1):
            if bsa_bounds_[i] <= bsa < bsa_bounds_[i+1]:
                point_fingerprint[i] = 1
                break
        fingerprint[point_idx] = point_fingerprint

    return fingerprint




def get_neighbors(cav_coords_, radius_=[1.51, 3.01, 4.51]):
    neighbors = []
    for r in radius_:
        neigh = NearestNeighbors(n_neighbors=len(cav_coords_), radius=r,
                               algorithm='ball_tree').fit(cav_coords_)
        _, indices = neigh.radius_neighbors(cav_coords_)
        neighbors.append(indices)

    #### xtract neighbors in concentric spheres:
    concentric_neighbors = []
    for r in range(len(radius_)):
        r_neigh = []
        if r == 0:
            #pass
            for i, n in enumerate(neighbors[r]):
                r_neigh.append(set(n).difference([i]))
            r_neigh = [np.array(list(n)) for n in r_neigh]    
            concentric_neighbors.append(r_neigh)
        else:
            for i, n in enumerate(neighbors[r]):
                r_neigh.append(set(n).difference(neighbors[r-1][i]))
            r_neigh = [np.array(list(n)) for n in r_neigh]    
            concentric_neighbors.append(r_neigh)

    # validation: OK
    #for c in concentric_neighbors:
    #    print(c[0])

            
    return concentric_neighbors




def calc_neigh_fingerprints(indexes_, concentric_neighbors_, ph4_features_, buriedness_,
                                                                bsa_bounds_):
    
    positions = {'CA': 0,
                 'CZ': 1,
                 'N': 2,
                 'NZ': 3,
                 'O': 4,
                 'OD1': 5,
                 'OG': 6,
                 'DU': 7,
                 }

    l_bsa = len(bsa_bounds_) - 1
    l_ph4 = len(positions)
    l_conc = len(concentric_neighbors_)
    bsa_fingerprint = {}
    count_fingerprint = {}
    for point_idx in indexes_:
        point_bsa_fingerprint = np.zeros(l_ph4*l_conc*l_bsa)
        point_count_fingerprint = np.zeros(l_ph4*l_conc)
        
        for c, neighbors in enumerate(concentric_neighbors_):
            for neigh_idx in neighbors[point_idx]:
                neigh_ph4 = ph4_features_[neigh_idx]
                neigh_bsa = buriedness_[neigh_idx]
                #print('neibhors', 'circle', c, neigh_idx, neigh_ph4, neigh_bsa)
                pos_in_fp = c*l_ph4 + positions[neigh_ph4]
                point_count_fingerprint[pos_in_fp] += 1

                for pos_bsa in range(l_bsa):
                    if bsa_bounds_[pos_bsa] <= neigh_bsa < bsa_bounds_[pos_bsa+1]:
                        pos_in_fp = c*l_ph4*l_bsa + l_ph4*positions[neigh_ph4] + pos_bsa
                        point_bsa_fingerprint[pos_in_fp] += 1/len(neighbors[point_idx])
                        break # once one interval is found, others become irrelevant

        #print(point_idx, point_count_fingerprint, np.sum(point_bsa_fingerprint[2*96:3*96]))
        bsa_fingerprint[point_idx] = point_bsa_fingerprint
        count_fingerprint[point_idx] = point_count_fingerprint

    return bsa_fingerprint, count_fingerprint

        
                            

def calc_dist_centroids(indexes_, cav_coords_, ph4_features_):
    fingerprint = {}
    centroid = [np.mean(cav_coords_[:, 0]),
                np.mean(cav_coords_[:, 1]),
                np.mean(cav_coords_[:, 2]),]

    idx_CA = [i for i, ph4 in enumerate(ph4_features) if ph4 == 'CA']
    centroid_CA = [np.mean(cav_coords_[idx_CA, 0]),
                   np.mean(cav_coords_[idx_CA, 1]),
                   np.mean(cav_coords_[idx_CA, 2]),]

    idx_CZ = [i for i, ph4 in enumerate(ph4_features)  if ph4 == 'CZ']
    centroid_CZ = [np.mean(cav_coords_[idx_CZ, 0]),
                   np.mean(cav_coords_[idx_CZ, 1]),
                   np.mean(cav_coords_[idx_CZ, 2]),]               


    idx_O = [i for i, ph4 in enumerate(ph4_features) if ph4 == 'O']
    centroid_O = [np.mean(cav_coords_[idx_O, 0]),
                   np.mean(cav_coords_[idx_O, 1]),
                   np.mean(cav_coords_[idx_O, 2]),]

    idx_OG = [i for i, ph4 in enumerate(ph4_features) if ph4 == 'OG']
    centroid_OG = [np.mean(cav_coords_[idx_OG, 0]),
                   np.mean(cav_coords_[idx_OG, 1]),
                   np.mean(cav_coords_[idx_OG, 2]),]


    idx_OD1 = [i for i, ph4 in enumerate(ph4_features) if ph4 == 'OD1']
    centroid_OD1 = [np.mean(cav_coords_[idx_OD1, 0]),
                   np.mean(cav_coords_[idx_OD1, 1]),
                   np.mean(cav_coords_[idx_OD1, 2]),]

    idx_N = [i for i, ph4 in enumerate(ph4_features) if ph4 == 'N']
    centroid_N = [np.mean(cav_coords_[idx_N, 0]),
                   np.mean(cav_coords_[idx_N, 1]),
                   np.mean(cav_coords_[idx_N, 2]),]
    
    idx_NZ = [i for i, ph4 in enumerate(ph4_features) if ph4 == 'NZ']
    centroid_NZ = [np.mean(cav_coords_[idx_NZ, 0]),
                   np.mean(cav_coords_[idx_NZ, 1]),
                   np.mean(cav_coords_[idx_NZ, 2]),]
    
    idx_DU = [i for i, ph4 in enumerate(ph4_features) if ph4 == 'DU']
    centroid_DU = [np.mean(cav_coords_[idx_DU, 0]),
                   np.mean(cav_coords_[idx_DU, 1]),
                   np.mean(cav_coords_[idx_DU, 2]),]


    for point_idx in indexes_:
        vect = cav_coords_[point_idx] - centroid
        #print(centroid, cav_coords_[point_idx], vect)
        dist = np.linalg.norm(vect)

        dist_CA = np.linalg.norm(cav_coords_[point_idx] - centroid_CA)
        dist_CZ = np.linalg.norm(cav_coords_[point_idx] - centroid_CZ)
        dist_O = np.linalg.norm(cav_coords_[point_idx] - centroid_O)
        dist_OG = np.linalg.norm(cav_coords_[point_idx] - centroid_OG)
        dist_OD1 = np.linalg.norm(cav_coords_[point_idx] - centroid_OD1)
        dist_N = np.linalg.norm(cav_coords_[point_idx] - centroid_N)
        dist_NZ = np.linalg.norm(cav_coords_[point_idx] - centroid_NZ)
        dist_DU = np.linalg.norm(cav_coords_[point_idx] - centroid_DU)


        # positions: 'CA': 0, 'CZ': 1, 'N': 2, 'NZ': 3, 'O': 4, 'OD1': 5, 'OG': 6, 'DU': 7,


        fingerprint[point_idx] = [dist, dist_CA, dist_CZ, dist_N, dist_NZ,
                                        dist_O, dist_OD1, dist_OG, dist_DU]

        # in precedent version 07/2022: fingerprint[point_idx] = [dist]

    return fingerprint




def towards_point(root_, exit_, point_, max_base_dist_=2.6,
                    max_height_dist_=4.0, angle_factor_=1,):
    vec_1 = np.array(exit_)-np.array(root_)
    vec_2 = np.array(point_)-np.array(root_)
    cos_point_angle = np.dot(vec_1, vec_2)/(np.linalg.norm(vec_1)*np.linalg.norm(vec_2))

    max_base_dist = 2.6 #
    max_height_dist = 4.0 #
    generatrice_cone = np.sqrt(max_height_dist**2 + (max_base_dist/2)**2)
    cos_semi_angle = (max_base_dist/2)/generatrice_cone
    semi_angle = np.arccos(cos_semi_angle)
    semi_angle *= angle_factor_
    angle = semi_angle * 2
    cos_angle = np.cos(angle)

    if cos_angle <= cos_point_angle <=  np.cos(0):
        towards = True
    else:
        towards = False
    
    return towards




def translate_point_by_vector(root_, exit_):

    exit_vec = np.array(exit_) - np.array(root_)
    new_exit = exit_ + exit_vec

    return new_exit





def in_cone(atom_1_, exit_point_, atom_2_, semi_angle_, h_=5):
    inside = False
    # tan = sin /cos, angles are between 0 and 90 degrees, tan is positive
    radius = np.tan(semi_angle_) * h_
    vect_axis = np.array(exit_point_) - np.array(atom_1_)
    #print(np.tan(semi_angle_), radius, exit_point_)

    axis = vect_axis / np.linalg.norm(vect_axis) # norm of axis should be 1
    #print(np.linalg.norm(axis))

    proj_axis_dist = np.dot((np.array(atom_2_)-np.array(atom_1_)), axis)
    #print(proj_axis_dist)
    #print(atom_2_)
    if h_ >= proj_axis_dist >= 0: # means atom found in the right exit direction and cone heigth limited to h_                
        vect_proj_axis = axis * proj_axis_dist
        vect_proj = (np.array(atom_2_) - np.array(atom_1_)) - vect_proj_axis

        proj_dist = np.linalg.norm(vect_proj)
        cone_radius_proj = proj_axis_dist * radius / h_

        if proj_dist <= cone_radius_proj:
            inside = True

    return inside




def interactions_desc(cav_coords_, ph4_features_, protein_heavy_coords_,
                        protein_heavy_features_, prot_coords_dict_, radius_=4.0):

        # need dict ('RES', 'NUM', 'ATOM') : coords

        # prot_coords_dict

    #print(protein_heavy_coords_)
    #print(prot_coords_dict_)
    neigh = NearestNeighbors(n_neighbors=len(protein_heavy_coords_), radius=radius_,
                               algorithm='ball_tree').fit(protein_heavy_coords_)
    _, indices = neigh.radius_neighbors(cav_coords_)

    interactions_desc = []
    interactions_type_desc = []
    interactions_type_all_desc = []

    co_angle = np.pi/3
    nh_angle = np.pi/6


    #['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
    #'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
    #'TYR', 'VAL', 'HOH', 'MG', 'FE', 'MN', 'ZN', 'CO', 'CA', 'NA']

    for i, prot_indices in enumerate(indices):
        force = ['attractive', 'neutral', 'repulsive']
        point_inter_count_desc =  {f:0 for f in force}

        # HYD, AR, HBD, HBA, HBDA, +, -, MET,
        inter_type = ['HYD', 'AR', 'HBD', 'HBA', 'HBDA', '+', '-', 'MET',]
        point_inter_fp = {inter:0 for inter in inter_type}
        point_all_resi_inter_fp = {inter:0 for inter in inter_type}
        
       

        for prot_idx in prot_indices:

            res_atm = protein_heavy_features_[prot_idx]
            
            

            # backbone C=O
            if res_atm[2] == 'O' and res_atm[0] != 'HOH':
                point_all_resi_inter_fp['HBA'] += 1
                root_tmp = prot_coords_dict_[(res_atm[0], res_atm[1], 'C')]
                exit_tmp = protein_heavy_coords_[prot_idx]
                exit = translate_point_by_vector(root_tmp, exit_tmp)
                root = exit_tmp
                point = cav_coords_[i]
                if in_cone(root, exit, point, semi_angle_=co_angle):
                    #print(i, res_atm)
                    if ph4_features_[i] in ['N', 'NZ', 'OG']:
                        point_inter_count_desc['attractive'] += 1
                    elif ph4_features_[i] in ['O', 'OD1']:
                        point_inter_count_desc['repulsive'] += 1
                    else:
                        point_inter_count_desc['neutral'] += 1
                    point_inter_fp['HBA'] += 1

            


            # backbone NH
            elif res_atm[2] == 'N':
                point_all_resi_inter_fp['HBD'] += 1
                if (res_atm[0], res_atm[1], 'H') in prot_coords_dict_:
                    root = protein_heavy_coords_[prot_idx]
                    exit = prot_coords_dict_[(res_atm[0], res_atm[1], 'H')]
                    point = cav_coords_[i]
                    if in_cone(root, exit, point, semi_angle_=nh_angle):
                        #print(i, res_atm)
                        if ph4_features_[i] in ['O', 'OD1', 'OG']:
                            point_inter_count_desc['attractive'] += 1
                        elif ph4_features_[i] in ['N', 'NZ']:
                            point_inter_count_desc['repulsive'] += 1
                        else:
                            point_inter_count_desc['neutral'] += 1
                        point_inter_fp['HBD'] += 1

            


            # ASP side chain COO
            elif (res_atm[0], res_atm[2]) == ('ASP', 'OD1'):
                #point_all_resi_inter_fp['HBA'] += 1
                point_all_resi_inter_fp['-'] += 1
                if ph4_features_[i] in ['NZ']:
                    #print(i, res_atm)
                    point_inter_count_desc['attractive'] += 1
                else:
                    root_tmp = prot_coords_dict_[(res_atm[0], res_atm[1], 'CG')]
                    exit_tmp = protein_heavy_coords_[prot_idx]
                    exit = translate_point_by_vector(root_tmp, exit_tmp)
                    root = exit_tmp
                    point = cav_coords_[i]
                    if in_cone(root, exit, point, semi_angle_=co_angle):
                        #print(i, res_atm)
                        if ph4_features_[i] in ['N', 'OG']:
                            point_inter_count_desc['attractive'] += 1
                        elif ph4_features_[i] in ['O', 'OD1',]:
                            point_inter_count_desc['repulsive'] += 1
                        else:
                            point_inter_count_desc['neutral'] += 1
                        point_inter_fp['-'] += 1
                        #point_inter_fp['HBA'] += 1



            elif (res_atm[0], res_atm[2]) == ('ASP', 'OD2'):
                point_all_resi_inter_fp['-'] += 1
                #point_all_resi_inter_fp['HBA'] += 1
                if ph4_features_[i] in ['NZ']:
                    #print(i, res_atm)
                    point_inter_count_desc['attractive'] += 1
                else:
                    root_tmp = prot_coords_dict_[(res_atm[0], res_atm[1], 'CG')]
                    exit_tmp = protein_heavy_coords_[prot_idx]
                    exit = translate_point_by_vector(root_tmp, exit_tmp)
                    root = exit_tmp
                    point = cav_coords_[i]
                    if in_cone(root, exit, point, semi_angle_=co_angle):
                        #print(i, res_atm)
                        if ph4_features_[i] in ['N', 'OG']:
                            point_inter_count_desc['attractive'] += 1
                        elif ph4_features_[i] in ['O', 'OD1',]:
                            point_inter_count_desc['repulsive'] += 1
                        else:
                            point_inter_count_desc['neutral'] += 1
                        point_inter_fp['-'] += 1
                        #point_inter_fp['HBA'] += 1



            # GLU side chain COO
            elif (res_atm[0], res_atm[2]) == ('GLU', 'OE1'):
                point_all_resi_inter_fp['-'] += 1
                #point_all_resi_inter_fp['HBA'] += 1
                if ph4_features_[i] in ['NZ']:
                    #print(i, res_atm)
                    point_inter_count_desc['attractive'] += 1
                else:
                    root_tmp = prot_coords_dict_[(res_atm[0], res_atm[1], 'CD')]
                    exit_tmp = protein_heavy_coords_[prot_idx]
                    exit = translate_point_by_vector(root_tmp, exit_tmp)
                    root = exit_tmp
                    point = cav_coords_[i]
                    if in_cone(root, exit, point, semi_angle_=co_angle):
                        #print(i, res_atm)
                        if ph4_features_[i] in ['N', 'OG']:
                            point_inter_count_desc['attractive'] += 1
                        elif ph4_features_[i] in ['O', 'OD1',]:
                            point_inter_count_desc['repulsive'] += 1
                        else:
                            point_inter_count_desc['neutral'] += 1
                        point_inter_fp['-'] += 1
                        #point_inter_fp['HBA'] += 1


            elif (res_atm[0], res_atm[2]) == ('GLU', 'OE2'):
                point_all_resi_inter_fp['-'] += 1
                #point_all_resi_inter_fp['HBA'] += 1
                if ph4_features_[i] in ['NZ']:
                    #print(i, res_atm)
                    point_inter_count_desc['attractive'] += 1
                else:
                    root_tmp = prot_coords_dict_[(res_atm[0], res_atm[1], 'CD')]
                    exit_tmp = protein_heavy_coords_[prot_idx]
                    exit = translate_point_by_vector(root_tmp, exit_tmp)
                    root = exit_tmp
                    point = cav_coords_[i]
                    if in_cone(root, exit, point, semi_angle_=co_angle):
                        #print(i, res_atm)
                        if ph4_features_[i] in ['N', 'OG']:
                            point_inter_count_desc['attractive'] += 1
                        elif ph4_features_[i] in ['O', 'OD1',]:
                            point_inter_count_desc['repulsive'] += 1
                        else:
                            point_inter_count_desc['neutral'] += 1
                        point_inter_fp['-'] += 1   
                        #point_inter_fp['HBA'] += 1


            # flexible aa, no direction check
            
            # ARG side chain
            elif (res_atm[0], res_atm[2]) == ('ARG', 'CZ') or\
                    (res_atm[0], res_atm[2]) == ('ARG', 'NE') or\
                    (res_atm[0], res_atm[2]) == ('ARG', 'NH1') or\
                    (res_atm[0], res_atm[2]) == ('ARG', 'NH2'):
                point_all_resi_inter_fp['+'] += 1
                point_inter_fp['+'] += 1
                #point_all_resi_inter_fp['HBD'] += 1
                #point_inter_fp['HBD'] += 1
                
                if ph4_features_[i] in ['OD1', 'O', 'OG', 'CZ']:
                    #print(i, res_atm)
                    point_inter_count_desc['attractive'] += 1
                elif ph4_features_[i] in ['NZ', 'N']:
                    #print(i, res_atm)
                    point_inter_count_desc['repulsive'] += 1
                else:
                    point_inter_count_desc['neutral'] += 1
            
            # LYS side chain
            elif (res_atm[0], res_atm[2]) == ('ARG', 'NZ'):
                point_all_resi_inter_fp['+'] += 1
                point_inter_fp['+'] += 1
                #point_all_resi_inter_fp['HBD'] += 1
                #point_inter_fp['HBD'] += 1
                if ph4_features_[i] in ['OD1', 'O', 'OG', 'CZ']:
                    #print(i, res_atm)
                    point_inter_count_desc['attractive'] += 1
                elif ph4_features_[i] in ['NZ', 'N']:
                    #print(i, res_atm)
                    point_inter_count_desc['repulsive'] += 1
                else:
                    point_inter_count_desc['neutral'] += 1


            # HIS side chain N
            elif (res_atm[0], res_atm[2]) == ('HIS', 'ND1'):
                if (res_atm[0], res_atm[1], 'HD1') not in prot_coords_dict_:
                    point_all_resi_inter_fp['HBA'] += 1
                    root_tmp = (prot_coords_dict_[(res_atm[0], res_atm[1], 'CG')] +
                                prot_coords_dict_[(res_atm[0], res_atm[1], 'CE1')])/2
                    exit_tmp = protein_heavy_coords_[prot_idx]
                    exit = translate_point_by_vector(root_tmp, exit_tmp)
                    root = exit_tmp
                    point = cav_coords_[i]
                    if in_cone(root, exit, point, semi_angle_=nh_angle):
                        #print(i, res_atm)
                        if ph4_features_[i] in ['N', 'NZ', 'OG']:
                            point_inter_count_desc['attractive'] += 1
                        elif ph4_features_[i] in ['O', 'OD1']:
                            point_inter_count_desc['repulsive'] += 1
                        else:
                            point_inter_count_desc['neutral'] += 1
                        point_inter_fp['HBA'] += 1

                elif (res_atm[0], res_atm[1], 'HD1') in prot_coords_dict_:
                    point_all_resi_inter_fp['HBD'] += 1
                    root = protein_heavy_coords_[prot_idx]
                    exit = prot_coords_dict_[(res_atm[0], res_atm[1], 'HD1')]
                    point = cav_coords_[i]
                    if in_cone(root, exit, point, semi_angle_=nh_angle):
                        #print(i, res_atm)
                        if ph4_features_[i] in ['N', 'NZ']:
                            point_inter_count_desc['repulsive'] += 1
                        elif ph4_features_[i] in ['O', 'OD1', 'OG']:
                            point_inter_count_desc['attractive'] += 1
                        else:
                            point_inter_count_desc['neutral'] += 1
                        point_inter_fp['HBD'] += 1



            elif (res_atm[0], res_atm[2]) == ('HIS', 'NE2'):
                if (res_atm[0], res_atm[1], 'HE2') not in prot_coords_dict_:
                    point_all_resi_inter_fp['HBA'] += 1
                    root_tmp = (prot_coords_dict_[(res_atm[0], res_atm[1], 'CD2')] +
                                prot_coords_dict_[(res_atm[0], res_atm[1], 'CE1')])/2
                    exit_tmp = protein_heavy_coords_[prot_idx]
                    exit = translate_point_by_vector(root_tmp, exit_tmp)
                    root = exit_tmp
                    point = cav_coords_[i]
                    if in_cone(root, exit, point, semi_angle_=nh_angle):
                        #print(i, res_atm)
                        if ph4_features_[i] in ['N', 'NZ', 'OG']:
                            point_inter_count_desc['attractive'] += 1
                        elif ph4_features_[i] in ['O', 'OD1']:
                            point_inter_count_desc['repulsive'] += 1
                        else:
                            point_inter_count_desc['neutral'] += 1
                        point_inter_fp['HBA'] += 1

                elif (res_atm[0], res_atm[1], 'HE2') in prot_coords_dict_:
                    point_all_resi_inter_fp['HBD'] += 1
                    root = protein_heavy_coords_[prot_idx]
                    exit = prot_coords_dict_[(res_atm[0], res_atm[1], 'HE2')]
                    point = cav_coords_[i]
                    if in_cone(root, exit, point, semi_angle_=nh_angle):
                        #print(i, res_atm)
                        if ph4_features_[i] in ['N', 'NZ']:
                            point_inter_count_desc['repulsive'] += 1
                        elif ph4_features_[i] in ['O', 'OD1', 'OG']:
                            point_inter_count_desc['attractive'] += 1
                        else:
                            point_inter_count_desc['neutral'] += 1
                        point_inter_fp['HBD'] += 1

            # same for other HIS: HIE, HIP, HID



            # ASN side chain
            elif (res_atm[0], res_atm[2]) == ('ASN', 'OD1'):
                point_all_resi_inter_fp['HBA'] += 1
                root_tmp = prot_coords_dict_[(res_atm[0], res_atm[1], 'CG')]
                exit_tmp = protein_heavy_coords_[prot_idx]
                exit = translate_point_by_vector(root_tmp, exit_tmp)
                root = exit_tmp
                point = cav_coords_[i]
                if in_cone(root, exit, point, semi_angle_=co_angle):
                    #print(i, res_atm)
                    if ph4_features_[i] in ['N', 'NZ', 'OG']:
                        point_inter_count_desc['attractive'] += 1
                    elif ph4_features_[i] in ['O', 'OD1']:
                        point_inter_count_desc['repulsive'] += 1
                    else:
                        point_inter_count_desc['neutral'] += 1
                    point_inter_fp['HBA'] += 1

           
            elif (res_atm[0], res_atm[2]) == ('ASN', 'ND2'):
                point_all_resi_inter_fp['HBD'] += 1
                root = protein_heavy_coords_[prot_idx]
                exit = prot_coords_dict_[(res_atm[0], res_atm[1], 'HD21')]
                point = cav_coords_[i]
                if in_cone(root, exit, point, semi_angle_=nh_angle):
                    #print(i, res_atm)
                    if ph4_features_[i] in ['O', 'OD1', 'OG']:
                        point_inter_count_desc['attractive'] += 1
                    elif ph4_features_[i] in ['N', 'NZ']:
                        point_inter_count_desc['repulsive'] += 1
                    else:
                        point_inter_count_desc['neutral'] += 1
                    point_inter_fp['HBD'] += 1

                exit = prot_coords_dict_[(res_atm[0], res_atm[1], 'HD22')]
                point = cav_coords_[i]
                if in_cone(root, exit, point, semi_angle_=nh_angle):
                    #print(i, res_atm)
                    if ph4_features_[i] in ['O', 'OD1', 'OG']:
                        point_inter_count_desc['attractive'] += 1
                    elif ph4_features_[i] in ['N', 'NZ']:
                        point_inter_count_desc['repulsive'] += 1
                    else:
                        point_inter_count_desc['neutral'] += 1
                    point_inter_fp['HBD'] += 1


            # GLN side chain
            elif (res_atm[0], res_atm[2]) == ('GLN', 'OE1'):
                point_all_resi_inter_fp['HBA'] += 1
                root_tmp = prot_coords_dict_[(res_atm[0], res_atm[1], 'CD')]
                exit_tmp = protein_heavy_coords_[prot_idx]
                exit = translate_point_by_vector(root_tmp, exit_tmp)
                root = exit_tmp
                point = cav_coords_[i]
                if in_cone(root, exit, point, semi_angle_=co_angle):
                    #print(i, res_atm)
                    if ph4_features_[i] in ['N', 'NZ', 'OG']:
                        point_inter_count_desc['attractive'] += 1
                    elif ph4_features_[i] in ['O', 'OD1']:
                        point_inter_count_desc['repulsive'] += 1
                    else:
                        point_inter_count_desc['neutral'] += 1
                    point_inter_fp['HBA'] += 1

            
            elif (res_atm[0], res_atm[2]) == ('GLN', 'NE2'):
                point_all_resi_inter_fp['HBD'] += 1
                root = protein_heavy_coords_[prot_idx]
                exit = prot_coords_dict_[(res_atm[0], res_atm[1], 'HE21')]
                point = cav_coords_[i]
                if in_cone(root, exit, point, semi_angle_=nh_angle):
                    #print(i, res_atm)
                    if ph4_features_[i] in ['O', 'OD1', 'OG']:
                        point_inter_count_desc['attractive'] += 1
                    elif ph4_features_[i] in ['N', 'NZ']:
                        point_inter_count_desc['repulsive'] += 1
                    else:
                        point_inter_count_desc['neutral'] += 1
                    point_inter_fp['HBD'] += 1

                exit = prot_coords_dict_[(res_atm[0], res_atm[1], 'HE22')]
                point = cav_coords_[i]
                if in_cone(root, exit, point, semi_angle_=nh_angle):
                    #print(i, res_atm)
                    if ph4_features_[i] in ['O', 'OD1', 'OG']:
                        point_inter_count_desc['attractive'] += 1
                    elif ph4_features_[i] in ['N', 'NZ']:
                        point_inter_count_desc['repulsive'] += 1
                    else:
                        point_inter_count_desc['neutral'] += 1
                    point_inter_fp['HBD'] += 1


            # SER
            elif (res_atm[0], res_atm[2]) == ('SER', 'OG'): # no need for direction
                #print(i, res_atm)
                point_all_resi_inter_fp['HBDA'] += 1
                point_inter_fp['HBDA'] += 1
                if ph4_features_[i] in ['O', 'OD1', 'OG', 'N', 'NZ']:
                    point_inter_count_desc['attractive'] += 1
                else:
                    point_inter_count_desc['neutral'] += 1
                

            # THR
            elif (res_atm[0], res_atm[2]) == ('THR', 'OG1'): # no need for direction
                #print(i, res_atm)
                point_all_resi_inter_fp['HBDA'] += 1
                point_inter_fp['HBDA'] += 1
                if ph4_features_[i] in ['O', 'OD1', 'OG', 'N', 'NZ']:
                    point_inter_count_desc['attractive'] += 1
                else:
                    point_inter_count_desc['neutral'] += 1


            elif res_atm[0] == 'HOH': # no need for direction
                #print(i, res_atm)
                point_all_resi_inter_fp['HBDA'] += 1
                point_inter_fp['HBDA'] += 1
                if ph4_features_[i] in ['O', 'OD1', 'OG', 'N', 'NZ']:
                    point_inter_count_desc['attractive'] += 1
                else:
                    point_inter_count_desc['neutral'] += 1
                

            # TYR
            elif (res_atm[0], res_atm[2]) == ('TYR', 'OH'): # no need for direction
                #print(i, res_atm)
                point_all_resi_inter_fp['HBDA'] += 1
                point_inter_fp['HBDA'] += 1
                if ph4_features_[i] in ['O', 'OD1', 'OG', 'N', 'NZ', 'CZ', 'CA']:
                    point_inter_count_desc['attractive'] += 1
                else:
                    point_inter_count_desc['neutral'] += 1                


            elif (res_atm[0], res_atm[2]) == ('TYR', 'CG') or\
                    (res_atm[0], res_atm[2]) == ('TYR', 'CD1') or\
                    (res_atm[0], res_atm[2]) == ('TYR', 'CD2') or\
                    (res_atm[0], res_atm[2]) == ('TYR', 'CE1') or\
                    (res_atm[0], res_atm[2]) == ('TYR', 'CE2') or\
                    (res_atm[0], res_atm[2]) == ('TYR', 'CZ') :
                #print(i, res_atm)
                point_all_resi_inter_fp['AR'] += 1
                point_inter_fp['AR'] += 1
                #point_all_resi_inter_fp['HYD'] += 1
                #point_inter_fp['HYD'] += 1
                if ph4_features_[i] in ['CZ', 'CA']:
                    point_inter_count_desc['attractive'] += 1
                else:
                    point_inter_count_desc['neutral'] += 1
                
                #point_all_resi_inter_fp['HYD'] += 1
                #point_inter_fp['HYD'] += 1

            
            # TRP
            elif (res_atm[0], res_atm[2]) == ('TRP', 'CG') or\
                    (res_atm[0], res_atm[2]) == ('TRP', 'CD1') or\
                    (res_atm[0], res_atm[2]) == ('TRP', 'CD2') or\
                    (res_atm[0], res_atm[2]) == ('TRP', 'CE2') or\
                    (res_atm[0], res_atm[2]) == ('TRP', 'CE3') or\
                    (res_atm[0], res_atm[2]) == ('TRP', 'CZ2') or\
                    (res_atm[0], res_atm[2]) == ('TRP', 'CZ3') or\
                    (res_atm[0], res_atm[2]) == ('TRP', 'CH2'):
                #print(i, res_atm)
                point_all_resi_inter_fp['AR'] += 1
                point_inter_fp['AR'] += 1
                #point_all_resi_inter_fp['HYD'] += 1
                #point_inter_fp['HYD'] += 1
                if ph4_features_[i] in ['CZ', 'CA']:
                    point_inter_count_desc['attractive'] += 1
                else:
                    point_inter_count_desc['neutral'] += 1
                
                #point_all_resi_inter_fp['HYD'] += 1
                #point_inter_fp['HYD'] += 1

            elif (res_atm[0], res_atm[2]) == ('TRP', 'NE1'):
                point_all_resi_inter_fp['AR'] += 1
                point_inter_fp['AR'] += 1
                #point_all_resi_inter_fp['HYD'] += 1
                #point_inter_fp['HYD'] += 1
                point_all_resi_inter_fp['HBD'] += 1
                if ph4_features_[i] in ['CZ', 'CA']:
                    point_inter_count_desc['attractive'] += 1
                else:
                    root = protein_heavy_coords_[prot_idx]
                    exit = prot_coords_dict_[(res_atm[0], res_atm[1], 'HE1')]
                    point = cav_coords_[i]
                    if in_cone(root, exit, point, semi_angle_=nh_angle):
                        #print(i, res_atm)
                        if ph4_features_[i] in ['O', 'OD1', 'OG']:
                            point_inter_count_desc['attractive'] += 1
                        elif ph4_features_[i] in ['N', 'NZ']:
                            point_inter_count_desc['repulsive'] += 1
                        point_inter_fp['HBD'] += 1




            # PHE
            elif (res_atm[0], res_atm[2]) == ('PHE', 'CG') or\
                    (res_atm[0], res_atm[2]) == ('PHE', 'CD1') or\
                    (res_atm[0], res_atm[2]) == ('PHE', 'CD2') or\
                    (res_atm[0], res_atm[2]) == ('PHE', 'CE1') or\
                    (res_atm[0], res_atm[2]) == ('PHE', 'CE2') or\
                    (res_atm[0], res_atm[2]) == ('PHE', 'CZ') :
                point_all_resi_inter_fp['AR'] += 1
                point_inter_fp['AR'] += 1
                #point_all_resi_inter_fp['HYD'] += 1
                #point_inter_fp['HYD'] += 1
                
                #print(i, res_atm)
                if ph4_features_[i] in ['CZ', 'CA']:
                    point_inter_count_desc['attractive'] += 1
                else:
                    point_inter_count_desc['neutral'] += 1

            # metals
            elif res_atm[0] in ['ZN', 'FE', 'MG', 'MN', 'CO', 'CA', 'NA']:
                point_all_resi_inter_fp['MET'] += 1
                point_inter_fp['MET'] += 1
                if ph4_features_[i] in ['O', 'OD1']:
                    point_inter_count_desc['attractive'] += 1
                elif ph4_features_[i] in ['N', 'NZ']:
                    point_inter_count_desc['repulsive'] += 1
                else:
                    point_inter_count_desc['neutral'] += 1


            # other
            elif res_atm[2] not in ['CA', 'CB']:
                point_all_resi_inter_fp['HYD'] += 1
                point_inter_fp['HYD'] += 1
                #print(i, res_atm)
                if ph4_features_[i] in ['CA']:
                    point_inter_count_desc['attractive'] += 1
                else:
                    point_inter_count_desc['neutral'] += 1



        interactions_desc.append([point_inter_count_desc[f] for f in force])
        interactions_type_desc.append([point_inter_fp[i] for i in inter_type])
        interactions_type_all_desc.append([point_all_resi_inter_fp[i] for i in inter_type])

        """
        interactions_desc.append([point_inter_count_desc])
        interactions_type_desc.append([point_inter_fp])
        interactions_type_all_desc.append([point_all_resi_inter_fp])
        """

    return interactions_desc, interactions_type_desc, interactions_type_all_desc


   

def relative_interaction_fingerprint(ph4_features_, inter_fingerprint_):

    features = set(ph4_features_)
    inter_quantiles = []
    for pos in range(len(inter_fingerprint_[0])):
        features_inter_quantiles = {}
        for ph4 in features:
            # quantiles per ph4 type
            ph4_inter = [inter_fingerprint_[i][pos]
                            for i in range(len(inter_fingerprint_))
                            if ph4_features_[i] == ph4]

            # quantiles all ph4 types merged
            #ph4_inter = [inter_fingerprint_[i][pos]
            #                for i in range(len(inter_fingerprint_))]

            quantiles = np.quantile(ph4_inter, [0.5, 0.75])

            features_inter_quantiles[ph4] = quantiles
            #print(pos, ph4, quantiles)

        inter_quantiles.append(features_inter_quantiles)

    # inter_quantiles content:
    #   [{ph4:[q1, q2, etc.]}, {ph4:[q1, q2, etc.]}, etc.]


    relative_fp = []
    for i in range(len(inter_fingerprint_)):
        relative_fp_point = []
        ph4 = ph4_features_[i]
        for bin_pos, count in enumerate(inter_fingerprint_[i]):
            if count < inter_quantiles[bin_pos][ph4][0]: # median
                relative_fp_point.append(0)
            elif inter_quantiles[bin_pos][ph4][0] <= count < inter_quantiles[bin_pos][ph4][1]:
                relative_fp_point.append(1)
            else:
                relative_fp_point.append(2)
        relative_fp.append(relative_fp_point)


    return relative_fp

        
        #print(pos, features_median)
        #print()










if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cavity', type=str, required=True,
                help='volsite cavity')
    parser.add_argument('-p', '--protein', type=str, required=True,
                help='protein file')
    parser.add_argument('-ibsa', '--inputbsa', type=str, required=False,
                help='buruiedness BSA file')
    parser.add_argument('-o', '--output', type=str, required=True,
                help='desc file')
    parser.add_argument('-obsa', '--outputbsa', type=str, required=False,
                help='output BSA file if calculated')
    args = parser.parse_args()



    cav_coords, ph4_features = iomol2.read_cav_mol2(args.cavity)
    
    prot_coords = iomol2.read_prot_mol2(args.protein)
    prot_coords_inter, prot_features_inter, prot_coords_dict_inter = iomol2.read_prot_mol2_features(args.protein)

    reduced_prot_indices = bsa.reduce_prot_coords(cav_coords, prot_coords)
    
    proj_points_origin = bsa.init_projections(False)
    bsa_bounds = init_bsa_bins(max_bsa_=len(proj_points_origin), step_=10)

    concentric_neighbors = get_neighbors(cav_coords)



    #### calc or read bsa | calculation is slow
    if args.inputbsa is not None:
        buriedness = iodesc.read_bsa(args.inputbsa)
    else:
        buriedness = bsa.calc_bsa(cav_coords, prot_coords, proj_points_origin, reduced_prot_indices)
    #print(buriedness)

   
    #### compute descriptors

    indexes = range(len(cav_coords)) # iterable, indexes of point whose descriptors are to be calculated
    #indexes = [0, 1]
    
    
    query_ph4_fingerprint = calc_query_ph4_fingerprint(indexes, ph4_features)
    #print(len(query_ph4_fingerprint[0]))

    query_bsa_fingerprint = calc_query_bsa_fingerprint(indexes, buriedness, bsa_bounds)
    #print(len(query_bsa_fingerprint[0]))
    
    query_dist_centroids = calc_dist_centroids(indexes, cav_coords, ph4_features)
    #print(query_dist_centroid)
    

    neigh_bsa_fingerprint, neigh_count_fingerprint = calc_neigh_fingerprints(
                    indexes, concentric_neighbors, ph4_features, buriedness, bsa_bounds)
    #print(len(neigh_count_fingerprint[0]))
    #print(len(neigh_bsa_fingerprint[0]))


    interactions_desc, interactions_type_desc, interactions_type_all_desc = interactions_desc(
                                    cav_coords, ph4_features, prot_coords_inter,
                                    prot_features_inter, prot_coords_dict_inter, radius_=4.0)

    relative_interactions_type_desc = relative_interaction_fingerprint(ph4_features, interactions_type_desc)
    #print(relative_interactions_type_desc)
    
    bsa_desc = [[buriedness[p]] for p in range(len(buriedness))]

    
    iodesc.write_desc(indexes, [query_ph4_fingerprint, query_bsa_fingerprint, query_dist_centroids,
                            neigh_count_fingerprint, neigh_bsa_fingerprint, bsa_desc,
                            interactions_desc, interactions_type_desc, interactions_type_all_desc],
                            args.output)


    if args.outputbsa is not None:
        iodesc.write_bsa(buriedness, args.outputbsa)


    #print(query_dist_centroids[0])
    indexes, descriptors = iodesc.read_desc(args.output)
    #print(indexes, descriptors[0])
    #print(indexes, descriptors[0][20:29])
    