


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



import networkx as nx
import os
import itertools
from time import strftime, localtime

import argparse
import numpy as np
from sklearn.neighbors import NearestNeighbors
from collections import defaultdict
from maximal_clique import *
import operator

import scipy
from scipy.spatial import distance
from scipy.spatial.transform import Rotation

from ph4_matches import ph4_matching, ph4_scores
from atomic_params import vdw_radius

import iomol2







def pairwise_distances_as_dict(coords_, atm_types_):
    pairs = {}
    for i in range(len(coords_)-1):
        for j in range(i+1, len(coords_)):
            dist = np.linalg.norm(np.array(coords_[i])-np.array(coords_[j]))
            pairs[(i, j)] = ((atm_types_[i], atm_types_[j]), round(dist, 4))

    return pairs



def associate_pairs(pairs_1_, pairs_2_, d_=1.5):
    product_edges = []
    scores = {}
    for p1, value in pairs_1_.items():
        #print(p1, value)
        ph4, d1 = value
        #print(ph4)
        if (ph4 == ('CA', 'CA') and d1 < 5): # to be modified
            continue
        #b_inf = max(0, d1 - d1 * d_)
        #b_sup = d1 + d1 * d_
        b_inf = d1 - d_
        b_sup = d1 + d_
        #print(b_inf, b_sup)
        ph4_rv = (ph4[1], ph4[0])
        for p2, v2 in pairs_2_.items():


            if b_inf <= v2[1] <= b_sup:
                #print(d1, v2[1], b_inf, b_sup)
                #if (p1[0] == 18 and p2[0] == 5) or (p1[1] == 18 and p2[1] == 5):
                #    print("found")


                if ((ph4[0], v2[0][0]), (ph4[1], v2[0][1])) in ph4_matching:
                    product_edges.append(((p1[0], p2[0]), (p1[1], p2[1])))
                    scores[((p1[0], p2[0]), (p1[1], p2[1]))] =\
                            ph4_matching[((ph4[0], v2[0][0]), (ph4[1], v2[0][1]))]

                if ((ph4_rv[0], v2[0][0]), (ph4_rv[1], v2[0][1])) in ph4_matching:
                    product_edges.append(((p1[0], p2[1]), (p1[1], p2[0])))
                    scores[((p1[0], p2[1]), (p1[1], p2[0]))] =\
                            ph4_matching[((ph4_rv[0], v2[0][0]), (ph4_rv[1], v2[0][1]))]
    return product_edges, scores




def nx_graph_from_edges(edges_):
    G = nx.Graph()
    G.add_edges_from(edges_)
    return G



   
def nx_to_adj_graph(nx_graph_):
    graph_adjacency = {}
    for n, d in nx_graph_.adjacency():
        graph_adjacency[n] = list(d.keys())
    return graph_adjacency




def combinatorial_nodes(nodes_):
    # return ascending sorted pairs
    edges = list(itertools.combinations(list(nodes_), 2))
    #print(nodes_)
    #print(edges)
    sorted_edges = [tuple(sorted(e)) for e in edges]
    #print(sorted_edges)
    return list(sorted_edges)



def reform_graph_product(product_graph_, graph_1_, graph_2_):
    G = nx.Graph()
    G.add_edges_from(graph_1_.edges())
    G.add_edges_from(graph_2_.edges())
    G.add_edges_from(product_graph_.nodes())
    return G



def score_ph4(scores_, pairs_, ph4_lig_, ph4_prot_):
    score = 0
    ph4_associations = []
    for idx_1, idx_2 in pairs_:
        score += scores_[(ph4_lig_[idx_1], ph4_prot_[idx_2])]
        ph4_associations.append((ph4_lig_[idx_1], ph4_prot_[idx_2]))
    # score = round(score/len(pairs), 2)
    return score, ph4_associations





def score_pl_clashes(coords_protein_, protein_atoms_, trans_coords_ligand_,
                                                ligand_atoms_, allow_=0.5):
    neigh = NearestNeighbors(n_neighbors=len(coords_protein_), radius=(2+2), # max atom vdw radii
                               algorithm='ball_tree').fit(coords_protein_)
    distances, indices = neigh.radius_neighbors(trans_coords_ligand_)

    clashes = []
    for i in range(len(trans_coords_ligand_)):
        if len(indices[i]) > 0:
            for j, prot_idx in enumerate(indices[i]):
                try:
                    total_vdw_radius = vdw_radius[protein_atoms_[prot_idx]] + vdw_radius[ligand_atoms_[i]]
                except KeyError:
                    total_vdw_radius = 1.2 + 1.2 # hydrogen
                if distances[i][j] + allow_ < total_vdw_radius:
                    #print(distances[i][j])
                    clashes.append(i)
    count = len(set(clashes))
    return count

            



def score_ligph4_coverage(trans_coords_lig_, coords_cav_):

    neigh = NearestNeighbors(n_neighbors=len(coords_cav_), 
                        radius=2, algorithm='ball_tree').fit(coords_cav_)
    distances, indices = neigh.radius_neighbors(trans_coords_lig_)
    count = 0
    for n in indices:
        if np.array(n).size > 0:
            count += 1
    return round(count/np.shape(trans_coords_lig_)[0]*100, 2)




def score_hetero_satisfaction(trans_coords_ligph4_, coords_cav_, features_lig_, features_cav_, radius_=2):

    neigh = NearestNeighbors(n_neighbors=len(coords_cav_), radius=radius_, 
                               algorithm='ball_tree').fit(coords_cav_)
    distances, indices = neigh.radius_neighbors(trans_coords_ligph4_)

    hetero = ['OG', 'O', 'OD1', 'N', 'NZ']

    hetero_total = 0
    hetero_satis = 0

    for i, cav_idx in enumerate(indices):
        if features_lig_[i] not in hetero:
            continue
        else:
            hetero_total += 1

        if features_lig_[i] in [features_cav_[j] for j in cav_idx]:
            hetero_satis += 1

    score = hetero_satis/hetero_total
    return score




def reject_pairing_dist(pairs_, coords_lig_, coords_cav_, distance_threshold_=1.5):
    pairing_distances = [np.linalg.norm(np.array(coords_lig_[i])-np.array(coords_cav_[j]))
                            for i, j in pairs_]
    reject = False
    if any([d > distance_threshold_ for d in pairing_distances]):
        reject = True
    return reject




def assign_weights(pairs_, ph4_cav_):
    weights = []
    for l, p in pairs_:
        if ph4_cav_[p] in ['O', 'N', 'OD1', 'OG', 'NZ', 'CZ']:
            weights.append(1)
        else:
            weights.append(0.3)

    return np.array(weights)




def estimate_transformation(pairs_, coords_lig_, coords_cav_, weights_,):

    # align #1 to #2
    corres_1 = np.array([coords_lig_[p[0]] for p in pairs_])
    corres_2 = np.array([coords_cav_[p[1]] for p in pairs_])

    center_1 = [np.mean(corres_1[:, i]) for i in range(3)]
    center_2 = [np.mean(corres_2[:, i]) for i in range(3)]
    translation = np.array(center_2) - np.array(center_1)
    #print(translation)
    trans_corres_1 = corres_1 + translation
    trans_coords_lig = coords_lig_ + translation

    if len(pairs_) > 2:
        rotation, rmse = Rotation.align_vectors(corres_2, trans_corres_1, weights=weights)
        trans_coords_lig = rotation.apply(trans_coords_lig)
        #print(rotation.as_matrix())
    else:
        # by default, a non-identity matrix. There is an infinity of solutions
        # the cases should be treated separately
        rotation = Rotation.identity()
        rmse = np.nan

    return trans_coords_lig, rotation, translation, rmse





def apply_transformation(rotation_, translation_, coords_):
    coords = coords_ + translation_
    trans_coords = rotation_.apply(coords)
    return trans_coords





def align(rotation_, translation_, mol2_, ofile_, sample_=None):
    with open(mol2_, 'r') as f:
        mol2 = f.read()
        molecule_area = mol2.split('@<TRIPOS>MOLECULE\n')[1].split('@<TRIPOS>')[0]
        molecule_area = "@<TRIPOS>MOLECULE\n" + molecule_area
        atom_area = mol2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0]
        atoms = atom_area.split('\n')
        #print(atoms)
        del atoms[-1]
        new_atom_area = "@<TRIPOS>ATOM"
        idx = 0
        for i, l in enumerate(atoms):
            if sample_ is not None and i not in sample_:
                continue
            idx += 1
            cols = l.split()
            coords = np.array([float(cols[2]), float(cols[3]), float(cols[4])])
            #print(coords)
            trans_coords = coords + translation_
            trans_coords = rotation_.apply(trans_coords)
            #print(trans_coords)
            new_atom_area += ("\n{:>7} {:<8} {:>9.4f} {:>9.4f} {:>9.4f} "
                              "{:<5} {:>5} {:<8} {:>9}".format(idx, #cols[0]
                                                     cols[1],
                                                     trans_coords[0],
                                                     trans_coords[1],
                                                     trans_coords[2],
                                                     cols[5],
                                                     cols[6],
                                                     cols[7],
                                                     cols[8]
                                                    ))
        new_atom_area += "\n"
        try:
            final_area = mol2.split('@<TRIPOS>BOND\n')[1]
            final_area = "@<TRIPOS>BOND\n" + final_area
        except IndexError:
            final_area = ""


    if ofile_[-5:] != '.mol2':
            ofile_ += '.mol2'
    basename = os.path.basename(ofile_)
    ofile_ = basename
    name = os.path.splitext(basename)[0]

    header = ""
    header += "# Modified by ProCare\n"
    header += "# Modification time: {}\n".format(
                            strftime("%a %d %b %Y %H:%M:%S", localtime()))
    header += "# Name: {}.mol2\n\n".format(name)

    with open(ofile_, 'w') as of:
        of.write('{}{}{}{}'.format(header, molecule_area, new_atom_area, final_area))



def write_report(ofile_, score_type_, poses):
    with open(ofile_, 'w') as of:
        of.write('pose\t{}\n'.format(score_type_))
        for name, score in poses:
            of.write('{}\t{}\n'.format(*poses))


if __name__ == '__main__':



    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cavity', type=str, required=True,
                help='cavity mol2')
    parser.add_argument('-p', '--protein', type=str, required=True,
                help='protein mol2')
    parser.add_argument('-l', '--ligph4', type=str, required=True,
                help='ligph4 mol2')
    parser.add_argument('-lig', '--ligand', type=str, required=True,
                help='ligand mol2')
    parser.add_argument('-r', '--ranking', type=str, required=False,
                help='ranking method, size: clique size or rmse: RMSE, both: combination of both',
                choices=['size', 'rmse', 'both'], default='rmse')
    parser.add_argument('-n', '--npose', type=int, required=False,
                help='number of poses', default=20)
    args = parser.parse_args()





    # extract coordinates and atom types:
    coords_lig, features_lig = iomol2.read_cav_mol2(args.ligph4)
    coords_cav, features_cav = iomol2.read_cav_mol2(args.cavity)
    coords_protein, features_protein = iomol2.read_mol2_features(args.protein)
    coords_ligand, features_ligand = iomol2.read_mol2_features(args.ligand)
    #print(coords_ligand.shape)


    # compute pairwise distances with Scipy --> pdist = upper triangular matrix
    pdist_lig = distance.pdist(coords_lig)
    pdist_cav = distance.pdist(coords_cav)


    # rewrite pairwise distances as dict 
    pairs_lig = pairwise_distances_as_dict(coords_lig, features_lig)
    pairs_cav = pairwise_distances_as_dict(coords_cav, features_cav)
    #print(pairs_lig)
    #print(pairs_cav)

    # create grah product
    product_edges, scores = associate_pairs(pairs_lig, pairs_cav)
    #print(product_edges)
    #print(scores)

    nx_subgraph_matching = nx_graph_from_edges(product_edges)
    adj_product = nx_to_adj_graph(nx_subgraph_matching)
    #print(adj_product)

    # search cliques
    cliques = list(find_cliques(adj_product))

    if len(cliques) > 0 and max([len(c) for c in cliques]) == 2:
        cliques = [tuple(c) for c in cliques if len(c) >= 2]
    else:
        cliques = [tuple(c) for c in cliques if len(c) >= 4]

   
    print(len(cliques))
    #print(cliques[0])





    solutions = []
    rejected = 0

    if len(cliques) > 0:
        #print("N sol = ", len(cliques))

        for s in cliques:
            weights = assign_weights(s, features_cav)
            trans_coords_ligph4, rotation, translation, clique_rmse = estimate_transformation(s, coords_lig, coords_cav, weights)
            if clique_rmse > 3:
                continue
            if reject_pairing_dist(s, trans_coords_ligph4, coords_cav, 2):
                continue
            cov = score_ligph4_coverage(trans_coords_ligph4, coords_cav)
            if cov < 50:
                continue

            #hetero_score = score_hetero_satisfaction(trans_coords_ligph4, coords_cav, features_lig, features_cav, radius_=2)
            #print(hetero_score)
            #if hetero_score == 0:
            #    continue

            trans_coords_ligand = apply_transformation(rotation, translation, coords_ligand)
            #print(trans_coords_ligand)
            #print(trans_coords_ligph4)
            clashes = score_pl_clashes(coords_protein, features_protein, trans_coords_ligand,
                                                features_ligand, allow_=0.5)

            if clashes > 30:
                continue


            score, ph4_associations = score_ph4(ph4_scores, s, features_lig, features_cav)

            solutions.append((s, clique_rmse, score, len(s)/clique_rmse, ph4_associations, 
                                    trans_coords_ligph4, cov, rotation, translation, len(s)))

           


    if args.ranking == 'rmse':
        sorted_solutions = sorted(solutions, key=operator.itemgetter(1))
    elif args.ranking == 'size':
        sorted_solutions = sorted(solutions, key=operator.itemgetter(-1), reverse=True)
    elif args.ranking == 'both':
        sorted_solutions = sorted(solutions, key=operator.itemgetter(3), reverse=True)



    i = 1
    selected_poses = []
    for rank, n in enumerate(sorted_solutions[:args.npose]):

        ligand_name = os.path.basename(args.ligand).split('.')[0]
        ligph4_name = os.path.basename(args.ligph4).split('.')[0]
        oligand = 'rot_{}_{}.mol2'.format(ligand_name, i)
        oligph4 = 'rot_{}_{}.mol2'.format(ligph4_name, i)



        align(n[-3], n[-2], args.ligand, oligand)

        align(n[-3], n[-2], args.ligph4, oligph4)

        i += 1

        if args.ranking == 'rmse':
            selected_poses.append((oligand, round(n[1], 4)))
        elif selected_poses.ranking == 'size':
            selected_poses.append((oligand, round(n[-1, 4])))
        elif args.ranking == 'both':
            selected_poses.append((oligand, round(n[3], 4)))


    write_report(args.ligand.replace('mol2', 'score'), args.ranking, selected_poses)



    print('file:', args.ligand, ' | # cliques processed:', len(cliques),
             ' | # solutions ranked:', len(solutions), ' | # poses:', len(selected_poses))
    
