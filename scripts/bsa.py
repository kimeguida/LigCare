

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
import copy
from sklearn.neighbors import NearestNeighbors
import iomol2




def init_projections(output_mol2_=False, proj_r_=8, angle_step_=np.pi/8):

    """
    initialize points at the surface of a sphere of radius proj_r
    projections are rays subdividing the 
    """

    proj_r = proj_r_

    pi = np.pi
    angle_step = angle_step_
    tetha_angles = np.arange(0, 2*pi, angle_step)
    phi_angles = np.arange(-pi/2, pi/2, angle_step)

    exluded_axis_angles = {0, pi/2, pi, 3*pi/2}
    exluded_phi_angles = {-pi/2, pi/2,}

    #### x = r.cos(phi).cos(tetha)
    #### y = r.cos(phi).sin(tetha)
    #### z = r.sin(phi)

    proj_points_origin = []    
    #set_proj = []

    for tetha in tetha_angles:
        for phi in phi_angles:
            if phi in exluded_phi_angles:
                continue
            
            x = proj_r * np.cos(phi) * np.cos(tetha)
            y = proj_r * np.cos(phi) * np.sin(tetha)
            z = proj_r * np.sin(phi)
            proj_points_origin.append(np.array([x, y, z]))
            #set_proj.append(( str(round(x, 4)), str(round(y, 4)), str(round(z, 4)) ))

    ##### add missing axis ####
    # X
    #proj_points_origin.append(proj_r * np.array([1, 0, 0]))
    #proj_points_origin.append(proj_r * np.array([-1, 0, 0]))
    
    # Y
    #proj_points_origin.append(proj_r * np.array([0, 1, 0]))
    #proj_points_origin.append(proj_r * np.array([0, -1, 0]))
    
    # Z
    proj_points_origin.append(proj_r * np.array([0, 0, 1]))
    proj_points_origin.append(proj_r * np.array([0, 0, -1]))


    # check duplicates of projections
    #print(len(set(set_proj)))
    #print(len(set_proj))
    #print('proj', len(proj_points_origin))


    if output_mol2_:
        out_proj_points_origin = copy.deepcopy(proj_points_origin)
        #out_proj_points_origin.append(proj_r * np.array([0, 0, 0]))
        iomol2.write_mol2('proj_points.mol2', out_proj_points_origin)
    
    return proj_points_origin





def reduce_prot_coords(cav_coords_, prot_coords_, d_=9):

    """ 
    nearest proteins atoms within radius d_
    that can be reached by projections
    d_max = 8.14
    """
    
    neigh = NearestNeighbors(n_neighbors=len(prot_coords_), radius=d_,
                               algorithm='ball_tree').fit(prot_coords_)
    distances, indices = neigh.radius_neighbors(cav_coords_)
    return indices





def calc_bsa(cav_coords_, prot_coords_, proj_points_origin_,
                            reduced_prot_indices_, resolution_=1.51):
    origin = np.array([0, 0, 0])
    buriedness = {}
    for i, point in enumerate(cav_coords_):

        #### align the projection sphere on the point in focus
        translation = np.subtract(point, origin)
        proj_points = np.add(proj_points_origin_, translation)

        #### speed up by considering only the protein atoms whose voxels 
        #### can be reached by projections 
        prot_coords = [prot_coords_[j] for j in reduced_prot_indices_[i]]
        
        #### bsa = number of projection vectors reaching at least one
        #### protein voxel : the height of the projection of the protein
        #### atom onto the projection vector is less than resolution_
        bsa = 0
        for j, proj in enumerate(proj_points):
            vect_proj = np.subtract(proj, point)
            norm_vect_proj = np.linalg.norm(vect_proj)
            for atm in prot_coords:              
                vect_cav_prot = np.subtract(atm, point)
                scalar_prot_proj = np.dot(vect_cav_prot, vect_proj)
                #### protein atom is in same direction as projection
                if scalar_prot_proj >= 0:
                    #### check whether scalar projection is whithin radius proj_r
                    norm_scalar_prot_proj = scalar_prot_proj/norm_vect_proj
                    #print(norm_scalar_prot_proj)
                    if norm_scalar_prot_proj <= 8 :
                        norm_vect_cav_prot = np.linalg.norm(vect_cav_prot)
                        sin_angle = np.linalg.norm(np.cross(vect_proj, vect_cav_prot))/\
                                    (norm_vect_proj*norm_vect_cav_prot)
                        #### norm of rejection vector
                        h = norm_vect_cav_prot * sin_angle
                        if h <= resolution_:
                            bsa += 1
                            break
        #print(i, bsa)

        buriedness[i] = bsa

    return buriedness






if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cavity', type=str, required=True,
                help='volsite cavity')
    parser.add_argument('-p', '--protein', type=str, required=True,
                help='protein file')
    parser.add_argument('-o', '--output', type=str, required=False,
                help='bsa file')
    args = parser.parse_args()

    proj_points_origin = init_projections(True)
    cav_coords, _ = iomol2.read_cav_mol2(args.cavity)
    prot_coords = iomol2.read_prot_mol2(args.protein)
    reduced_prot_indices = reduce_prot_coords(cav_coords, prot_coords)

    #print(len(prot_coords))

    buriedness = calc_bsa(cav_coords, prot_coords, proj_points_origin, reduced_prot_indices)
    print(buriedness)

    if args.output is not None:
        with open(args.output, 'w') as of:
            of.write('mol2_index\tbsa\n')
            for i, bsa in buriedness.items():
                of.write(f'{i+1}\t{bsa}\n')




