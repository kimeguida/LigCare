import numpy as np
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, required=True, help="Input mol2")
parser.add_argument('-o', '--output', type=str, required=True, help="Output mol2")
args = parser.parse_args()


# Rotation along x axis = move y towards z axis, here tetha=pi
# [1     0               0         ]
# [0     cos(tetha)     -sin(tetha)]
# [1     sin(tetha)      cos(tetha)]

"""
matrix_Rx = np.array([[1, 0, 0],
                      [0, 0.5, -np.sqrt(3)/2],
                      [0, np.sqrt(3)/2, 1/2]])
"""
matrix_Rx = np.array([[1, 0, 0],
                      [0, -1, 0],
                      [0, 0, -1]])


vector_T = np.array([20,
                     20,
                     20])


"""
vector_T = np.array([0,
                     0,
                     10])
"""
vector_T.reshape(1, 3)
#print(vector_T)

def transform(coords_, matrix_R_=matrix_Rx, vector_t_=vector_T):

    RT_coord = []


    centroid = np.array([np.mean([p[0] for p in coords_]),
                        np.mean([p[1] for p in coords_]),
                        np.mean([p[2] for p in coords_]),])
                
    #print(centroid)


    rotated_coords = np.matmul(matrix_R_, np.transpose(coords_))
    rotated_coords = np.transpose(rotated_coords)
    rot_centroid = np.array([np.mean([p[0] for p in rotated_coords]),
                            np.mean([p[1] for p in rotated_coords]),
                            np.mean([p[2] for p in rotated_coords]),])

    rotated_coords = rotated_coords - (rot_centroid-centroid)
    trans_coords = rotated_coords + vector_t_
    #print(trans_coords)


    """
    for point in coord_:
        rotated_point = np.matmul(matrix_R_, point)
        rotated_point_centered

        translated_point = np.array(point) + vector_t_
        RT_coord.append(rotated_point)
    """


    """
    R_point = np.matmul(matrix_R_, point-centroid)
    print(point, R_point)
    RT_point = R_point + centroid + np.array(vector_T)
    RT_coord.append(RT_point)
    """

    return trans_coords


## get coordinates
with open(args.input, "r") as crdfile:
    mol2 = crdfile.read().split("\n")

start = mol2.index("@<TRIPOS>ATOM") + 1
end = mol2.index("@<TRIPOS>BOND") - 1

crd = []
ALL_Atom = []
for i in range(start, end+1):
    cols = mol2[i].split()
    x = float(cols[2])
    y = float(cols[3])
    z = float(cols[4])

    crd.append([x,y,z])
    ALL_Atom.append(cols)

RT_coord = transform(crd)

with open(args.output, "w") as flip_mol2:
    for i in range(start):
        flip_mol2.write(mol2[i]+"\n")

    for i, point in enumerate(RT_coord):
        flip_mol2.write("{:>7} {:<8} {:>9.4f} {:>9.4f} {:>9.4f} {:<5} {:>5} {:<8} {:>9}\n".format(ALL_Atom[i][0],
                                                                                                  ALL_Atom[i][1],
                                                                                                  point[0],
                                                                                                  point[1],
                                                                                                  point[2],
                                                                                                  ALL_Atom[i][5],
                                                                                                  ALL_Atom[i][6],
                                                                                                  ALL_Atom[i][7],
                                                                                                  ALL_Atom[i][8]))
    for i in range(end+1, len(mol2)):
        flip_mol2.write(mol2[i]+"\n")
