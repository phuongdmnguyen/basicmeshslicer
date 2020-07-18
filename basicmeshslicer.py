import numpy as np
import argparse
from stl import mesh

# Construct the argument parser
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True,
                help="path to the input image")
ap.add_argument("-p", "--planes", required=True, help="path to .txt files containing plane info")
# ap.add_argument("-o", "--output", required=True, nargs=2,
#     help="path to the output image")
args = vars(ap.parse_args())

# load the .stl file
prism = mesh.Mesh.from_file(args["input"])

# prism = mesh.Mesh.from_file("Prism.stl")

# maximum vertices of prism
x_upper = prism.max_[0]
y_upper = prism.max_[1]
z_upper = prism.max_[2]

# minimum vertices of prism
x_lower = prism.min_[0]
y_lower = prism.min_[1]
z_lower = prism.min_[2]

# These are the vertices of the front and back faces (the faces that are parallel to the y-z plane)
front_vertices = np.array([[x_upper, y_upper, z_upper],
                           [x_upper, y_upper, z_lower],
                           [x_lower, y_upper, z_upper],
                           [x_lower, y_upper, z_lower]])
back_vertices = np.array([[x_upper, y_lower, z_upper],
                          [x_upper, y_lower, z_lower],
                          [x_lower, y_lower, z_upper],
                          [x_lower, y_lower, z_lower]])


# arguments: a point on the plane, and a normal vector passed each as a tuple of 3
# returns: the corresponding equation of the plane Ax + By + Cz = D as a tuple: (A, B, C, D)
def calculate_plane_equation(point_on_plane, normal):
    x = (-1 * point_on_plane[0]) * normal[0]
    y = (-1 * point_on_plane[1]) * normal[1]
    z = (-1 * point_on_plane[2]) * normal[2]
    sum_all = x + y + z
    equation = (normal[0], normal[1], normal[2], (sum_all * -1))
    return equation


# load plane information from .txt
# initialize plane information as a list with plane[0] = points and plane[1] = normal
planes = []

with open(args["planes"], 'r') as f:
    normal = []
    points = []
    temp_normal = f.readline().split()
    for x in temp_normal:
        y = float(x)
        normal.append(y)
    temp_points = f.readline().split()
    for x in temp_points:
        y = float(x)
        points.append(y)
    planes.append((normal, points))


# Given in .txt files
# plane1_normal = (0.094168, -0.993819, 0.0587916)
# plane1_point = (-3.35276126862e-08, -0.00140627007931, -0.0237716548145)
# plane1_equation = calculate_plane_equation(plane1_point, plane1_normal)

# plane2_normal = (0.203912, 0.97879, -0.0197536)
# plane2_point = (0.00557016581297, -8.00095510483, 0.0101777203381)
# plane2_equation = calculate_plane_equation(plane2_point, plane2_normal)


# arguments: equation of the plane, max and min vertex of prism (the a 3 element tuple - (x, y, z))
# returns: 4 vertices at intersection between the plane and prism (2 front and 2 back)
#       1. (x_max, y_a, z_max)  -> front_top
#       2. (x_max, y_b, z_min)  -> front_bottom
#       3. (x_min, y_c, z_max)  -> back_top
#       4. (x_min, y_d, z_min)  -> back_bottom

def get_vertices_at_intersection(plane_equation, max_vertex, min_vertex):
    # plane equation : Ax_0 + By_0 + Cz_0 = D
    # since we know x and z and A, B, C, D -> we can find y we can rearrange to find y
    # y_0 = (D - Ax_0 - Cz_0) / B
    # here we are calculating (D - Ax_0 - Cz_0)
    y_a_numerator = plane_equation[3] - (plane_equation[0] * max_vertex[0]) - (plane_equation[2] * max_vertex[2])
    # then divide known by the numerator by B
    y_a = y_a_numerator // plane_equation[1]

    # Similarly for y_b
    y_b_numerator = plane_equation[3] - (plane_equation[0] * max_vertex[0]) - (plane_equation[2] * min_vertex[2])
    y_b = y_b_numerator // plane_equation[1]

    y_c_numerator = plane_equation[3] - (plane_equation[0] * min_vertex[0]) - (plane_equation[2] * max_vertex[2])
    y_c = y_c_numerator // plane_equation[1]

    y_d_numerator = plane_equation[3] - (plane_equation[0] * min_vertex[0]) - (plane_equation[2] * min_vertex[2])
    y_d = y_d_numerator // plane_equation[1]

    front_top = np.array([max_vertex[0], y_a, max_vertex[2]])
    front_bottom = np.array([max_vertex[0], y_b, min_vertex[2]])
    back_top = np.array([min_vertex[0], y_c, max_vertex[2]])
    back_bottom = np.array([min_vertex[0], y_d, min_vertex[2]])
    return front_top, front_bottom, back_top, back_bottom


# get plane equations
plane_eqs = []
for plane in planes:
    eq = calculate_plane_equation(plane[1], plane[0])
    plane_eqs.append(eq)

# the vertices where plane intersects the prism
ftop_p1, fbot_p1, btop_p1, bbot_p1 = get_vertices_at_intersection(plane_eqs[0], prism.max_, prism.min_)
plane1_vertices = np.array([ftop_p1, fbot_p1, btop_p1, bbot_p1])

# # the vertices where plane 2 intersects the prism
# ftop_p2, fbot_p2, btop_p2, bbot_p2 = get_vertices_at_intersection(plane2_equation, prism.max_, prism.min_)
# plane2_vertices = np.array([ftop_p2, fbot_p2, btop_p2, bbot_p2])

# Note: y_a, y_b, y_c, y_d of plane 2 is -9.0, -9.0, -8.0 -8.0 (obtained from printing out the above values)
#            "                plane 1 is 1, 0, -1, -1
# Therefore, plane 2 is closer to the back face than plane 1


# arrays of vertices for the 3 sliced pieces
# since the plane 2 is close to the back face, the back piece will contain vertices where they intersect along with the
# vertices of the back face
# back_piece_vertices = np.concatenate((back_vertices, plane2_vertices), 0)
# print(back_piece_vertices)

# since the plane 1 is close to the front face, the front piece will contain vertices where they intersect along with
# the vertices of the front face
front_piece_vertices = np.concatenate((plane1_vertices, front_vertices), 0)

back_piece_vertices = np.concatenate((back_vertices, plane1_vertices), 0)

# The middle piece will contain vertices at the points where the two planes intersect the prism
# mid_piece_vertices = np.concatenate((plane2_vertices, plane1_vertices), 0)


# these denote the indices from the vertices that make up the 12 triangles of the pieces
# The standard convention was as follows:
#   2------0        6------4
#   |      |        |      |
#   |      |        |      |
#   3------1        7------5
# vertices numbers 0, 1, 2, 3 is for faces that are more negative in y (back piece)
# vertices numbers 5, 6, 7, 8 is for faces that are more positive in y (front piece)
# triangles are built so that each triangle share 2 vertices with one another
faces = np.array([
    [2, 0, 1],
    [2, 3, 1],
    [0, 4, 1],
    [4, 1, 5],
    [6, 7, 4],
    [4, 5, 7],
    [2, 3, 6],
    [6, 7, 3],
    [6, 4, 2],
    [2, 0, 4],
    [3, 1, 5],
    [3, 7, 5]])


# arguments: using the vertices info for each piece and the faces and the mesh obj with no data
# returns:   a mesh object containing the points of the triangle
def build_mesh_data(vertices_info, faces_info, mesh_obj):
    for i, f in enumerate(faces_info):
        for j in range(3):
            mesh_obj.vectors[i][j] = vertices_info[f[j], :]


# initializing empty mesh objects with mesh data (where the vertices are initialized as 0)
back_piece = mesh.Mesh(np.zeros(faces.shape[0], mesh.Mesh.dtype))
# mid_piece = mesh.Mesh(np.zeros(faces.shape[0], mesh.Mesh.dtype))
front_piece = mesh.Mesh(np.zeros(faces.shape[0], mesh.Mesh.dtype))

# build data of each piece using their information
build_mesh_data(back_piece_vertices, faces, back_piece)
build_mesh_data(front_piece_vertices, faces, front_piece)
# build_mesh_data(mid_piece_vertices, faces, mid_piece)

# updating normal vectors for each vertices in each piece
back_piece.update_normals()
front_piece.update_normals()
# mid_piece.update_normals()

# save each piece as an stl
back_piece.save("backpiece.stl")
# mid_piece.save("midpiece.stl")
front_piece.save("frontpiece.stl")

