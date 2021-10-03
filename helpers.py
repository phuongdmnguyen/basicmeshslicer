import numpy as np


def get_plane_from_text(filename):
    # load plane information from .txt
    # initialize plane information as a list with planes[0] = normal and planes[1] = points
    with open(filename, 'r') as f:
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
    return [normal, points]


# arguments: a point on the plane, and a normal vector passed each as a tuple of 3
# returns: the corresponding equation of the plane Ax + By + Cz = D as a tuple: (A, B, C, D)
def calculate_plane_equation(point_on_plane, normal):
    x = (-1 * point_on_plane[0]) * normal[0]
    y = (-1 * point_on_plane[1]) * normal[1]
    z = (-1 * point_on_plane[2]) * normal[2]
    sum_all = x + y + z
    equation = (normal[0], normal[1], normal[2], (sum_all * -1))
    return equation


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


# arguments: using the vertices info for each piece and the faces and the mesh obj with no data
# returns:   a mesh object containing the points of the triangle
def build_mesh_data(vertices_info, faces_info, mesh_obj):
    for i, f in enumerate(faces_info):
        for j in range(3):
            mesh_obj.vectors[i][j] = vertices_info[f[j], :]
