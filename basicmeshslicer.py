import numpy as np
import argparse
from stl import mesh
import helpers

# CONSTANTS
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


# load the .stl file
def main(args):
    prism = mesh.Mesh.from_file(args["input"])

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

    # load plane information from .txt
    plane = helpers.get_plane_from_text(args["planes"])
    plane_eq = helpers.calculate_plane_equation(plane[1], plane[0])


    # the vertices where plane intersects the prism
    ftop_p1, fbot_p1, btop_p1, bbot_p1 = helpers.get_vertices_at_intersection(plane_eq, prism.max_, prism.min_)
    plane1_vertices = np.array([ftop_p1, fbot_p1, btop_p1, bbot_p1])

    front_piece_vertices = np.concatenate((plane1_vertices, front_vertices), 0)
    back_piece_vertices = np.concatenate((back_vertices, plane1_vertices), 0)

    # initializing empty mesh objects with mesh data (where the vertices are initialized as 0)
    back_piece = mesh.Mesh(np.zeros(faces.shape[0], mesh.Mesh.dtype))
    front_piece = mesh.Mesh(np.zeros(faces.shape[0], mesh.Mesh.dtype))

    # build data of each piece using their information
    helpers.build_mesh_data(back_piece_vertices, faces, back_piece)
    helpers.build_mesh_data(front_piece_vertices, faces, front_piece)

    # updating normal vectors for each vertices in each piece
    back_piece.update_normals()
    front_piece.update_normals()

    # save each piece as an stl
    back_piece.save("backpiece.stl")
    front_piece.save("frontpiece.stl")


if __name__ == '__main__':
    # Construct the argument parser
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True,
                    help="path to the input image")
    ap.add_argument("-p", "--planes", required=True, help="path to .txt files containing plane info")
    # ap.add_argument("-o", "--output", required=True, nargs=2,
    #     help="path to the output image")
    args = vars(ap.parse_args())
    main(args)
