import meshio

mesh_filename = './cube_fine_128.msh'

import os

print(os.path.basename(mesh_filename))

mesh_filename_base= os.path.basename(mesh_filename).split('.')[0]

mesh = meshio.read(
    mesh_filename,  # string, os.PathLike, or a buffer/open file
    file_format="gmsh",  # optional if filename is a path; inferred from extension
)

print(mesh)

meshio.write(
    mesh_filename_base+'.stl',  # string, os.PathLike, or a buffer/open file
    mesh,
    file_format="stl",  # optional if filename is a path; inferred from extension
)


