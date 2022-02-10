import meshio

mesh_filename = './tetra.msh'

mesh = meshio.read(
    mesh_filename,  # string, os.PathLike, or a buffer/open file
    file_format="gmsh",  # optional if filename is a path; inferred from extension
)

print(mesh)

meshio.write(
    'tetra.stl',  # string, os.PathLike, or a buffer/open file
    mesh,
    file_format="stl",  # optional if filename is a path; inferred from extension
)


