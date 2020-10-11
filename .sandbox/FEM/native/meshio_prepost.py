import meshio

filename = 'beam'
def create_Meshcpp_for_siconos(mesh):

    f = open(filename+".cpp", "w")
    f.write("#include <Mesh.hpp>\n");

    f.write("static Mesh * createMesh()\n{\n");

    p_cnt =0
    f.write("std::vector<MVertex *> vertices;\n");
    
    for p in mesh.points:
        f.write("vertices.push_back(new MVertex({0}, {1}, {2}, {3}));\n".format(p_cnt,p[0], p[1], p[2]))
        p_cnt = p_cnt+1


    f.write("std::vector<MElement *> elements;\n");
    e_cnt=0
    dim=0
    for mc in mesh.cells:
        print("mc.type", mc.type)
        for me in mc[1]:
            print("me nodes:" ,me)
            str_v ="{"
            for p in me:
                str_v += "vertices[{0}],".format(p) 
                
                p_cnt= p_cnt+1
            str_v += "}" 
            print(str_v)
            f.write("std::vector<MVertex *> vertices{0} = {1};\n".format(e_cnt, str_v));
            
            if mc.type == 'triangle':
                type = 2
                dim= max(dim,2)
            elif mc.type == 'line':
                type = 1
                dim= max(dim,1)
            elif mc.type == 'vertex':
                type=15
                dim= max(dim,0)
            elif  mc.type == 'tetra':
                type = 4
                dim= max(dim,3)
            f.write("elements.push_back(new MElement({0}, {1}, vertices{2}));\t".format(e_cnt, type, e_cnt))

            e_cnt=e_cnt+1

    f.write("Mesh * m =    new Mesh({0}, vertices, elements);\n".format(dim));

    f.write("return m;\n }; \n ");
    f.close()


mesh = meshio.read(
    filename+'.msh',  # string, os.PathLike, or a buffer/open file
    file_format="gmsh",  # optional if filename is a path; inferred from extension
)

import numpy as np
print(mesh.points, mesh.cells, mesh.cells_dict)

# ------------------------------------- vtk output #
meshio.write(
    "cube.vtk",  # str, os.PathLike, or buffer/ open file
    mesh,
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)

# ------------------------------------- gmsh v2 output #
meshio.write(
    "cube.msh2",  # str, os.PathLike, or buffer/ open file
    mesh,
    file_format="gmsh22",  # optional if first argument is a path; inferred from extension
    binary=False,
    float_fmt='.16e'
)

#------------------------------------- siconos output #
create_Meshcpp_for_siconos(mesh)

#------------------------------------- post processing for paraview"
from displacement import *
import os
import os.path
from os import path

if  not path.exists('vtk'):
    os.mkdir('vtk')
n_samples = len(x)

for i in range(n_samples):
    point_data={}
    point_data['u_x']= x[i]
    point_data['u_y']= y[i]
    point_data['u_z']= z[i]
    point_data['u']= np.column_stack((x[i], y[i], z[i]))

    #print('point_data[u]', point_data['u'])

    foutput = '{0}{1:03d}.vtk'.format(filename,i)
    print(foutput)
    meshio.write_points_cells(
        './vtk/'+foutput,
        points=mesh.points,
        cells=mesh.cells,
        # Optionally provide extra data on points, cells, etc.
        point_data=point_data
        # cell_data=cell_data,
        # field_data=field_data
    )

