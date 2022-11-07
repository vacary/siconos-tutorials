import numpy as np
import meshio

def export_vtk_mesh_from_manifold(manifold):


    #print(vertices)
    #print(triangles)

    ### get triangles from original mesh
    #print(manifold['m_body0']['data']['mesh_array']['triangle'])
    triangle = np.array(manifold['m_body0']['data']['mesh_array']['triangle'])
    vertices =  np.array(manifold['m_body0']['data']['mesh_array']['vertices'])

    idx = triangle[:, 1:4]
    points = vertices[:,1:4]
    cells =[("triangle", idx)]

    #print('points', points)
    #print('cells', cells)

    mesh = meshio.Mesh(
        points,
        cells,
    )

    mesh.write(
        'colliding_debug/colliding_results_mesh.vtk',  # str, os.PathLike, or buffer/open file
        # file_format="vtk",  # optional if first argument is a path; inferred from extension
    )

def export_vtk_contact_points_from_manifold(manifold, file_number):
    #print(manifold['contact_points'])
    points = []
    idx= []
    k=0
    point_data={}
    point_data['body']= []
    for cp in manifold['contact_points']:
        points.append(cp['pt_A'])
        point_data['body'].append(0)



        idx.append([k])
        k=k+1
        
    for cp in manifold['contact_points']:
        points.append(cp['pt_B'])
        point_data['body'].append(1)
        idx.append([k])
        k=k+1

    cells = [("vertex",idx)]    
    
    #print('points', points)
    #print('cells', cells)

    
    mesh = meshio.Mesh(
        points,
        cells,
    )
    
    meshio.write_points_cells(
        'colliding_debug/colliding_results_contact_points_'+file_number+'.vtk',  # str, os.PathLike, or buffer/open file
        points=mesh.points,
        cells=mesh.cells,
        # Optionally provide extra data on points, cells, etc.
        point_data=point_data
        # cell_data=cell_data,
        # field_data=field_data
    )


    
def export_vtk_colliding_results(vertices, triangles, contact_points, file_number):
    
    points = []
    idx= []

    k=0
    point_data={}
    point_data['triangle']= []
    point_data['body']= []
    
    #print(contact_points)
    points = []
    idx= []
    k=0
    point_data={}
    point_data['body']= []
    for cp in contact_points:
        points.append(cp['pt_A'])
        point_data['body'].append(0)
        idx.append([k])
        k=k+1
        
    for cp in contact_points:
        points.append(cp['pt_B'])
        point_data['body'].append(1)
        idx.append([k])
        k=k+1
    cells = [("vertex",idx)]    
    
    #print('points', points)
    #print('point_data', point_data)
    #print('cells', cells)
    

    
    mesh = meshio.Mesh(
        points,
        cells,
    )
    
    
    meshio.write_points_cells(
        'colliding_debug/colliding_results_all_contact_points_'+file_number+'.vtk',  # str, os.PathLike, or buffer/open file
        points=mesh.points,
        cells=mesh.cells,
        # Optionally provide extra data on points, cells, etc.
        point_data=point_data
        # cell_data=cell_data,
        # field_data=field_data
    )




    
    pass


file_number='00498'
file_number='00246'
filename = 'colliding_debug/colliding_results_' + file_number + '.py'
exec(open(filename).read())
#print(manifold)
#export_vtk_mesh_from_manifold(manifold)
#export_vtk_contact_points_from_manifold(manifold,file_number)
          
#export_vtk_colliding_results(vertices, triangles, contact_points)




import glob
import os
cnt_file =0
for name in glob.glob('colliding_debug/colliding_results_*.py'):
    print('\n ######################################## ')
    print('vtk output of ', name)
    filename = os.path.basename(name).split('.')[0]
#    print(filename)
    file_number = filename.split('_')[-1]
#    print(file_number)
    exec(open(name).read())
    export_vtk_contact_points_from_manifold(manifold,file_number)
    export_vtk_colliding_results(vertices, triangles, contact_points, file_number)
    cnt_file = cnt_file+1
    #if cnt_file > 10:
    #    break
    
    
