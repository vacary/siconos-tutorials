import numpy


file_number='00498'
file_number='00246'


idx = numpy.loadtxt('colliding_debug/output_triangle_i_'+ file_number+'.py', dtype ='int')
print(idx)

v = numpy.loadtxt('colliding_debug/output_triangle_v_'+file_number+'.py')
#print(v)

n_vertex = int(numpy.max(v[:,0]))+1
print(n_vertex)
v_sorted = numpy.zeros((n_vertex,3))

for i in range(v.shape[0]):
    v_sorted[int(v[i,0]), :] = v [i, 1:]  

print('v_sorted', v_sorted)

import meshio

print('idx', idx)

points = v_sorted
cells =[("triangle", idx)]

print('points', points)
print('cells', cells)

mesh = meshio.Mesh(
    points,
    cells,
)
mesh.write(
    'foo'+file_number+'.vtk',  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)

