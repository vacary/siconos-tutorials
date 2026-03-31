import sys, os
import numpy as np


from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
from siconos.mechanics.collision.convexhull import ConvexHull
#from convexhull_modif import ConvexHull

import siconos.numerics as sn

def normalize_shape(vertices,y_aspect_ratio,z_aspect_ratio,dest_vol):
    vertices= np.array(vertices) 
    siconos_convex_hull = ConvexHull(vertices)
    phch = siconos_convex_hull.hull # "PyHull Convex Hull"
    
    #1- Cleanup vertices (= remove concave points):
    #NOTE: phch.vertices are not really vertices, but rather a list of index of points (into self.__base_vertices__) that forms the faces of the convex shell : [[i_face1_pt1, i_face1_pt2, i_face3_pt3], [i_face2_pt1, i_face2_pt2, i_face2_pt3], ...]
    #So using numpy "flatten" then "unique" is a way to select only the points into self.__base_vertices__ that are involved in the convex shell (concave points are filtered out).

    vertices = phch.points[np.sort(np.unique(np.array(phch.vertices).flatten()))]
    
    #2- transform the shape to fit y_aspect_ratio and z_aspect_ratio:
    
    #        2.1- transform x coordinates to exactly fit [0, 1] interval.
    vertices[:,0] -= vertices[:,0].min()
    vertices[:,0] /= vertices[:,0].max()
    
    #        2.2- transform y coordinates to exactly fit [0, y_aspect_ratio] interval.
    vertices[:,1] -= vertices[:,1].min()
    vertices[:,1] /= vertices[:,1].max()
    vertices[:,1] *= y_aspect_ratio
    
    #        2.3- transform z coordinates to exactly fit [0, z_aspect_ratio] interval.
    vertices[:,2] -= vertices[:,2].min()
    vertices[:,2] /= vertices[:,2].max()
    vertices[:,2] *= z_aspect_ratio
    siconos_convex_hull = ConvexHull(vertices)#new convex hull, cleaned and with correct aspect ratio.
    
    #3- Move the shape center of mass to origin [0, 0, 0]:
    vertices -= siconos_convex_hull.centroid()
    siconos_convex_hull = ConvexHull(vertices) #center the shape
    
    ##4- Set volume
    base_inertia , start_volume = siconos_convex_hull.inertia([0., 0., 0.]) #assume density:=1 for now, so inertia will be recomputed later as it is proportionnal to density
    vol_factor     = dest_vol/start_volume
    coord_factor   = vol_factor**(1/3)
    vertices *= coord_factor # homothetic transform -> now vol=dest_vol
    
    return vertices


def generate_shape(nbPts,y_aspect_ratio,z_aspect_ratio,dest_vol):
       
    #np.random.seed(0)

    #NOTE: generate way too much points as pyhull (qhull) will have to ignore a lot of them to get a CONVEX hull. So at first we don't have a good control on nbPts...
    base_vertices = np.random.rand(nbPts*100,3)
    norm_vertices = normalize_shape(base_vertices,y_aspect_ratio,z_aspect_ratio,dest_vol)
    
    #NOTE: after normalize_shape, the shape is now convex, but we don't respect nbPts vertices, we hope to have more so we just have to remove some of them.
    nb_points_to_remove = len(norm_vertices) - nbPts
    
    '''
    print('nb_points_to_remove',nb_points_to_remove)
    print('norm_vertices.shape',norm_vertices.shape[0])
    '''
    if nb_points_to_remove < 0:
        Debug.error('In random shape generation, nb_points_to_remove was negative:',nb_points_to_remove)
        return
    
    del_index_vector = np.random.randint(norm_vertices.shape[0],size=nb_points_to_remove)   
    norm_vertices = np.delete(norm_vertices, del_index_vector, axis=0)

    # renormalize shape as we likely changed its volume and size above.
    norm_vertices = normalize_shape(norm_vertices,y_aspect_ratio,z_aspect_ratio,dest_vol)
    
    return norm_vertices

'''

    MAIN

'''

    
nbPts = 20
y_aspect_ratio = 2
z_aspect_ratio = 1
dest_vol       = 2.5 

final_vertices = generate_shape(nbPts,y_aspect_ratio,z_aspect_ratio,dest_vol)

print(np.shape(final_vertices))

ch = ConvexHull(final_vertices)
cm = ch.centroid()
print('orginal centroid', cm)
          
         
# computation of inertia and volume
inertia,volume=ch.inertia(cm)
print(inertia,volume)
