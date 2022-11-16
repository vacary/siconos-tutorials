"""
TODO:
"""

#
# Minimal working example for coordinate system deviation between simulated terrain and input height map grid
#

import siconos.numerics as sn
import siconos.kernel as sk

from siconos.mechanics.collision.convexhull import ConvexHull
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner

import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)

  
######################################################################
########################## FUNCTIONS #################################
######################################################################
 
def sphere_setup(io_hdf5, mZ, resolution):

    ####
    # Create spheres at specified (below) position and add them to siconos
    
    # -IN- #    
    # io_hdf5; <siconos -> MechanicsHdf5Runner>:    Siconos simulation
    # mZ; <np.array(float)>:                        Heightmap of simulation to add blocks to
    # resolution; <float>:                          Length of one heightmap's cell side / horizontal/vertical distance between grid points

    # -OUT- #    
    ####
    
    
    # Get all grid points
    positions = []
    for j in range(mZ.shape[0]):
        for k in range(mZ.shape[1]):
            positions.append([j, k])
    
    io_hdf5.add_primitive_shape('Ball', 'Sphere', [0.25])
    
    for j, start_index in enumerate(positions):
    
        # Get the corresponding 'real' coordinates
        ## x -> rows ; y -> cols (x ^ y ->) !!
        ## Center of coordiante system  (0, 0, 0) in the middle!
        ## Corners of heightmap are grid points, e.g. lower right = [-1, -1], top left = [0, 0]
        
        Xpos = resolution* (start_index[0] - 0.5*mZ.shape[0] )
        Ypos = resolution* (start_index[1] - 0.5*mZ.shape[1] )
        Zpos = mZ[start_index[0], start_index[1]]
        
        
        io_hdf5.add_object('ball_%d'%j, [Contactor('Ball', collision_group=100)], translation=[Xpos, Ypos, Zpos], velocity=[0,0,0,0,0,0], mass=1)  
    
    return 0


  
   
def create_incline(length, width, inclination, resolution):

    ####
    # Create and return heightmap-array

    # -IN- #    
    # length; <float>:      Length of inclined plane (m)
    # width; <float>:       Width of inclined plane (m)
    # inclination; <float>: Inclination of plane (degrees)
    # resolution; <float>:  Length of the cell's sides (m)
    
    # -OUT- #
    # mX; <np.array(<float>)>: X-axis grid (matrix), shape adjusted for Y-axis
    # mY; <np.array(<float>)>: Y-axis grid (matrix), shape adjusted for X-axis
    # mZ; <np.array(<float>)>: Z-axis grid / height grid (matrix), same shape as X and Y
    ####
    

    # Number of coordinate pairs i.e. segment resolution per axis
    nX = int(np.ceil(length/resolution))
    nY = int(np.ceil(width/resolution))

    # Coordinate axis
    x_axis = np.linspace(0, length, nX)    
    y_axis = np.linspace(-width/2, width/2, nY)
        
    # Coordinate grid (adjust shapes)
    ## adjusted for siconos' coordinate system definition (x ^ y ->)
    mX, mY = np.meshgrid(x_axis, y_axis, indexing='ij')


    # Compute z for given inclination (measured from x-axis -> lower for higher x) 
    mZ = np.tan(inclination*np.pi/180) *np.flip(mX, axis=0)
    
    return mX, mY, mZ
    
    
    
######################################################################
############################ MAIN ####################################
######################################################################

    # Create inclined plane

# Inclination of plane in degrees
inclination = 60
# Resolution of mesh-grid (used for heightmap) in m: Determines length of one quadratic-cell
resolution = 5

heightmap = create_incline(100, 20, inclination, resolution)[2]
# Soil properties
e = 0
mu = 0.1
mu_r = 0

# Max time to simulate
T = 50.


with MechanicsHdf5Runner() as io:


         # Setup soil
         
            # Add map and soil params to simulation
    io.add_height_map('MyTerrain_1', heightmap, (resolution*heightmap.shape[0], resolution*heightmap.shape[1]), outsideMargin=0.01, insideMargin=0.01)
    io.add_object('ground_1', [Contactor('MyTerrain_1', collision_group=1)], translation=[0, 0, 0])        
    
    io.add_Newton_impact_rolling_friction_nsl('contact_1', e=e, mu=mu, mu_r=mu_r, collision_group1=100, collision_group2=1)
         
         
        # Setup spheres
    sphere_setup(io, heightmap, resolution)
                      

        # Setup siconos solver
            
    # Create solver options
    options = sk.solver_options_create(sn.SICONOS_ROLLING_FRICTION_3D_NSGS)
    options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
    options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-4

    
    # Start the simulation
    io.run(
        with_timer=False,               # ?
        multipoints_iterations=False,   # ?
        gravity_scale=1,                # Divisor to 9.81 , i.e. if e.g. g = 1.62 is desired: gravity_scale=9.81/1.62 
        t0=0,                           # Starting time
        T=T,                            # End time
        h=1e-3,                         # Size of time step
        theta=0.50001,                  # ?
        Newton_max_iter=1,              # ?
        set_external_forces=None,       # ?
        solver_options=options,         # Specify the solver ?
        violation_verbose=False,        # ?
        numerics_verbose=False,         # ?
        output_frequency=None           # What datapoints (time points) to save to .hdf5 file
    )
        
        
        
