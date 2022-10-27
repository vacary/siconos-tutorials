# Some examples and uniyt tests of czm use

## Binary CZM

BinaryCohesiveNSL.hpp  Implementation in Mechanics.

+ For this NSLAW, a rule to update beta is coded direclty in the NSLAW and is based on the internal variable of the interaction.
  + it could be better to perform the explicit integration in MoreauJeanOSI or in CohesiveFrictionalContact (in the latter the structure of the problem will change.)

+ To take into account the tangential part:
   - in 2D, we need either from the LCP
   - in 3D, we need to construct a new fdc3d problem with the LCP variable for \beta
  
+ Assembly type
  + due to the structure of the problem that needs to compute V. the assembly type REDUCED_BLOCK is really slow.
  
+ Initialization of the Interaction through Bullet
  + we need to be careful on the initilization of beta if bullet creates new contact afterwards.

+ output cohesion forces in hdf5?
  + how to do that?
  + how to visualize cohesion forces?
  
+ add tangential part in explicit ?  
  


## BouncingBallTS-BinaryCZM

In this example, a ball is glued on the floor with a binary CZM

We apply piecewise constant external force.
