# EpsilonThetaWallFunctions

Boundary conditions for epsilonTheta. To be used with the 2 equation models for the thermal variance and its dissipation rate. The difference between EpsilonThetaWallFunction and EpsilonTheta2WallFunction is that the latter sets the production of the turbulent kinetic energy in the first cell off the wall to the appropriate value (zero in the case of low-Re turbulence models). The current implementation does not set the production of the thermal variance and the dissipation rate to appropriate
values in the case of y<sup>+</sup> > 30 ( high Re cases). 
