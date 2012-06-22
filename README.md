aum_monte_carlo3
================
this simulation is similar to the previous with the exception of 
boundary conditions. periodic boundary conditions still apply to the
x and y axes, but the percentage to fill the domain only applies to 
the bottom layer of the z axis. if the particle tries to move below
the z axis then it will remain in the same location. if it moves up
then another particle will take its place keeping the fill percentage
of the bottom layer constant. if it move past the top of the z-axis,
the particle will escape. diffusion.
