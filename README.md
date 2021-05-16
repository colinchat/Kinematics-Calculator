# Kinematics-Calculator
Iterative suspension kinematics calculator in MATLAB (original idea from Tim.Wright on FSAE forums)

## How it works
This calculator takes in a set of inboard (on car) and outboard (on upright) points which determine one quarter of a double a-arm suspension. 
It then converts the points into vectors to represent the linkages and wheel. 
Then it sets up velocities at each outboard point in terms of 5 kinematic variables Toe, Camber, Spin, Track, and Wheelbase.
Then it solves a system of linear equations created by taking the dot product of each outboard points' velocity with one of the 5 linkages it connects to.
Then it moves the points slightly, records data, draws a visual representation using quiver3, then repeats the process. 

## Current state
Currently, the calculator can compute all the kinematic variables in simple bump or droop movements. 
Simulating roll accurately as well as finding roll center accurately is currently being developed.
