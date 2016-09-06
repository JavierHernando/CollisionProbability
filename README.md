# CollisionProbability
Methods for calculating Collision probability between two bodies in the short-term encounter scenario 


Files:
    CollisionProbability_Methods.f90    This file collects the other 6 modules
    CollisionProbability_GPHA.f90       GPHA method
    CollisionProbability_Serra.f90      Serra's method
    CollisionProbability_Chan.f90       Chan's method
    CollisionProbability_Alfano.f90     Alfano's method
    CollisionProbability_Foster.f90     Foster's method
    factorial.f90                       Factorial function
    test.f90                            Sample test/example program
    CompileProbFiles.sh                 compile script for Linux using gfortran


All methods assume main axis of the covariance matrix.
Common inputs:
    covariance matrix eigenvalues: sigmax, sigmay
    miss distance: x0, y0 (in main axis)
    combined radius of the bodies: R
Inputs for sonly some methods:
    maximum summation index: N


 Function syntax:
    Garcia-Pelayo Hernando-Ayuso (GPHA) method
        P = GPHA(x0, y0, sigmax, sigmay, R, N)

    Serra's method
        P = Serra(x0, y0, sigmax, sigmay, R, N)
 
    Chan's method
        P = Chan(x0, y0, sigmax, sigmay, R, N)

    Alfano's method
        P = Alfano(x0, y0, sigmax, sigmay, R, N)

    Foster's method
        P = Foster(x0, y0, sigmax, sigmay, R)


