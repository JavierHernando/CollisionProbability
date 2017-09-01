! 2016 / 09 / 06
! Javier Hernando-Ayuso
! 

program test


    use CollisionProbability_GPHA
    use CollisionProbability_Serra
    use CollisionProbability_Chan
    use CollisionProbability_Alfano
    use CollisionProbability_Foster



implicit none

    real(kind=8) :: sigmax, sigmay, x0, y0, R
    real(kind=8) :: PGPHA, PSerra, PChan, PAlfano, PFoster

    R = 1
    sigmax = 100
    sigmay = 20
    x0 = 50
    y0 = 10

    PGPHA = GPHA(x0, y0, sigmax, sigmay, R, 10) 
    PSerra = Serra(x0, y0, sigmax, sigmay, R, 10) 
    PChan = Chan(x0, y0, sigmax, sigmay, R, 10) 
    PAlfano = Alfano(x0, y0, sigmax, sigmay, R, 50) 
    PFoster = Foster(x0, y0, sigmax, sigmay, R) 

    print *, 'GPHA   ', PGPHA
    print *, 'Serra  ', PSerra
    print *, 'Chan   ', PChan
    print *, 'Alfano ', PAlfano
    print *, 'Foster ', PFoster

   
end program test
