! 2016 / 09 / 06
! Javier Hernando-Ayuso
!
! Assuming main axis of the covariance matrix.
! Iinputs:
!   covariance matrix eigenvalues: sigmax, sigmay
!   miss distance: x0, y0 (in main axis)
!   combined radius of the bodies: R
!   maximum summation index: N
!
! Function syntax:
!   Chan's method
!       P = Chan(x0, y0, sigmax, sigmay, R, N)


module CollisionProbability_Chan
    implicit none


    contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Chan's method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Chan(x0, y0, sigmax, sigmay, R, imax) result (P)
        real(kind=8), intent(in) :: x0, y0, sigmax, sigmay, R
        integer, intent(in) :: imax
        real(kind=8) :: P

        P = RicianIntegral(R**2 / sigmax / sigmay, (x0 / sigmax)**2 + (y0 / sigmay)**2, imax)
    end

! Rician Integral function using the iterative method by Chan 
    function RicianIntegral(u, v, m) result(P)
        implicit none

        real(kind=8) :: t, s, SS
        real(kind=8), intent(in) :: u, v
        real(kind=8) :: P
        integer, intent(in) :: m
        integer :: i

        t = 1.0D0
        s = 1.0D0
        SS = 1.0D0
        !first iteration
        P = exp( -v/2.0 ) * t - exp( -(u+v)/2.0 ) * t * SS
        !iterative expression
        do i = 1, m
            t = (v/2.0) /i * t
            s = (u/2.0) /i * s
            SS = SS + s
            P = P + exp( -v/2.0 )*t - exp( -(u+v)/2.0 ) * t * SS
        end do
    end function RicianIntegral





end module CollisionProbability_Chan
