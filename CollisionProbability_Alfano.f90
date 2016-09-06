! 2016 / 09 / 06
! Javier Hernando-Ayuso
!
! Assuming main axis of the covariance matrix.
! Inputs:
!   covariance matrix eigenvalues: sigmax, sigmay
!   miss distance: x0, y0 (in main axis)
!   combined radius of the bodies: R
!   maximum summation index: N
!
! Function syntax:
!   Alfano's method
!       P = Alfano(x0, y0, sigmax, sigmay, R, N)


module CollisionProbability_Alfano
    implicit none

    real(kind=8), parameter :: PI = 4D0 * atan(1D0)

    contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Alfano Method: 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Alfano (xm, ym, a, b, R, m) result (P)
    ! Alfano recommends:
    ! m = floor(5*sA/ min([sigmax, sigmaz, sqrt(xe^2 + ze^2)]));
    ! m = max(m, 10);
    ! m = min(m, 50);
        implicit none
        real(kind=8), intent(in) :: xm, ym, a, b, R
        integer, intent(in) :: m

        real(kind=8) :: P, dx, x, y, mOdd, mEven, m0
        integer :: i

        dx = R / (2.0D0 * m)

        !Odd terms
        mOdd = 0D0
        do i = 1, m
            x = (2D0*real(i, 8) - 1D0) * dx - R
            y = sqrt(R**2 - x**2)
            mOdd = mOdd + ( erf( (ym + y)/(b * sqrt(2D0)) ) - erf( (ym - y)/(b * sqrt(2D0)) ) ) * ( exp( -(xm + x)**2 / (2D0 * a**2) ) + exp( -(xm - x)**2 /(2D0 * a**2) ) )
        end do
        mOdd = 4D0 * mOdd

        !Even terms
        mEven = ( erf( (ym + R)/(b * sqrt(2D0)) ) - erf( (ym - R)/(b * sqrt(2D0)) ) ) * exp( -xm**2 / (2D0 * a**2) )
        do i = 1, (m-1)
            x = 2D0 * real(i, 8) * dx - R
            y = sqrt(R**2 - x**2)
            mEven = mEven + ( erf( (ym + y)/(b * sqrt(2D0)) ) - erf( (ym - y)/(b * sqrt(2D0)) ) ) * ( exp( -(xm + x)**2 / (2D0 * a**2) ) + exp( -(xm - x)**2 /(2D0 * a**2) ) ) 
        end do
        mEven = 2D0 * mEven

        ! Zero order term
        x = 0.015D0 * dx - R
        y = sqrt(R**2 - x**2)
        m0 = 2D0 * ( erf( (ym + y)/(b * sqrt(2D0)) ) - erf( (ym - y)/(b * sqrt(2D0)) ) ) * ( exp( -(xm + x)**2 / (2D0 * a**2) ) + exp( -(xm - x)**2 /(2D0 * a**2) ) )

        !Combining the previous results
        P = dx * ( m0 + mEven + mOdd) / ( 3D0 * sqrt(8 * PI) * a)
    end function Alfano








end module CollisionProbability_Alfano
