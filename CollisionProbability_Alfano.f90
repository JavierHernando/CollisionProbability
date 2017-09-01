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
    use kind_constants
    implicit none


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
        real(dp), intent(in) :: xm, ym, a, b, R
        integer, intent(in) :: m

        real(dp) :: P, dx, x, y, mOdd, mEven, m0
        integer :: i

        dx = R / (2.0_dp * m)

        !Odd terms
        mOdd = 0.0_dp
        do i = 1, m
            x = (2.0_dp*real(i, dp) - 1.0_dp) * dx - R
            y = sqrt(R**2 - x**2)
            mOdd = mOdd + ( erf( (ym + y)/(b * sqrt(2.0_dp)) ) - erf( (ym - y)/(b * sqrt(2.0_dp)) ) ) * ( exp( -(xm + x)**2 / (2.0_dp * a**2) ) + exp( -(xm - x)**2 /(2.0_dp * a**2) ) )
        end do
        mOdd = 4D0 * mOdd

        !Even terms
        mEven = ( erf( (ym + R)/(b * sqrt(2.0_dp)) ) - erf( (ym - R)/(b * sqrt(2.0_dp)) ) ) * exp( -xm**2 / (2.0_dp * a**2) )
        do i = 1, (m-1)
            x = 2.0_dp * real(i, 8) * dx - R
            y = sqrt(R**2 - x**2)
            mEven = mEven + ( erf( (ym + y)/(b * sqrt(2.0_dp)) ) - erf( (ym - y)/(b * sqrt(2.0_dp)) ) ) * ( exp( -(xm + x)**2 / (2.0_dp * a**2) ) + exp( -(xm - x)**2 /(2.0_dp * a**2) ) ) 
        end do
        mEven = 2.0_dp * mEven

        ! Zero order term
        x = 0.015_dp * dx - R
        y = sqrt(R**2 - x**2)
        m0 = 2.0_dp * ( erf( (ym + y)/(b * sqrt(2.0_dp)) ) - erf( (ym - y)/(b * sqrt(2.0_dp)) ) ) * ( exp( -(xm + x)**2 / (2.0_dp * a**2) ) + exp( -(xm - x)**2 /(2.0_dp * a**2) ) )

        !Combining the previous results
        P = dx * ( m0 + mEven + mOdd) / ( 3.0_dp * sqrt(8.0_dp * PI) * a)
    end function Alfano








end module CollisionProbability_Alfano
