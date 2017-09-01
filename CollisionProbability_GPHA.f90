! 2016 / 09 / 06
! Javier Hernando-Ayuso
!
! Assuming main axis of the covariance matrix.
! Inputs:
!   covariance matrix eigenvalues: sigmax, sigmay
!   miss distance: x0, y0 (in main axis)
!   combined radius of the bodies: R
!   maximum summation index: imax
!
! Function syntax:
!   Garcia-Pelayo Hernando-Ayuso (GPHA) method
!       P = GPHA(x0, y0, sigmax, sigmay, R, N)
!

module CollisionProbability_GPHA
    use kind_constants
    implicit none


    interface
        function factorial(n)
            integer :: n, factorial
        end function
    end interface

    contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GPHA method (Garcia-Pelayo Hernando-Ayuso)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function GPHA(x0, y0, a, b, R, imax) result (P)
        implicit none
      
        real(dp), intent(in) :: x0, y0, a, b, R
        integer, intent(in) :: imax
        real(dp) :: P
        integer :: i, j
        real(dp) :: auxSum, AR, R2a
        real(dp), dimension(2*imax+1) :: HermiteX, HermiteY

        !Hermite Polynomials
        !H_i is stored in the i+1 memory position of the array
        HermiteX(1) = 1.0_dp    !H0(x) = 1
        HermiteX(2) = x0 / a    !H1(x) = x
        HermiteY(1) = 1.0_dp    !H0(x) = 1
        HermiteY(2) = y0 / b    !H1(x) = x
        do i = 2, 2*imax
            HermiteX(i+1) = (x0/a) * HermiteX(i) - (i-1) * HermiteX(i-1)
            HermiteY(i+1) = (y0/b) * HermiteY(i) - (i-1) * HermiteY(i-1)
        end do

        !Double summation
        AR = a/b
        R2a = R/2.0_dp/a 
        P = 0
        do i = 0, imax
            auxSum = 0
            do j = 0, i
                auxSum = auxSum + HermiteX(2*(i-j)+1) * AR**(2*j) * HermiteY(2*j+1) / factorial(j) / factorial(i-j)
            end do
            P = P + auxSum * R2a**(2*i) / factorial(i+1)
            ! print *, P
        end do

        ! Constant term left outside the double summation
        P = P * exp(-((x0/a)**2 + (y0/b)**2)/2.0_dp) * (R / a) * (R / b) / 2.0_dp
    end function GPHA












end module CollisionProbability_GPHA
