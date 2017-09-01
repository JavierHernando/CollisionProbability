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
!   Serra's method
!       P = Serra(x0, y0, sigmax, sigmay, R, N) 
!
!

module CollisionProbability_Serra
    use kind_constants
    implicit none

    interface
        function factorial(n)
            integer :: n, factorial
        end function
    end interface

    contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Serra et al algorithm
! 'Serra' calculates the collision probability using the N+1 first terms of Serra's Method
! Algorithm 2 calculates the collision probability with a guaranteed error bound (note: it's bases on exponential functions so it can be very conservative)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function Serra (x, y, sigmax, sigmay, R, N) result (Pc)    
        ! bounds and positivity of the terms are only guaranteed if sigmax > sigmay
        real(dp), intent(in) :: R
        integer, intent(in) :: N
        real(dp) :: x, y, sigmax, sigmay
        real(dp) :: Pc, p, phi, wx, wy, alpha0 
        integer ::  i
        real(dp), dimension(0:max(N-1, 0)) :: c

        if (sigmay > sigmax) then
! Note: bounds and positivity of the terms are only guaranteed if sigmax > sigmay
! if this is not true, we swap the axes using p as dummy variable
            p = sigmax;
            sigmax = sigmay;
            sigmax = p;

            p = x;
            x = y;
            y = p;
        end if

        ! Auxiliary Variables
        p = 1.0_dp / ( 2.0_dp * sigmay**2 )
        phi = 1.0_dp - ( sigmay / sigmax )**2
        wx = x**2 / ( 4.0_dp * sigmax**4 )
        wy = y**2 / ( 4.0_dp * sigmay**4 )
        alpha0 = exp( -( (x/sigmax)**2 + (y/sigmay)**2 )/ 2.0_dp ) / ( 2.0_dp * sigmax * sigmay )

        !Calculating the terms in the series
        c(0) = R**2
        if (N > 1) then
            c(1) = R**4 / 2.0_dp * ( p*(phi/2.0_dp + 1.0_dp) + wx + wy )
            if (N > 2) then 
                c(2) = R**6 / 12.0_dp * ( ( p*(phi/2.0_dp+1.0_dp) + wx + wy )**2 +  p**2*(phi**2/2.0_dp + 1.0_dp) + 2.0_dp*p*phi*wx )
                if (N > 3) then
                    c(3) = R**8 / 144.0_dp *  ( ( p*(phi/2+1) + wx + wy)**3 + 3*(p*(phi/2+1) + wx + wy) * ( p**2*(phi**2/2.0_dp+1.0_dp) + 2.0_dp*p*phi*wx ) + 2.0_dp*( p**3*(phi**3/2.0_dp+1.0_dp) + 3.0_dp*p**2*phi**2*wx ) );
                endif
            endif
        endif
        
        do i = 0, N-5
            c(i+4) = -(R**8 * p**3 * phi**2 * wy) / ((i+2) * (i+3) * (i+4)**2 * (i+5)) * c(i) + ( R**6*p**2*phi*( p*phi*(i+2.5_dp)+2*wy*(phi/2+1) ) )/((i+3) * (i+4)**2 * (i+5)) * c(i+1) - ( R**4*p*( p*phi*(phi/2.0_dp+1.0_dp)*(2*i+5) + phi*(2.0_dp*wy+1.5_dp*p) + wx + wy) )/((i+4)**2 * (i+5)) * c(i+2) + (R**2*( p*(2*phi+1)*(i+3) + p*(phi/2+1) + wx + wy ))/((i+4) * (i+5))*c(i+3)
        enddo

        !Adding it up
        Pc = sum(c) * alpha0 * exp( -p*R**2)
    end function Serra





    function SerraAlg2 (x, y, sigmax, sigmay, R, delta) result (Pc)
        real(dp), intent(in) :: R, delta
        real(dp) :: x, y, sigmax, sigmay
        real(dp) :: p, phi, wx, wy, alpha0, upper, lower
        real(dp), dimension(3) :: Pc 
        integer :: n, N1, N2


        if (sigmay > sigmax) then
! Note: bounds and positivity of the terms are only guaranteed if sigmax > sigmay
! if this is not true, we swap the axes using p as dummy variable
            p = sigmax;
            sigmax = sigmay;
            sigmax = p;

            p = x;
            x = y;
            y = p;
        end if

        ! Auxiliary Variables
        p = 1.0_dp / ( 2.0_dp * sigmay**2 )
        phi = 1.0_dp - ( sigmay / sigmax )**2
        wx = x**2 / ( 4.0_dp * sigmax**4 )
        wy = y**2 / ( 4.0_dp * sigmay**4 )
        alpha0 = exp( -( (x/sigmax)**2 + (y/sigmay)**2 )/ 2.0_dp ) / ( 2.0_dp * sigmax * sigmay )

        lower = alpha0 * ( 1 - exp(-p*R**2) ) / p
        upper = alpha0 * ( exp( p*R**2*( phi/2.0_dp + (wx+wy)/p ) ) - exp(-p*R**2) ) / ( p*(1 + phi/2) + wx + wy ) 
        if (upper - lower .le. delta) then
            Pc = (/ (upper+lower)/2.0_dp  ,  lower  ,  upper /)
        else
            N1 = 2.0_dp * ceiling( exp(1.0_dp) * p * R**2 * ( 1.0_dp + phi/2.0_dp + (wx+wy)/p ) )
            N2 = ceiling( log( ( alpha0 * exp( p*R**2*(phi/2.0_dp + (wx+wy)/p) ) ) / ( log(2.0_dp) * delta * p * sqrt(2.0_dp*PI*N1) * ( 1.0_dp + phi/2.0_dp + (wx+wy)/p ) ) ) )
            n = max(N1, N2) - 1
            Pc(1) = Serra (x, y, sigmax, sigmay, R, n)
            upper = alpha0 * exp( -p*R**2) * (p*R**2)**(n+1) / ( p * factorial(n+1) )
            lower = alpha0 * exp( p*R**2*(phi/2.0_dp + (wx+wy)/p) ) * ( p*R**2 * ( 1 + phi/2.0_dp + (wx+wy)/p ) )**(n+1) / ( p*(1.0_dp + phi/2.0_dp + (wx+wy)/p ) * factorial(n+1)  )

            Pc(2) = Pc(1) - lower
            Pc(3) = Pc(1) + upper

        endif
    end function SerraAlg2
  







end module CollisionProbability_Serra
