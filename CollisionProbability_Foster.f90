! 2016 / 09 / 06
! Javier Hernando-Ayuso
!
! Assuming main axis of the covariance matrix.
! Inputs:
!   covariance matrix eigenvalues: sigmax, sigmay
!   miss distance: x0, y0 (in main axis)
!   combined radius of the bodies: R
!
! Function syntax:
!   Foster's method
!       P = Foster(x0, y0, sigmax, sigmay, R)
!

module CollisionProbability_Foster
    use kind_constants
    implicit none


    contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Foster Method: 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Foster (x0, y0, sigmax, sigmay, R) result (P)
        implicit none
        real(dp), intent(in) :: x0, y0, sigmax, sigmay, R
        real(dp) :: P
        real(dp) :: xe, phi

        xe = sqrt( x0**2 + y0**2 )
        phi = atan2(x0, y0)

        P = FosterAngleDist (xe, phi, sigmax, sigmay, R)

    end function Foster
    
    function FosterAngleDist (xe, phi, sigmax, sigmay, R) result (P)
        implicit none
        real(dp), intent(in) :: xe, phi, sigmax, sigmay, R
        real(dp) :: P, theta, dTheta, rho, dR, sinPhi, cosPhi, sinTheta, cosTheta

        integer :: nR, nTheta, iR, iTheta 

        !Foster sets dTheta = 0.5 deg, this means 720 points in theta
        nTheta = 720
        !Foster sets dr = R/12, this means 13 points in r
        nR = 13;

        dR = R / (nR - 1.0_dp);
        dTheta = 2.0_dp * PI / nTheta
        sinPhi = sin(phi)
        cosPhi = cos(phi)

        P = 0.0_dp
        theta = 0.0_dp
        do iTheta = 1, nTheta
            sinTheta = sin(theta)
            cosTheta = cos(theta)
            rho = 0.0_dp
            do iR = 1, nR
                if (iR == 1 .or. iR == nR) then
                     P = P + rho * exp( -xe**2/2*( ( sinPhi/sigmax )**2 + ( cosPhi/sigmay )**2 ) ) * exp( -rho**2/2.0_dp*( ( sinTheta/sigmax )**2 + ( cosTheta/sigmay )**2 ) +  rho*xe*( ( sinTheta*sinPhi )/sigmax**2 + ( cosTheta*cosPhi )/sigmay**2 ) )   / ( 2.0_dp * PI * sigmax * sigmay ) / 2.0_dp
                else
                     P = P + rho * exp( -xe**2/2*( ( sinPhi/sigmax )**2 + ( cosPhi/sigmay )**2 ) ) * exp( -rho**2/2.0_dp*( ( sinTheta/sigmax )**2 + ( cosTheta/sigmay )**2 ) +  rho*xe*( ( sinTheta*sinPhi )/sigmax**2 + ( cosTheta*cosPhi )/sigmay**2 ) )   / ( 2.0_dp * PI * sigmax * sigmay )
                end if
                rho = rho + dR
            end do
            theta = theta + dTheta
        end do
        P = P * dTheta * dR
    end function FosterAngleDist







end module CollisionProbability_Foster
