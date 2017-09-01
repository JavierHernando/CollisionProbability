module kind_constants
!this module defines the kind for floating numbers we use
    implicit none

    integer, parameter  :: dp=kind(0.d0)
    real(dp), parameter :: infty = huge(1.0_dp)

end module kind_constants
