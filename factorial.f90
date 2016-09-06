! 2016 / 09 / 06
! Javier Hernando-Ayuso
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Factorial function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function factorial(n)
    implicit none
    integer :: i, n, factorial
    integer, pointer :: iP
    if (n == 0) then
        factorial = 1
    else
        factorial = 1
        do i = 2, n
            factorial = factorial*i
        end do
    end if
end function factorial

