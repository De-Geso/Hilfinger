! Generic utilities and functions.
module utilities
use kind_parameters
implicit none
public

contains

pure function percent_difference (a, b) result (x)
! Calculates the percent difference between two real(8) numbers.
	real(dp), intent(in) :: a, b
	real(dp) :: x
	x = 100._dp * abs(a-b) / ((a+b)/2)
end function


!function chi2(t, obs, func) result(chi)
!	real(dp) :: chi
!	real(dp), intent(in) :: t, obs(:)
!	real(dp), external :: func
!	integer :: i
	
!	chi = sum((obs(i) - func(t(i)))**2) / len(obs)
!	do i = 1, len(obs)
!		chi = chi + (obs(i) - func(t(i))
!	end do
!end function

end module
