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
	x = 100._dp * abs(a-b) / ((a+b) / 2)
end function

end module
