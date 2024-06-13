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
end function percent_difference

subroutine get_command_line_arg(x)
	real(dp), intent(out) :: x
	character(len=32) :: arg
	
	! Get first command line argument
	call get_command_argument(1, arg)
	! Convert argument from string to real
	read(arg, *) x
end subroutine get_command_line_arg

end module
