! Generic utilities and functions.
module utilities
use ieee_arithmetic
use kind_parameters
implicit none
public

real(dp), parameter :: utilities_eps = 1E-14


contains


function relative_change(x, y, method) result(change)
! Calculate the relative change between two numbers using specified
! measure. y-x
	real(dp), intent(in) :: x, y
	character(len=*), intent(in) :: method
	real(dp) :: change, f
	
	! Unequal x and y
	select case (trim(method))
	case ("classical")
		f = 1._dp * x
	case ("logarithmic")
		f = 1._dp * (y-x) / log(y/x)
	case ("arithmetic")
		f = 0.5_dp * (x+y)
	case default
		! Unrecognized case
		write(*,*) "ERROR: unknown method. Options are 'classical', 'arithmetic', 'logarithmic'"
		change = ieee_value(change, ieee_quiet_nan)  ! Return NaN
		return
	end select
	
	! Equal x and y
	if (y-x .eq. 0._dp) then
		change = 0._dp
	! Throw error if we are dividing by zero
	else if (f .eq. 0._dp) then
		write(*,*) "ERROR: division by zero in utilities/relative_change!"
		change = ieee_value(change, ieee_quiet_nan)  ! Return NaN
	! No exceptions
	else
		change = 1._dp * (y-x) / f
	end if
end function relative_change


! Functions below this line might be shitty ============================


pure function percent_difference (a, b) result (x)
! Calculates the percent difference between two real(8) numbers.
	real(dp), intent(in) :: a, b
	real(dp) :: x
	x = 100._dp * abs(a-b) / ((a+b)/2)
end function percent_difference

subroutine get_command_line_arg(x, i)
	integer, intent(in) :: i
	real(dp), intent(out) :: x
	character(len=32) :: arg
	
	! Get first command line argument
	call get_command_argument(i, arg)
	! Convert argument from string to real
	read(arg, *) x
end subroutine get_command_line_arg

subroutine generate_ISO_filename(path, prefix, suffix, filename)
! Generates an ISO 8601 compliant filename suffix. Used to create
! unique filenames up to milliseconds.
	character(len=*), intent(in) :: path, prefix, suffix
	character(len=*), intent(out) :: filename
	character(len=32) :: datetime
	integer, dimension(8) :: values  ! Array to store date and time components
	
	call date_and_time(VALUES=values)
	! Format date and time to ISO 8601 standards.
	write(datetime, '(I4.4, I2.2, I2.2, "T", I2.2, I2.2, I2.2, ".", I3.3)') &
		values(1), values(2), values(3), values(5), values(6), values(7), values(8)
	! Create filename
	filename = trim(adjustl(path)) // trim(adjustl(prefix)) // trim(adjustl(datetime)) // trim(adjustl(suffix))
end subroutine

function linspace(start,end,num,endpoint,step) result(samples)
	! PARAMETERS
	real(dp), intent(in) :: start 
	!! The starting value of the sequence.
	real(dp), intent(in) :: end
	!! The end value of the sequence, unless `endpoint` is set to `.false.`. 
	!! In that case, the sequence consists of all but the last of `num + 1` 
	!! evenly spaced samples, so that `end` is excluded. Note that the 
	!! step size changes when `endpoint` is `.false.`.
	integer, intent(in), optional :: num
	!! Number of samples to generate. Default value is 50.
	logical, intent(in), optional :: endpoint
	!! If `.true.`, `end` is the last sample. Otherwise, it is not included. Default is `.true.`.
	real(dp), intent(out), optional :: step
	!! If present, `step` is the size of spacing between samples.

	! RETURNS
	real(dp), allocatable :: samples(:)
	!! There are `num` equally spaced samples in the closed interval `[start, stop]` or 
	!! the half-open interval `[start, stop)` (depending on whether `endpoint` is `.true.` or `.false.`).
	integer :: num_, i
	logical :: endpoint_
	real(dp) :: step_

	num_ = 50
	if (present(num)) num_ = num
	
	endpoint_ = .true.
	if (present(endpoint)) endpoint_ = endpoint
	
	! find step size
	if (endpoint_) then
		step_ = (end - start)/real(num_-1,dp)
	else
		step_ = (end - start)/real(num_,dp)
	end if
	
	if (present(step)) step = step_
	
	allocate(samples(num_))
	do i = 1, num_
		samples(i) = start + (i-1)*step_
	end do
end function linspace

end module
