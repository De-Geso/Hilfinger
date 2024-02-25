module randf
use kind_parameters
implicit none
public

contains

subroutine random_exp(L, u)
	! Converts continuous random variable (0,1] to exponential random
	! variable with mean = L.
	real(dp), intent(out) :: u
	real(dp) L, x
	call random_number(x)
	! 1.0-x saves us if we manage to roll a 0.
	u = -L*log(1.0-x)
end subroutine


subroutine random_uniform(u, a, b)
	! Generate a uniformly distributed random variable a <= u < b.
	real(dp), intent(inout) :: u
	real(dp), intent(in) :: a, b
	call random_number(u)
	u = (b-a)*u + a
end subroutine

end module
