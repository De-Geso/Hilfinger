program gillespie
! Simulate the system given in Hilfinger 2016 Eq. 10 using the Gillespie
! algorithm. Imagine a feedback control system. mRNA is produced, which
! triggers the production of a Protein. Both decay at some rate which
! depends on their abundance.
! x0 == R >> x0 + a			x1 == lmbda x0 >> x1 + b
! x0 == x0/tau1 >> x0 - 1	x1 == x1/tau1 >> x1 - 1
use kind_parameters
implicit none

! == System parameters ==
integer, dimension(2), parameter :: burst = [1, 1]	! Size of burst for [x0, x1]
real(dp), parameter :: alpha = 1.0_dp	! x0 production rate
real(dp), parameter :: lmbda = 1.0_dp	! x1 production rate constant
real(dp), parameter, dimension(2) :: r = [1.0_dp, 1.0_dp]	! Decay rates
integer, parameter, dimension(4,2) :: abund_update = &	! Abundance update matrix
	reshape((/burst(1), 0, &
			-1, 0, &
			0, burst(2), &
			0, -1/), shape(abund_update))

! == Program hyperparameters ==
real(dp), parameter :: tmax = 100._dp	! When we want to stop
character(*), parameter :: f_abund = "abundances.dat"	! File for abundances

! == Other variables ==
integer, dimension(2) :: abund = [0, 0]	!	(Initial) abundances
real(dp), dimension(4) :: propensity = 0._dp	! The rates, but with a fancy name
real(dp) :: time = 0._dp
real(dp), allocatable :: t(:)	! Array for time
integer, allocatable :: x(:,:)	! Array for abundances
integer :: steps, i, io


call random_seed()	! Initialize random seed. Later it might be useful to use the Numerical Recipes version?

steps = 0
! Open file to hold abundances
open(newunit=io, file=f_abund, action="write")
do while (time < tmax)
	write(io,*) time, abund
	! Do gillespie algorithm
	call gillespie_iter(abund, time)
	steps = steps + 1
end do
! Take care of last loop iteration where we go over tmax.
write(io,*) time, abund
close(io)
steps = steps + 1

! Fill allocatable arrays with data we just generated
call get_abundance(steps)
! Make the covariance matrix
call covariance_matrix(steps, x(:,1), x(:,2))



contains !==============================================================



subroutine get_abundance(n)
	integer, intent(in) :: n
	integer :: i, io
	
	allocate(x(n,2))
	allocate(t(n))
	open(newunit=io, file=f_abund, action="read")
	do i = 1,n
		read(io, *) t(i), x(i,1), x(i,2)
	end do
	close(io)
		
end subroutine


subroutine covariance_matrix(n, x, y)
	integer, intent(in) :: n
	integer, dimension(n), intent(in) :: x, y
	real(dp) :: cov(2,2)
	integer :: i, j
	
	! Naive algorithm
	cov(1,1) = 1._dp*(sum(x*x) - sum(x)*sum(x)/n) / n
	cov(1,2) = 1._dp*(sum(x*y) - sum(x)*sum(y)/n) / n
	cov(2,1) = 1._dp*(sum(y*x) - sum(y)*sum(x)/n) / n
	cov(2,2) = 1._dp*(sum(y*y) - sum(y)*sum(y)/n) / n
	
	write(*,*) "Covariance matrix:"
	do i = 1,2
		write(*,*) cov(i,1), cov(i,2)
	end do
end subroutine


! Run one step of the Gillespie algorithm: update the propensities, take
! a time step, get the reaction, update abundances. 
subroutine gillespie_iter(x, t)
	real(dp), intent(inout) :: t
	integer, intent(inout) :: x(2)
	integer :: i
	real(dp) :: psum, propsum, roll, dt

	! Update propensity
	call update_propensity(x)
	propsum = sum(propensity)
	
	! Update time
	call random_exp(1./propsum, dt)
	t = t + dt
	
	! Get reaction
	call random_number(roll)
	i = 0
	psum = 0._dp
	do while (psum < roll*propsum)
		i = i + 1
		psum = psum + propensity(i)
	end do
	
	! Update abundances
	x = x + abund_update(i,:)
end subroutine


! Updates propensities depending on the state of the system
subroutine update_propensity(x)
	integer, intent(in) :: x(2)
	
	propensity(1) = alpha		! Make x0 (mRNA)
	propensity(2) = x(1)*r(1)	! Degrade x0 (mRNA)
	propensity(3) = lmbda*x(1)	! Make x1 (Protein)
	propensity(4) = x(2)*r(2)	! Degrade x1 (Protein)
end subroutine


! Converts continuous random variable (0,1] to exponential random
! variable with mean = l. Set up as subroutine to match random_number
subroutine random_exp(l,u)
	real(dp), intent(out) :: u
	real(dp) l, x
1	call random_number(x)
	if (x == 0) goto 1	! save us if we manage to roll a 0
	u = -l*log(x)
end subroutine

end program
