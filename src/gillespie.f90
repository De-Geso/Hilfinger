program gillespie
! Simulate the system given in Hilfinger 2016 Eq. 10 using the Gillespie
! algorithm. Imagine a feedback control system. mRNA is produced, which
! triggers the production of a Protein. Both decay at some rate which
! depends on their abundance.
! x0 == alpha >> x0 + a			x1 == beta x0 >> x1 + b
! x0 == x0/tau1 >> x0 - 1	x1 == x1/tau1 >> x1 - 1
use kind_parameters
implicit none

! System parameters ====================================================
integer, dimension(2), parameter :: burst = [1, 1]	! Size of burst for [x0, x1]
real(dp), parameter :: alpha = 1.0_dp	! x0 production rate
real(dp), parameter :: beta = 1.0_dp	! x1 production rate constant
real(dp), parameter, dimension(2) :: decay = [1.0_dp, 1.0_dp]	! Decay rates
integer, parameter, dimension(4,2) :: abund_update = &	! Abundance update matrix
	reshape((/burst(1), 0, &
			-1, 0, &
			0, burst(2), &
			0, -1/), shape(abund_update), order=[2,1])

! Program hyperparameters ==============================================
integer, parameter :: dmax = 10**3	! When we want to stop
integer, parameter :: abund_max = 10**2
character(*), parameter :: f_abund = "abundance.dat"	! File for abundances
character(*), parameter :: f_prop = "propensity.dat"	! File for propensities

! Other variables ======================================================
integer, dimension(2) :: abund = [0, 0]	!	(Initial) abundances
real(dp), dimension(4) :: propensity = 0._dp	! The rates, but with a fancy name
real(dp) :: time = 0._dp, dt
real(dp), dimension(2, abund_max) :: probx = 0._dp
real(dp), allocatable :: t(:), r(:,:)	! Array for time and rates
integer, allocatable :: x(:,:)	! Arrays for abundances
integer :: steps=0, ndecay(2)=0, i, io1, io2, event


! ======================================================================
call random_seed()	! Initialize random seed. Later it might be useful to use the Numerical Recipes version?

! Open file to hold abundances and rates
open(newunit=io1, file=f_abund, action="write")
open(newunit=io2, file=f_prop, action="write")

do while (minval(ndecay) < dmax)
	write(io1,*) time, abund
	write(io2,*) time, propensity
	! Do gillespie algorithm
	call gillespie_iter(abund, dt, event)
	time = time + dt
	! If we performed a decay step, add it to the counter.
	if (event == 2 .or. event == 4) then
		ndecay(event/2) = ndecay(event/2) + 1
	end if
	! Add the time step to whatever the current abundances are.
	do i = 1,2
		probx(i,abund(i)+1) = probx(i,abund(i)+1) + dt
	end do
	abund = abund + abund_update(event,:)
	steps = steps + 1
end do
! Take care of last loop iteration where we go over tmax.
write(io1,*) time, abund;	close(io1)
write(io2,*) time, propensity;	close(io2)
steps = steps + 1

probx = probx / time
write(*,*) sum(probx)

do i = 1, abund_max
	write(1,*) i-1, probx(1,i), probx(2,i)
end do

call checks(probx)

! Fill allocatable arrays with data we just generated.
! call get_allocatable(steps)

! Make the covariance matrix.
! call covariance_matrix(x(:,1), x(:,2), steps)


contains !==============================================================


subroutine checks(p)
	real(dp), dimension(2, abund_max), intent(in) :: p
	real(dp) :: mean(2), covar(2,2)=0.
	integer i, j, k, l
	
	do i = 1, abund_max
		mean(:) = 1._dp * (mean(:) + (i-1.) * p(:,i))
	end do
	write(*,*) 'Theory mean: ', alpha/decay(1), mean(1)*beta/decay(2)
	write(*,*) 'Sim mean: ', mean
	
	do i = 1, 2
	do j = 1, 2
		do k = 1, abund_max
			if (i /= j) then
				do l = 1, abund_max
					covar(i,j) = covar(i,j) + p(i,k)*p(j,l)*((k-1.)-mean(i))*((l-1.)-mean(j))
				end do
			else
				covar(i,j) = 1._dp * covar(i,j) + p(i,k) * ((k-1.)-mean(i))**2
			end if
		end do
		! covar(i,j) = 1._dp * covar(i,j) / (mean(1)*mean(2))
	end do
	end do
	
	write(*,*) 'Theory covar: ', &
		1./mean(1) * decay(2)/sum(decay), &
		1./mean(1) * decay(2)/sum(decay), &
		1./mean(2) + 1./mean(1) * decay(2)/sum(decay)
	write(*,*) 'Sim covar: ', covar(1,2), covar(2,1), covar(2,2)
	
	
end subroutine


function covariance(x, y, n) result(cov)
	real(dp), dimension(n), intent(in) :: x, y
	integer, intent(in) :: n
	real(dp) :: cov
	
	cov = 1._dp * (sum(x*y) - sum(x)*sum(y)/n) / n
end function

subroutine get_allocatable(n)
	integer, intent(in) :: n
	integer :: i, io1, io2
	real :: dum
	
	allocate(x(n,2))
	allocate(r(n,4))
	allocate(t(n))
	
	open(newunit=io1, file=f_abund, action="read")
	open(newunit=io2, file=f_prop, action="read")
	do i = 1,n
		read(io1, *) t(i), x(i,1), x(i,2)
		read(io2, *) dum, r(i,1), r(i,2), r(i,3), r(i,4)
	end do
	close(io1)
	close(io2)
end subroutine


subroutine covariance_matrix(x, y, n)
	integer, intent(in) :: n
	integer, dimension(n), intent(in) :: x, y
	real(dp) :: cov(2,2)
	integer :: i
	
	! Naive algorithm
	cov(1,1) = (sum(x*x) - 1._dp*sum(x)*sum(x)/n) / n
	cov(1,2) = (sum(x*y) - 1._dp*sum(x)*sum(y)/n) / n
	cov(2,1) = (sum(y*x) - 1._dp*sum(y)*sum(x)/n) / n
	cov(2,2) = (sum(y*y) - 1._dp*sum(y)*sum(y)/n) / n
	
	write(*,*) "Covariance matrix:"
	do i = 1,2
		write(*,*) cov(i,1), cov(i,2)
	end do
end subroutine


! Run one step of the Gillespie algorithm: update the propensities, take
! a time step, get the reaction, update abundances. 
subroutine gillespie_iter(x, tstep, i)
	real(dp), intent(out) :: tstep
	integer, intent(inout) :: x(2)
	integer, intent(out) :: i
	real(dp) :: psum, propsum, roll

	! Update propensity
	call update_propensity(x)
	propsum = sum(propensity)
	
	! Update time
	call random_exp(1./propsum, tstep)
	
	! Get reaction
	call random_number(roll)
	i = 0
	psum = 0._dp
	do while (psum < roll*propsum)
		i = i + 1
		psum = psum + propensity(i)
	end do
	
	! Update abundances
	! x = x + abund_update(i,:)
end subroutine


! Updates propensities depending on the state of the system
subroutine update_propensity(x)
	integer, intent(in) :: x(2)
	
	propensity(1) = alpha		! Make x0 (mRNA)
	propensity(2) = x(1)*decay(1)	! Degrade x0 (mRNA)
	propensity(3) = beta*x(1)	! Make x1 (Protein)
	propensity(4) = x(2)*decay(2)	! Degrade x1 (Protein)
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
