program mrna_gean
use kind_parameters
implicit none

! Hyperparameters
! ======================================================================
! Number of decay events before stopping
integer, parameter :: decay_min = 10**6
! Maximum abundances, pad this.
integer, parameter :: abund_max = 2**6

! Parameters
! ======================================================================
! Size of burst for [x0, x1]
integer, parameter, dimension(2) :: burst = [1, 1]
! x1 production rate
real(dp) :: alpha = 1.0
! x2 production rate
real(dp) :: beta = 1.0
! Decay rates
real(dp), dimension(2) :: decay = [1.0, 1.0]
! Abundance update matrix.
integer, parameter, dimension(2,4) :: abund_update = &
	reshape((/burst(1), 0, &
			-1, 0, &
			0, burst(2), &
			0, -1/), shape(abund_update))

! Variables
! ======================================================================
! Propensity of each event
real(dp), dimension(4) :: propensity = 0.0
! (initial) Abundances of each species, number of decay events for each species
integer, dimension(2) :: x = [0, 0], ndecay = [0, 0]
! Probability matrix
real(dp) :: prob_cond(abund_max, abund_max)
integer :: i, event, io
real(dp) :: t
character(*), parameter :: fout = "mrna_data.dat"

! Program start
! ======================================================================

call random_seed()

call random_uniform(alpha, 1._dp, 10._dp)
call random_uniform(beta, 1._dp, 10._dp)
call random_uniform(decay(1), 1._dp, 10._dp)
call random_uniform(decay(2), 1._dp, 10._dp)

write(*,*) alpha, beta, decay

open(newunit=io, file=fout, position="append", action="write")
write(io, "(9f20.14)", advance="no") alpha, beta, decay

prob_cond = 0._dp

do while (minval(ndecay) < decay_min)
	! If we go over what we expect to be the maximum, crash.
	if (maxval(x) >= abund_max) then
		write(*,*) "Maximum abundance exceeded."
		call exit()
	end if

	! Get time step and event from Gillespie algorithm
	call gillespie_iter(x, t, event)

	! Add time step to probability matrix
	prob_cond(x(1)+1, x(2)+1) = prob_cond(x(1)+1, x(2)+1) + t

	! Update abundances
	x = x + abund_update(:,event)

	! Add decay events to counter
	if (event .eq. 2 .or. event .eq. 4) then
		ndecay(event/2) = ndecay(event/2) + 1
	end if
end do

! Normalize probability
prob_cond = prob_cond / sum(prob_cond)

call checks(prob_cond)
close(io)


contains


subroutine checks(pij)
	real(dp), dimension(abund_max, abund_max), intent(in) :: pij
	real(dp), dimension(2, abund_max) :: p
	real(dp) :: mean(2)=0., covar(2,2)=0.
	integer labels(abund_max), i, j
	
	! Create probability distributions for x1 and x2 from joint
	! probability distribution, create labels while we're at it.
	do i = 1, abund_max
		labels(i) = i-1
		p(1,i) = sum(pij(i,:))
		p(2,i) = sum(pij(:,i))
	end do
	
	! Check means
	mean = matmul(p, labels)
	write(*,*) "Theoretical mean: ", alpha/decay(1), alpha*beta/(decay(1)*decay(2))
	write(*,*) "Simulation mean: ", mean(1), mean(2)
	
	write(*,*) 'Theory covar: ', &
		1./mean(1), &	! 11
		1./mean(1) * decay(2)/sum(decay), & ! 12
		1./mean(1) * decay(2)/sum(decay), & ! 21
		1./mean(2) + 1./mean(1) * decay(2)/sum(decay) ! 22
	! Check covariances
	covar(1,1) = dot_product(p(1,:), (labels-mean(1))**2) / mean(1)**2
	covar(2,2) = dot_product(p(2,:), (labels-mean(2))**2) / mean(2)**2
	do i = 1, abund_max
	do j = 1, abund_max
		covar(1,2) = covar(1,2) + pij(i,j)*(labels(i)-mean(1))*(labels(j)-mean(2))
	end do
	end do
	covar(1,2) = covar(1,2) / (mean(1)*mean(2))
	covar(2,1) = covar(1,2)
	
	write(*,*) "Simulation covar: ", covar(1,1), &
		covar(1,2), &
		covar(2,1), &
		covar(2,2)
		
	write(io, "(9f20.14)", advance="no") mean(1), mean(2), covar(1,1), covar(1,2), covar(2,1), covar(2,2)
end subroutine


subroutine gillespie_iter(xx, tstep, event)
	real(dp), intent(out) :: tstep
	integer, intent(inout) :: xx(2)
	integer, intent(out) :: event
	real(dp) :: psum, propsum, roll
	! Run one step of the Gillespie algorithm given the current state of the
	! system. Update the propensities, return time step and reaction

	! Update propensity
	call update_propensity(xx)
	propsum = sum(propensity)
	
	! Update time
	call random_exp(1./propsum, tstep)
	
	! Get reaction
	call random_number(roll)
	event = 0
	psum = 0.
	do while (psum < roll*propsum)
		event = event + 1
		psum = psum + propensity(event)
	end do
	
end subroutine


subroutine update_propensity(x)
	integer, intent(in) :: x(2)
	! Updates propensities depending on the state of the system
	propensity(1) = alpha		! Make x1 (mRNA)
	propensity(2) = x(1)*decay(1)	! Degrade x1 (mRNA)
	propensity(3) = beta*x(1)	! Make x2 (Protein)
	propensity(4) = x(2)*decay(2)	! Degrade x2 (Protein)
end subroutine


subroutine random_exp(l,u)
	real(dp), intent(out) :: u
	real(dp) l, x
	! Converts continuous random variable (0,1] to exponential random
	! variable with mean = l. Set up as a subroutine to match
	! intrinsic random_number.
	call random_number(x)
	! 1.0-x saves us if we manage to roll a 0.
	u = -l*log(1.0-x)
end subroutine


subroutine random_uniform(u,a,b)
	! Generate a uniformly distributed random variable a <= u < b.
	real(dp), intent(inout) :: u
	real(dp), intent(in) :: a, b
	
	call random_number(u)
	u = (b-a)*u + a
end subroutine

end program
