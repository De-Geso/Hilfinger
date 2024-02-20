program mrna_gean
use kind_parameters
use init_mrna_gene
implicit none

!! Hyperparameters
!! ======================================================================
!! Number of decay events before stopping
!integer, parameter :: decay_min = 10**6
!! Maximum abundances, pad this.
!integer, parameter :: abund_max = 2**4

!! Parameters
!! ======================================================================
!! Size of burst for [x0, x1]
!integer, parameter, dimension(2) :: burst = [1, 1]
!! x1 production rate
!real(dp) :: alpha = 1._dp
!! x2 production rate
!real(dp) :: beta = 1._dp
!! Decay rates
!real(dp), dimension(2) :: decay = [1._dp, 1._dp]
!real(dp) :: k = 1._dp
!real(dp) :: n = 1._dp
!! Abundance update matrix.
!integer, parameter, dimension(2,4) :: abund_update = &
!	reshape((/burst(1), 0, &
!			-1, 0, &
!			0, burst(2), &
!			0, -1/), shape(abund_update))

!! Variables
!! ======================================================================
!! Propensity of each event
!real(dp), dimension(4) :: propensity = 0.0
!! (initial) Abundances of each species, number of decay events for each species
!integer, dimension(2) :: x = [0, 0], ndecay = [0, 0]
!! Probability matrix
!real(dp) :: prob_abund(abund_max, abund_max)
!real(dp) :: prob_rate(abund_max)
!integer :: i, event, io
!real(dp) :: t
!character(*), parameter :: fout = "mrna_data.dat"

!! Program start
!! ======================================================================

call random_seed()

! call random_uniform(alpha, 1._dp, 10._dp)
! call random_uniform(beta, 1._dp, 10._dp)
! call random_uniform(decay(1), 1._dp, 10._dp)
! call random_uniform(decay(2), 1._dp, 10._dp)

write(*,*) alpha, beta, tau

open(newunit=io, file=fout, position="append", action="write")
write(io, "(9f20.14)", advance="no") alpha, beta, tau

prob_abund = 0._dp
prob_rate = 0._dp

do while (minval(ndecay) < decay_min)
	! If we go over what we expect to be the maximum, crash.
	if (maxval(x) >= abund_max) then
		write(*,*) "Maximum abundance exceeded."
		call exit()
	end if

	! Get time step and event from Gillespie algorithm
	call gillespie_iter(x, t, event)

	! Add time step to probability matrix
	prob_abund(x(1)+1, x(2)+1) = prob_abund(x(1)+1, x(2)+1) + t
	prob_rate(x(2)+1) = prob_rate(x(2)+1) + t
	
	! Update abundances
	x = x + abund_update(:,event)

	! Add decay events to counter
	if (event .eq. 2 .or. event .eq. 4) then
		ndecay(event/2) = ndecay(event/2) + 1
	end if
end do

! Normalize probability
write(*,*) sum(prob_rate), sum(prob_abund)
prob_rate = prob_rate / sum(prob_rate)
prob_abund = prob_abund / sum(prob_abund)

call checks(prob_abund)
close(io)


contains


subroutine checks(pij)
	real(dp), dimension(abund_max, abund_max), intent(in) :: pij
	real(dp), dimension(2, abund_max) :: p
	real(dp) :: theory_mean(2), mean(2)=0._dp, covar(2,2)=0._dp, mean_rate, theory_covar(2,2)=0._dp
	real(dp), dimension(abund_max) :: hill_values
	integer labels(abund_max), i, j
	
	! Create probability distributions for x1 and x2 from joint
	! probability distribution, create labels while we're at it.
	do i = 1, abund_max
		hill_values(i) = hill(i-1, k, n)
		labels(i) = i-1
		p(1,i) = sum(pij(i,:))
		p(2,i) = sum(pij(:,i))
	end do
	
	! Check means
	mean_rate = dot_product(hill_values, prob_rate)
	theory_mean(1) = mean_rate*tau(1)
	theory_mean(2) = theory_mean(1)*beta*tau(2)
	mean = matmul(p, labels)
	write(*,*) 'Theory mean: ', theory_mean
	write(*,*) 'Simulation mean: ', mean
	
	! Check covariances
	covar(1,1) = dot_product(p(1,:), (labels-mean(1))**2) / mean(1)**2
	covar(2,2) = dot_product(p(2,:), (labels-mean(2))**2) / mean(2)**2
	do i = 1, abund_max
	do j = 1, abund_max
		covar(1,2) = covar(1,2) + pij(i,j)*(labels(i)-mean(1))*(labels(j)-mean(2))
		theory_covar(1,1) = theory_covar(1,1) + pij(i,j)*(labels(i)-theory_mean(1))*(hill_values(j)-mean_rate)
		! Notice we use labels(j), because this is for x1.
		theory_covar(1,2) = theory_covar(1,2) + pij(i,j)*(labels(j)-theory_mean(2))*(hill_values(j)-mean_rate)
	end do
	end do
	theory_covar(1,1) = theory_covar(1,1)/(theory_mean(1)*mean_rate) + 1./theory_mean(1)
	theory_covar(1,2) = theory_covar(1,2)/(theory_mean(2)*mean_rate) * tau(2)/sum(tau) &
		+ theory_covar(1,1) * tau(1)/sum(tau)
	theory_covar(2,1) = theory_covar(1,2)
	theory_covar(2,2) = theory_covar(1,2) + 1./theory_mean(2)
	
	covar(1,2) = covar(1,2) / (mean(1)*mean(2))
	covar(2,1) = covar(1,2)
	
	write(*,*) 'Theory covar: ', &
		theory_covar(1,1), &
		theory_covar(1,2), &
		theory_covar(2,1), &
		theory_covar(2,2)
	
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
!	propensity(1) = alpha
	propensity(1) = alpha * hill(x(2), k, n)	! Make x1 (mRNA)
	propensity(2) = x(1)/tau(1)	! Degrade x1 (mRNA)
	propensity(3) = beta*x(1)	! Make x2 (Protein)
	propensity(4) = x(2)/tau(2)	! Degrade x2 (Protein)
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


pure function hill(x, k, n)
	real(dp), intent(in) :: k, n
	integer, intent(in) :: x
	real(dp) :: hill
	hill = 1._dp / (1. + (1._dp*x/k)**n)
end function


subroutine random_uniform(u,a,b)
	! Generate a uniformly distributed random variable a <= u < b.
	real(dp), intent(inout) :: u
	real(dp), intent(in) :: a, b
	
	call random_number(u)
	u = (b-a)*u + a
end subroutine

end program
