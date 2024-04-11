program mrna_protein_bench
! Compile with:
! gfortran -O3 -fdefault-real-8 -o gillespie_speed_test.out gillespie_speed_test.f90
implicit none

! Number of events before stopping
integer, parameter :: event_min = 10**7
! Maximum abundances. Program will exit if this is exceeded
integer, parameter :: abund_max = 2**6

! Parameters
! ======================================================================
! Size of burst for [x0, x1]
integer, parameter, dimension(2) :: burst = [1, 1]
! x1 production rate
real :: alpha = 7.2
! x2 production rate
real :: beta = 10.
! Decay rates. Can always leave tau_1=1
real, dimension(2) :: tau = [1., 0.1]
! Abundance update matrix.
integer, parameter, dimension(2,4) :: abund_update = &
	reshape((/burst(1), 0, &
			-1, 0, &
			0, burst(2), &
			0, -1/), shape(abund_update))
			
! Variables
! ======================================================================
! Propensity of each event
real, dimension(4) :: propensity = 0.0
! Abundances of each species, number of decay events for each species
integer, dimension(2) :: x = [0, 0], nevents(4)=0
! Probability matrices
real :: prob_cond(abund_max, abund_max), prob(2, abund_max), prob_rate(abund_max)
! Moments
real :: mean(2), cov(2,2), mean_thry(2), cov_thry(2,2)
! Timers
real :: t, tstep
integer :: i, event, io



call random_seed()

prob_cond = 0.
prob_rate = 0.

t = 0.

do while (minval(nevents) < event_min)
	! If we go over maximum abundance, stop.
	if (maxval(x) >= abund_max) then
		write(*,*) "Maximum abundance exceeded."
		call exit()
	end if

	! Get time step and event from Gillespie algorithm. Update propensity
	call gillespie_iter(x, tstep, event, propensity)
	t = t + tstep

	! Add time step to probability matrices
	prob_cond(x(1)+1, x(2)+1) = prob_cond(x(1)+1, x(2)+1) + tstep
	prob_rate(x(2)+1) = prob_rate(x(2)+1) + tstep
	
	! Update abundances
	x = x + abund_update(:,event)

	! Add decay events to counter
	nevents(event) = nevents(event) + 1
end do

! Simulation ends. =====================================================
! Now do calculations.

! Normalize probability
prob_rate = prob_rate / sum(prob_rate)
prob_cond = prob_cond / sum(prob_cond)


! Get the simulation and theoretical moments
mean = mean_sim(prob_cond)
cov = covariance_sim(mean, prob_cond)
call moments_theory(mean_thry, cov_thry)

write(*,*) "alpha =", alpha, "beta =", beta, "tau_m =", tau(1), "tau_p =", tau(2)
write(*,*) "Minimum number of reactions = ", event_min
write(*,*) '<m>', mean(1), 'eta-mm', cov(1,1)
write(*,*) '<p>', mean(2), 'eta-pp', cov(2,2)



contains



function mean_sim(p_cond) result(mean)
! Get the mean abundances.
	real :: mean(2)
	real, intent(in) :: p_cond(abund_max, abund_max)
	real :: p(2, abund_max)
	integer :: xval(abund_max)
	
	do i = 1, abund_max
		xval(i) = i-1
		p(1,i) = sum(p_cond(i,:))
		p(2,i) = sum(p_cond(:,i))
	end do
	mean = matmul(p, xval)
end function


function covariance_sim(mean, p_cond) result(cov)
! Get the coviance matrix.
	real :: cov(2,2)
	real, intent(in) :: p_cond(:,:), mean(2)
	real :: p(2,abund_max)
	integer :: xval(abund_max), i, j
	
	! Create probability distributions for x1 and x2 from joint probability distribution.
	do i = 1, abund_max
		xval(i) = i-1
		p(1,i) = sum(p_cond(i,:))
		p(2,i) = sum(p_cond(:,i))
	end do
	
	! Diagonal covariances are just variances.
	cov(1,1) = dot_product(p(1,:), (xval-mean(1))**2) / mean(1)**2
	cov(2,2) = dot_product(p(2,:), (xval-mean(2))**2) / mean(2)**2
	
	! Off diagonal covariances are symmetric.
	do i = 1, abund_max
	do j = 1, abund_max
		cov(1,2) = cov(1,2) + (p_cond(i,j)*(xval(i)-mean(1))*(xval(j)-mean(2))) / (mean(1)*mean(2))
	end do
	end do
	cov(2,1) = cov(1,2)
end function


subroutine moments_theory(mean, cov)
	real, intent(inout) :: mean(2), cov(2,2)
	
	mean(1) = alpha*tau(1)
	mean(2) = mean(1)*beta*tau(2)
	
	cov(1,1) = 1./mean(1)
	cov(1,2) = cov(1,1) * tau(1)/sum(tau)
	cov(2,1) = cov(1,2)
	cov(2,2) = cov(1,2) + 1./mean(2)
end subroutine


subroutine gillespie_iter(xx, tstep, event, propensity)
! Run one step of the Gillespie algorithm given the current state of the
! system. Update the propensities, return time step and reaction
	real, intent(out) :: tstep
	integer, intent(inout) :: xx(2)
	real, intent(inout) :: propensity(4)
	integer, intent(out) :: event
	real :: psum, propsum, roll

	! Update propensity
	propensity = update_propensity(xx)
	propsum = sum(propensity)
	
	! Get time step time
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


pure function update_propensity(x) result(propensity)
! Updates propensities depending on the state of the system
	integer, intent(in) :: x(2)
	real :: propensity(4)
	propensity(1) = alpha	! Make x1 (mRNA)
	propensity(2) = x(1)/tau(1)	! Degrade x1 (mRNA)
	propensity(3) = beta*x(1)	! Make x2 (Protein)
	propensity(4) = x(2)/tau(2)	! Degrade x2 (Protein)
end function


subroutine random_exp(L, u)
	! Converts continuous random variable (0,1] to exponential random
	! variable with mean = L.
	real, intent(out) :: u
	real :: L, x
	call random_number(x)
	! 1.0-x saves us if we manage to roll a 0.
	u = -L*log(1.0-x)
end subroutine


end program
