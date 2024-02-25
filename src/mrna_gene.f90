program mrna_gean
use kind_parameters
! use init_mrna_gene
use randf
implicit none
include "fftw3.f"

! Hyperparameters
! ======================================================================
real(dp), parameter :: pi = 4.D0*DATAN(1.D0)
! Number of decay events before stopping.
integer, parameter :: decay_min = 10**6
! Maximum abundances. Program will exit if this is exceeded.
integer, parameter :: abund_max = 2**8
! Number of abundance updates to remember for autocovariance.
integer, parameter :: nrec = 2**5
!! Length of correlation in memory.
! integer, parameter :: corr = 10**2

! Parameters
! ======================================================================
! Size of burst for [x0, x1]
integer, parameter, dimension(2) :: burst = [1, 1]
! x1 production rate
real(dp) :: alpha = 1._dp
! x2 production rate
real(dp) :: beta = 1._dp
! Decay rates
real(dp), dimension(2) :: tau = [1._dp, 1._dp]
! Hill function parameters
real(dp) :: k = 1._dp
real(dp) :: n = 1._dp
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
! Abundances of each species, number of decay events for each species
integer, dimension(2) :: x = [0, 0], ndecay = [0, 0]
! Probability matrix
real(dp), dimension(abund_max, abund_max) :: prob_cond
real(dp), dimension(2, abund_max) :: prob
real(dp) :: prob_rate(abund_max)
character(*), parameter :: fout = "mrna_accum.dat"
integer, dimension(2, nrec) :: xrec = 0
real(dp), dimension(nrec) :: trec = 0._dp
integer :: i, event, io
real(dp) :: t, dtmin=1.0

! Program begins
! ======================================================================
call random_seed()

! call random_uniform(alpha, 1._dp, 10._dp)
! call random_uniform(beta, 1._dp, 10._dp)
! call random_uniform(tau(1), 0.1_dp, 1._dp)
! call random_uniform(tau(2), 0.1_dp, 1._dp)

! write(*,*) alpha, beta, tau

open(newunit=io, file=fout, position="append", action="write")
write(io, "(4f20.14)", advance="no") alpha, beta, tau

prob_cond = 0._dp
prob_rate = 0._dp

do while (minval(ndecay) < decay_min)
	! If we go over maximum abundance, crash.
	if (maxval(x) >= abund_max) then
		write(*,*) "Maximum abundance exceeded."
		call exit()
	end if

	! Get time step and event from Gillespie algorithm
	call gillespie_iter(x, t, event)

	! Add time step to probability matrices
	prob_cond(x(1)+1, x(2)+1) = prob_cond(x(1)+1, x(2)+1) + t
	prob_rate(x(2)+1) = prob_rate(x(2)+1) + t
	
	! Update abundances
	x = x + abund_update(:,event)

	! Add decay events to counter
	if (event .eq. 2 .or. event .eq. 4) then
		ndecay(event/2) = ndecay(event/2) + 1
	end if
	
	! If we're getting close to the end of the simulation start tracking
	! the abundances for autocorrelation.
	if (minval(ndecay) > decay_min - nrec) then
		if (t .lt. dtmin) dtmin = t
		xrec(:,:nrec-1) = xrec(:,2:nrec)
		xrec(:,nrec) = x(:)
		trec(:nrec-1) = trec(2:nrec)
		trec(nrec) = trec(nrec-1) + t
	end if
end do
! Shift our time series so that our record is at t=0
trec = trec - trec(1)

! do i = 1, nrec
!	write(*,*) trec(i), xrec(:,i)
! end do

! Normalize probability
prob_rate = prob_rate / sum(prob_rate)
prob_cond = prob_cond / sum(prob_cond)

call autocorr_wiener_khinchin(trec, xrec(1,:), ceiling(trec(nrec)/dtmin))
! call checks(prob_cond)
close(io)


contains


subroutine autocorr_wiener_khinchin(tin, xin, n)
	real(dp), dimension(:), intent(in) :: tin
	integer, dimension(:), intent(in) :: xin
	integer, intent(in) :: n
	real(dp) :: x(n), Sxx(n), corr(n)
	complex(8), dimension(n) :: xout
	integer :: i, tindex
	real(dp) :: t, dt
	integer*8 :: plan
	
	dt = maxval(tin)/(n-1.)
	
	do i = 1, n
		t = (i-1.)*dt
		tindex = maxloc(tin, dim=1, mask=tin-t .lt. 1E-12)
		x(i) = xin(tindex)! * sin(pi*(i-1)/(n-1))**2
!		write(*,*) t, tin(tindex), x(i)
	end do
	
	call dfftw_plan_dft_r2c_1d(plan, n, x, xout, FFTW_ESTIMATE)
	call dfftw_execute_dft_r2c(plan, x, xout)
	call dfftw_destroy_plan(plan)
	Sxx = abs(xout)
	call dfftw_plan_dft_r2c_1d(plan, n, Sxx, xout, FFTW_ESTIMATE)
	call dfftw_execute_dft_r2c(plan, Sxx, xout)
	call dfftw_destroy_plan(plan)
	
	corr = abs(xout)
	
	do i = 1, n/2
		write(1,*) i/dt/n, Sxx(i), corr(i)
	end do
	
end subroutine


subroutine checks(pij)
	real(dp), dimension(abund_max, abund_max), intent(in) :: pij
	real(dp), dimension(2, abund_max) :: p
	real(dp) :: theory_mean(2), mean(2)=0._dp, covar(2,2)=0._dp, mean_rate, theory_covar(2,2)=0._dp
	real(dp), dimension(abund_max) :: rate_values
	integer labels(abund_max), i, j
	
	! Create probability distributions for x1 and x2 from joint
	! probability distribution, create labels while we're at it.
	do i = 1, abund_max
		rate_values(i) = alpha*R([0,i-1])
		labels(i) = i-1
		p(1,i) = sum(pij(i,:))
		p(2,i) = sum(pij(:,i))
	end do
	
	! Check means
	mean_rate = dot_product(rate_values, prob_rate)
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
		theory_covar(1,1) = theory_covar(1,1) + pij(i,j)*(labels(i)-theory_mean(1))*(rate_values(j)-mean_rate)
		! Notice we use labels(j), because this is for x1.
		theory_covar(1,2) = theory_covar(1,2) + pij(i,j)*(labels(j)-theory_mean(2))*(rate_values(j)-mean_rate)
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
		
	write(io, "(10f20.14)", advance="no") mean, theory_mean, &
		covar(1,1), theory_covar(1,1), &
		covar(1,2), theory_covar(1,2), &
		covar(2,2), theory_covar(2,2)
end subroutine


subroutine gillespie_iter(xx, tstep, event)
! Run one step of the Gillespie algorithm given the current state of the
! system. Update the propensities, return time step and reaction
	real(dp), intent(out) :: tstep
	integer, intent(inout) :: xx(2)
	integer, intent(out) :: event
	real(dp) :: psum, propsum, roll

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


pure function R(x)
! Rate function for x1 production
	integer, dimension(2), intent(in) :: x
	real(dp) :: R
	! R = hill(x(2), k, n)
	R = 1._dp
end function


subroutine update_propensity(x)
! Updates propensities depending on the state of the system
	integer, intent(in) :: x(2)
	propensity(1) = alpha * R(x)	! Make x1 (mRNA)
	propensity(2) = x(1)/tau(1)	! Degrade x1 (mRNA)
	propensity(3) = beta*x(1)	! Make x2 (Protein)
	propensity(4) = x(2)/tau(2)	! Degrade x2 (Protein)
end subroutine


pure function hill(x, k, n)
! Hill function
	real(dp), intent(in) :: k, n
	integer, intent(in) :: x
	real(dp) :: hill
	hill = 1._dp / (1. + (1._dp*x/k)**n)
end function


end program
