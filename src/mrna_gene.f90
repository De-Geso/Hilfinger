program mrna_gean
use kind_parameters
! use init_mrna_gene
use randf
use init_mrna_gene
implicit none
include "fftw3.f"

! Program begins
! ======================================================================
call random_seed()

! call random_uniform(alpha, 1._dp, 10._dp)
! call random_uniform(beta, 1._dp, 10._dp)
! call random_uniform(tau(1), 0.1_dp, 1._dp)
! call random_uniform(tau(2), 0.1_dp, 1._dp)

open(newunit=io, file=fout, position="append", action="write")
write(io, "(4f20.14)", advance="no") alpha, beta, tau

prob_cond = 0._dp
prob_rate = 0._dp

t = 0._dp
acorr_mean = 0._dp

do while (minval(ndecay) < decay_min)
	! If we go over maximum abundance, stop.
	if (maxval(x) >= abund_max) then
		write(*,*) "Maximum abundance exceeded."
		call exit()
	end if

	! Get time step and event from Gillespie algorithm
	call gillespie_iter(x, tstep, event)
	t = t + tstep
	
	! Track the abundances and time in a window for correlations.
	xtail(:,:ntail-1) = xtail(:,2:ntail)
	xtail(:,ntail) = x(:)
	ttail(:ntail-1) = ttail(2:ntail)
	ttail(ntail) = t
	
	call update_autocorr(acorr, xtail, ttail, acorr_mean, tstep, acorr_tstep)

	! Add time step to probability matrices
	prob_cond(x(1)+1, x(2)+1) = prob_cond(x(1)+1, x(2)+1) + tstep
	prob_rate(x(2)+1) = prob_rate(x(2)+1) + tstep
	
	! Update abundances
	x = x + abund_update(:,event)

	! Add decay events to counter
	if (event .eq. 2 .or. event .eq. 4) then
		ndecay(event/2) = ndecay(event/2) + 1
	end if
end do

! Don't know total time until end.
acorr = acorr / t
acorr_mean = acorr_mean / t
write(*,*) acorr(2), acorr_mean(2), acorr_mean(1)

do i = 1, acorr_n
	acorr(i) = (acorr(i) - acorr_mean(i)*acorr_mean(1))
end do

! This LOOKS like cheating, but it's actaully not. acoor(1) is exactly
! the normalization
acorr = acorr / acorr(1)

do i = 1, acorr_n
	t = (i-1)*acorr_tstep
	write(2,*) t, acorr(i), exp(-t/tau(1))
end do



! Normalize probability
prob_rate = prob_rate / sum(prob_rate)
prob_cond = prob_cond / sum(prob_cond)

call checks(prob_cond)
close(io)


! Functions and Subroutines ============================================
contains


subroutine update_autocorr(acorr, xvec, tvec, mean, tint, dt)
	integer, intent(in) :: xvec(2, ntail)
	real(dp), intent(in) :: tvec(ntail), tint, dt
	real(dp), intent(inout) :: acorr(acorr_n), mean(acorr_n)
	real(dp) :: t, ta, tb
	integer :: i, j, ti, tai, tbi
	
	! Crash if we aren't carrying enough information around
	if (tvec(ntail)-tvec(1) < lag_max .and. tvec(1) /= 0.) then
		write(*,*) tvec(ntail)-tvec(1)
		write(*,*) "Error in update_autocorr: Maximum lag larger than recorded time."
		call exit()
	end if
	
	! For the leading edge (no lag) we always take exactly one step.
	mean(1) = mean(1) + xvec(1,ntail) * tint
	acorr(1) = acorr(1) + xvec(1,ntail)**2 * tint
	
	do i = 2, acorr_n
		! For points with lag, we have to check some things.
		tb = tvec(ntail) - i*dt
		ta = tb - tint
		t = ta
		! Most recent points smaller than ta and tb
		tai = maxloc(tvec, dim=1, mask=ta-tvec .gt. 0.)
		tbi = maxloc(tvec, dim=1, mask=tb-tvec .gt. 0.)
		if (tai == tbi) then
			acorr(i) = acorr(i) + xvec(1,tai)*xvec(1,ntail) * tint
			mean(i) = mean(i) + xvec(1,tai) * tint
		else
			! ta side
			acorr(i) = acorr(i) + xvec(1,tai)*xvec(1,ntail) * (tvec(tai+1)-ta)
			mean(i) = mean(i) + xvec(1,tai) * (tvec(tai+1)-ta)
			! tb side
			acorr(i) = acorr(i) + xvec(1,tbi+1)*xvec(1,ntail) * (tb-tvec(tbi))
			mean(i) = mean(i) + xvec(1,tbi+1) * (tb-tvec(tbi))
			
			! integrate from ta to tb
			do j = tai+1, tbi-1
				acorr(i) = acorr(i) + xvec(1,j+1)*xvec(1,ntail) * (tvec(j+1)-tvec(j))
				mean(i) = mean(i) + xvec(1,j+1) * (tvec(j+1)-tvec(j))
			end do
		end if	
	end do
end subroutine


!subroutine autocorr_wiener_khinchin(tin, xin, n)
!	real(dp), dimension(:), intent(in) :: tin
!	integer, dimension(:), intent(in) :: xin
!	integer, intent(in) :: n
!	real(dp) :: x(n), Sxx(n), corr(n)
!	complex(8), dimension(n) :: xout
!	integer :: i, tindex
!	real(dp) :: t, dt
!	integer*8 :: plan
	
!	dt = maxval(tin)/(n-1.)
	
!	do i = 1, n
!		t = (i-1.)*dt
!		tindex = maxloc(tin, dim=1, mask=tin-t .lt. 1E-12)
!		! x(i) = xin(tindex)
!		x(i) = xin(tindex) * sin(pi*(i-1)/(n-1))**2
!		! write(*,*) t, tin(tindex), x(i), xin(tindex)
!	end do
	
!	call dfftw_plan_dft_r2c_1d(plan, n, x, xout, FFTW_ESTIMATE)
!	call dfftw_execute_dft_r2c(plan, x, xout)
!	call dfftw_destroy_plan(plan)
!	Sxx = abs(xout)
	
!	call dfftw_plan_dft_1d(plan, n/2, Sxx, corr, FFTW_BACKWARD, FFTW_ESTIMATE)
!	call dfftw_execute_dft(plan, Sxx, corr)
!	call dfftw_destroy_plan(plan)
	
!	corr = corr/maxval(corr)
	
!	do i = 1, n/2
!		write(1,*) i*dt, Sxx(i), corr(i), exp(-i*dt/tau(1))
!	end do
	
!end subroutine


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
