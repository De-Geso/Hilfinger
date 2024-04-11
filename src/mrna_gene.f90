program mrna_gean
use kind_parameters
use randf
use init_mrna_gene
implicit none

call random_seed(put=seed)
! call random_seed()

! Randomize variables when testing, if we so choose.
! call random_uniform(alpha, 1._dp, 5._dp)
! call random_uniform(beta, 1._dp, 5._dp)
! call random_uniform(tau(2), 0.2_dp, 1._dp)

! Start abundances near their averages. Runs slightly faster, less starting artifacts.
x(1) = alpha*tau(1); x(2) = alpha*tau(1)*beta*tau(2)

open(newunit=io, file=fout, position="append", action="write")
write(io, "(4f20.14)", advance="no") alpha, beta, tau
write(*,*) "alpha=", alpha, "beta=", beta, "Tau=", tau

prob_cond = 0._dp 
prob_rate = 0._dp

t = 0._dp
corr_mean = 0._dp

do while (minval(nevents) < event_min)
	! If we go over maximum abundance, stop.
	if (maxval(x) >= abund_max) then
		write(*,*) "Maximum abundance exceeded."
		call exit()
	end if

	! Get time step and event from Gillespie algorithm. Update propensity
	call gillespie_iter(x, tstep, event, propensity)
	t = t + tstep
	
	! Track the abundances and time in a window for correlations.
	xtail(:,:ntail-1) = xtail(:,2:ntail)
	xtail(:,ntail) = x(:)
	ttail(:ntail-1) = ttail(2:ntail)
	ttail(ntail) = t
	
	! Start calculating correlations once window is big enough.
	if (ttail(1) /= 0.0) then
 		call update_corr(corr_mean2, corr_mean(1,:), corr_mean(2,:), xtail(2,:), xtail(2,:), ttail, tstep)
	end if

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

! Get the mean abundances
mean = get_mean(prob_cond)

! Get the covariances
covar = get_covar(mean, prob_cond)

! Don't know total time until end.
corr_mean2 = corr_mean2 / t
corr_mean = corr_mean / t

do i = 1, corr_n
	! Combine the variances and means into the correlation (normalized by variance)
	corr(i) = (corr_mean2(i) - corr_mean(1,i)*corr_mean(2,1)) / &
			(corr_mean2(1)-corr_mean(1,i)*corr_mean(2,1))
!	This works if you normalize by the mean as well. This is because our covariance is normalized by the mean
!	but in the definition of the correlation, it is not.
!	corr(i) = (corr_mean2(i) - corr_mean(1,i)*corr_mean(2,1)) / (covar(1,2)*mean(1)*mean(2))
end do

call checks(prob_cond)
call dump()


! Functions and Subroutines ============================================
contains


subroutine update_corr(mean2, meanx, meany, y, x, tvec, dt)
! Iteratively updates correlation variables (covariance and mean) every time step.
	integer, intent(in) :: x(ntail), y(ntail)
	real(dp), intent(in) :: tvec(ntail), dt
	real(dp), intent(inout) :: mean2(corr_n), meanx(corr_n), meany(corr_n)
	real(dp) :: t, ta, tb
	integer :: i, j, ti, ita, itb
	
	! Crash if we aren't carrying enough information around
	if (tvec(ntail)-tvec(1) < lag_max .and. tvec(1) /= 0.) then
		write(*,*) "Error in update_corr: Required lag larger than recorded lag. Increase ntail."
		write(*,*) "Recorded time lag: ", tvec(ntail)-tvec(1)
		write(*,*) "Required time lag: ", lag_max
		call exit()
	end if
	
	! For the leading edge (no lag) we always take exactly one step.
	! Integrate a rectangle.
	meanx(1) = meanx(1) + x(ntail) * dt
	meany(1) = meany(1) + y(ntail) * dt
	mean2(1) = mean2(1) + x(ntail)*y(ntail) * dt
	
	do i = 2, corr_n
		! For points with lag, we have to check some things.
		! Range of time we're interested in.
		tb = tvec(ntail) - (i-1.)*corr_tstep
		ta = tb - dt
		t = ta
		! Most recent points smaller than ta and tb
		ita = maxloc(tvec, dim=1, mask=ta-tvec .gt. 0.)
		itb = maxloc(tvec, dim=1, mask=tb-tvec .gt. 0.)
		! Special case where no step occured in time range. Integrate a rectangle
		if (ita == itb) then
			mean2(i) = mean2(i) + x(ita+1) * y(ntail) * dt
			meanx(i) = meanx(i) + x(ita+1) * dt
			meany(i) = meany(i) + y(ntail) * dt
		else
			! ta side
			mean2(i) = mean2(i) + x(ita+1) * y(ntail) * (tvec(ita+1)-ta)
			meanx(i) = meanx(i) + x(ita+1) * (tvec(ita+1)-ta)
			meany(i) = meany(i) + y(ntail) * (tvec(ita+1)-ta)
			! tb side
			mean2(i) = mean2(i) + x(itb+1) * y(ntail) * (tb-tvec(itb))
			meanx(i) = meanx(i) + x(itb+1) * (tb-tvec(itb))
			meany(i) = meany(i) + y(ntail) * (tb-tvec(itb))
			
			! integrate from ta to tb
			do j = ita+1, itb-1
				mean2(i) = mean2(i) + x(j+1) * y(ntail) * (tvec(j+1)-tvec(j))
				meanx(i) = meanx(i) + x(j+1) * (tvec(j+1)-tvec(j))
				meany(i) = meany(i) + y(ntail) * (tvec(j+1)-tvec(j))
			end do
		end if	
	end do
end subroutine


function get_mean(p_cond) result(mean)
! Get the mean abundances.
	real(dp) :: mean(2)
	real(dp), intent(in) :: p_cond(abund_max, abund_max)
	real(dp) :: p(2, abund_max)
	integer :: xval(abund_max)
	
	do i = 1, abund_max
		xval(i) = i-1
		p(1,i) = sum(p_cond(i,:))
		p(2,i) = sum(p_cond(:,i))
	end do
	mean = matmul(p, xval)
end function

function get_covar(mean, p_cond) result(cov)
! Get the covariance matrix.
	real(dp) :: cov(2,2)
	real(dp), intent(in) :: p_cond(:,:), mean(2)
	real(dp) :: p(2,abund_max)
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


subroutine checks(pij)
	real(dp), dimension(abund_max, abund_max), intent(in) :: pij
	real(dp), dimension(2, abund_max) :: p
	real(dp) :: theory_mean(2), mean(2)=0._dp, covar(2,2)=0._dp, mean_rate, theory_covar(2,2)=0._dp
	real(dp), dimension(abund_max) :: rate_values
	integer labels(abund_max), i, j
	
	! Create probability distributions for x1 and x2 from joint
	! probability distribution, create labels while we're at it.
	do i = 1, abund_max
		rate_values(i) = alpha
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


subroutine gillespie_iter(xx, tstep, event, propensity)
! Run one step of the Gillespie algorithm given the current state of the
! system. Update the propensities, return time step and reaction
	real(dp), intent(out) :: tstep
	integer, intent(inout) :: xx(2)
	real(dp), intent(inout) :: propensity(4)
	integer, intent(out) :: event
	real(dp) :: psum, propsum, roll

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
	real(dp) :: propensity(4)
	propensity(1) = alpha	! Make x1 (mRNA)
	propensity(2) = x(1)/tau(1)	! Degrade x1 (mRNA)
	propensity(3) = beta*x(1)	! Make x2 (Protein)
	propensity(4) = x(2)/tau(2)	! Degrade x2 (Protein)
end function


subroutine dump()
	real(dp) :: t, tmax
	integer :: io, nt
	
	tmax = 10.
	nt = 100*tmax
	
	write(*,*) "Simulation mean: ", mean
	write(*,*) "Covariance matrix: "
	write(*,*) covar(1,1), covar(1,2)
	write(*,*) covar(2,1), covar(2,2)
	
	! Correlations from simulations
	open(newunit=io, file='corr_sim.dat', action='write')
	do i = 1, corr_n
		t = (i-1)*corr_tstep
		write(io,*) t, corr(i)
	end do
	close(io)
	
	! Correlations and their derivative from theory
	open(1, file='corr_thry.dat', action='write')
	open(2, file='dcorr_thry.dat', action='write')
	do i = 1, nt
		t = (i-1)*(tmax/(nt-1))
		write(1,*) t, corr_thry(t)
		write(2,*) t, dcorr_thry(t)
	end do
	close(1)
	close(2)
end subroutine


pure function corr_thry(t) result (corr)
	real(dp) :: corr(4)
	real(dp), intent(in) :: t
	! column major
	
	! Amm
	corr(1) = exp(-t/tau(1))
	! Amp
	if (tau(1) /= tau(2)) then
		corr(2) = exp(-t/tau(2)) + (exp(-t/tau(1))-exp(-t/tau(2)))*beta*tau(1)*tau(2) & 
				* covar(1,1)/covar(1,2) * mean(1)/mean(2) / (tau(1)-tau(2))
	else 
		corr(2) = exp(-t/tau(2)) + t*exp(-t/tau(2))*beta*tau(1) & 
				* covar(1,1)/covar(1,2) * mean(1)/mean(2) / tau(2)
	end if
	! Apm
	corr(3) = exp(-t/tau(1))
	! App
	if (tau(1) /= tau(2)) then
		corr(4) = exp(-t/tau(2)) + (exp(-t/tau(1))-exp(-t/tau(2)))*beta*tau(1)*tau(2) & 
				* covar(2,1)/covar(2,2) * mean(1)/mean(2) / (tau(1)-tau(2))
	else
		corr(4) = exp(-t/tau(2)) + t*exp(-t/tau(2))*beta*tau(1) & 
				* covar(2,1)/covar(2,2) * mean(1)/mean(2) / tau(2)
	end if
end function


pure function dcorr_thry(t) result (dcorr)
	real(dp) :: dcorr(4)
	real(dp), intent(in) :: t
	! column major
	
	! Amm
	dcorr(1) = -1/tau(1)*exp(-t/tau(1))
	! Apm
	if (tau(1) /= tau(2)) then
		dcorr(2) = -1./tau(2)*exp(-t/tau(2)) + (exp(-t/tau(2))/tau(2) - exp(-t/tau(1))/tau(1)) &
				*beta*tau(1)*tau(2) * covar(1,1)/covar(1,2) * mean(1)/mean(2) / (tau(1)-tau(2))
	else 
		dcorr(2) = exp(-t/tau(2)) * ( -1./tau(1) + &
				(1.-t/tau(1)) * beta * covar(1,1)/covar(1,2) * mean(1)/mean(2) )
	end if
	! Amp
	dcorr(3) = -1/tau(1)*exp(-t/tau(1))
	! App
	if (tau(1) /= tau(2)) then
		dcorr(4) = -1./tau(2)*exp(-t/tau(2)) + (exp(-t/tau(2))/tau(2)-exp(-t/tau(1))/tau(1)) &
				*beta*tau(1)*tau(2) * covar(2,1)/covar(2,2) * mean(1)/mean(2) / (tau(1)-tau(2))
	else
		dcorr(4) = exp(-t/tau(2)) * ( -1./tau(1) + &
				(1.-t/tau(1)) * beta * covar(2,1)/covar(2,2) * mean(1)/mean(2))
	end if
end function


!pure function R(x)
!! Rate function for x1 production
!	integer, dimension(2), intent(in) :: x
!	real(dp) :: R
!	! R = hill(x(2), k, n)
!	R = 1._dp
!end function


!pure function hill(x, k, n)
!! Hill function
!	real(dp), intent(in) :: k, n
!	integer, intent(in) :: x
!	real(dp) :: hill
!	hill = 1._dp / (1. + (1._dp*x/k)**n)
!end function


end program
