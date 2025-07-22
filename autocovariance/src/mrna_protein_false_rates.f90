program mrna_protein_false_rates
! mRNA-protein system with feedback. Here, we want to try incorrect
! predictions. For example, the system obeys p == alpha*m => p+1, but
! we can calculate a correlation for an incorrect rate such as alpha*m^2.
! Then we can check our derivative rule to see if it still holds.

use kind_parameters
use randf
use mrna_protein_system_parameters
use utilities
implicit none

! Program Hyperparameters ==============================================
real(dp), parameter :: eps = 1E-12_dp
! Number of events before stopping
integer, parameter :: event_min = 10**4
! Maximum abundance. Program will exit if exceeded.
integer, parameter :: abund_max = 2**9
! Number of abundance updates to remember for correlation. Reducing this gives big time savings.
integer, parameter :: nwindow = 500
! Number of points in correlation vector
integer, parameter :: ncorr = 2**7
! Maximum time lag for correlation vector
real(dp), parameter :: maxlag = 4._dp
! Time step for correlation
real(dp), parameter :: corr_tstep = 1._dp*maxlag/(ncorr-1)

! Fake parameters
real(dp) :: param(2) = [0._dp, 0._dp]
! Fake power of m and p dependence. This should be received as a command line argument.
! real(dp) :: l(2) = [ell(1), ell(2)]
! Fake hill function parameters. This should be received as a command line argument.
! real(dp) :: k_guess = k
! real(dp) :: n_guess = n


! Variables ============================================================
! Time
real(dp) :: t=0._dp, tcorr=0._dp, tstep
! Probabilities
real(dp) :: pcond(abund_max, abund_max)=0._dp
! Moments
real(dp) :: mean(3), meanR(6), eta(2,2), thry_mean(2), thry_eta(2,2)
! Correlation pieces. Amm, Apm, Amp, App, ARp+p, ARp-p
real(dp) :: corr(ncorr), corr_mean(2,ncorr)=0._dp, corr_mean2(ncorr)=0._dp, &
	twindow(nwindow)=0._dp
real(dp) :: xwindow(nwindow)=0
! Fake moments
real(dp) :: fake_mean(2)
! Fake correlation pieces. ARp+p, ARp-p
real(dp) :: fake_corr(ncorr,2), fake_corr_mean(2,ncorr,2)=0._dp, fake_corr_mean2(ncorr,2)=0._dp
real(dp) :: fake_Rwindow(2,nwindow)=0
! Miscellaneous 
real(dp) :: propensity(6), roll
integer :: i, j, mp(3)=0, nevents(6)=0, event, nseed
integer, allocatable :: rseed(:)

integer :: mp_max(3)=0

! ======================================================================

! nevents(1:2) = event_min+1

call get_command_line_arg(param(1), 1)
call get_command_line_arg(param(2), 2)

write(*,*) "parameter guess:", param

! Set a seed
call random_seed(put=seed)
! Generate a seed
! call random_seed()

! Get random seed for output in metadata
call random_seed(size=nseed)
allocate(rseed(nseed))
call random_seed(get=rseed)

! Randomize variables when testing, if we so choose.
! call random_uniform(roll, -1._dp, 1._dp)
! lmbda = 10._dp**roll
! call random_uniform(roll, -1._dp, 1._dp)
! alpha = 10._dp**roll
! call random_uniform(roll, -1._dp, 1._dp)
! beta(2) = 10._dp**roll

do while (minval(nevents) < event_min)
! 	write(*,*) nevents
	! Exit the program if we exceed maximum abundance.
	if (maxval(mp) >= abund_max) then
		write(*,*) "Maximum abundance exceeded."
		write(*,*) "lmbda=", lmbda, "alpha=", alpha, "Beta=", beta
		call exit()
	end if
	
	do i = 1, size(mp)
		if (mp(i) > mp_max(i)) then
			mp_max(i) = mp(i)
		end if
	end do
	
	if (event_min <= 10**4) then
		write(1,*) t, mp
	end if
	
	! Get timestep and event from Gillespie algorithm. Update propensity
	call gillespie_iter(mp, tstep, event, propensity)
	
	! Update time first to record the end times for an abundance, rather than the start times.
	t = t + tstep

	! Track the abundances and time in a window for correlations.
	xwindow(:nwindow-1) = xwindow(2:nwindow)
	xwindow(nwindow) = mp(3)
	twindow(:nwindow-1) = twindow(2:nwindow)
	twindow(nwindow) = t
		
	! Track the fake rates, m^l1, p^l2*beta_p
	fake_Rwindow(:, :nwindow-1) = fake_Rwindow(:, 2:nwindow)
	! fake_Rwindow(:, nwindow) = [Pbirth(real(mp,dp), param), Pdecay(real(mp,dp), param)]
	fake_Rwindow(:, nwindow) = [x3_production(mp, param), x3_decay(mp, param)]
	
	if (twindow(1) > eps) then
		! Record a special time to divide correlation means by. Otherwise we're off by ~20 seconds
		tcorr = tcorr + tstep
		!$omp parallel
		!$omp sections
			!$omp section
			! Autocovariances of components
			call update_correlation(corr_mean2, corr_mean(1,:), corr_mean(2,:), &
				xwindow, xwindow, twindow, tstep)
			!$omp section
			call update_correlation(fake_corr_mean2(:,1), fake_corr_mean(1,:,1), fake_corr_mean(2,:,1), &
				fake_Rwindow(1,:), xwindow, twindow, tstep)
			!$omp section
			! ARpdecay,p for the fake rate.
			call update_correlation(fake_corr_mean2(:,2), fake_corr_mean(1,:,2), fake_corr_mean(2,:,2), &
				fake_Rwindow(2,:), xwindow, twindow, tstep)
		!$omp end sections
		!$omp end parallel
	end if

	! Add time step to probability matrices
	pcond(mp(2)+1, mp(3)+1) = pcond(mp(2)+1, mp(3)+1) + tstep
	
	! Update abundances LAST
	mp = mp + abund_update(:,event)
	nevents(event) = nevents(event) + 1
end do

! Normalize probabilities.
write(*,*) t, sum(pcond)
pcond = pcond / sum(pcond)

! Get moments
call simulation_moments(pcond, nevents, mean, meanR, eta)
write(*,*) mean
call theory_moments(pcond, nevents, thry_mean, thry_eta)

! Normalize correlation means.
corr_mean2 = corr_mean2 / tcorr
corr_mean = corr_mean / tcorr

fake_corr_mean2 = fake_corr_mean2 / tcorr
fake_corr_mean = fake_corr_mean / tcorr


do i = 1, ncorr
	! Combine the variances and means into the correlation
	! These are not normalized for this program because I don't calculate sigma_Rp(m) explicitly.
	corr(i) = (corr_mean2(i) - corr_mean(1,i)*corr_mean(2,1))
end do

! Assemble guessed correlation
do i = 1, ncorr
	fake_corr(i,:) = fake_corr_mean2(i,:) - fake_corr_mean(1,i,:)*fake_corr_mean(2,1,:)
end do

call dump()



contains ! =============================================================


! Calculate simulation mean abundances and covariance matrix
! (eta--normalized) from conditional probability distribution.
subroutine simulation_moments(pcond, nevents, mean, meanR, eta)
	real(dp), intent(in) :: pcond(abund_max, abund_max)
	integer, intent(in) :: nevents(6)
	real(dp), intent(out) :: mean(3), meanR(6), eta(2,2)
	real(dp) :: p(abund_max, 3), x(abund_max)
	integer :: i, j
	
	! First moments
	! Abundance means
	do i = 1, abund_max
		x(i) = i-1
		p(i,1) = sum(pcond(i,:))
		p(i,2) = sum(pcond(i,:))
		p(i,3) = sum(pcond(:,i))
	end do
	mean = matmul(x,p)
	write(*,*) mean
	
	meanR = nevents/t
	
	! Second moments
!	eta = 0._dp
!	! Diagonal covariances are just variances.
!	eta(1,1) = dot_product(p(:,1), (x-mean(1))**2) / mean(1)**2
!	eta(2,2) = dot_product(p(:,2), (x-mean(2))**2) / mean(2)**2
!	! Off diagonal covariances are symmetric.
!	do i = 1, abund_max
!	do j = 1, abund_max
!		! These are etas. They are normalized by the means.
!		eta(1,2) = eta(1,2) + (pcond(i,j)*(x(i)-mean(1))*(x(j)-mean(2))) / (mean(1)*mean(2))
!	end do
!	end do
!	eta(2,1) = eta(1,2)
end subroutine


! Calculate the theoretical means and etas.
subroutine theory_moments(pcond, nevents, mean, eta)
	real(dp), intent(in) :: pcond(abund_max, abund_max)
	integer, intent(in) :: nevents(4)
	real(dp), intent(out) :: mean(2), eta(2,2)
	real(dp) :: meanR(3), p(abund_max, 2), s(2), x(abund_max)
	integer :: i, j
	
	do i = 1, abund_max
		x(i) = i-1
		p(i,1) = sum(pcond(i,:))
		p(i,2) = sum(pcond(:,i))
	end do
	mean = matmul(x,p)
	
	! Step sizes, if we ever want them
	s(1) = 1._dp*(burst(1)+1) / 2
	s(2) = 1._dp*(burst(2)+1) / 2
	
	! First moments
	
	! Mean rate
	! There isn't really an analytic mean rate except for the mrna, but we need something for the etas
	meanR = nevents(2:)/t
	meanR(1) = mean(1)*beta(1)

	! Mean abundances
	! Depending on the ells, there may not be an analytic solution here
	! mean(1) = 1._dp * burst(1)*meanR(1)
	! if (abs(ell(1) - 1._dp) .lt. eps .and. abs(ell(2) - 1._dp) .lt. eps) then
	! 	mean(2) = 1._dp * burst(2)*alpha*mean(1)/beta(2)
	! end if

	! Second moments
	eta = 0._dp
	
	do i = 1, abund_max
	do j = 1, abund_max
	!	eta(1,1) = eta(1,1) + (pcond(i,j) * (x(i) - mean(1)) * (R([x(i), x(j)]) - meanR(1)))
	!	eta(1,2) = eta(1,2) + (pcond(i,j) * (x(j) - mean(2)) * (R([x(i), x(j)]) - meanR(1)))
	end do
	end do
	
	eta(1,1) = eta(1,1)/(mean(1)*meanR(1)) + s(1)/mean(1)
	eta(1,2) = beta(2)/sum(beta) * eta(1,1) + beta(1)/sum(beta) * eta(1,2)/(mean(2)*meanR(1))
	eta(2,1) = eta(1,2)
	eta(2,2) = s(2)/mean(2) + eta(1,2)
end subroutine


subroutine update_correlation(mean2, meanx, meany, y, x, tvec, dt)
! Iteratively updates correlation variables (covariance and mean) every time step.
!	real(dp), parameter :: corr_tstep = 1._dp*maxlag/(ncorr-1)
	real(dp), intent(in) :: x(nwindow), y(nwindow)
	real(dp), intent(in) :: tvec(nwindow), dt
	real(dp), intent(inout) :: mean2(ncorr), meanx(ncorr), meany(ncorr)
	real(dp) :: t, ta, tb
	integer :: i, j, ti, ita, itb
	
	! Crash if we aren't carrying enough information around
	if (tvec(nwindow)-tvec(1) < maxlag .and. tvec(1) /= 0.) then
		write(*,*) "Error in update_corr: Required lag larger than recorded lag. Increase nwindow."
		write(*,*) "Recorded time lag: ", tvec(nwindow)-tvec(1)
		write(*,*) "Required time lag: ", maxlag
		write(*,*) "lmbda=", lmbda, "alpha=", alpha, "Beta=", beta
		call exit()
	end if
	
	! For the leading edge (no lag) we always take exactly one step.
	! Integrate a rectangle.
	meanx(1) = meanx(1) + x(nwindow) * dt
	meany(1) = meany(1) + y(nwindow) * dt
	mean2(1) = mean2(1) + x(nwindow)*y(nwindow) * dt
	
	do i = 2, ncorr
		! For points with lag, we have to check some things.
		! Range of time we're interested in.
		tb = tvec(nwindow) - (i-1.)*corr_tstep
		ta = tb - dt
		t = ta
		! Most recent points smaller than ta and tb
		ita = maxloc(tvec, dim=1, mask=ta-tvec .gt. 0.)
		itb = maxloc(tvec, dim=1, mask=tb-tvec .gt. 0.)
		! Special case where no step occured in time range. Integrate a rectangle
		if (ita == itb) then
			mean2(i) = mean2(i) + x(ita+1) * y(nwindow) * dt
			meanx(i) = meanx(i) + x(ita+1) * dt
			meany(i) = meany(i) + y(nwindow) * dt
		else
			! ta side
			mean2(i) = mean2(i) + x(ita+1) * y(nwindow) * (tvec(ita+1)-ta)
			meanx(i) = meanx(i) + x(ita+1) * (tvec(ita+1)-ta)
			meany(i) = meany(i) + y(nwindow) * (tvec(ita+1)-ta)
			! tb side
			mean2(i) = mean2(i) + x(itb+1) * y(nwindow) * (tb-tvec(itb))
			meanx(i) = meanx(i) + x(itb+1) * (tb-tvec(itb))
			meany(i) = meany(i) + y(nwindow) * (tb-tvec(itb))
			
			! integrate from ta to tb
			do j = ita+1, itb-1
				mean2(i) = mean2(i) + x(j+1) * y(nwindow) * (tvec(j+1)-tvec(j))
				meanx(i) = meanx(i) + x(j+1) * (tvec(j+1)-tvec(j))
				meany(i) = meany(i) + y(nwindow) * (tvec(j+1)-tvec(j))
			end do
		end if	
	end do
end subroutine


subroutine gillespie_iter(abund, tstep, event, propensity)
! Run one step of the Gillespie algorithm given the current state of the
! system. Update the propensities, return time step and reaction
	real(dp), intent(out) :: tstep
	integer, intent(inout) :: abund(3)
	real(dp), intent(inout) :: propensity(6)
	integer, intent(out) :: event
	real(dp) :: psum, propsum, roll

	! Update propensity
	propensity = update_propensity(abund)
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


pure function update_propensity(x) result(prop)
! Updates propensities depending on the state of the system
! LAMBDA AND ALPHA LIVE IN R AND Pbirth RESPECTIVELY.
	integer, intent(in) :: x(3)
	real(dp) :: prop(6)
	! Make x1
	! prop(1) = 1._dp * lmbda(1) * hill_neg(real(x(3),dp), k(1), n(1))
	prop(1) = 1._dp * lmbda(1)
	! Degrade x1
	prop(2) = 1._dp * x(1)
	! Make x2
	! prop(3) = 1._dp * lmbda(2) * hill_pos(real(x(1),dp), k(2), n(2))
	prop(3) = 1._dp * lmbda(2) * x(1)
	! Degrade x2
	prop(4) = 1._dp * x(2)
	! Degrade x3
	prop(5) = x3_production(x, [k(3), n(3)])
	! Degrade x3
	prop(6) = x3_decay(x, ell)
end function

pure function x3_production(x, params) result(f)
	integer, intent(in) :: x(3)
	real(dp), intent(in) :: params(2)
	real(dp) :: f
	! f = 1._dp * lmbda(3) * hill_pos(real(x(2),dp), params(1), params(2))
	f = 1._dp * lmbda(3) * x(2)!**params(1)
end function

pure function x3_decay(x, params) result(f)
	integer, intent(in) :: x(3)
	real(dp), intent(in) :: params(2)
	real(dp) :: f
	f = 1._dp * x(3)
	! f = 1._dp * x(3)**params(2)
end function

pure function hill_pos(x, k, n) result(f)
! Hill function. k and n controlled in system parameters.
	real(dp), intent(in) :: x, k, n
	real(dp) :: f
	f = 1._dp * x**n / (x**n + k**n)	
end function

pure function hill_neg(x, k, n) result(f)
! Hill function. k and n controlled in system parameters.
	real(dp), intent(in) :: x, k, n
	real(dp) :: f
	f = 1._dp / (1. + (x/k)**n)
end function


subroutine dump()
! Output results.
	real(dp) t, x
	integer :: i, io
	character(len=256) :: fname, prefix="infer_hill_"
	character(len=64) :: directory = "data/", suffix = ".dat"
	
	real(dp) :: pcond_23(abund_max, abund_max)
	
	write(prefix, '(A, F0.3, A, F0.3, A)') "infer_power_law_l1_", param(1), "_l2_", param(2), "_"
	! Generate a unique filename
	call generate_ISO_filename(directory, prefix, suffix, fname)
	

	! Console output
	write(*,*) "Events: ", event_min
	write(*,*) "lmbda=", lmbda, "alpha=", alpha, "Beta=", beta, "k=", k, "n=", n, "ell=", ell
	write(*,*) "Mean rates (m-, p+, p-): ", meanR
	write(*,*) "Theory mean abundances: ", thry_mean
	write(*,*) "Simulation mean abundances: ", mean
	write(*,*) "Abundance percent difference: ", &
			percent_difference(thry_mean(1), mean(1)), &
			percent_difference(thry_mean(2), mean(2))
	write(*,*) "Theory etas: "
	write(*,*) thry_eta
	write(*,*) "Simulation etas: "
	write(*,*) eta
	write(*,*) "Percent difference: "
	write(*,*) percent_difference(thry_eta(1,1), eta(1,1)), &
			percent_difference(thry_eta(2,1), eta(2,1)), &
			percent_difference(thry_eta(1,2), eta(1,2)), &
			percent_difference(thry_eta(2,2), eta(2,2))

	write(*,*) "File output at: ", fname
	
	open(newunit=io, file=fname, status='new', action='write')
	! Write metadata
	write(io,*) "# Program Metadata"
	write(io,*) "# Program: mrna_protein_false_rates.f90"
	write(io,*) "# Description: &
		mRNA-protein system simulation unnormalized correlations and derivative of App from correlation relations."
	write(io,*) "# Creation date: ", fdate()
	write(io,*) "# Seed: ", rseed
	write(io,*) "# Min events: ", event_min
	write(io,*) "# Maximum abundance: ", abund_max
	write(io,*) "# Correlation window size: ", nwindow
	write(io,*) "# Correlation vector length", ncorr
	write(io,*) "# Correlation max lag", maxlag
	write(io,*) ""
	write(io,*) "# Parameter Metadata"
	write(io,*) "# lmbda1: ", lmbda(1)
	write(io,*) "# lmbda2: ", lmbda(2)
	write(io,*) "# lmbda3: ", lmbda(3)
	write(io,*) "# alpha: ", alpha
	write(io,*) "# beta_m: ", beta(1)
	write(io,*) "# beta_p: ", beta(2)
	write(io,*) "# k1: ", k(1)
	write(io,*) "# k2: ", k(2)
	write(io,*) "# k3: ", k(3)
	write(io,*) "# n1: ", n(1)
	write(io,*) "# n2: ", n(2)
	write(io,*) "# n3: ", n(3)
	write(io,*) "# k_guess: ", param(1)
	write(io,*) "# n_guess: ", param(2)
	write(io,*) "# ell1: ", ell(1)
	write(io,*) "# ell2: ", ell(2)
	write(io,*) "# l1: ", param(1)
	write(io,*) "# l2: ", param(2)
	write(io,*) "# x1avg: ", mean(1)
	write(io,*) "# x2avg: ", mean(2)
	write(io,*) "# x3avg: ", mean(3)
	write(io,*) "# Rx1+avg: ", meanR(1)
	write(io,*) "# Rx1-avg: ", meanR(2)
	write(io,*) "# Rx2+avg: ", meanR(3)
	write(io,*) "# Rx2-avg: ", meanR(4)
	write(io,*) "# Rx3+avg: ", meanR(5)
	write(io,*) "# Rx3-avg: ", meanR(6)
	write(io,*) ""
	write(io,*) "# Data"
	write(io,*) "# Time, Amm, Apm, Amp, App, fake_ARpbirthp, fake_ARpdecayp"
		
	! Write simulation data
	do i = 1, ncorr
		t = (i-1)*corr_tstep
		write(io,*) t, corr(i), fake_corr(i,:), fake_corr(i,2) - fake_corr(i,1)
	end do
!	! Write probability matrix, but only the parts that aren't empty
!	write(io,*) ""
!	write(io,*) "# Probability Matrix"
!	pcond_23 = sum(pcond, DIM=3)
!	do i = 1, mp_max(2)+1
!		write(io,*) pcond_23(:mp_max(1)+1, i)
!	end do
!	close(io)
end subroutine
end program
