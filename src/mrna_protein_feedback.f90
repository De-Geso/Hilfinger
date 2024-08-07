program mrna_protein_feedback
use kind_parameters
use randf
use mrna_protein_system_parameters
use utilities
implicit none

! Program Hyperparameters ==============================================
real(dp), parameter :: eps = 1E-16_dp
! Number of events before stopping
integer, parameter :: event_min = 10**6
! Maximum abundance. Program will exit if exceeded.
integer, parameter :: abund_max = 2**9

! Number of abundance updates to remember for correlation. Reducing this gives big time savings.
integer, parameter :: nwindow = 2**9
! Number of points in correlation vector
integer, parameter :: ncorr = 2**6
! Maximum time lag for correlation vector
real(dp), parameter :: maxlag = 5._dp
! Time step for correlation
real(dp), parameter :: corr_tstep = 1._dp*maxlag/(ncorr-1)

! Variables ============================================================
! Time
real(dp) :: t=0._dp, tcorr=0._dp, tstep
! Probabilities
real(dp) :: pcond(abund_max, abund_max)=0._dp
! Moments
real(dp) :: mean(2), meanR, cov(2,2), thry_mean(2), thry_cov(2,2)
! Correlation pieces
real(dp) :: corr(ncorr,4), corr_mean(2,ncorr,4)=0._dp, corr_mean2(ncorr,4)=0._dp, &
	twindow(nwindow)=0._dp
real(dp) :: dcorr(ncorr)
integer :: xwindow(2, nwindow)=0
! Miscellaneous 
real(dp) :: propensity(4), roll
integer :: i, j, mp(2)=0, nevents(4)=0, event, nseed
integer, allocatable :: rseed(:)


! call random_seed(put=seed)
call random_seed()

! Get random seed for output in metadata
call random_seed(size=nseed)
allocate(rseed(nseed))
call random_seed(get=rseed)

! Randomize variables when testing, if we so choose.
call random_uniform(roll, -1._dp, 1._dp)
lmbda = 10._dp**roll
call random_uniform(roll, -1._dp, 1._dp)
alpha = 10._dp**roll
call random_uniform(roll, -1._dp, 1._dp)
tau(2) = 10._dp**roll

do while (minval(nevents) < event_min)
	! Exit the program if we exceed maximum abundance.
	if (maxval(mp) >= abund_max) then
		write(*,*) "Maximum abundance exceeded."
		write(*,*) "lmbda=", lmbda, "alpha=", alpha, "Tau=", tau
		call exit()
	end if
	
	! Get timestep and event from Gillespie algorithm. Update propensity
	call gillespie_iter(mp, tstep, event, propensity)
	! Update time first to record the end times for an abundance, rather than the start times.
	t = t + tstep

	! Track the abundances and time in a window for correlations.
	xwindow(:,:nwindow-1) = xwindow(:,2:nwindow)
	xwindow(:,nwindow) = mp(:)
	twindow(:nwindow-1) = twindow(2:nwindow)
	twindow(nwindow) = t
	
	if (twindow(1) > eps) then
		tcorr = tcorr + tstep
		!$omp parallel
		!$omp sections
		!$omp section
		call update_correlation(corr_mean2(:,1), corr_mean(1,:,1), corr_mean(2,:,1), &
			xwindow(1,:), xwindow(1,:), twindow, tstep)
		!$omp section
		call update_correlation(corr_mean2(:,2), corr_mean(1,:,2), corr_mean(2,:,2), &
			xwindow(2,:), xwindow(1,:), twindow, tstep)
		!$omp section
		call update_correlation(corr_mean2(:,3), corr_mean(1,:,3), corr_mean(2,:,3), &
			xwindow(1,:), xwindow(2,:), twindow, tstep)
		!$omp section
		call update_correlation(corr_mean2(:,4), corr_mean(1,:,4), corr_mean(2,:,4), &
			xwindow(2,:), xwindow(2,:), twindow, tstep)
		!$omp end sections
		!$omp end parallel
	end if

	! Add time step to probability matrices
	pcond(mp(1)+1, mp(2)+1) = pcond(mp(1)+1, mp(2)+1) + tstep
	
	! Update abundances LAST
	mp = mp + abund_update(:,event)
	nevents(event) = nevents(event) + 1
end do

! Normalize probabilities.
pcond = pcond / sum(pcond)

! Get moments
call simulation_moments(pcond, mean, meanR, cov)
call theory_moments(pcond, thry_mean, thry_cov)

! Normalize correlation means.
corr_mean2 = corr_mean2 / tcorr
corr_mean = corr_mean / tcorr

do j = 1, 4
	do i = 1, ncorr
		! Combine the variances and means into the correlation (no longer normalized)
		corr(i,j) = (corr_mean2(i,j) - corr_mean(1,i,j)*corr_mean(2,1,j))
		!		/ (corr_mean2(1,j)-corr_mean(1,i,j)*corr_mean(2,1,j))
	end do
end do
dcorr = -1._dp/tau(2) * (corr(:,4) - corr(:,3)*cov(2,1)/cov(2,2))

call dump()


contains


subroutine simulation_moments(pcond, mean, meanR, cov)
! Calculate mean abundances and covariance matrix (eta--normalized)
! from conditional probability distribution.
	real(dp), intent(out) :: mean(2), meanR, cov(2,2)
	real(dp), intent(in) :: pcond(abund_max, abund_max)
	real(dp) :: p(abund_max, 2), x(abund_max)
	integer :: i, j
	
	! First moments
	! Abundance means
	do i = 1, abund_max
		x(i) = i-1
		p(i,1) = sum(pcond(i,:))
		p(i,2) = sum(pcond(:,i))
	end do
	mean = matmul(x,p)
	
	! Rate mean
	meanR = 0._dp
	do i = 1, abund_max
	do j = 1, abund_max
		meanR = meanR + pcond(i,j)*R([x(i),x(j)])
	end do
	end do
	
	! Second moments
	cov = 0._dp
	! Diagonal covariances are just variances.
	cov(1,1) = dot_product(p(:,1), (x-mean(1))**2) / mean(1)**2
	cov(2,2) = dot_product(p(:,2), (x-mean(2))**2) / mean(2)**2
	! Off diagonal covariances are symmetric.
	do i = 1, abund_max
	do j = 1, abund_max
		! These are etas. They are normalized by the means.
		cov(1,2) = cov(1,2) + (pcond(i,j)*(x(i)-mean(1))*(x(j)-mean(2))) / (mean(1)*mean(2))
	end do
	end do
	cov(2,1) = cov(1,2)
end subroutine


subroutine theory_moments(pcond, mean, cov)
! Calculate the theoretical means and covariances.
	real(dp), intent(in) :: pcond(abund_max, abund_max)
	real(dp), intent(out) :: mean(2), cov(2,2)
	real(dp) :: meanR, p(abund_max, 2), s(2), x(abund_max)
	integer :: i, j
	
	do i = 1, abund_max
		x(i) = i-1
		p(i,1) = sum(pcond(i,:))
		p(i,2) = sum(pcond(:,i))
	end do
	
	! Step sizes, if we ever want them
	s(1) = 1._dp*(burst(1)+1) / 2
	s(2) = 1._dp*(burst(2)+1) / 2
	
	! First moment
	! Mean rate
	meanR = 0._dp
	do i = 1, abund_max
	do j = 1, abund_max
		meanR = meanR + pcond(i,j)*R([x(i),x(j)])
	end do
	end do
	
	! Mean abundances
	mean(1) = 1._dp * burst(1)*tau(1)*meanR
	mean(2) = 1._dp * burst(2)*alpha*tau(2)*mean(1)
	
	! Second moments
	cov = 0._dp
	
	do i = 1, abund_max
	do j = 1, abund_max
		cov(1,1) = cov(1,1) + (pcond(i,j) * (x(i) - mean(1)) * (R([x(i), x(j)]) - meanR))
		cov(1,2) = cov(1,2) + (pcond(i,j) * (x(j) - mean(2)) * (R([x(i), x(j)]) - meanR))
	end do
	end do
	cov(1,1) = cov(1,1)/(mean(1)*meanR) + s(1)/mean(1)
	cov(1,2) = tau(1)/sum(tau) * cov(1,1) + tau(2)/sum(tau) * cov(1,2)/(mean(2)*meanR)
	cov(2,1) = cov(1,2)
	cov(2,2) = s(2)/mean(2) + cov(1,2)
end subroutine


subroutine update_correlation(mean2, meanx, meany, y, x, tvec, dt)
! Iteratively updates correlation variables (covariance and mean) every time step.
!	real(dp), parameter :: corr_tstep = 1._dp*maxlag/(ncorr-1)
	integer, intent(in) :: x(nwindow), y(nwindow)
	real(dp), intent(in) :: tvec(nwindow), dt
	real(dp), intent(inout) :: mean2(ncorr), meanx(ncorr), meany(ncorr)
	real(dp) :: t, ta, tb
	integer :: i, j, ti, ita, itb
	
	! Crash if we aren't carrying enough information around
	if (tvec(nwindow)-tvec(1) < maxlag .and. tvec(1) /= 0.) then
		write(*,*) "Error in update_corr: Required lag larger than recorded lag. Increase nwindow."
		write(*,*) "Recorded time lag: ", tvec(nwindow)-tvec(1)
		write(*,*) "Required time lag: ", maxlag
		write(*,*) "lmbda=", lmbda, "alpha=", alpha, "Tau=", tau
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
	integer, intent(inout) :: abund(2)
	real(dp), intent(inout) :: propensity(4)
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
! LAMBDA LIVES IN R.
	integer, intent(in) :: x(2)
	real(dp) :: prop(4)
	prop(1) = 1._dp * R(real(x,dp))	! Make mRNA
	prop(2) = 1._dp * x(1)/tau(1)	! Degrade mRNA
	prop(3) = 1._dp * alpha * x(1)	! Make protein
	prop(4) = 1._dp * x(2)/tau(2)	! Degrade protein
end function


pure function R(x) result(f)
! Production function for mRNA.
	real(dp), intent(in) :: x(2)
	real(dp) :: f
	associate(m => x(1), p => x(2))
	f = 1._dp * hill(p)
	! Multiply by lambda
	f = lmbda * f
	end associate
end function


pure function hill(x) result(f)
! Hill function. k and n controlled in system parameters.
	real(dp), intent(in) :: x
	real(dp) :: f
	
	f = 1._dp / (1. + (x/k)**n)
end function


subroutine dump()
! Output results.
	real(dp) t, x
	integer :: i, io, fnum
	character(len=*), parameter :: path='data/'
	character(len=:), allocatable :: prefix, suffix
	character(len=256) :: fname
	logical :: ex
	
	! Check file existances, and relabel with correct numbering
	prefix = "feedback_correlation_"
	suffix = ".dat"
	fnum = 0
	! Make filename. Check if it already exists. Increase number on filename.
	write (fname, '(a, a, i0, a)') path, prefix, fnum, suffix
	inquire(file=fname, exist=ex)
	do while (ex)
		fnum = fnum + 1
		write (fname, '(a, a, i0, a)') path, prefix, fnum, suffix
		inquire(file=fname, exist=ex)
	end do
	write(*,*) "File output at: ", fname

	! Console output
	write(*,*) "Events: ", event_min
	write(*,*) "lmbda=", lmbda, "alpha=", alpha, "Tau=", tau
	write(*,*) "Mean rate: ", meanR
	write(*,*) "Theoretical mean: ", thry_mean
	write(*,*) "Simulation mean: ", mean
	write(*,*) "Percent difference: ", &
			percent_difference(thry_mean(1), mean(1)), &
			percent_difference(thry_mean(2), mean(2))
	write(*,*) "Theoretical covariance matrix: "
	write(*,*) thry_cov
	write(*,*) "Simulation covariance matrix (eta): "
	write(*,*) cov
	write(*,*) "Percent difference: "
	write(*,*) percent_difference(thry_cov(1,1), cov(1,1)), &
			percent_difference(thry_cov(2,1), cov(2,1)), &
			percent_difference(thry_cov(1,2), cov(1,2)), &
			percent_difference(thry_cov(2,2), cov(2,2))

	! Correlations from simulation
	open(newunit=io, file=fname, action='write')
	call write_metadata_sim(io, &
		"Simulation unnormalized correlations and derivative of App from correlation relations.", &
		"Time, Amm, Apm, Amp, App, dApp")

	do i = 1, ncorr
		t = (i-1)*corr_tstep
		write(io,*) t, corr(i,:), dcorr(i)
	end do
	close(io)
end subroutine


subroutine write_metadata_sim(io, desc, headers)
! Write metadata to file
	integer, intent(in) :: io
	character(len=*), intent(in) :: desc, headers
	write(io,*) "# Program Metadata"
	write(io,*) "# Program: mrna_protein_feedback.f90"
	write(io,*) "# Description: ", desc
	write(io,*) "# Creation date: ", fdate()
	write(io,*) "# Seed: ", rseed
	write(io,*) "# Min events: ", event_min
	write(io,*) "# Maximum abundance: ", abund_max
	write(io,*) "# Correlation window size: ", nwindow
	write(io,*) "# Correlation vector length", ncorr
	write(io,*) "# Correlation max lag", maxlag
	write(io,*) ""
	write(io,*) "# Parameter Metadata"
	write(io,*) "# lmbda: ", lmbda
	write(io,*) "# alpha: ", alpha
	write(io,*) "# tau_m: ", tau(1)
	write(io,*) "# tau_p: ", tau(2)
	write(io,*) "# k: ", k
	write(io,*) "# n: ", n
	write(io,*) "# mavg: ", mean(1)
	write(io,*) "# pavg: ", mean(2)
	write(io,*) "# Ravg: ", meanR
	write(io,*) "# etamm: ", cov(1,1)
	write(io,*) "# etamp: ", cov(1,2)
	write(io,*) "# etapm: ", cov(2,1)
	write(io,*) "# etapp: ", cov(2,2)
	write(io,*) ""
	write(io,*) "# Data"
	write(io,*) "# ", headers
end subroutine

end program
