program cascade_2D_discrete

use kind_parameters
use stochastics
use randf
use utilities

implicit none

! Program hyper parameters =============================================
! Output filename
character(len=*), parameter :: filename = "hill_2D"
character(len=*), parameter :: directory = "data/"
! Epsilon
real(dp), parameter :: eps = 1E-14
! Number of events of each reaction before stopping.
integer, parameter :: event_min = 10**7
! Maximum abundance in each direction
integer, parameter, dimension(2) :: abund_max = [256, 256]
! Maximum time lag for autocovariance.
real(dp), parameter :: umax = 3._dp
! Discretization of time
real(dp), parameter :: tdisc = 0.1_dp
! Number of points in autocovariance.
integer, parameter :: nacov = ceiling(umax/tdisc) + 1


! System parameters ====================================================
! Production rates
real(dp), parameter, dimension(2) :: lmbda = [1._dp, 20._dp]
! Decay rates
real(dp), parameter, dimension(2) :: beta = [1._dp, 1._dp]
! Hill function parameters
real(dp), parameter, dimension(2) :: k = [1._dp, 10._dp]
real(dp), parameter, dimension(2) :: n = [0._dp, -2._dp]
! Constant offset
real(dp), parameter, dimension(2) :: c = [0._dp, 0._dp]
! Abundance update matrix.
integer, parameter, dimension(2, 4) :: abund_update = &
	reshape((/1, 0, &
			-1, 0, &
			0, 1, &
			0, -1 &
			/), shape(abund_update))

! Autocovariance variables =============================================
! Windows
real(dp) :: guess_Rwindow(2, nacov) = 0._dp
integer(dp) :: xwindow(nacov) = 0
! Covariance arrays
real(dp) :: acov(3, nacov)=0._dp, acov_mean(3, 2, nacov)=0._dp, acov_mean2(3, nacov)=0._dp

! Variables ============================================================
! Timers
real(dp) :: t = 0._dp, tstep, tlast=0._dp, tstart=-1._dp
! Probability distribution
real(dp), dimension(abund_max(1), abund_max(2)) :: pcond = 0._dp
! Moments
type(OnlineCovariance) :: stats(2,2), covxx
real(dp) :: mean_x(2), mean_r(4), cov(2,2)
! Gillespie
real(dp) :: propensity(4), roll
integer :: x(2)=0, xlast(2), event_count(4) = 0, event
! Counters
integer :: i, j, nseed, tcount=0
integer, allocatable :: rseed(:)
! Guessed k and n
real(dp) :: guess_params(2) = [k(2), n(2)]


! Here the program begins ==============================================

call get_command_line_arg(guess_params(1), 1)
call get_command_line_arg(guess_params(2), 2)

call random_seed(put=seed)
! Get random seed for output in metadata
call random_seed(size=nseed)
allocate(rseed(nseed))
call random_seed(get=rseed)

call dump_params()

xlast = x

do while (minval(event_count) < event_min)
	do i = 1, size(abund_max)
		if (x(i) .eq. abund_max(i)) then
			write(*,"(A, I0, A, I0)") "Element ", i, " exceeded maximum abundance ", abund_max(i)
			call exit()
		end if
	end do
	
	! Update the propensity before taking a Gillespie step
	call update_propensity(propensity, x, lmbda, k, n, c, beta)
	! Get timestep and event from Gillespie algorithm
	call gillespie_iter(tstep, event, propensity)
	! Update time by adding how long we were in the previous state
	t = t + tstep
	
	do i = 1, 2
	do j = 1, 2
		call stats(i,j)%update(real(x(i),dp), propensity(2*j) - propensity(2*j-1), tstep)
	end do
	end do
	call covxx%update(real(x(2),dp), real(x(2),dp), tstep)

	if (t .ge. (tlast + tdisc)) then
		! Online correlation calculation
		! Track the abundances and time in a window for correlations.
		do i = 1, floor((t-tlast)/tdisc)
			! Record a special time to divide autocovariance means by.
			tcount = tcount + 1
			tlast = 1._dp * tdisc * tcount

			xwindow(:nacov-1) = xwindow(2:nacov)
			xwindow(nacov) = x(2)

			! Track the guessed rates
			guess_Rwindow(:, :nacov-1) = guess_Rwindow(:, 2:nacov)
			guess_Rwindow(:, nacov) = guess_R(x, lmbda(2), guess_params(1), guess_params(2), c(2), beta(2))
			
			! Online calculation of autocovariance
			if (t .gt. umax) then
				if (tstart < 0._dp) tstart = t
				! We don't start the autocovariance timer at the start of the simulation.
				!$omp parallel
				!$omp sections
					!$omp section
					! Autocovariances of component of interest
					call update_discretized_autocovariance(acov_mean2(1,:), acov_mean(1,1,:), acov_mean(1,2,:), &
						real(xwindow, dp), real(xwindow, dp), nacov, tdisc)
					!$omp section
					! Autocovariances of Rin and component
					call update_discretized_autocovariance(acov_mean2(2,:), acov_mean(2,1,:), acov_mean(2,2,:), &
						guess_Rwindow(1,:), real(xwindow, dp), nacov, tdisc)
					!$omp section
					! Autocovariances of Rout and component
					call update_discretized_autocovariance(acov_mean2(3,:), acov_mean(3,1,:), acov_mean(3,2,:), &
						guess_Rwindow(2,:), real(xwindow, dp), nacov, tdisc)
				!$omp end sections
				!$omp end parallel
			end if
		end do
	end if
	! Update state of system in preparation for next step.
	! Update probability
	pcond(x(1)+1, x(2)+1) = pcond(x(1)+1, x(2)+1) + tstep
	! Update counter
	event_count(event) = event_count(event) + 1
	! Update abundances according to what happened
	xlast = x
	x = x + abund_update(:, event)
end do

! Normalize probabilities.
write(*,*) "Recorded time:", t
write(*,*) "Summed probability:", sum(pcond)
pcond = pcond / sum(pcond)

! Assign means
mean_x(1) = stats(1,1)%mean_x
mean_x(2) = stats(2,2)%mean_x
mean_r = event_count / t
write(*,*) mean_x

do i = 1,2
do j = 1,2
	cov(i,j) = stats(i,j)%get_cov()
end do
end do

! Finalize autocovariances
do i = 1, 3
	call finalize_autocovariance(acov_mean2(i,:), acov_mean(i,1,:), acov_mean(i,2,:), tlast-tstart, acov(i,:))
end do

write(*,*) covxx%get_cov()
write(*,*) acov(1,1)

call check_moments(cov, mean_x, mean_r)

! Output results
call dump_data()


contains


pure function guess_R(x, lmbda, k, n, c, beta) result(rate)
	integer, intent(in) :: x(2)
	real(dp), intent(in) :: lmbda, beta, k, n, c
	real(dp) :: rate(2)
	
	! Rate in
	rate(1) = 1._dp * lmbda * hill(x(1), k, n, c)
	! Rate out
	rate(2) = 1._dp * beta * x(2)
end function


pure function production_rates(x, lmbda, k, n, c) result(rate)
	integer, intent(in) :: x(2)
	real(dp), intent(in) :: lmbda(2), k(2), n(2), c(2)
	real(dp) :: rate(2)
	
	! Write it like this so we kind of flip x1 and x2.
	! R1 is proportional to x2. R2 is R(U, x1, x2)
	rate(1) = 1._dp * lmbda(1) * x(2)
	rate(2) = 1._dp * lmbda(2) * hill(x(1), k(2), n(2), c(2))
	! rate(2) = 1._dp * lmbda(2)
end function


pure function decay_rates(x, beta) result(rate)
	integer, intent(in), dimension(2) :: x
	real(dp), intent(in), dimension(2) :: beta
	real(dp) :: rate(2)
	
	rate(1) = 1._dp * x(1) * beta(1)
	rate(2) = 1._dp * x(2) * beta(2)
end function


subroutine update_propensity(prop, x, lmbda, k, n, c, beta)
	integer, intent(in), dimension(2) :: x
	real(dp), intent(in), dimension(2) :: lmbda, k, n, c, beta
	real(dp), dimension(4), intent(out) :: prop
	real(dp), dimension(2) :: Rin, Rout
	
	Rin = production_rates(x, lmbda, k, n, c)
	Rout = decay_rates(x, beta)
	
	! Make and degrade x1.
	prop(1) = Rin(1);	prop(2) = Rout(1)
	! Make and degrade x2
	prop(3) = Rin(2);	prop(4) = Rout(2)
end subroutine


subroutine check_moments(cov, mean_x, mean_r)
	character(len=*), parameter :: change_method = "logarithmic"
	real(dp), intent(in) :: cov(2,2), mean_x(2), mean_r(4)
	real(dp) :: tau(2), sumA, sumB, sumAB, LHS, RHS, rel_change
	integer :: i, j
		
	! Calculate lifetimes using Little's Law
	do i = 1, 2
		tau(i) = mean_x(i) / mean_r(2*i)
	end do
		
	write(*,*) "Covariance balance checks: "
	write(*,*) "Relative change between LHS and RHS, or term1 and -term2 if <S>'s are zero: "
	do i = 1, 2
	do j = 1, 2
		if (i .eq. j) then
			LHS = 2._dp * cov(i,j)
			RHS = 2._dp * mean_r(2*i)
		else
			LHS = 1._dp * ( 1/tau(j) * cov(i,j) / mean_r(2*j) / mean_x(i))
			RHS = -1._dp * ( 1/tau(i) * cov(j,i) / mean_r(2*i) / mean_x(j))
		end if		
		rel_change = relative_change(LHS, RHS, change_method)
		write(*,*) i, j, rel_change
	end do
	end do
end subroutine


subroutine dump_params()
! Output system information to console
    integer :: values(8)
    character(19) :: formatted_time

	call date_and_time(VALUES=values)
	write(formatted_time, '(I2.2, ":", I2.2, ":", I2.2)') values(5), values(6), values(7)
	write(*,*) formatted_time
	write(*,*) "==== System parameters ===="
	write(*,*) "Minimum events: ", event_min
	write(*,*) "Lambda: ", lmbda
	write(*,*) "k: ", k
	write(*,*) "n: ", n
	write(*,*) "c: ", c
	write(*,*) "Beta: ", beta
end subroutine


subroutine dump_data()
! Output data to file
	real(dp) t, x
	integer :: i, j, io
	character(len=256) :: fname, prefix
	character(len=*), parameter :: suffix = ".dat"
	
	write(prefix, '(A, A, F0.3, A, F0.3, A)') filename, "_var1_", guess_params(1), "_var2_", guess_params(2), "_"

	call generate_ISO_filename(directory, prefix, suffix, fname)
	
	write(*,*) "File output at: ", fname

	open(newunit=io, file=fname, status='new', action='write')
	! Write metadata
	write(io,*) "# Program Metadata"
	write(io,*) "# Program: cascade_2D_discrete.f90"
	write(io,*) "# Creation date: ", fdate()
	write(io,*) "# Seed: ", rseed
	write(io,*) "# Min events: ", event_min
	write(io,*) "# Maximum abundance: ", abund_max
	write(io,*) "# Time discretization: ", tdisc
	write(io,*) "# Autocovariance max lag", umax
	write(io,*) ""
	write(io,*) ""
	write(io,*) "# Parameter Metadata"
	write(io,*) "# k_guess: ", guess_params(1)
	write(io,*) "# n_guess: ", guess_params(2)
	write(io,*) "# lmbda1: ", lmbda(1)
	write(io,*) "# lmbda2: ", lmbda(2)
	write(io,*) "# k1: ", k(1)
	write(io,*) "# k2: ", k(2)
	write(io,*) "# n1: ", n(1)
	write(io,*) "# n2: ", n(2)
	write(io,*) "# c1: ", c(1)
	write(io,*) "# c2: ", c(2)
	write(io,*) "# beta1: ", beta(1)
	write(io,*) "# beta2: ", beta(2)
	write(io,*) "# x1avg: ", mean_x(1)
	write(io,*) "# x2avg: ", mean_x(2)
	write(io,*) "# Rx1+avg: ", mean_r(1)
	write(io,*) "# Rx1-avg: ", mean_r(2)
	write(io,*) "# Rx2+avg: ", mean_r(3)
	write(io,*) "# Rx2-avg: ", mean_r(4)
	write(io,*) "# covxx: ", covxx%get_cov()
	write(io,*) ""
	write(io,*) ""
	write(io,*) "# Data"
	write(io,*) "# Time, A22, guess_Ap2, guess_Ad2"
	
	! Write simulation data
	do i = 1, nacov
		t = (i-1)*tdisc
		write(io,*) t, acov(:, i)
	end do

	write(io,*) ""
	write(io,*) ""
	write(io,*) "# Probability Matrix"
	do i = 1, abund_max(1)
	do j = 1, abund_max(2)
		if (pcond(i,j) > eps) then
			write(io, *) i-1, j-1, pcond(i,j)
		end if
	end do
	end do
	
	close(io)
end subroutine

end program
