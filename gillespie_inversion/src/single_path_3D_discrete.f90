program single_path_3D

use kind_parameters
use stochastics
use randf
use utilities

implicit none

! Program hyper parameters =============================================
! Output filename
character(len=*), parameter :: filename = "3D_discrete_oscillating"
! Epsilon
real(dp), parameter :: eps = 1E-15
! Number of events of each reaction before stopping.
integer, parameter :: event_min = 10**6
! Maximum abundance in each direction
integer, parameter, dimension(3) :: abund_max = [2**9, 2**7, 2**7]
! Maximum time lag for autocovariance.
real(dp), parameter :: maxlag = 15._dp
! Discretization of time
real(dp), parameter :: tdisc = 0.1_dp
! Number of points in autocovariance.
integer, parameter :: nacov = ceiling(maxlag/tdisc) + 1


! System parameters ====================================================
! Production rates
real(dp), parameter, dimension(3) :: lmbda = [500._dp, 80._dp, 80._dp]
! Decay rates
real(dp), parameter, dimension(3) :: beta = [1._dp, 1._dp, 1._dp]
! Hill function parameters
real(dp), parameter, dimension(3) :: k = [0.1_dp, 100._dp, 40._dp]
real(dp), parameter, dimension(3) :: n = [10._dp, 1._dp, 2._dp]
! Constant offset
real(dp), parameter, dimension(3) :: c = [0._dp, 0._dp, 0._dp]
! Abundance update matrix.
integer, parameter, dimension(3,6) :: abund_update = &
	reshape((/1, 0, 0, &
			-1, 0, 0, &
			0, 1, 0, &
			0, -1, 0, &
			0, 0, 1, &
			0, 0, -1/), shape(abund_update))

! Autocovariance variables =============================================
! Windows
real(dp) :: guess_Rwindow(2, nacov) = 0._dp
integer(dp) :: xwindow(nacov) = 0
! Covariance arrays
real(dp) :: acov(3, nacov)=0._dp, acov_mean(3, 2, nacov)=0._dp, acov_mean2(3, nacov)=0._dp

! Variables ============================================================
! Timers
real(dp) :: t = 0._dp, tstep, tlast = 0._dp
! Probability distribution
real(dp), dimension(abund_max(1), abund_max(2), abund_max(3)) :: pcond = 0._dp
! Moments
real(dp) :: meanX(3), meanR(6), eta(3,3)
! Gillespie
real(dp) :: propensity(6), roll
integer :: x(3)=0, xlast(3), maxX(3) = 0, event_count(6) = 0, event
! Counters
integer :: i, nseed, tcount=0
integer, allocatable :: rseed(:)
! Guessed k and n
real(dp) :: guess_params(2) = [k(3), n(3)]


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
	do i = 1, 3
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

	if (t .ge. (tlast + tdisc)) then
		! Online correlation calculation
		! Track the abundances and time in a window for correlations.
		do i = 1, floor((t-tlast)/tdisc)
			tcount = tcount + 1
			tlast = 1._dp * tdisc * tcount

			xwindow(:nacov-1) = xwindow(2:nacov)
			xwindow(nacov) = x(3)

			! Track the guessed rates
			guess_Rwindow(:, :nacov-1) = guess_Rwindow(:, 2:nacov)
			guess_Rwindow(:, nacov) = guess_R(x, lmbda(3), guess_params(1), guess_params(2), beta(3))
			
			! Online calculation of autocovariance
			if (t .gt. maxlag) then
				! Record a special time to divide autocovariance means by.
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
	pcond(x(1)+1, x(2)+1, x(3)+1) = pcond(x(1)+1, x(2)+1, x(3)+1) + tstep
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

! Get moments
call get_moments(pcond, event_count, t, meanX, meanR, eta)
! Check theory
call check_moments(pcond, event_count, t, meanX, meanR)

do i = 1, 3
	call finalize_autocovariance(acov_mean2(i,:), acov_mean(i,1,:), acov_mean(i,2,:), tlast, acov(i,:))
end do

! Output results
call dump_data()


contains


pure function guess_R(x, lmbda, k, n, beta) result(rate)
	integer, intent(in) :: x(3)
	real(dp), intent(in) :: lmbda, beta, k, n
	real(dp) :: rate(2)
	
	! Rate in
	rate(1) = 1._dp * lmbda * hill_pos(x(2), k, n)
	! Rate out
	rate(2) = 1._dp * beta * x(3)
end function


pure function production_rates(x, lmbda, k, n, c) result(rate)
	integer, intent(in) :: x(3)
	real(dp), intent(in) :: lmbda(3), k(3), n(3), c(3)
	real(dp) :: rate(3)
	
	! PLOS R3+
	rate(3) = 1._dp * lmbda(3) * hill_pos(x(2), k(3), n(3))

	! Simple testing
	! rate(1) = 1._dp * lmbda(1)
	! rate(2) = 1._dp * lmbda(2) * x(1)
	! rate(3) = 1._dp * lmbda(3) * x(2)
	
	! Oscillating
	rate(1) = 1._dp * lmbda(1) * hill_neg(x(3), k(1), n(1))
	rate(2) = 1._dp * lmbda(2) * hill_pos(x(1), k(2), n(2))
end function


pure function decay_rates(x, beta) result(rate)
	integer, intent(in), dimension(3) :: x
	real(dp), intent(in), dimension(3) :: beta
	real(dp) :: rate(3)
	
	rate(1) = 1._dp * x(1) * beta(1)
	rate(2) = 1._dp * x(2) * beta(2)
	rate(3) = 1._dp * x(3) * beta(3)
end function


subroutine update_propensity(prop, x, lmbda, k, n, c, beta)
	integer, intent(in), dimension(3) :: x
	real(dp), intent(in), dimension(3) :: lmbda, k, n, c, beta
	real(dp), dimension(6), intent(inout) :: prop
	real(dp), dimension(3) :: Rin, Rout
	
	Rin = production_rates(x, lmbda, k, n, c)
	Rout = decay_rates(x, beta)
	
	! Make and degrade x1.
	prop(1) = Rin(1);	prop(2) = Rout(1)
	! Make and degrade x2
	prop(3) = Rin(2);	prop(4) = Rout(2)
	! Make and degrade x3
	prop(5) = Rin(3);	prop(6) = Rout(3)
end subroutine


pure function hill_pos(x, k, n) result(f)
! Hill function. k and n controlled in system parameters.
	integer, intent(in) :: x
	real(dp), intent(in) :: k, n
	real(dp) :: f
	f = 1._dp * x**n / (x**n + k**n)	
end function


pure function hill_neg(x, k, n) result(f)
! Hill function. k and n controlled in system parameters.
	integer, intent(in) :: x
	real(dp), intent(in) :: k, n
	real(dp) :: f
	f = 1._dp / (1. + (x/k)**n)
end function

! Calculate mean abundances, mean rates, and etas from probability distribution.
subroutine get_moments(pcond, event_count, t, meanX, meanR, eta)
	real(dp), intent(in) :: pcond(:,:,:), t
	integer, intent(in), dimension(6) :: event_count
	real(dp), intent(out) :: meanX(3), meanR(6), eta(3,3)
	
	real(dp), dimension(maxval(shape(pcond))) :: x
	real(dp), dimension(maxval(shape(pcond)), maxval(shape(pcond))) :: p2d
	real(dp), dimension(maxval(shape(pcond)), size(meanX)) :: prob
	integer :: i, j, k, l
	
	! Calculate mean rates
	! Abuse knowledge of event numbers and time to get <R> easily
	meanR = event_count/t
	write(*,*) "Simulation mean rates (in/out): "
	! Rates in
	write(*,*) meanR(1), meanR(3), meanR(5)
	! Rates out
	write(*,*) meanR(2), meanR(4), meanR(6)
	
	
	! Calculate mean abundances
	! Make the abundance vectors
	do i = 1, size(x)
		x(i) = i-1
	end do
	
	! Make marginal probability distributions
	prob(:,1) = sum(sum(pcond, dim=3), dim=2)
	prob(:,2) = sum(sum(pcond, dim=3), dim=1)
	prob(:,3) = sum(sum(pcond, dim=1), dim=1)
	
	! Calculate means
	meanX = matmul(x, prob)
	write(*,*) "Simulation mean abundances: ", meanX

	! Calculate etas
	eta = 0
	do i = 1, 3
		do j = 1, i
			if (i .eq. j) then
				eta(i,i) = dot_product(prob(:,i), (x-meanX(i))**2) / meanX(i)**2
			else
				p2d = 0
				p2d = sum(pcond, dim=3-mod((i+j),3))
				eta(i,j) = prob2covariance(p2d, x, x) / meanX(i) / meanX(j)
				eta(j,i) = eta(i,j)
			end if
		end do
	end do	
	write(*,*) "Simulation etas:"
	do i = 1, 3
		write(*,*) eta(i,:)
	end do
end subroutine


! AFTER calculating the means, check flux balance and covariance balance
subroutine check_moments(pcond, event_count, t, meanX, meanR)
	real(dp), intent(in) :: pcond(:,:,:), t, meanX(3), meanR(6)
	integer, intent(in) :: event_count(6)
	real(dp) :: tau(3), cov(3,3) = 0._dp, Rin(3), Rout(3), sumA, sumB, sumAB, LHS, RHS, rel_change
	integer :: vec(3), i, j, x, y, z
	character(len=*), parameter :: change_method = "logarithmic"
	
	! Check flux balance
	write(*,*) "Flux balance checks: "
	write(*,*) "Relative change between Ri+, Ri-: "
	do i = 1, 3
		rel_change = relative_change(meanR(2*i-1), meanR(2*i), change_method)
		write(*,*) i, rel_change
	end do
	
	
	! Calculate lifetimes using Little's Law
	do i = 1, 3
		tau(i) = meanX(i) / meanR(2*i)
	end do
	
	do i = 1, 3
	do j = 1, 3
		sumA = 0._dp
		sumB = 0._dp
		sumAB = 0._dp
		do x = 1, size(pcond, 1)
		do y = 1, size(pcond, 2)
		do z = 1, size(pcond, 3)
			vec = [x-1, y-1, z-1]
			Rin = production_rates(vec, lmbda, k, n, c)
			Rout = decay_rates(vec, beta)
			sumA = sumA + pcond(x,y,z) * (Rout(i)-Rin(i))
			sumB = sumB + pcond(x,y,z) * vec(j)
			sumAB = sumAB + pcond(x,y,z) * (Rout(i)-Rin(i))*vec(j)
		end do
		end do
		end do
		cov(i,j) = sumAB - sumA*sumB
	end do
	end do
	
	write(*,*) "Covariance balance checks: "
	write(*,*) "Relative change between LHS and RHS, or term1 and -term2 if <S>'s are zero: "
	do i = 1, 3
	do j = 1, 3
		if (i .eq. j) then
			LHS = 1._dp * ( 1/tau(j) * cov(j,i)/meanR(2*j)/meanX(i) + (1/tau(i) * cov(i,j)/meanR(2*i)/meanX(j)) )
			RHS = 1._dp * ( (1/tau(i) * 1/meanX(j)) + 1/tau(i) * 1/meanX(j) )
		else
			LHS = 1._dp * ( 1/tau(j) * cov(j,i)/meanR(2*j)/meanX(i) )
			RHS = -1._dp * ( 1/tau(i) * cov(i,j)/meanR(2*i)/meanX(j) )
		end if
		
		rel_change = relative_change(LHS, RHS, change_method)
		write(*,*) i, j, rel_change
	end do
	end do
	
end subroutine


! Calculate covariance given 2D conditional probability distribution and values
function prob2covariance(prob, x, y) result(cov)
	real(dp), intent(in) :: prob(:,:), x(:), y(:)
	real(dp) :: cov, xavg, yavg
	integer :: nx, ny, i, j
	
	nx = size(x)
	ny = size(y)

	xavg = dot_product(sum(prob, dim=1) , x)
	yavg = dot_product(sum(prob, dim=2) , y)
	
	cov = 0
	do i = 1, nx
		do j = 1, ny
			cov = cov + (prob(i,j) * (x(i)-xavg)*(x(j)-xavg))
		end do
	end do
end function


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
	integer :: i, io
	character(len=256) :: fname, prefix
	character(len=*), parameter :: directory = "data/", suffix = ".dat"
	
	write(prefix, '(A, A, F0.3, A, F0.3, A)') filename, "_kguess_", guess_params(1), "_nguess_", guess_params(2), "_"

	call generate_ISO_filename(directory, prefix, suffix, fname)
	
	write(*,*) "File output at: ", fname

	open(newunit=io, file=fname, status='new', action='write')
	! Write metadata
	write(io,*) "# Program Metadata"
	write(io,*) "# Program: single_path_3D_discrete.f90"
	write(io,*) "# Creation date: ", fdate()
	write(io,*) "# Seed: ", rseed
	write(io,*) "# Min events: ", event_min
	write(io,*) "# Maximum abundance: ", abund_max
	write(io,*) "# Time discretization: ", tdisc
	write(io,*) "# Autocovariance max lag", maxlag
	write(io,*) ""
	write(io,*) "# Parameter Metadata"
	write(io,*) "# k_guess: ", guess_params(1)
	write(io,*) "# n_guess: ", guess_params(2)
	write(io,*) "# lmbda1: ", lmbda(1)
	write(io,*) "# lmbda2: ", lmbda(2)
	write(io,*) "# lmbda3: ", lmbda(3)
	write(io,*) "# k1: ", k(1)
	write(io,*) "# k2: ", k(2)
	write(io,*) "# k3: ", k(3)
	write(io,*) "# n1: ", n(1)
	write(io,*) "# n2: ", n(2)
	write(io,*) "# n3: ", n(3)
	write(io,*) "# c1: ", c(1)
	write(io,*) "# c2: ", c(2)
	write(io,*) "# c3: ", c(3)
	write(io,*) "# beta1: ", beta(1)
	write(io,*) "# beta2: ", beta(2)
	write(io,*) "# beta3: ", beta(3)
	write(io,*) "# x1avg: ", meanX(1)
	write(io,*) "# x2avg: ", meanX(2)
	write(io,*) "# x3avg: ", meanX(3)
	write(io,*) "# Rx1+avg: ", meanR(1)
	write(io,*) "# Rx1-avg: ", meanR(2)
	write(io,*) "# Rx2+avg: ", meanR(3)
	write(io,*) "# Rx2-avg: ", meanR(4)
	write(io,*) "# Rx3+avg: ", meanR(5)
	write(io,*) "# Rx3-avg: ", meanR(6)
	write(io,*) ""
	write(io,*) "# Data"
	write(io,*) "# Time, A33, guess_Ap3, guess_Ad3"
	
	! Write simulation data
	do i = 1, nacov
		t = (i-1)*tdisc
		write(io,*) t, acov(:, i)
	end do
	
	close(io)
end subroutine

end program
