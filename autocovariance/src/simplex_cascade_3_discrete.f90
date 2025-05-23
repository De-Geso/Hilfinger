program simplex_cascade_3_discrete
use nr_minmax
use kind_parameters
use stochastics
use randf
use utilities
implicit none

! Program hyper parameters =============================================
! Number of events of each reaction before stopping.
integer, parameter :: event_min = 10**6
! Maximum abundance in each direction
! integer, parameter, dimension(3) :: abund_max = [2**9, 2**7, 2**7]
! Maximum time lag for autocovariance.
real(dp), parameter :: maxlag = 3._dp
! Discretization of time
real(dp), parameter :: tdisc = 0.1_dp
! Number of points in autocovariance.
integer, parameter :: nacov = ceiling(maxlag/tdisc) + 1

! System parameters ====================================================
! Production rates
real(dp), parameter, dimension(3) :: lmbda = [25._dp, 25._dp, 80._dp]
! Decay rates
real(dp), parameter, dimension(3) :: beta = [1._dp, 1._dp, 1._dp]
! Hill function parameters
real(dp), parameter, dimension(3) :: k = [0._dp, 50._dp, 40._dp]
real(dp), parameter, dimension(3) :: n = [0._dp, 4._dp, 2._dp]
! Constant offset
real(dp), parameter, dimension(3) :: c = [0._dp, 8._dp, 0._dp]

! Abundance update matrix.
integer, parameter, dimension(3,6) :: abund_update = &
	reshape((/1, 0, 0, &
			-1, 0, 0, &
			0, 1, 0, &
			0, -1, 0, &
			0, 0, 1, &
			0, 0, -1/), shape(abund_update))

! Variables ============================================================
! real(dp), dimension(nacov) :: t, Axx, Ainx, Aoutx, dAxx
real(dp), dimension(3, 2) :: xx = &
	reshape([15., 17., 14., &
			-2.9, -2.8, -2.7], shape(xx))
real(dp), dimension(3) :: yy
integer :: i, iter


! Here the program begins. =============================================


call random_seed()
! do i = 1, 3
!	call random_number(xx(i,1))
!	call random_number(xx(i,2))
!	xx(i,1) = 100._dp * xx(i,1)
!	xx(i,2) = -10._dp * xx(i,2)
! end do

! Initial guesses
do i = 1, 3
	yy(i) = acov_SSE(xx(i,:))
end do

write(*,*) yy

call amoeba(xx, yy, 3, 2, 2, 1.E-6_dp, acov_SSE, iter)

write(*,*) iter
write(*,*) "x, y, f(x)"
write(*,*) xx(1,:), yy(1)


contains


function acov_SSE(xx) result(SSE)
! Return the SSE of the autocovariance for a given set of guessed
! parameters xx. This gets fed to the amoeba
	real(dp), intent(in) :: xx(:)
	real(dp), dimension(nacov) :: t, Axx, Ainx, Aoutx, dAxx
	real(dp), parameter :: large=huge(1._dp)
	real(dp) :: SSE
	integer :: i

	! If we have negative parameters, don't bother running the simulation
	! just return a large number.
	if (xx(1) <= 0._dp) then
		SSE = large
	else
		call simulate_system(Axx, Ainx, Aoutx, xx)
		
		! Initialize time
		do i = 1, nacov
			t(i) = (i-1)*tdisc
		end do
		
		dAxx = gradient(Axx, tdisc, nacov)
		SSE = sum((dAxx-(Ainx-Aoutx))**2)
	end if
end function


subroutine simulate_system(Axx, Ainx, Aoutx, guess_params)
! Simulate the system and calculate the autocovariances for a given set
! of guessed parameters. Will be the workhorse when doing the search.
	real(dp), intent(in) :: guess_params(2)
	real(dp), intent(out), dimension(nacov) :: Axx, Ainx, Aoutx

	! Autocovariance variables =============================================
	! Windows
	real(dp) :: guess_Rwindow(2, nacov)=0._dp
	integer(dp) :: xwindow(nacov)=0._dp
	! Covariance arrays
	real(dp) :: acov(3,nacov)=0._dp, acov_mean(3,2,nacov)=0._dp, acov_mean2(3,nacov)=0._dp

	! Variables ============================================================
	! Timers
	real(dp) :: t=0._dp, tstep, tlast=0._dp
	! Probability distribution
!	real(dp), dimension(abund_max(1), abund_max(2), abund_max(3)) :: pcond=0._dp
	! Moments
	real(dp) :: meanX(3), meanR(6), eta(3,3)
	! Gillespie
	real(dp) :: propensity(6), roll
	integer :: x(3)=0, xlast(3)=0, maxX(3)=0, event_count(6)=0, event
	! Other
	integer :: i, nseed, tcount=0
	integer, allocatable :: rseed(:)
	
	guess_Rwindow=0._dp
	xwindow=0._dp
	acov=0._dp; acov_mean=0._dp; acov_mean2=0._dp
	t=0._dp; tlast=0._dp
!	pcond=0._dp
	x=0; xlast=0; maxX=0; event_count=0
	tcount=0
	
	! Here the program begins ==============================================
	
	! call random_seed(put=seed)
	call random_seed()
	! Get random seed for output in metadata
	call random_seed(size=nseed)
	allocate(rseed(nseed))
	call random_seed(get=rseed)
		
	xlast = x
	
	do while (minval(event_count) < event_min)
!		do i = 1, 3
!			if (x(i) .eq. abund_max(i)) then
!				write(*,"(A, I0, A, I0)") "Element ", i, " exceeded maximum abundance ", abund_max(i)
!				call exit()
!			end if
!		end do
	
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
				guess_Rwindow(:, nacov) = guess_R(x, guess_params(1), guess_params(2))
				
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
!		pcond(x(1)+1, x(2)+1, x(3)+1) = pcond(x(1)+1, x(2)+1, x(3)+1) + tstep
		! Update counter
		event_count(event) = event_count(event) + 1
		! Update abundances according to what happened
		xlast = x
		x = x + abund_update(:, event)
	end do
		
	! Normalize probabilities.
!	pcond = pcond / sum(pcond)
	
	do i = 1, 3
		call finalize_autocovariance(acov_mean2(i,:), acov_mean(i,1,:), acov_mean(i,2,:), tlast, acov(i,:))
	end do
	
	Axx = acov(1,:)
	Ainx = acov(2,:)
	Aoutx = acov(3,:)
end subroutine


pure function guess_R(x, var1, var2) result(rate)
	integer, intent(in) :: x(3)
	real(dp), intent(in) :: var1, var2
	real(dp) :: rate(2)
	
	! Rate in
	rate(1) = 1._dp * lmbda(3) * hill(x(2), var1, var2, c(3))
	! Rate out
	rate(2) = 1._dp * beta(3) * x(3)
end function


pure function production_rates(x, lmbda, k, n, c) result(rate)
	integer, intent(in) :: x(3)
	real(dp), intent(in) :: lmbda(3), k(3), n(3), c(3)
	real(dp) :: rate(3)
	
	! Testing
!	rate(1) = 1._dp * lmbda(1)
!	rate(2) = 1._dp * lmbda(2) * x(1)
!	rate(3) = 1._dp * lmbda(3) * x(2)
	
	! PLOS R3+
	rate(3) = 1._dp * lmbda(3) * hill(x(2), k(3), n(3), c(3))
	
	! Bistable
!	rate(1) = 1._dp * lmbda(1) * hill(x(2), k(1), n(1), c(1))
!	rate(2) = 1._dp * x(1)
	
	! Oscillating
!	rate(1) = 1._dp * lmbda(1) * hill_neg(x(3), k(1), n(1))
!	rate(2) = 1._dp * lmbda(2) * hill_pos(x(1), k(2), n(2))

	! Noise Enhancing
	rate(1) = 5._dp
	rate(2) = 1._dp * lmbda(2) * hill(x(3), k(2), n(2), c(2)*x(1))
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


pure function gradient(f, h, n) result (fx)
	integer, intent(in) :: n
	real(dp), intent(in) :: f(n), h
	real(dp) :: fx(n)
	integer :: i
	
	! First point
	fx(1) = (-3._dp*f(1) + 4._dp*f(2) - f(3)) / (2._dp*h)
	! Last point
	fx(n) = (3._dp*f(n) - 4._dp*f(n-1) + f(n-2)) / (2._dp*h)
	
	do i = 2, size(fx)-1
		fx(i) = (f(i+1) - f(i-1)) / (2._dp*h)
	end do
end function


!pure function hill_pos(x, k, n) result(f)
!! Hill function. k and n controlled in system parameters.
!	integer, intent(in) :: x
!	real(dp), intent(in) :: k, n
!	real(dp) :: f
!	f = 1._dp * x**n / (x**n + k**n)	
!end function


!pure function hill_neg(x, k, n) result(f)
!! Hill function. k and n controlled in system parameters.
!	integer, intent(in) :: x
!	real(dp), intent(in) :: k, n
!	real(dp) :: f
!	f = 1._dp / (1. + (x/k)**n)
!end function


end program
