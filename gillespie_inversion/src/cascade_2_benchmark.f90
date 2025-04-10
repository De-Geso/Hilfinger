program gillespie
! Simulate the system given in Hilfinger 2016 Eq. 10 using the Gillespie
! algorithm. Imagine a feedback control system. mRNA is produced, which
! triggers the production of a Protein. Both decay at some rate which
! depends on their abundance.
! x0 == R >> x0 + a			x1 == lmbda x0 >> x1 + b
! x0 == x0/tau1 >> x0 - 1	x1 == x1/tau1 >> x1 - 1
use kind_parameters
use stochastics
use randf
implicit none

! == Program hyperparameters ==
logical, parameter :: output = .true.
integer, parameter :: event_min = 10**7

! == System parameters ==
real(dp), parameter, dimension(2) :: lmbda=[7.2_dp, 10._dp], beta=[1._dp, 0.1_dp]
integer, parameter, dimension(2,4) :: abund_update = &
	reshape((/1, 0, &
			-1, 0, &
			0, 1, &
			0, -1/), shape(abund_update))
			
! Variables ============================================================
! Timers
real(dp) :: t = 0._dp, tstep
! Gillespie
real(dp) :: propensity(4), roll
integer :: x(2)=0, maxX(2)=0, event_count(4)=0, event
! Counters
integer :: i, io

integer, allocatable :: rseed(:)

! Here the program begins ==============================================
	
call random_seed(put=seed)		
	
if (output) then
	! Write to file in binary
	open(unit=io, file='binary_data.dat', form='unformatted', access='stream', status='replace')
end if

do while (minval(event_count) < event_min)
	if (output) write(io) t, x
		
	! Update the propensity before taking a Gillespie step
	call update_propensity(propensity, x, lmbda, beta)
	! Get timestep and event from Gillespie algorithm
	call gillespie_iter(tstep, event, propensity)
	! Update time by adding how long we were in the previous state
	t = t + tstep
	
	! Update state of system in preparation for next step.
	! Update counter
	event_count(event) = event_count(event) + 1
	! Update abundances according to what happened
	x = x + abund_update(:, event)
end do

if (output) write(io) t, x


contains


pure function production_rates(x, lmbda) result(rate)
	integer, intent(in) :: x(2)
	real(dp), intent(in), dimension(2) :: lmbda
	real(dp) :: rate(2)
	
	! Testing
	rate(1) = 1._dp * lmbda(1)
	rate(2) = 1._dp * lmbda(2) * x(1)
end function


pure function decay_rates(x, beta) result(rate)
	integer, intent(in), dimension(2) :: x
	real(dp), intent(in), dimension(2) :: beta
	real(dp) :: rate(2)
	
	rate(1) = 1._dp * beta(1) * x(1)
	rate(2) = 1._dp * beta(2) * x(2)
end function


subroutine update_propensity(prop, x, lmbda, beta)
	integer, intent(in), dimension(2) :: x
	real(dp), intent(in), dimension(2) :: lmbda, beta
	real(dp), dimension(4), intent(inout) :: prop
	real(dp), dimension(2) :: Rin, Rout
	
	Rin = production_rates(x, lmbda)
	Rout = decay_rates(x, beta)
	
	! Make and degrade x1.
	prop(1) = Rin(1);	prop(2) = Rout(1)
	! Make and degrade x2
	prop(3) = Rin(2);	prop(4) = Rout(2)
end subroutine
	
end program
