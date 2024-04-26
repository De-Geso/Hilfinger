program mrna_protein_feedback
use kind_parameters
use mrna_protein_system_parameters
use randf
implicit none

! Program Hyperparameters ==============================================
! Number of events before stopping
integer, parameter :: event_min = 10**6
! Maximum abundance. Program will exit if exceeded.
integer, parameter :: abund_max = 2**6

! Variables ============================================================
! Time
real(dp) :: t=0._dp, tstep
! Probabilities
real(dp) :: pcond(abund_max, abund_max)=0._dp, prate(abund_max)=0._dp, &
	p(abund_max,2)=0._dp
! Miscellaneous 
real(dp) :: propensity(4), mean(2)
integer :: i, mp(2)=0, nevents(4)=0, event

call random_seed(put=seed)
! call random_seed()

do while (minval(nevents) < event_min)
	! Exit the program if we exceed maximum abundance.
	if (maxval(mp) >= abund_max) then
		write(*,*) "Maximum abundance exceeded."
		write(*,*) "alpha=", alpha, "beta=", beta, "Tau=", tau
		call exit()
	end if
	
	! Get timestep and event from Gillespie algorithm. Update propensity
	call gillespie_iter(mp, tstep, event, propensity)
	! Update time first to record the end times for an abundance, rather than the start times.
	t = t + tstep

	! Add time step to probability matrices
	pcond(mp(1)+1, mp(2)+1) = pcond(mp(1)+1, mp(2)+1) + tstep
	prate(mp(2)+1) = prate(mp(2)+1) + tstep
	p(mp(1)+1,1) = p(mp(1)+1,1) + tstep
	p(mp(2)+1,2) = p(mp(2)+1,2) + tstep

	! Update abundances LAST
	mp = mp + abund_update(:,event)
	nevents(event) = nevents(event) + 1
end do

! Normalize probabilities.
pcond = pcond / sum(pcond)
prate = prate / sum(prate)
p = p / (sum(p)/2)

mean = simulation_mean(pcond)
write(*,*) mean

contains


pure function simulation_mean(pcond) result(mean)
! Calculate mean abundances from joint(conditional) probability dist.
	real(dp) :: mean(2)
	real(dp), intent(in) :: pcond(abund_max, abund_max)
	real(dp) :: p(abund_max, 2)
	integer :: x(abund_max), i
	
	do i = 1, abund_max
		x(i) = i-1
		p(i,1) = sum(pcond(i,:))
		p(i,2) = sum(pcond(:,i))
	end do
	mean = matmul(x,p)
end function


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
	integer, intent(in) :: x(2)
	real(dp) :: prop(4)
	prop(1) = 1._dp * R(x)	! Make x1 (mRNA)
	prop(2) = 1._dp * x(1)/tau(1)	! Degrade x1 (mRNA)
	prop(3) = 1._dp * beta*x(1)	! Make x2 (Protein)
	prop(4) = 1._dp * x(2)/tau(2)	! Degrade x2 (Protein)
end function


pure function R(x) result(f)
	integer, intent(in) :: x(2)
	real(dp) :: f
	f = alpha
end function

end program
