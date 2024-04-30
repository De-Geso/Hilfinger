program mrna_protein_feedback
use kind_parameters
use randf
use mrna_protein_system_parameters
use utilities
implicit none

! Program Hyperparameters ==============================================
! Number of events before stopping
integer, parameter :: event_min = 10**6
! Maximum abundance. Program will exit if exceeded.
integer, parameter :: abund_max = 2**4

! Variables ============================================================
! Time
real(dp) :: t=0._dp, tstep
! Probabilities
real(dp) :: pcond(abund_max, abund_max)=0._dp
! Moments
real(dp) :: mean(2), mean_rate, cov(2,2), thry_mean(2), thry_cov(2,2)
! Miscellaneous 
real(dp) :: propensity(4)
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
	
	! Update abundances LAST
	mp = mp + abund_update(:,event)
	nevents(event) = nevents(event) + 1
end do

! Normalize probabilities.
pcond = pcond / sum(pcond)

call simulation_moments(pcond, mean, mean_rate, cov)
call theory_moments(pcond, thry_mean, thry_cov)
call dump()



contains



subroutine simulation_moments(pcond, mean, mean_rate, cov)
! Calculate mean abundances and covariance matrix (eta--normalized)
! from conditional probability distribution.
	real(dp), intent(out) :: mean(2), mean_rate, cov(2,2)
	real(dp), intent(in) :: pcond(abund_max, abund_max)
	real(dp) :: p(abund_max, 2)
	integer :: i, j, x(abund_max)
	
	! First moments
	! Abundance means
	do i = 1, abund_max
		x(i) = i-1
		p(i,1) = sum(pcond(i,:))
		p(i,2) = sum(pcond(:,i))
	end do
	mean = matmul(x,p)
	
	! Rate mean
	mean_rate = 0._dp
	do i = 1, abund_max
	do j = 1, abund_max
		mean_rate = mean_rate + pcond(i,j)*R([x(i),x(j)])
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
	real(dp) :: mean_rate, p(abund_max, 2), s(2)
	integer :: i, j, x(abund_max)
	
	do i = 1, abund_max
		x(i) = i-1
		p(i,1) = sum(pcond(i,:))
		p(i,2) = sum(pcond(:,i))
	end do
	
	! First moment
	! Mean rate
	mean_rate = 0._dp
	do i = 1, abund_max
	do j = 1, abund_max
		mean_rate = mean_rate + pcond(i,j)*R([x(i),x(j)])
	end do
	end do
	
	! Mean abundances
	mean(1) = 1._dp * tau(1)*mean_rate
	mean(2) = 1._dp * beta*tau(2)*mean(1)
	
	! Second moments
	cov = 0._dp
	s(1) = (burst(1)+1) / 2
	s(2) = (burst(2)+1) / 2
	
	do i = 1, abund_max
	do j = 1, abund_max
		cov(1,1) = cov(1,1) + (pcond(i,j) * (x(i) - mean(1)) * (R([x(i), x(j)]) - mean_rate))
		cov(1,2) = cov(1,2) + (pcond(i,j) * (x(j) - mean(2)) * (R([x(i), x(j)]) - mean_rate))
	end do
	end do
	cov(1,1) = cov(1,1)/(mean(1)*mean_rate) + s(1)/mean(1)
	cov(1,2) = tau(1)/sum(tau) * cov(1,1) + tau(2)/sum(tau) * cov(1,2)/(mean(2)*mean_rate)
	cov(2,1) = cov(1,2)
	cov(2,2) = s(2)/mean(2) + cov(1,2)
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
	integer, intent(in) :: x(2)
	real(dp) :: prop(4)
	prop(1) = 1._dp * R(x)	! Make x1 (mRNA)
	prop(2) = 1._dp * x(1)/tau(1)	! Degrade x1 (mRNA)
	prop(3) = 1._dp * beta*x(1)	! Make x2 (Protein)
	prop(4) = 1._dp * x(2)/tau(2)	! Degrade x2 (Protein)
end function


pure function R(x) result(f)
! Production function for mRNA.
	integer, intent(in) :: x(2)
	real(dp) :: f
	associate(m => x(1), p => x(2))
	
	f = alpha * 1./(p+1)
	f = 1._dp * alpha
	
	end associate
end function


subroutine dump()
! Output results to console or file.

! Console output =======================================================
	write(*,*) "Events: ", event_min
	write(*,*) "alpha=", alpha, "beta=", beta, "Tau=", tau
	write(*,*) "Mean rate: ", mean_rate
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
end subroutine

end program
