program inversion_cascade_1D
use kind_parameters
use stochastics
use randf
use utilities
implicit none

character(len=*), parameter :: fpath = "data/"
character(len=*), parameter :: fname = "dimer_1D"

real(dp), parameter :: eps=tiny(eps)

integer, parameter :: event_min = 10**6
integer, parameter :: abund_max = 128

! System parameters ====================================================
real(dp), parameter :: lmbda = 1._dp
real(dp), parameter :: beta = 1._dp
real(dp), parameter :: k = 10._dp
real(dp), parameter :: n = 2._dp
real(dp), parameter :: c = 0._dp
integer, parameter, dimension(2) :: burst = [1, -1]
integer, parameter, dimension(2) :: abund_update = [burst(1), burst(2)]

real(dp) :: t=0._dp, tstep
real(dp) :: pcond(abund_max)=0._dp, r_infer(abund_max,2)
type(OnlineCovariance) :: cov_balance, x_cov
integer :: visits_exits(abund_max,3)=0
real(dp) :: x_mean, r_mean(2)
! Gillespie
real(dp) :: propensity(2), roll
integer :: x=0, x_last, x_max=0, event_count(2)=0, event
! Other
integer :: i, nseed
integer, allocatable :: rseed(:)


! Here the program begins ==============================================


call random_seed(put=seed)
! Get random seed for output in metadata
call random_seed(size=nseed)
allocate(rseed(nseed))
call random_seed(get=rseed)
		
do while (minval(event_count) < event_min)
	if (x .eq. abund_max) then
		write(*,"(A, I0)") "Element 1 exceeded maximum abundance ", abund_max
		call exit()
	end if
	
	if (x .gt. x_max) x_max = x
	
	! Update the propensity before taking a Gillespie step
	call update_propensity(propensity, x, lmbda, beta)
	! Get timestep and event from Gillespie algorithm
	call gillespie_iter(tstep, event, propensity)
	
	! Online calculation of mean and covariance
	call cov_balance%update(real(x,dp), abs(burst(2))*propensity(2) - abs(burst(1))*propensity(1), tstep)
	call x_cov%update(real(x,dp), real(x,dp), tstep)
		
	visits_exits(x+1, 1) = visits_exits(x+1, 1) + 1
	visits_exits(x+1, event+1) = visits_exits(x+1, event+1) + 1

	! Update state of system in preparation for next step.
	t = t + tstep
	pcond(x+1) = pcond(x+1) + tstep
	x = x + abund_update(event)
	event_count(event) = event_count(event) + 1
end do

pcond = pcond / sum(pcond)

r_infer(:,1) = visits_exits(:,2) / (pcond(:) * t + eps)
r_infer(:,2) = visits_exits(:,3) / (pcond(:) * t + eps)

x_mean = x_cov%mean_x
r_mean = event_count/t

call dump()
write(*,*) "x mean:", x_mean
write(*,*) "x covariance:", x_cov%get_cov()
write(*,*) "r mean:", r_mean
call check_covariance_balance(cov_balance%get_cov(), x_mean, r_mean, burst, "arithmetic")


contains


subroutine check_covariance_balance(cov, x_avg, r_avg, burst, change_method)
	character(len=*), intent(in) :: change_method
	real(dp), intent(in) :: cov, x_avg, r_avg(2)
	integer, intent(in) :: burst(2)
	real(dp) :: rel_change, U, D, tau, s, flux_avg(2)
	
	flux_avg = abs(r_avg * burst)
	
	tau = 1._dp * x_avg / flux_avg(2)
	s = 0.5_dp * (burst(1) - burst(2))
		
	U = 1._dp/tau * cov / x_avg / (sum(flux_avg)/size(flux_avg))
	D = 2._dp/tau * s / x_avg
	rel_change = relative_change(2._dp*U, D, change_method)
	
	write(*,*) change_method, " relative change between 2U and D:", rel_change
	write(*,*) "2U:", 2._dp*U, "D:", D
end subroutine

pure function Rin(x) result(f)
	integer, intent(in) :: x
	real(dp) :: f
	
	! Linear system
	f = 1._dp * lmbda
end function


pure function Rout(x) result(f)
	integer, intent(in) :: x
	real(dp) :: f
	
	f = 1._dp * x * beta
end function


subroutine update_propensity(prop, x, lmbda, beta)
	integer, intent(in) :: x
	real(dp), intent(in) :: lmbda, beta
	real(dp), dimension(2), intent(out) :: prop
	
	prop(1) = Rin(x)
	prop(2) = Rout(x)
end subroutine


subroutine dump()
	character(len=64) :: filename
	character(len=20) :: headers(9)
	integer :: i, fnum, io, ios
	
	headers(1) = "x"
	headers(2) = "probability"
	headers(3) = "visits"
	headers(4) = "births"
	headers(5) = "deaths"
	headers(6) = "Rin_infer"
	headers(7) = "Rin_true"
	headers(8) = "Rout_infer"
	headers(9) = "Rout_true"
	
	call generate_ISO_filename(fpath, fname, ".dat", filename)
	
	open(io, file=trim(filename), status="replace", action="write")
	write(*,*) "Output at: ", filename
	
	! Write metadata
	write(io,*) "# Program Metadata"
	write(io,*) "# Program: inversion_cascade_1.f90"
	write(io,*) "# Creation date: ", fdate()
	write(io,*) "# Seed: ", rseed
	write(io,*) "# Min events: ", event_min
	write(io,*) "# Maximum abundance: ", abund_max
	write(io,*) ""
	write(io,*) ""
	write(io,*) "# Parameter Metadata"
	write(io,*) "# lmbda: ", lmbda
	write(io,*) "# k: ", k
	write(io,*) "# n: ", n
	write(io,*) "# c: ", c
	write(io,*) "# beta: ", beta
	write(io,*) "# x_mean: ", x_mean
	write(io,*) "# rx+_mean: ", r_mean(1)
	write(io,*) "# rx-_avg: ", r_mean(2)
	write(io,*) ""
	write(io,*) ""
	write(io,*) "# Data"
	
	! Write headers	
	write(io, '(9(A20))') trim(headers(1)), &
		trim(headers(2)), &
		trim(headers(3)), &
		trim(headers(4)), &
		trim(headers(5)), &
		trim(headers(6)), &
		trim(headers(7)), &
		trim(headers(8)), &
		trim(headers(9))
		
	do i = 1, x_max+1
		write(io,'(I20, G20.12, I20, I20, I20, G20.12, G20.12, G20.12, G20.12)') &
			i-1, &
			pcond(i), &
			visits_exits(i,1), &
			visits_exits(i,2), &
			visits_exits(i,3), &
			r_infer(i,1), &
			Rin(i-1), &
			r_infer(i,2), &
			Rout(i-1)
	end do
end subroutine


end program
