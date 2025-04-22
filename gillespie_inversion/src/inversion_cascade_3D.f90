program inversion_cascade_3D
use kind_parameters
use stochastics
use randf
use utilities
use stdlib_hashmaps, only: chaining_hashmap_type, int_calls
use stdlib_hashmap_wrappers, only: key_type, set, get, &
	seeded_water_hasher
implicit none

character(len=*), parameter :: fpath = "data/"
character(len=*), parameter :: fname = "dimer_3D"

real(dp), parameter :: eps=tiny(eps)

integer, parameter :: event_min = 10**6

! System parameters ====================================================
real(dp), parameter, dimension(3) :: lmbda = 2._dp
real(dp), parameter, dimension(3) :: beta = 1._dp
real(dp), parameter, dimension(3) :: k = 10._dp
real(dp), parameter, dimension(3) :: n = 2._dp
real(dp), parameter, dimension(3) :: c = 0._dp
integer, parameter, dimension(3,6) :: burst = &
	reshape((/1, 0, 0, &
			-1, 0, 0, &
			0, 1, 0, &
			0, -1, 0, &
			0, 0, 1, &
			0, 0, -1/), shape(burst))
integer, parameter, dimension(3,6) :: abund_update = burst

! Variables ============================================================
! Timers and counters
real(dp) :: t=0._dp, tstep
integer, allocatable :: visits(:)
integer, allocatable :: exits(:,:)
! Stats/calculated
real(dp), allocatable :: r_infer(:,:)
type(OnlineCovariance) :: online_cov_balance(3,3), online_x_cov(3)
real(dp) :: cov_balance(3,3), x_cov(3)
real(dp) :: x_mean(3), r_mean(6)
! Gillespie
real(dp) :: propensity(6), roll
integer :: x(3)=0, event_count(6)=0, event
! Other
integer :: i, j, nseed
integer, allocatable :: rseed(:)

! Hashmap ==============================================================
type :: state_data
	real(dp) :: time_spent
	integer :: visit_count
	integer :: exit_count(6)
end type state_data
type(chaining_hashmap_type) :: map
type(key_type) :: key
type(key_type), allocatable :: keys(:)
type(state_data) :: values
class(*), allocatable :: retrieved
logical :: conflict, key_exists


! Initialize hashmap
call map%init(seeded_water_hasher)

! Initialize random seeding, and get random seed for output in metadata
call random_seed(put=seed)
call random_seed(size=nseed)
allocate(rseed(nseed))
call random_seed(get=rseed)
		
do while (minval(event_count) < event_min)
	! Update the propensity before taking a Gillespie step
	! call update_propensity(propensity, x)
	propensity = rates(x)
	! Get timestep and event from Gillespie algorithm
	call gillespie_iter(tstep, event, propensity)
	
	! Online calculation of mean and covariance
	do i = 1,3
		call online_x_cov(i)%update(real(x(i),dp), real(x(i),dp), tstep)
		do j = 1,3
			call online_cov_balance(i,j)%update(&
				real(x(i),dp), &
				abs(burst(j,2*j))*propensity(2*j) - abs(burst(j,j))*propensity(j), tstep)
		end do
	end do
	
	! Update hashmap
	call set(key, [x])
	call update_hashmap(map, key, tstep, event)
	
	! Update state of system for next step
	t = t + tstep
	x = x + abund_update(:,event)
	event_count(event) = event_count(event) + 1
end do

! Get all the keys (states visited)
call map%get_all_keys(keys)

! Use retrieved keys and values to infer rates
allocate(r_infer(size(keys),size(burst,dim=2)))
do i = 1, size(keys)
	call map%get_other_data(keys(i), retrieved)
	select type (retrieved)
	type is (state_data)	
		r_infer(i,:) = retrieved%exit_count(:) / retrieved%time_spent
	end select
end do

! Stats stuff for validity of simulation, i.e. check covariance balance
x_mean = online_x_cov%mean_x
r_mean = event_count/t
do i = 1,3
	x_cov(i) = online_x_cov(i)%get_cov()
do j = 1,3
	cov_balance(i,j) = online_cov_balance(i,j)%get_cov()
end do
end do

write(*,*) "x mean:", x_mean
write(*,*) "x covariance:", x_cov
write(*,*) "r mean:", r_mean

call check_covariance_balance(3, cov_balance, x_mean, r_mean, burst, "arithmetic")
call dump()


contains


subroutine update_hashmap(this, key, dt, r_channel)
! Put a new entry into the hashmap, or update an entry if it already exists
	type(chaining_hashmap_type), intent(inout) :: this
	type(key_type), intent(in) :: key
	real(dp), intent(in) :: dt
	integer, intent(in) :: r_channel
	type(state_data) :: new
	class(*), allocatable :: other
	logical :: key_exists
	
	! Check key existence
	call this%key_test(key, key_exists)
	! If exists, update entry
	if (key_exists) then
	! other is polymorphic, so need to set type to work with it
		call this%get_other_data(key, other)
		select type (other)
		type is (state_data)
			other%time_spent = other%time_spent + dt
			other%visit_count = other%visit_count + 1
			other%exit_count(event) = other%exit_count(r_channel) + 1
			call this%set_other_data(key, other)
		end select		
	! If doesn't exists, initialize new entry
	else
		new%time_spent = tstep
		new%visit_count = 1
		new%exit_count = 0
		new%exit_count(event) = 1
		call this%map_entry(key, new, conflict)
	end if
end subroutine


subroutine check_covariance_balance(n, cov, x_avg, r_avg, burst, change_method)
	character(len=*), intent(in) :: change_method
	real(dp), intent(in) :: cov(n,n), x_avg(n), r_avg(:)
	integer, intent(in) :: n, burst(n,2*n)
	real(dp) :: rel_change, U(n,n), D(n,n), tau(n), s(n,n), flux_avg(n,2)
	integer :: i, j, k
	
	
	! Get fluxes
	flux_avg = 0
	do i = 1, n
	do k = 1, size(r_avg)
		if (burst(i,k) > 0) then
			flux_avg(i,1) = flux_avg(i,1) + r_avg(k)*abs(burst(i,k))
		elseif (burst(i,k) < 0) then
			flux_avg(i,2) = flux_avg(i,2) + r_avg(k)*abs(burst(i,k))
		end if
	end do
	end do
	write(*,*) "flux: ", flux_avg
	
	! Get lifetimes
	tau = 1._dp * x_avg / flux_avg(:,2)
	write(*,*) "tau: ", tau
	
	! Get step sizes
	s = 0._dp
	do i = 1, n
	do j = 1, n
	do k = 1, size(r_avg)
		s(i,j) = s(i,j) + (r_avg(k)*abs(burst(i,k))/sum(r_avg(:)*abs(burst(i,:))) &
			* abs(burst(j,k)) * sign(1, burst(i,k)*burst(j,k)))
	end do
	end do
	end do
	print *, s
	
	write(*,*) change_method, " relative change between U(i,j)+U(j,i) and D(i,j) (D!=0):"
	do i = 1, n
	do j = 1, n
		if (s(i,j) > eps) then
			U(i,j) = 1._dp/tau(j) * cov(i,j)/(x_avg(i)*0.5*(flux_avg(j,1)+flux_avg(j,2)))
			D(i,j) = s(i,j)/tau(i)/x_avg(j) + s(j,i)/tau(j)/x_avg(i)
			rel_change = relative_change(U(i,j)+U(j,i), D(i,j), change_method)
			write(*,*) "i:", i, "j:", j, "Relative change:", rel_change, "U+U^T: ", U(i,j)+U(j,i), "D:", D(i,j)
		end if
	end do
	end do
end subroutine


pure function rates(x) result(r)
	integer, intent(in) :: x(3)
	real(dp) :: r(6)
	
	r(1) = 1._dp * lmbda(1)
	r(2) = 1._dp * x(1) * beta(1)
	
	r(3) = 1._dp * lmbda(2)
	r(4) = 1._dp * x(2) * beta(2)
	
	r(5) = 1._dp * lmbda(3)
	r(6) = 1._dp * x(3) * beta(3)
end function


subroutine dump()
	character(len=64) :: filename
	character(len=20) :: headers(9)
	integer :: i, fnum, io, ios
	integer, allocatable :: state(:)
	real(dp) :: r(6)
	
	print '("Number of keys in the hashmap = ", I0)', size(keys)
	
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
	! write(filename, "(A,A,I0,A)") fname, ".dat"
	
	open(io, file=trim(filename), status="replace", action="write")
	write(*,*) "Output at: ", filename
	
	! Write metadata
	write(io,*) "# Program Metadata"
	write(io,*) "# Program: inversion_cascade_1.f90"
	write(io,*) "# Creation date: ", fdate()
	write(io,*) "# Seed: ", rseed
	write(io,*) "# Min events: ", event_min
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
		
!	do i = 1, size(keys)
!!		call set(key, keys(i))
!		call map%get_other_data(keys(i), retrieved)
!		select type (retrieved)
!		type is (state_data)
!			write(io,'(I20, G20.12, I20, I20, I20, G20.12, G20.12, G20.12, G20.12)') &
	do i = 1, size(keys)
	    call get(keys(i), state)
	    r = rates(state)
		call map%get_other_data(keys(i), retrieved)
		select type (retrieved)
		type is (state_data)
			write(io,*) & 
				state, &
				retrieved%time_spent/t, &
				retrieved%visit_count, &
				r_infer(i,1), &
				r(1), &
				r_infer(i,2), &
				r(2), &
				r_infer(i,3), &
				r(3), &
				r_infer(i,4), &
				r(4), &
				r_infer(i,5), &
				r(5), &
				r_infer(i,6), &
				r(6)
		end select		
	end do

end subroutine


end program
