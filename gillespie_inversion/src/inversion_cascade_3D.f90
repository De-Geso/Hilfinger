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
integer, parameter :: n_species = 3, n_reactions = 6
real(dp), parameter, dimension(n_species) :: lmbda = 2._dp
real(dp), parameter, dimension(n_species) :: beta = 1._dp
real(dp), parameter, dimension(n_species) :: k = 10._dp
real(dp), parameter, dimension(n_species) :: n = 2._dp
real(dp), parameter, dimension(n_species) :: c = 0._dp
! integer, parameter, dimension(n_species, n_reactions) :: burst = &
integer, parameter, dimension(n_species, n_reactions) :: burst = &
	reshape((/1, 0, 0, &
			-1, 0, 0, &
			0, 1, 0, &
			0, -1, 0, &
			0, 0, 1, &
			0, 0, -1/), shape(burst))

! Variables ============================================================
! Timers and counters
real(dp) :: t=0._dp, tstep
integer, allocatable :: visits(:)
integer, allocatable :: exits(:,:)
! Stats/calculated
real(dp), allocatable :: r_infer(:,:)
type(OnlineCovariance) :: online_cov_balance(n_species, n_species), online_x_cov(n_species)
real(dp) :: cov_balance(n_species, n_species), x_cov(n_species)
real(dp) :: x_mean(n_species), r_mean(n_reactions)
! Gillespie
real(dp) :: propensity(n_reactions), roll
integer :: x(n_species)=0, event_count(n_reactions)=0, event
! Other
integer :: i, j, nseed
integer, allocatable :: rseed(:)

! Hashmap ==============================================================
type :: state_data
	real(dp) :: time_spent
	integer :: visit_count
	integer :: exit_count(n_reactions)
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
	do i = 1, n_species
		call online_x_cov(i)%update(real(x(i),dp), real(x(i),dp), tstep)
		do j = 1, n_species
			call online_cov_balance(i,j) & 
				%update(real(x(i),dp), &
				abs(burst(j,2*j))*propensity(2*j) - abs(burst(j,2*j-1))*propensity(2*j-1), tstep)
		end do
	end do
	
	! Update hashmap
	call set(key, [x])
	call update_hashmap(map, key, tstep, event)
	
	! Update state of system for next step
	t = t + tstep
	x = x + burst(:,event)
	event_count(event) = event_count(event) + 1
end do

! Get all the keys (states visited)
call map%get_all_keys(keys)

! Use retrieved keys and values to infer rates
allocate(r_infer(size(keys), n_reactions))
do i = 1, size(keys)
	call map%get_other_data(keys(i), retrieved)
	select type (retrieved)
	type is (state_data)	
		r_infer(i,:) = retrieved%exit_count(:) / retrieved%time_spent
	end select
end do

! Stats for simulation sampling, i.e. check covariance balance
x_mean = online_x_cov%mean_x
r_mean = event_count/t
do i = 1, n_species
	x_cov(i) = online_x_cov(i)%get_cov()
do j = 1, n_species
	cov_balance(i,j) = online_cov_balance(i,j)%get_cov()
end do
end do

write(*,*) "x mean:", x_mean
write(*,*) "x covariance:", x_cov
write(*,*) "r mean:", r_mean

call check_covariance_balance(n_species, n_reactions, cov_balance, x_mean, r_mean, burst, "arithmetic")
call dump()


contains


subroutine update_hashmap(this, key, dt, r_channel)
! Put a new entry into the hashmap, or update an entry if it already exists
! Would love to make this general, but the derived type makes this difficult
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


pure function rates(x) result(r)
	integer, intent(in) :: x(n_species)
	real(dp) :: r(n_reactions)
	
	r(1) = 1._dp * lmbda(1)
	r(2) = 1._dp * x(1) * beta(1)
	
	r(3) = 1._dp * lmbda(2) * x(1)
	r(4) = 1._dp * x(2) * beta(2)
	
	r(5) = 1._dp * lmbda(3) * x(2)
	r(6) = 1._dp * x(3) * beta(3)
end function


subroutine dump()
	character(len=64) :: filename
	character(len=20) :: headers(17)
	integer :: i, fnum, io, ios
	integer, allocatable :: state(:)
	real(dp) :: r(6)
	
	print '("Number of keys in the hashmap = ", I0)', size(keys)
	
	headers(1) = "x1"
	headers(2) = "x2"
	headers(3) = "x3"
	headers(4) = "visits"
	headers(5) = "probability"
	headers(6) = "r1_inf"
	headers(7) = "r1_true"
	headers(8) = "r2_inf"
	headers(9) = "r2_true"
	headers(10) = "r3_inf"
	headers(11) = "r3_true"
	headers(12) = "r4_inf"
	headers(13) = "r4_true"
	headers(14) = "r5_inf"
	headers(15) = "r5_true"
	headers(16) = "r6_inf"
	headers(17) = "r6_true"
	
	call generate_ISO_filename(fpath, fname, ".dat", filename)
	
	open(io, file=trim(filename), status="replace", action="write")
	write(*,*) "Output at: ", filename
	
	! Write metadata
	write(io,*) "# Program Metadata"
	write(io,*) "# Program: inversion_cascade_3D.f90"
	write(io,*) "# Creation date: ", fdate()
	write(io,*) "# Seed: ", rseed
	write(io,*) "# Min events: ", event_min
	write(io,*) "# n_species: ", n_species
	write(io,*) "# n_reactions: ", n_reactions
	write(io,*) ""
	write(io,*) ""
	write(io,*) "# Parameter Metadata"
	write(io,*) "# Burst: ", burst
	write(io,*) "# lmbda: ", lmbda
	write(io,*) "# k: ", k
	write(io,*) "# n: ", n
	write(io,*) "# c: ", c
	write(io,*) "# beta: ", beta
	write(io,*) "# x_mean: ", x_mean
	write(io,*) "# r_mean: ", r_mean
	write(io,*) ""
	write(io,*) ""
	write(io,*) "# Data"
	
	! Write headers	
	write(io, '(17(A20))') trim(headers(1)), &
		trim(headers(2)), &
		trim(headers(3)), &
		trim(headers(4)), &
		trim(headers(5)), &
		trim(headers(6)), &
		trim(headers(7)), &
		trim(headers(8)), &
		trim(headers(9)), &
		trim(headers(10)), &
		trim(headers(11)), &
		trim(headers(12)), &
		trim(headers(13)), &
		trim(headers(14)), &
		trim(headers(15)), &
		trim(headers(16)), &
		trim(headers(17))
		
	do i = 1, size(keys)
	    call get(keys(i), state)
	    r = rates(state)
		call map%get_other_data(keys(i), retrieved)
		select type (retrieved)
		type is (state_data)
			write(io,'(4(I20), 13(G20.12))') & 
				state, &
				retrieved%visit_count, &
				retrieved%time_spent/t, &
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
