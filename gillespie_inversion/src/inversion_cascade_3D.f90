program inversion_cascade_3D
use kind_parameters
use stochastics
use randf
use utilities
use hashmap_utilities
use stdlib_hashmaps, only: chaining_hashmap_type, int_calls
use stdlib_hashmap_wrappers, only: key_type, set, get, &
	seeded_water_hasher
implicit none

character(len=*), parameter :: fpath = "data/"
character(len=*), parameter :: fname = "3D_oscillating"

real(dp), parameter :: eps=tiny(eps)

integer, parameter :: event_min = 10**6

! System parameters ====================================================
integer, parameter :: n_species = 3, n_reactions = 6
real(dp), parameter, dimension(n_species) :: lmbda = [500._dp, 80._dp, 80._dp]
real(dp), parameter, dimension(n_species) :: beta = [1._dp, 1._dp, 1._dp]
real(dp), parameter, dimension(n_species) :: k = [0.1_dp, 100._dp, 40._dp]
real(dp), parameter, dimension(n_species) :: n = [-10._dp, 1._dp, 2._dp]
real(dp), parameter, dimension(n_species) :: c = [0._dp, 0._dp, 0._dp]
integer, parameter, dimension(n_species, n_reactions) :: burst = reshape( &
	(/1, 0, 0, & ! Reaction 1
	-1, 0, 0, &
	0, 1, 0, &
	0, -1, 0, &
	0, 0, 1, &
	0, 0, -1/), & ! Reaction n
	shape(burst))

! Variables ============================================================
! Timers and counters
real(dp) :: t=0._dp, tstep
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
!type :: state_exits
!	real(dp) :: time_spent
!	integer :: visit_count
!	integer :: exit_count(10*n_reactions)
!end type state_exits
type(chaining_hashmap_type) :: map
type(key_type) :: key
type(key_type), allocatable :: keys(:)
type(state_exits) :: values
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
	call set(key, x)
	call update_hashmap(map, key, tstep, event)
	
	! Update state of system for next step
	t = t + tstep
	x = x + burst(:,event)
	event_count(event) = event_count(event) + 1
end do

! Get all the keys (states visited)
call map%get_all_keys(keys)
call sort_keys(keys)

! Use retrieved keys and values to infer rates
allocate(r_infer(size(keys), n_reactions))
do i = 1, size(keys)
	call map%get_other_data(keys(i), retrieved)
	select type (retrieved)
	type is (state_exits)	
!		r_infer(i,:) = retrieved%exit_count(:) / retrieved%time_spent
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


pure function rates(x) result(r)
	integer, intent(in) :: x(n_species)
	real(dp) :: r(n_reactions)
	
	! Oscillating
	r(1) = 1._dp * lmbda(1) * hill(real(x(3), dp), k(1), n(1), c(1))
	r(3) = 1._dp * lmbda(2) * hill(real(x(1), dp), k(2), n(2), c(2))

	! f_3(x_2) doesn't change
	r(5) = 1._dp * lmbda(3) * hill(real(x(2), dp), k(3), n(3), c(3))
	! Linear decay doesn't change
	r(2) = 1._dp * x(1) * beta(1)
	r(4) = 1._dp * x(2) * beta(2)
	r(6) = 1._dp * x(3) * beta(3)
end function


!subroutine update_hashmap(this, key, dt, r_channel)
!! Put a new entry into the hashmap, or update an entry if it already exists
!! Would love to make this general, but the derived type makes this difficult
!	type(chaining_hashmap_type), intent(inout) :: this
!	type(key_type), intent(in) :: key
!	real(dp), intent(in) :: dt
!	integer, intent(in) :: r_channel
!	type(state_exits) :: new
!	class(*), allocatable :: other
!	logical :: key_exists
	
!	! Check key existence
!	call this%key_test(key, key_exists)
!	! If exists, update entry
!	if (key_exists) then
!	! other is polymorphic, so need to set type to work with it
!		call this%get_other_data(key, other)
!		select type (other)
!		type is (state_exits)
!			other%time_spent = other%time_spent + dt
!			other%visit_count = other%visit_count + 1
!			other%exit_count(event) = other%exit_count(r_channel) + 1
!			call this%set_other_data(key, other)
!		end select		
!	! If doesn't exists, initialize new entry
!	else
!		new%time_spent = tstep
!		new%visit_count = 1
!		new%exit_count = 0
!		new%exit_count(event) = 1
!		call this%map_entry(key, new, conflict)
!	end if
!end subroutine


!subroutine sort_keys(keys)
!	type(key_type), intent(inout) :: keys(:)
!	if (size(keys) > 1) then
!		call quicksort_keys(keys, 1, size(keys))
!	end if
!end subroutine


!subroutine quicksort_keys(keys, left, right)
!	type(key_type), intent(inout) :: keys(:)
!	integer, intent(in) :: left, right
!	integer :: i, j
!	type(key_type) :: temp
!	integer, allocatable :: pivot(:), key_val(:)
		
!	if (left >= right) return
	
!	call get(keys((left + right) / 2), pivot)
!	i = left
!	j = right

!	do
!		call get(keys(i), key_val)
!		do while (lex_less_than(key_val, pivot))
!			i = i + 1
!			call get(keys(i), key_val)
!		end do
!		call get(keys(j), key_val)
!		do while (lex_less_than(pivot, key_val))
!			j = j - 1
!			call get(keys(j), key_val)
!		end do
!		if (i <= j) then
!			temp = keys(i)
!			keys(i) = keys(j)
!			keys(j) = temp
!			i = i + 1
!			j = j - 1
!		end if
!		if (i > j) exit
!	end do

!	call quicksort_keys(keys, left, j)
!	call quicksort_keys(keys, i, right)
!end subroutine


subroutine dump()
	character(len=64) :: fname_rates, fname_exits
	character(len=20) :: headers_rates(17), headers_exits(11)
	integer :: i, fnum, io, ios
	integer, allocatable :: state(:)
	real(dp) :: r(6)
	
	print '("Number of keys in the hashmap = ", I0)', size(keys)
	
	headers_rates(1) = "x1"
	headers_rates(2) = "x2"
	headers_rates(3) = "x3"
	headers_rates(4) = "visits"
	headers_rates(5) = "probability"
	headers_rates(6) = "r1_inf"
	headers_rates(7) = "r1_true"
	headers_rates(8) = "r2_inf"
	headers_rates(9) = "r2_true"
	headers_rates(10) = "r3_inf"
	headers_rates(11) = "r3_true"
	headers_rates(12) = "r4_inf"
	headers_rates(13) = "r4_true"
	headers_rates(14) = "r5_inf"
	headers_rates(15) = "r5_true"
	headers_rates(16) = "r6_inf"
	headers_rates(17) = "r6_true"
	
	headers_exits(1) = "x1"
	headers_exits(2) = "x2"
	headers_exits(3) = "x3"
	headers_exits(4) = "visits"
	headers_exits(5) = "probability"
	headers_exits(6) = "exits1"
	headers_exits(7) = "exits2"
	headers_exits(8) = "exits3"
	headers_exits(9) = "exits4"
	headers_exits(10) = "exits5"
	headers_exits(11) = "exits6"
	
	call generate_ISO_filename(fpath, trim(fname)//"_rates", ".dat", fname_rates)
	call generate_ISO_filename(fpath, trim(fname)//"_exits", ".dat", fname_exits)
	write(*,*) "Output at: ", fname_rates, fname_exits

	! Write rate file
	open(io, file=fname_rates, status="replace", action="write")
	! Write metadata
	write(io,*) "# Program Metadata"
	write(io,*) "# Program: inversion_cascade_3D.f90"
	write(io,*) "# Creation date: ", fdate()
	write(io,*) "# Seed: ", rseed
	write(io,*) "# Min events: ", event_min
	write(io,*) "# n_species: ", n_species
	write(io,*) "# n_reactions: ", n_reactions
	write(io,*) "# time: ", t
	write(io,*) ""
	write(io,*) ""
	write(io,*) "# Parameter Metadata"
	write(io,*) "# burst: ", burst
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
	
	! Write rates headers
	write(io, '(17(A20))') trim(headers_rates(1)), &
		trim(headers_rates(2)), &
		trim(headers_rates(3)), &
		trim(headers_rates(4)), &
		trim(headers_rates(5)), &
		trim(headers_rates(6)), &
		trim(headers_rates(7)), &
		trim(headers_rates(8)), &
		trim(headers_rates(9)), &
		trim(headers_rates(10)), &
		trim(headers_rates(11)), &
		trim(headers_rates(12)), &
		trim(headers_rates(13)), &
		trim(headers_rates(14)), &
		trim(headers_rates(15)), &
		trim(headers_rates(16)), &
		trim(headers_rates(17))
		
	do i = 1, size(keys)
	    call get(keys(i), state)
	    r = rates(state)
		call map%get_other_data(keys(i), retrieved)
		select type (retrieved)
		type is (state_exits)
			write(io,'(4(I20), 13(ES20.12))') & 
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
	close(io)
	
	! Write exits file
	open(io, file=fname_exits, status="replace", action="write")
	! Write metadata
	write(io,*) "# Program Metadata"
	write(io,*) "# Program: inversion_cascade_3D.f90"
	write(io,*) "# Creation date: ", fdate()
	write(io,*) "# Seed: ", rseed
	write(io,*) "# Min events: ", event_min
	write(io,*) "# n_species: ", n_species
	write(io,*) "# n_reactions: ", n_reactions
	write(io,*) "# time: ", t
	write(io,*) ""
	write(io,*) ""
	write(io,*) "# Parameter Metadata"
	write(io,*) "# burst: ", burst
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
	
	! Write exits headers
	write(io, '(17(A20))') trim(headers_exits(1)), &
		trim(headers_exits(2)), &
		trim(headers_exits(3)), &
		trim(headers_exits(4)), &
		trim(headers_exits(5)), &
		trim(headers_exits(6)), &
		trim(headers_exits(7)), &
		trim(headers_exits(8)), &
		trim(headers_exits(9)), &
		trim(headers_exits(10)), &
		trim(headers_exits(11))
				
	do i = 1, size(keys)
	    call get(keys(i), state)
	    r = rates(state)
		call map%get_other_data(keys(i), retrieved)
		select type (retrieved)
		type is (state_exits)
			write(io,'(4(I20), ES20.12, 6(I20))') & 
				state, &
				retrieved%visit_count, &
				retrieved%time_spent/t, &
				retrieved%exit_count(1), &
				retrieved%exit_count(2), &
				retrieved%exit_count(3), &
				retrieved%exit_count(4), &
				retrieved%exit_count(5), &
				retrieved%exit_count(6)
		end select		
	end do
	close(io)

end subroutine


end program
