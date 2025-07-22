program inversion_cascade_2D
use kind_parameters
use stochastics
use randf
use utilities
use hashmap_utilities
use stdlib_hashmaps, only: chaining_hashmap_type, int_calls
use stdlib_hashmap_wrappers, only: key_type, set, get, &
	seeded_water_hasher, fnv_1_hasher
implicit none

character(len=*), parameter :: fpath = "data/"
character(len=*), parameter :: fname = "2D_simple"

real(dp), parameter :: eps=tiny(eps)

integer, parameter :: event_min = 10**6

! System parameters ====================================================
integer, parameter :: n_species = 2, n_reactions = 4
real(dp), parameter :: lmbda(2) = [1._dp, 1._dp]
real(dp), parameter :: beta(2) = [100._dp, 100._dp]
real(dp), parameter :: k(2) = [0._dp, 0._dp]
real(dp), parameter :: n(2) = [0._dp, 0._dp]
real(dp), parameter :: c(2) = [0._dp, 0._dp]
integer, parameter, dimension(n_species, n_reactions) :: burst = reshape( & 
	(/1, 0, & ! Reaction 1
	-1, 0, &
	0, 1, &
	0, -1/), & ! Reaction n
	shape(burst))

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
real(dp) :: fluxes(n_species,2)
integer :: i, j, nseed
integer, allocatable :: rseed(:)
! Hashmap
type(chaining_hashmap_type) :: state_map
type(key_type) :: key
type(key_type), allocatable :: keys(:)
type(state_exits) :: values
class(*), allocatable :: retrieved
logical :: conflict, exists


! Here the program begins ==============================================


! Initialize hashmap
call state_map%init(seeded_water_hasher, slots_bits=10)

! Initialize seed and get seed for metadata output
call random_seed(put=seed)
call random_seed(size=nseed)
allocate(rseed(nseed))
call random_seed(get=rseed)

! Simulate system
do while (minval(event_count) < event_min)	
	! print *, event_count, x
	
	! Update the propensity before taking a Gillespie step
	propensity = rates(x)
	! Get timestep and event from Gillespie algorithm
	call gillespie_iter(tstep, event, propensity)
	
	! Online calculation of mean and covariance
	fluxes = flux(propensity, burst, n_reactions, n_species)
	do i = 1, n_species
		call online_x_cov(i)%update(real(x(i),dp), real(x(i),dp), tstep)
		do j = 1, n_species
			call online_cov_balance(i,j)%update &
				(real(x(i),dp), fluxes(j,2) - fluxes(j,1), tstep)
		end do
	end do

	! Update hashmap
	call set(key, x)
	call update_hashmap(state_map, key, tstep, event)
			
	! Update state of system in preparation for next step.
	t = t + tstep
	event_count(event) = event_count(event) + 1
	x = x + burst(:,event)
end do

! SIMULATION DONE

! Get and sort all the keys (states visited)
call state_map%get_all_keys(keys)
call sort_keys(keys)

! Use retrieved keys and values to infer rates
allocate(r_infer(size(keys), n_reactions))
do i = 1, size(keys)
	call state_map%get_other_data(keys(i), retrieved)
	select type (retrieved)
	type is (state_exits)	
		r_infer(i,:) = retrieved%exit_count(:n_reactions) / retrieved%time_spent
	end select
end do

! Stats to check simulation sampling, i.e. check covariance balance
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


pure function rates(x) result (r)
	integer, intent(in) :: x(n_species)
	real(dp) :: r(n_reactions)
	
	! Production rates
	r(1) = 1._dp * lmbda(1)
	r(3) = 1._dp * lmbda(2) * x(1)
	! r(3) = 1._dp * lmbda(2) * hill(real(x(1)*x(2), dp), k(2), n(2))
	
	! Decay rates
	r(2) = 1._dp * beta(1) * x(1)
	r(4) = 1._dp * beta(2) * x(2)
end function rates


subroutine dump()
	character(len=64) :: fname_rates, fname_exits
	character(len=20) :: headers_rates(12), headers_exits(8)
	integer :: i, fnum, io, ios
	integer, allocatable :: state(:)
	real(dp) :: r(n_reactions)
	
	print '("Number of keys in the hashmap = ", I0)', size(keys)
	
	headers_rates(1) = "x1"
	headers_rates(2) = "x2"
	headers_rates(3) = "visits"
	headers_rates(4) = "probability"
	headers_rates(5) = "r1_inf"
	headers_rates(6) = "r1_true"
	headers_rates(7) = "r2_inf"
	headers_rates(8) = "r2_true"
	headers_rates(9) = "r3_inf"
	headers_rates(10) = "r3_true"
	headers_rates(11) = "r4_inf"
	headers_rates(12) = "r4_true"
	
	headers_exits(1) = "x1"
	headers_exits(2) = "x2"
	headers_exits(3) = "visits"
	headers_exits(4) = "probability"
	headers_exits(5) = "exits1"
	headers_exits(6) = "exits2"
	headers_exits(7) = "exits3"
	headers_exits(8) = "exits4"
	
	call generate_ISO_filename(fpath, trim(fname)//"_rates", ".dat", fname_rates)
	call generate_ISO_filename(fpath, trim(fname)//"_exits", ".dat", fname_exits)
	write(*,*) "Output at: ", fname_rates, fname_exits

	! Write rate file
	open(io, file=fname_rates, status="replace", action="write")
	! Write metadata
	write(io,*) "# Program Metadata"
	write(io,*) "# Program: inversion_cascade_2D.f90"
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
	write(io, '(12(A20))') trim(headers_rates(1)), &
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
		trim(headers_rates(12))
		
	do i = 1, size(keys)
	    call get(keys(i), state)
	    r = rates(state)
		call state_map%get_other_data(keys(i), retrieved)
		select type (retrieved)
		type is (state_exits)
			write(io,'(3(I20), 9(ES20.12))') & 
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
				r(4)
		end select		
	end do
	close(io)
	
	! Write exits file
	open(io, file=fname_exits, status="replace", action="write")
	! Write metadata
	write(io,*) "# Program Metadata"
	write(io,*) "# Program: inversion_cascade_2D.f90"
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
		trim(headers_exits(8))
				
	do i = 1, size(keys)
	    call get(keys(i), state)
	    r = rates(state)
		call state_map%get_other_data(keys(i), retrieved)
		select type (retrieved)
		type is (state_exits)
			write(io,'(3(I20), ES20.12, 6(I20))') & 
				state, &
				retrieved%visit_count, &
				retrieved%time_spent/t, &
				retrieved%exit_count(1), &
				retrieved%exit_count(2), &
				retrieved%exit_count(3), &
				retrieved%exit_count(4)
		end select		
	end do
	close(io)

end subroutine

end program
