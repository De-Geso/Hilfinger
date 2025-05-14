module hashmap_utilities
! Hashmap operations and variables that may be useful across many programs
use kind_parameters
use utilities
use stdlib_hashmaps, only: chaining_hashmap_type, int_calls
use stdlib_hashmap_wrappers, only: key_type, set, get, &
	seeded_water_hasher
implicit none

! Maximum number of reactions in the derived type. This can be very large.
! As long as we aren't saving them, there is negligable performance difference.
! See gillespie_inversion/src/profiling_data_type.f90 
integer, parameter :: max_reactions = 20

type :: state_exits
	real(dp) :: time_spent
	integer :: visit_count
	integer :: exit_count(max_reactions)
end type state_exits


contains


subroutine update_hashmap(this, key, dt, channel)
! Put a new entry into the hashmap, or update an entry if it already exists
	type(chaining_hashmap_type), intent(inout) :: this
	type(key_type), intent(in) :: key
	real(dp), intent(in) :: dt
	integer, intent(in) :: channel
	type(state_exits) :: new
	class(*), allocatable :: other
	logical :: key_exists, conflict
	
	! Check key existence
	call this%key_test(key, key_exists)
	! If exists, get entry and update it
	if (key_exists) then
	! other is polymorphic, so need to set type to work with it
		call this%get_other_data(key, other)
		select type (other)
		type is (state_exits)
			other%time_spent = other%time_spent + dt
			other%visit_count = other%visit_count + 1
			other%exit_count(channel) = other%exit_count(channel) + 1
			call this%set_other_data(key, other)
		end select		
	! If doesn't exists, initialize new entry
	else
		new%time_spent = dt
		new%visit_count = 1
		new%exit_count = 0
		new%exit_count(channel) = 1
		call this%map_entry(key, new, conflict)
	end if
end subroutine


subroutine sort_keys(keys)
! Interface to sort hashmap keys lexigraphically
	type(key_type), intent(inout) :: keys(:)
	if (size(keys) > 1) then
		call quicksort_keys(keys, 1, size(keys))
	end if
end subroutine


subroutine quicksort_keys(keys, left, right)
! Sort hashmap keys lexigraphically using quicksort algorithm.
	type(key_type), intent(inout) :: keys(:)
	integer, intent(in) :: left, right
	integer :: i, j
	type(key_type) :: temp
	integer, allocatable :: pivot(:), key_val(:)
		
	if (left >= right) return
	
	call get(keys((left + right) / 2), pivot)
	i = left
	j = right

	do
		call get(keys(i), key_val)
		do while (lex_less_than(key_val, pivot))
			i = i + 1
			call get(keys(i), key_val)
		end do
		call get(keys(j), key_val)
		do while (lex_less_than(pivot, key_val))
			j = j - 1
			call get(keys(j), key_val)
		end do
		if (i <= j) then
			temp = keys(i)
			keys(i) = keys(j)
			keys(j) = temp
			i = i + 1
			j = j - 1
		end if
		if (i > j) exit
	end do

	call quicksort_keys(keys, left, j)
	call quicksort_keys(keys, i, right)
end subroutine


end module hashmap_utilities
