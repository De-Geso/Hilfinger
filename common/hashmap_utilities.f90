module hashmap_utilities
! Hashmap operations and variables that may be useful across many programs
use stdlib_hashmaps, only: chaining_hashmap_type, int_calls
use stdlib_hashmap_wrappers, only: key_type, set, get, &
	seeded_water_hasher
implicit none

! Maximum number of reactions in the derived type. This can be very large.
! As long as we aren't saving them, there is negligable performance difference.
! See gillespie_inversion/src/profiling_data_type.f90 
integer, parameter :: max_reactions = 100

type :: state_exits
	real(dp) :: time_spent
	integer :: visit_count
	integer :: exit_count(max_reactions)
end type state_exits


contains


subroutine update_hashmap(this, key, dt, r_channel)
! Put a new entry into the hashmap, or update an entry if it already exists
! Would love to make this general, but the derived type makes this difficult
	type(chaining_hashmap_type), intent(inout) :: this
	type(key_type), intent(in) :: key
	real(dp), intent(in) :: dt
	integer, intent(in) :: r_channel
	type(state_exits) :: new
	class(*), allocatable :: other
	logical :: key_exists
	
	! Check key existence
	call this%key_test(key, key_exists)
	! If exists, update entry
	if (key_exists) then
	! other is polymorphic, so need to set type to work with it
		call this%get_other_data(key, other)
		select type (other)
		type is (state_exits)
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


end module hashmap_utilities
