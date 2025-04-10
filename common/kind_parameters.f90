! Numerical storage size parameters for real and integer values.
module kind_parameters
implicit none
public

! Single precision real numbers, 6 digits, range 10^-37 to 10^37-1; 32 bits
integer, parameter :: sp = selected_real_kind(6,37)
! Double precision real numbers, 15 digits, range 10^-307 to 10^307-1; 64 bits
integer, parameter :: dp = selected_real_kind(15,307)
! Quad precision real numbers, 33 digits, range 10^-4931 to 10^4931-1; 128 bits
integer, parameter :: qp = selected_real_kind(33,4931)

! Char length for integers, range -2^7 to 2^7-1; 8 bits
integer, parameter :: i1 = selected_int_kind(2)
! Short length for integers, range -2^15 to 2^15-1; 16 bits
integer, parameter :: i2 = selected_int_kind(4)
! Length of default integers, range -2^31 to 2^31-1; 32 bits
integer, parameter :: i4 = selected_int_kind(9)
! Long length for integers, range -2^63 to 2^63-1; 64 bits
integer, parameter :: i8 = selected_int_kind(18)

end module kind_parameters
