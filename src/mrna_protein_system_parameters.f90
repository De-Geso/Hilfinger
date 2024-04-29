! Parameters for mRNA-protein system
module mrna_protein_system_parameters
use kind_parameters
implicit none
public

! Parameters
! ======================================================================
! Size of burst for [x0, x1]
integer, parameter, dimension(2) :: burst = [1, 1]
! x1 production rate
real(dp) :: alpha = 1._dp
! x2 production rate
real(dp) :: beta = 1._dp
! Decay rates. Can always leave tau_1=1
real(dp), dimension(2) :: tau = [1._dp, 1._dp]
! Abundance update matrix.
integer, parameter, dimension(2,4) :: abund_update = &
	reshape((/burst(1), 0, &
			-1, 0, &
			0, burst(2), &
			0, -1/), shape(abund_update))
! Hill function parameters
real(dp) :: k = 1._dp
real(dp) :: n = 1._dp

end module
