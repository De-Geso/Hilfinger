! Parameters for mRNA-protein system
module mrna_protein_system_parameters
use kind_parameters
implicit none
public

! Parameters
! ======================================================================
! Size of burst for [x0, x1]
integer, parameter, dimension(3) :: burst = [1, 1, 1]
! mRNA production rate
real(dp), dimension(3) :: lmbda = [1._dp, 2._dp, 3._dp]
! real(dp), dimension(3) :: lmbda = [0._dp, 30._dp, 1._dp]
! Protein production rate
real(dp) :: alpha = 1._dp
! Lifetimes. Can always leave tau_m=1
real(dp), dimension(2) :: tau = [1._dp, 1._dp]
! Decay rates. Can always leave beta_m=1
real(dp), dimension(2) :: beta = [1._dp, 1._dp]
! Protein birth and decay rate powers
real(dp), parameter, dimension(2) :: ell = [1._dp, 1._dp]
! Abundance update matrix.
integer, parameter, dimension(3,6) :: abund_update = &
	reshape((/burst(1), 0, 0, &
			-1, 0, 0, &
			0, burst(2), 0, &
			0, -1, 0, &
			0, 0, burst(3), &
			0, 0, -1/), shape(abund_update))
! Hill function parameters
! Coefficient
real(dp), parameter, dimension(3) :: k = [0.1_dp, 100._dp, 40._dp]
! Power
real(dp), parameter, dimension(3) :: n = [10._dp, 1._dp, 2._dp]

end module
