! Initialize parameters for mrna_gene.f90
module init_mrna_gene
use kind_parameters
implicit none
public

! Hyperparameters
! ======================================================================
! Number of decay events before stopping
integer, parameter :: decay_min = 10**6
! Maximum abundances, pad this.
integer, parameter :: abund_max = 2**5

! Parameters
! ======================================================================
! Size of burst for [x0, x1]
integer, parameter, dimension(2) :: burst = [1, 1]
! x1 production rate
real(dp) :: alpha = 1.0
! x2 production rate
real(dp) :: beta = 1.0
! Decay rates
real(dp), dimension(2) :: decay = [1.0, 1.0]
! Abundance update matrix.
integer, parameter, dimension(2,4) :: abund_update = &
	reshape((/burst(1), 0, &
			-1, 0, &
			0, burst(2), &
			0, -1/), shape(abund_update))

end module
