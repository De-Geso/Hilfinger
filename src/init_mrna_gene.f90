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
real(dp) :: alpha = 1._dp
! x2 production rate
real(dp) :: beta = 1._dp
! Decay rates
real(dp), dimension(2) :: tau = [1._dp, 1._dp]
! Hill function parameters
real(dp) :: k = 1._dp
real(dp) :: n = 1._dp

! Abundance update matrix.
integer, parameter, dimension(2,4) :: abund_update = &
	reshape((/burst(1), 0, &
			-1, 0, &
			0, burst(2), &
			0, -1/), shape(abund_update))

! Variables
! ======================================================================
! Propensity of each event
real(dp), dimension(4) :: propensity = 0.0
! (initial) Abundances of each species, number of decay events for each species
integer, dimension(2) :: x = [0, 0], ndecay = [0, 0]
! Probability matrix
real(dp) :: prob_abund(abund_max, abund_max)
real(dp) :: prob_rate(abund_max)
integer :: i, event, io
real(dp) :: t
character(*), parameter :: fout = "mrna_data.dat"

end module init_mrna_gene
