! Initialize parameters for mrna_gene.f90
module init_mrna_gene
use kind_parameters
implicit none
public

! Hyperparameters
! ======================================================================
real(dp), parameter :: pi = 4.D0*DATAN(1.D0)
! Number of events before stopping
integer, parameter :: event_min = 10**5
! Maximum abundances. Program will exit if this is exceeded
integer, parameter :: abund_max = 2**7
! Number of abundance updates to remember for autocorrelation
integer, parameter :: ntail = 2**8
! Length of autocorrelation vector
integer, parameter :: corr_n = 2**3
! Maximum time lag for autocorrelation
real(dp), parameter :: lag_max = 2._dp
! Time step for autocorrelation
real(dp), parameter :: corr_tstep = 1._dp*lag_max/corr_n



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
! Abundances of each species, number of decay events for each species
integer, dimension(2) :: x = [0, 0], nevents(4)=0
! Probability matrices
real(dp) :: prob_cond(abund_max, abund_max), prob(2, abund_max), prob_rate(abund_max)
! Autocorrelation
real(dp) :: corr(corr_n), corr_mean(2,corr_n), corr_mean2(corr_n)
! Covariance
real(dp) :: covar(2,2), mean(2)
! Timers
real(dp) :: ttail(ntail) = 0._dp, t, tstep
character(*), parameter :: fout = "mrna_accum.dat"
integer :: xtail(2, ntail) = 0
integer :: i, event, io


end module init_mrna_gene
