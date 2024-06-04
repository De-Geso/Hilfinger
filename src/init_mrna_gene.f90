! Initialize parameters for mrna_gene.f90
module init_mrna_gene
use kind_parameters
implicit none
public

! Hyperparameters
! ======================================================================
real(dp), parameter :: pi = 4.D0*DATAN(1.D0)
! Number of events before stopping
integer, parameter :: event_min = 10**6
! Maximum abundances. Program will exit if this is exceeded
integer, parameter :: abund_max = 10**3

! Number of abundance updates to remember for correlation.
! Reducing this gives big time savings.
integer, parameter :: ntail = 2**7
! Length of correlation vector
integer, parameter :: corr_n = 2**1
! Maximum time lag for correlation
real(dp), parameter :: lag_max = 0.1_dp
! Time step for correlation
real(dp), parameter :: corr_tstep = 1._dp*lag_max/(corr_n-1)



!! Parameters
!! ======================================================================
!! Size of burst for [x0, x1]
!integer, parameter, dimension(2) :: burst = [1, 1]
!! x1 production rate
!real(dp) :: alpha = 1._dp
!! x2 production rate
!real(dp) :: beta = 1._dp
!! Decay rates. Can always leave tau_1=1
!real(dp), dimension(2) :: tau = [1._dp, 1._dp]
!! Hill function parameters
!real(dp) :: k = 1._dp
!real(dp) :: n = 1._dp
!! Abundance update matrix.
!integer, parameter, dimension(2,4) :: abund_update = &
!	reshape((/burst(1), 0, &
!			-1, 0, &
!			0, burst(2), &
!			0, -1/), shape(abund_update))

! Variables
! ======================================================================
! Propensity of each event
real(dp), dimension(4) :: propensity = 0.0
! Abundances of each species, number of decay events for each species
integer, dimension(2) :: x = [0, 0], nevents(4)=0
! Probability matrices
real(dp) :: prob_cond(abund_max, abund_max), prob(2, abund_max), prob_rate(abund_max)
! Correlation
real(dp) :: corr(corr_n, 4), corr_mean(2, corr_n, 4), corr_mean2(corr_n, 4)
real(dp) :: dcorr(corr_n)
! Moments
real(dp) :: mean(2), cov(2,2), mean_thry(2), cov_thry(2,2)
! Timers
real(dp) :: ttail(ntail) = 0._dp, t, tstep, tcorr=0._dp

real(dp) :: roll
character(*), parameter :: fout = "mrna_accum.dat"
integer :: xtail(2, ntail) = 0
integer :: i, j, event, io


end module init_mrna_gene
