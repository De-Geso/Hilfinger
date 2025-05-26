program profile_visit_count
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: n = 10**9  ! Number of data structs
  integer, parameter :: m_small = 10
  integer, parameter :: m_large = 100000
  integer :: i, j
  real(dp) :: t1, t2

  type :: data_small
    real(dp) :: time_spent
    integer :: visit_count(m_small)
  end type

  type :: data_large
    real(dp) :: time_spent
    integer :: visit_count(m_large)
  end type

  type(data_small) :: a_small
  type(data_large) :: a_large

  ! Initialize
  a_small%visit_count = 0
  a_large%visit_count = 0

  ! ==== Test with m_small ====
  call cpu_time(t1)
  do i = 1, n
    do j = 1, m_small
      a_small%visit_count(j) = a_small%visit_count(j) + 1
    end do
  end do
  call cpu_time(t2)
  write(*,'(A,I0,A,1PE13.6,A)') 'Time for visit_count(', m_small, '):', t2 - t1, ' seconds'

  ! ==== Test with m_large (but still only use first 10 bins) ====
  call cpu_time(t1)
  do i = 1, n
    do j = 1, m_small   ! Still only using first 10
      a_large%visit_count(j) = a_large%visit_count(j) + 1
    end do
  end do
  call cpu_time(t2)
  write(*,'(A,I0,A,1PE13.6,A)') 'Time for visit_count(', m_large, '):', t2 - t1, ' seconds'

end program
