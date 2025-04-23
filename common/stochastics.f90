module stochastics
use kind_parameters
use randf
use utilities
implicit none
public

type :: OnlineCovariance
	real(dp) :: total_time = 0._dp
	real(dp) :: mean_x = 0._dp
	real(dp) :: mean_y = 0._dp
	real(dp) :: cov_xy = 0._dp
	
	contains
	
	procedure :: update => update_cov
	procedure :: get_cov => get_covariance
end type OnlineCovariance


contains


subroutine check_covariance_balance(nx, nr, cov, x_avg, r_avg, burst, change_method)
	real(dp), parameter :: eps=tiny(eps)
	character(len=*), intent(in) :: change_method
	real(dp), intent(in) :: cov(nx,nx), x_avg(nx), r_avg(nr)
	integer, intent(in) :: nx, nr, burst(nx, nr)
	real(dp) :: rel_change, U(nx, nx), D(nx, nx), tau(nx), s(nx,nx), flux_avg(nx,2)
	integer :: i, j, k
	
		
	! Get fluxes
	flux_avg = 0
	do i = 1, nx
	do k = 1, nr
		if (burst(i,k) > 0) then
			flux_avg(i,1) = flux_avg(i,1) + r_avg(k)*abs(burst(i,k))
		elseif (burst(i,k) < 0) then
			flux_avg(i,2) = flux_avg(i,2) + r_avg(k)*abs(burst(i,k))
		end if
	end do
	end do
	write(*,*) "flux: ", flux_avg
	
	! Get lifetimes
	tau = 1._dp * x_avg / flux_avg(:,2)
	write(*,*) "tau: ", tau
	
	! Get step sizes
	s = 0._dp
	do i = 1, nx
	do j = 1, nx
	do k = 1, nr
		s(i,j) = s(i,j) + (r_avg(k)*abs(burst(i,k))/sum(r_avg(:)*abs(burst(i,:))) &
			* abs(burst(j,k)) * sign(1, burst(i,k)*burst(j,k)))
	end do
	end do
	end do
	
	! Calculate diffusion and correlation matrices
	do i = 1, nx
	do j = 1, nx
		U(i,j) = 1._dp/tau(j) * cov(i,j)/x_avg(i)/(0.5*(flux_avg(j,1)+flux_avg(j,2)))
		D(i,j) = s(i,j)/tau(i)/x_avg(j) + s(j,i)/tau(j)/x_avg(i)
	end do
	end do
	
	! Output results
	write(*,*) change_method, " relative change between U(i,j)+U(j,i) and D(i,j), or U(i,j) and -U(i,j) if D(i,j)==0:"
	write(*,'(*(A20))') "D(ij)==0","i", "j", "Relative change", "U(ij)+U(ji)", "D(ij)"
	do i = 1, nx
	do j = 1, nx
		if (D(i,j) > eps) then
			rel_change = relative_change(U(i,j)+U(j,i), D(i,j), change_method)
			write(*,'(A20,2(I20),3(G20.12))') "FALSE", i, j, rel_change, U(i,j)+U(j,i), D(i,j)
		else
			rel_change = relative_change(U(i,j), -U(j,i), change_method)
			write(*,'(A20,2(I20),*(G20.12))') "TRUE", i, j, rel_change, U(i,j), -U(j,i)
		end if
	end do
	end do
end subroutine


subroutine gillespie_iter(tstep, event, propensity)
! Run one step of the Gillespie algorithm given the current state of the
! system. Update the propensities, return time step and reaction
	real(dp), intent(out) :: tstep
	real(dp), intent(in) :: propensity(:)
	integer, intent(out) :: event
	real(dp) :: psum, propsum, roll

	! Update propensity
	propsum = sum(propensity)
	
	! Get time step time
	call random_exp(1./propsum, tstep)
	
	! Get reaction
	call random_number(roll)
	event = 0
	psum = 0.
	do while (psum < roll*propsum)
		event = event + 1
		psum = psum + propensity(event)
	end do
end subroutine


subroutine finalize_autocovariance(mean2, meanx, meany, t, acov)
	real(dp), intent(in) :: t
	real(dp), intent(inout) :: mean2(:), meanx(:), meany(:)
	real(dp), intent(out) :: acov(:)
	integer :: i
	
	! Take the time average
	mean2 = mean2 / t
	meanx = meanx / t
	meany = meany / t
	
	do i = 1, size(mean2)
		acov(i) = mean2(i) - meanx(i)*meany(1)
	end do
end subroutine


subroutine update_autocovariance(mean2, meanx, meany, y, x, tvec, maxlag, nwindow, nacov, tstep, dt)
! Update online autocovariance variables (covariance and mean) every time step.
! Don't think too hard about the variable names. Second one is the one lagging.
! Ex. A_ab = update_online_correlation(mean2, meana, meanb, a, b, tvec, dt
	real(dp), intent(in) :: x(nwindow), y(nwindow), tvec(nwindow), maxlag, tstep, dt
	real(dp), intent(inout) :: mean2(nacov), meanx(nacov), meany(nacov)
	integer, intent(in) :: nwindow, nacov
	real(dp) :: t, ta, tb
	integer :: i, j, ti, ita, itb
	
	! Crash if we aren't carrying enough information around
	if (tvec(nwindow)-tvec(1) < maxlag .and. tvec(1) /= 0.) then
		write(*,*) "Error in update_corr: Required lag larger than recorded lag. Increase nwindow."
		write(*,*) "Recorded time lag: ", tvec(nwindow)-tvec(1)
		write(*,*) "Required time lag: ", maxlag
		call exit()
	end if
	
	! For the leading edge (no lag) we always take exactly one step.
	! Integrate a rectangle.
	meanx(1) = meanx(1) + x(nwindow) * dt
	meany(1) = meany(1) + y(nwindow) * dt
	mean2(1) = mean2(1) + x(nwindow)*y(nwindow) * dt
	
	do i = 2, nacov
		! For points with lag, we have to check some things.
		! Range of time we're interested in.
		tb = tvec(nwindow) - (i-1.)*tstep
		ta = tb - dt
		t = ta
		! Most recent points smaller than ta and tb
		ita = maxloc(tvec, dim=1, mask=ta-tvec .gt. 0.)
		itb = maxloc(tvec, dim=1, mask=tb-tvec .gt. 0.)
		! Special case where no step occured in time range. Integrate a rectangle
		if (ita == itb) then
			mean2(i) = mean2(i) + x(ita+1) * y(nwindow) * dt
			meanx(i) = meanx(i) + x(ita+1) * dt
			meany(i) = meany(i) + y(nwindow) * dt
		else
			! ta side
			mean2(i) = mean2(i) + x(ita+1) * y(nwindow) * (tvec(ita+1)-ta)
			meanx(i) = meanx(i) + x(ita+1) * (tvec(ita+1)-ta)
			meany(i) = meany(i) + y(nwindow) * (tvec(ita+1)-ta)
			! tb side
			mean2(i) = mean2(i) + x(itb+1) * y(nwindow) * (tb-tvec(itb))
			meanx(i) = meanx(i) + x(itb+1) * (tb-tvec(itb))
			meany(i) = meany(i) + y(nwindow) * (tb-tvec(itb))
			
			! integrate from ta to tb
			do j = ita+1, itb-1
				mean2(i) = mean2(i) + x(j+1) * y(nwindow) * (tvec(j+1)-tvec(j))
				meanx(i) = meanx(i) + x(j+1) * (tvec(j+1)-tvec(j))
				meany(i) = meany(i) + y(nwindow) * (tvec(j+1)-tvec(j))
			end do
		end if	
	end do
end subroutine


subroutine update_discretized_autocovariance(mean2, meany, meanx, x, y, n, tstep)
! Update online autocovariance variables (covariance and mean) at discrete times.
! Again, don't think too hard about the variable names. Second one is the one lagging.
! Ex. A_ab = update_online_correlation(mean2, meana, meanb, a, b, ...)
	real(dp), intent(in) :: x(n), y(n), tstep
	real(dp), intent(inout) :: mean2(n), meanx(n), meany(n)
	integer, intent(in) :: n
	integer :: i
	
	! Calculate the autocovariance
	do i = 1, n
		meanx(i) = 1._dp * (meanx(i) + x(i) * tstep)
		meany(i) = 1._dp * (meany(i) + y(i) * tstep)
		mean2(i) = 1._dp * (mean2(i) + x(n)*y(n-(i-1)) * tstep)
	end do
end subroutine


subroutine update_cov(this, x, y, dt)
	class(OnlineCovariance), intent(inout) :: this
	real(dp), intent(in) :: x, y, dt
	real(dp) :: dx, dy, w

	if (dt <= 0._dp) return
	
	! Update total time
	this%total_time = this%total_time + dt
	
	! Time-weighted mean updates
	w = dt / this%total_time
	dx = x - this%mean_x
	dy = y - this%mean_y
	this%mean_x = this%mean_x + w * dx
	this%mean_y = this%mean_y + w * dy
	
	! Time-weighted covariance updates
	this%cov_xy = this%cov_xy + dt * dx * (y - this%mean_y)
end subroutine update_cov


function get_covariance(this) result(cov)
! Finalize covariance
	class(OnlineCovariance), intent(in) :: this
	real(dp) :: cov
	
	if (this%total_time > 0._dp) then
		cov = this%cov_xy / (this%total_time)
	else
		cov = 0._dp
	end if
end function get_covariance


pure function hill(x, k, n, c) result(f)
! Hill function. n<0 for negative hill function, n>0 for positive.
	real, intent(in) :: x
	real(dp), intent(in) :: k, n, c
	real(dp) :: f
	
	if (n > 0._dp) then
		f = x**n / (k**n + x**n)
	else if (n < 0._dp) then
		f = k**(-n) / (k**(-n) + x**(-n))
	else
		f = 0.5
	end if	
end function


end module
