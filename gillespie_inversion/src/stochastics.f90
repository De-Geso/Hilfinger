module stochastics
use kind_parameters
use randf
implicit none
public


contains


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


end module
