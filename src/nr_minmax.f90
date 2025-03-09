module nr_minmax
! Subroutines and functions from minimizing and maximizing functions
! from Numerical Recipes
use kind_parameters
implicit none
public


abstract interface
	function funk_type(x)
		use kind_parameters
		real(dp), intent(in), dimension(:) :: x
		real(dp) :: funk_type
	end function funk_type
end interface


contains


! Numerical recipes Sec. 10.4
subroutine amoeba(p, y, mp, np, ndim, ftol, funk, iter)
	real(dp), intent(in) :: ftol
	integer, intent(in) ::  mp, np, ndim
	real(dp), intent(inout) :: p(mp,np), y(mp)
	integer, intent(inout) :: iter
	! Maximum allowed dimensions and function evaluations, and a small number.
	real(dp), parameter :: TINY=1.e-10
	integer, parameter :: ITMAX=5000
	procedure(funk_type) :: funk
!	real :: funk
!	external :: funk
!	USES amotry,funk
!	Multidimensional minimization of the function funk(x) where x(1:ndim) is a vector
!	in ndim dimensions, by the downhill simplex method of Nelder and Mead. The matrix
!	p(1:ndim+1,1:ndim) is input. Its ndim+1 rows are ndim-dimensional vectors which are
!	the vertices of the starting simplex. Also input is the vector y(1:ndim+1), whose compo-
!	nents must be pre-initialized to the values of funk evaluated at the ndim+1 vertices (rows)
!	of p; and ftol the fractional convergence tolerance to be achieved in the function value
!	(n.b.!). On output, p and y will have been reset to ndim+1 new points all within ftol of
!	a minimum function value, and iter gives the number of function evaluations taken.
	integer :: i, ihi, ilo, inhi, j, m, n
	real(dp) :: rtol, sum, swap, ysave, ytry, psum(ndim)
	
	iter = 0
!	Enter here when starting or have just overall contracted.
1	do n = 1, ndim
		! Recompute psum.
		sum = 0.
		do m = 1, ndim+1
			sum = sum + p(m,n)
		end do
		psum(n) = sum
	end do

	! Enter here when we have just changed a single point.
2	ilo = 1 
	! Determine which point is the highest (worst), next-highest, and lowest (best)
	if (y(1) .gt. y(2)) then
		ihi = 1
		inhi = 2
	else
		ihi = 2
		inhi = 1
	end if
	
	! by looping over the points in the simplex.
	do i = 1, ndim + 1
		if (y(i) .le. y(ilo)) ilo = i
		if (y(i) .gt. y(ihi)) then
			inhi = ihi
			ihi = i
		else if (y(i) .gt. y(inhi)) then
			if (i .ne. ihi) inhi = i
		end if
	end do

	write(*,*) iter, p(1,:), p(2,:), p(3,:), y(ihi)
	
!	Compute the fractional range from highest to lowest and return if satisfactory.
	rtol = 2.*abs(y(ihi)-y(ilo)) / (abs(y(ihi))+abs(y(ilo))+TINY)
!	If returning, put best point and value in slot 1.
	if (rtol .lt. ftol) then
		swap = y(1)
		y(1) = y(ilo)
		y(ilo) = swap
		do n = 1, ndim
			swap = p(1,n)
			p(1,n) = p(ilo,n)
			p(ilo,n) = swap
		end do
		return
	end if

	if (iter .ge. ITMAX) then
		write(*,*) p(1,:), y(1)
		stop 'ITMAX exceeded in amoeba'
	end if
		
	iter = iter + 2
!	Begin a new iteration. First extrapolate by a factor −1 through the face of the simplex across
!	from the high point, i.e., reflect the simplex from the high point.
	call amotry(ytry, p, y, psum, mp, np, ndim, funk, ihi, -1._dp)
	if (ytry .le. y(ilo)) then
!		Gives a result better than the best point, so try an additional extrapolation by a factor 2.
		call amotry(ytry, p, y, psum, mp, np, ndim, funk, ihi, 2._dp)
	else if (ytry .ge. y(inhi)) then
!		The reflected point is worse than the second-highest, so look for an intermediate lower point,
!		i.e., do a one-dimensional contraction.
		ysave = y(ihi)
		call amotry(ytry, p, y, psum, mp, np, ndim, funk, ihi, 0.5_dp)
		if (ytry .ge. ysave) then
!			Can’t seem to get rid of that high point. Better contract around the lowest (best) point
			do i = 1, ndim+1
				if(i .ne. ilo)then
					do j = 1, ndim
						psum(j) = 0.5*(p(i,j)+p(ilo,j))
						p(i,j) = psum(j)
					end do
					y(i) = funk(psum)
				end if
			end do
!			Keep track of function evaluations.
			iter = iter + ndim
!			Go back for the test of doneness and the next iteration.
			go to 1
		end if
	else
!	Correct the evaluation count.
		iter=iter-1
	end if
	goto 2
end subroutine


subroutine amotry(res, p, y, psum, mp, np, ndim, funk, ihi, fac)
	real(dp), intent(out) :: res
	real(dp), intent(in) :: fac
	integer, intent(in) :: ihi, mp, ndim, np
	real(dp), intent(inout) :: p(mp,np), psum(np), y(mp)
	real(dp), parameter :: TINY=1.e-8
	procedure(funk_type) :: funk
!	real :: funk
!	external :: funk
!	USES funk
!	Extrapolates by a factor fac through the face of the simplex across from the high point,
!	tries it, and replaces the high point if the new point is better.	
	integer :: j
	real(dp) :: fac1, fac2, ptry(ndim), ytry
	
	ptry = 0.
	fac1 = (1.-fac)/ndim
	fac2 = fac1-fac
	do j = 1, ndim
		ptry(j) = psum(j)*fac1 - p(ihi,j)*fac2
	end do
!	Evaluate the function at the trial point.
	ytry=funk(ptry)
!	If it’s better than the highest, then replace the highest.
	if (ytry .lt. y(ihi)) then
		y(ihi)=ytry
		do j=1, ndim
			psum(j) = psum(j) - p(ihi,j) + ptry(j)
			p(ihi,j) = ptry(j)
		end do
	end if
	res = ytry
end subroutine


end module
