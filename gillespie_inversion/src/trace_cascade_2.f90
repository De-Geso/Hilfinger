program cascade_2
use kind_parameters
use stochastics
use randf
implicit none

character(len=*), parameter :: prefix = "trace_cascade_2D_"

real(dp), parameter :: eps=tiny(eps)

integer, parameter :: event_min = 10**4
integer, parameter :: abund_max(2) = [64, 64]

! System parameters ====================================================
real(dp), parameter :: lmbda(2) = [10._dp, 1._dp]
real(dp), parameter :: tau(2) = [1._dp, 1._dp]
real(dp), parameter :: k(2) = 10._dp
real(dp), parameter :: n(2) = 2._dp
real(dp), parameter :: c(2) = 0._dp
integer, parameter, dimension(2, 4) :: abund_update = &
	reshape((/2, 0, &
			-1, 0, &
			0, 1, &
			0, -1 &
			/), shape(abund_update))


real(dp) :: t=0._dp, tstep
real(dp) :: pcond(abund_max(1), abund_max(2))=0._dp, rate_inf(abund_max(1), abund_max(2), 4)
real(dp) :: pcond_1d(abund_max(2))=0._dp, rate_inf_1d(abund_max(2), 2)
integer :: visits(abund_max(1), abund_max(2))=0, exits(abund_max(1), abund_max(2), 4)=0
integer :: visits_1d(abund_max(2))=0, exits_1d(abund_max(2), 2)=0
real(dp) :: meanX(2)=0._dp, meanR(4)=0._dp, eta(2,2)=0._dp, cov(2,2)=0._dp, deltaX(2)
! Gillespie
real(dp) :: propensity(4), roll
integer :: x(2)=0, maxX(2)=0, event_count(4)=0, event
! Other
integer :: i, nseed
integer, allocatable :: rseed(:)

real(dp), parameter :: tdisc = 0.001_dp
real(dp) :: tlast = 0._dp
integer :: xlast(2), tcount


! Here the program begins ==============================================


call random_seed(put=seed)
! Get random seed for output in metadata
call random_seed(size=nseed)
allocate(rseed(nseed))
call random_seed(get=rseed)

xlast = x

do while (minval(event_count) < event_min)
	do i = 1, size(x)
		if (x(i) .eq. abund_max(i)) then
			write(*,"(A, I0)") "Element 1 exceeded maximum abundance ", abund_max(i)
			call exit()
		end if
	end do
	
	do i = 1, size(x)
		if (x(i) .gt. maxX(i)) maxX(i) = x(i)
	end do
	
	! Update the propensity before taking a Gillespie step
	call update_propensity(propensity, x)
	! Get timestep and event from Gillespie algorithm
	call gillespie_iter(tstep, event, propensity)
	
	! Update time by adding how long we were in the previous state
	t = t + tstep
	do i = 1, floor((t-tlast)/tdisc)
		tcount = tcount + 1
		tlast = 1._dp * tdisc * tcount
		write(10,*) tlast, x-xlast
		xlast = x
	end do
	
	! Update online mean
	deltaX = x - meanX
	meanX = meanX + tstep/t * deltaX
	
	visits(x(1)+1, x(2)+1) = visits(x(1)+1, x(2)+1) + 1
	exits(x(1)+1, x(2)+1, event) = exits(x(1)+1, x(2)+1, event) + 1
	
	if (event <= 2) then
		visits_1d(x(2)+1) = visits_1d(x(2)+1) + 1
		exits_1d(x(2)+1, event) = exits_1d(x(2)+1, event) + 1
	end if
	

	! Update state of system in preparation for next step.
	! Update probability
	pcond(x(1)+1, x(2)+1) = pcond(x(1)+1, x(2)+1) + tstep
	pcond_1d(x(2)+1) = pcond_1d(x(2)+1) + tstep
	! Update counter
	event_count(event) = event_count(event) + 1
	! Update abundances according to what happened
	x = x + abund_update(:,event)
end do

pcond = pcond / sum(pcond)
pcond_1d = pcond_1d / sum(pcond_1d)

do i = 1,4
	rate_inf(:,:,i) = exits(:,:,i) / (pcond(:,:) * t + eps)
end do

do i = 1,2
	rate_inf_1d(:,i) = exits_1d(:,i) / (pcond_1d(:) * t + eps)
end do

call dump()
call dump_1D()

write(*,*) meanX


contains


pure function production_rates(x) result(rate)
	integer, intent(in) :: x(2)
	real(dp) :: rate(2)
	
	! Linear system
	rate(1) = 1._dp * lmbda(1)
	rate(2) = 1._dp * x(1) * lmbda(2)
	
	
	! Positive hill function
	! f = 1._dp * lmbda * x**n / (x**n + k**n) + c
	
	! Negative hill function
	! f = 1._dp * lmbda * k**n / (k**n + x**n)
end function


pure function decay_rates(x) result(rate)
	integer, intent(in) :: x(2)
	real(dp) :: rate(2)
	
	rate(1) = 1._dp * x(1) / tau(1)
	rate(2) = 1._dp * x(2) / tau(2)
end function


subroutine update_propensity(prop, x)
	integer, intent(in) :: x(2)
	real(dp), intent(out) :: prop(4)
	real(dp) :: Rin(2), Rout(2)
	
	Rin = production_rates(x)
	Rout = decay_rates(x)
	
	! X1 rates
	prop(1) = Rin(1)
	prop(2) = Rout(1)
	! X2 rates
	prop(3) = Rin(2)
	prop(4) = Rout(2)
end subroutine


subroutine dump()
	character(len=32) :: filename
	integer :: i, j, fnum, io, ios
	real(dp) :: rin(2), rout(2)
	
	fnum = int(log10(real(event_min)))
	write(filename, "(A,I0,A)") prefix, fnum, ".dat"
	
	open(io, file=trim(filename), status="replace", action="write")
	
	write(io,*) "x1,x2,pcond,visits,rin1_inf,rout1_inf,rin2_inf,rout2_inf,rin1,rout1,rin2,rout2"
			
	do i = 1, maxX(1)+1
		do j = 1, maxX(2)+1
			rin = production_rates([i-1,j-1])
			rout = decay_rates([i-1,j-1])
			write(io,*) i-1, j-1, &
				pcond(i, j), &
				visits(i,j), &
				rate_inf(i,j,1), &
				rate_inf(i,j,2), &
				rate_inf(i,j,3), &
				rate_inf(i,j,4), &
				rin(1), rout(1), &
				rin(2), rout(2)
			end do
	end do
end subroutine


subroutine dump_1D()
	integer :: io, i
	
	open(io, file="cloud.dat", status="replace", action="write")
	
	write(io,*) "x1,p,rin1,rout1"
			
	do i = 1, maxX(2)+1
		write(io,*) i-1, &
			pcond_1d(i), &
			rate_inf_1d(i,1), &
			rate_inf_1d(i,2)
	end do
end subroutine

end program
