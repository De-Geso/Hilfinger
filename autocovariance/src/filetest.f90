program filetest
implicit none

character(len=*), parameter :: path='data/'
character(len=:), allocatable :: prefix, suffix
character(len=256) :: fname
integer :: num, iostat, io
logical :: ex

prefix = "test_"
suffix = ".dat"
num = 0
write (fname, '(a, a, i0, a)') path, prefix, num, suffix
inquire(file=fname, exist=ex)
do while (ex)
	num = num + 1
	write (fname, '(a, a, i0, a)') path, prefix, num, suffix
	inquire(file=fname, exist=ex)
end do
open(newunit=io, file=fname, action='write')
write(io, *) "Hello world!"
write(io,*) trim(fname)


end program
