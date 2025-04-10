program read_variables_in_chunks
implicit none
character(len=256) :: filename
integer, parameter :: chunk_size = 100
real :: time_buffer(chunk_size), x1_buffer(chunk_size), x2_buffer(chunk_size)
integer :: n_read
logical :: eof

! Ensure a filename is provided
if (command_argument_count() < 1) then
	print *, "Error: Please provide a filename as an argument."
	stop
end if

! Get the filename from the command-line argument
call get_command_argument(1, filename)
print *, "Filename provided:", trim(filename)

! Open the file for binary read
open(unit=10, file=trim(filename), form='unformatted', access='stream', status='old')

eof = .false.

! Read chunks
do
	! Read the next chunk of all variables together
	read(10, iostat=n_read) time_buffer, x1_buffer, x2_buffer

	if (n_read /= 0) exit ! End of file or error

	! Process the chunk (example: print the first value of this chunk)
!	print *, "Time:", time_buffer(1)
!	print *, "x1:", x1_buffer(1)
!	print *, "x2:", x2_buffer(1)
end do

close(10)

print *, "Finished reading the file."

end program read_variables_in_chunks
