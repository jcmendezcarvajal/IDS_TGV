Subroutine Write_SolutionRestart
  use Variables
  implicit none

  ! Create a Restart file
  open(62,file = 'SolutionRestart.dat', form = 'Unformatted',access = 'stream', status = 'replace')

  do j = 1, Jmax
    do i = 1, Imax
      write(62) r_old(i,j),u_old(i,j),v_old(i,j),T_old(i,j)
    enddo
  enddo

  ! close file
  close(62)

End Subroutine Write_SolutionRestart

Subroutine Read_SolutionRestart
  use Variables
  implicit none

  ! Reading the Restart file
  open(62,file = 'SolutionRestart.dat', form = 'Unformatted',access = 'stream')

  do j = 1, Jmax
    do i = 1, Imax
      read(62) r_old(i,j),u_old(i,j),v_old(i,j),T_old(i,j)
    enddo
  enddo

  ! close file
  close(62)

End Subroutine Read_SolutionRestart
