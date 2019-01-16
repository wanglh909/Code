subroutine data_folder
  use kind
  use data, only: folder, NEM, NEM_alge, x_alge

  implicit none

  integer(kind=ik):: status, l!length of the argument
  character(LEN=40):: arg, temp, name
  logical:: e

  call get_command_argument(1,arg,l,status)
  if (status.ne.0) then
     print *, 'No command line arg'
     name = 'test'
     call system( 'rm ../test/*' )
     print *, 'rm files in ../test/'
  else
     read(arg,*) name   !the changing parameter
     print *, 'Using arg:',  name
  end if
  ! call get_command_argument(2,arg,l,status)
  ! if (status.ne.0) then
  !    print *, 'No command line arg'
  !    NEM_alge = 70!150
  ! else
  !    read(arg,*) NEM_alge   !the changing parameter
  !    print *, 'Using arg:',  NEM_alge
  ! end if
  ! ! call get_command_argument(3,arg,l,status)
  ! ! if (status.ne.0) then
  ! !    print *, 'No command line arg'
  ! !    x_alge = 0.8
  ! ! else
  ! !    read(arg,*) x_alge   !the changing parameter
  ! !    print *, 'Using arg:',  x_alge
  ! ! end if

  !Begins folder name
  write(folder,'(A)')  name
  ! write(folder,'(i4)')  NEM
  ! folder = ADJUSTL(folder)   !adjustl: Removes leading whitespace
  ! !Adds a second parameter
  ! write(temp,'(i4)') NEM_alge
  ! temp = ADJUSTL(temp)
  ! folder = trim(folder)//'_'//trim(temp) !concatenate     !trim: Removes trailing whitespace
  ! ! !Adds a third parameter
  ! ! write(temp,'(f6.3)') x_alge
  ! ! temp = ADJUSTL(temp)
  ! ! folder = trim(folder)//'_'//trim(temp) !concatenate     !trim: Removes trailing whitespace

  folder = '../'//trim(folder)//'/'       !Always end folder in /

  !!CREATING FOLDER
  !see if this folder already exists
  INQUIRE(FILE=trim(folder),EXIST=e)    !Inquire: inquire information of a file
  if (.not.e) then
     print *, 'Created directory '//trim(folder)
     call system('mkdir '//trim(folder))
  else
     print *, 'Acessing directory '//trim(folder)!//'?'
     !pause
  end if

end subroutine data_folder
