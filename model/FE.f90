! A dead-simple Fortran (2003) wrapper to instantiate instances of the
! Environment and FiniteElement classes, call 'init', 'step', and 'finialize',
! and then to delete the class instances again.

program main

  ! The neXtSIM module presents an interface to the C++ FiniteElement class
  use neXtSIM
  use, intrinsic :: ISO_C_Binding, only: C_double

  type(CXXClass) :: FE, env

  ! Local variables
  integer :: i, i_init, argc
  character(len=128), target, allocatable :: argv(:)
  real(C_double), pointer :: conc(:,:)

  ! Read the command line arguments
  ! argv is (and must be) contiguous
  argc = command_argument_count()
  allocate(argv(argc+1))
  do i=0, argc
    call get_command_argument(i, argv(i+1))
  end do

  ! Instantiate new instances of the classes
  write(*,*) "Fortran instantiate class instances"
  call new(env, FE, argc, argv)

  ! Initialise the model
  write(*,*) "Fortran call init"
  call init(FE, i_init)
  write(*,*) "Fortran call init done"

  ! Run the model
  !call run(FE)
  do i=i_init,i_init+10
    write(*,*) "Fortran call step with i=", i
    call step(FE, i)
    write(*,*) "Step completed"
  end do

  ! Get the concentration field
  call getConc(FE, conc)
  do i=1,size(conc,1)
    write(10,*) (conc(i,j), j=1,size(conc,2))
  enddo

  ! Finalise the model
  write(*,*) "Fortran call finalise"
  call finalise(FE)
  write(*,*) "Fortran call finalise done"

  ! Delete the class instance, now that we're done
  call delete(env, FE)

  write(*,*) 'Good bye from FORTRAN'

end program main
