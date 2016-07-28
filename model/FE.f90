! A dead-simple Fortran (2003) wrapper to instantiate instances of the
! Environment and FiniteElement classes, call 'init', 'step', and 'finialize',
! and then to delete the class instances again.

program main

  ! The neXtSIM module presents an interface to the C++ FiniteElement class
  use neXtSIM

  type(CXXClass) :: FE, env

  ! Local variables
  integer :: i, i_init

  ! Instantiate new instances of the classes
  call new(env, FE)

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

  ! Finalise the model
  write(*,*) "Fortran call finalise"
  call finalise(FE)
  write(*,*) "Fortran call finalise done"

  ! Delete the class instance, now that we're done
  call delete(env, FE)

  write(*,*) 'Good bye from FORTRAN'

end program main
