! A dead-simple fortran (2003) wrapper to initiate a FiniteElement class, call
! the run function and then delete the class.
! 
! An ESMF/CPL/etc compatable version needs to have 'init', 'step', and
! 'finialize' calls ... but this requires some rewriting of the C++ code

program main

  ! The CFE_module module presents an interface to the C++ FiniteElement class
  use CFE_module

  type(CFE_type) :: FE, env

  ! Local variables
  integer :: i, i_init

  ! Instantiate the class
  env = new_env()
  FE = new()

  ! write(*,*) "Hello from FORTRAN"
  ! call run(FE)

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

  write(*,*) 'Good bye from FORTRAN'

  ! Delete the class instance, now that we're done
  call delete(FE)
  call delete_env(env)

end program main
