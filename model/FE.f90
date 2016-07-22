! A dead-simple fortran (2003) wrapper to initiate a FiniteElement class, call
! the run function and then delete the class.
! 
! An ESMF/CPL/etc compatable version needs to have 'init', 'step', and
! 'finialize' calls ... but this requires some rewriting of the C++ code

program main

  ! The CFE_module module presents an interface to the C++ FiniteElement class
  use CFE_module

  type(CFE_type) :: FE

  ! Instantiate the class
  call new(FE)

  ! Run the model
  call run(FE)

  write(*,*) 'Good bye form FORTRAN'

  ! Delete the class instance, now that we're done
  call delete(FE)

end program main
