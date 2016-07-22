! A module file that represents an interface (Fortran 2003, standard)
! between the C++ code and a Fortran (2003) program.
!
! An ESMF/CPL/etc compatable version needs to have 'init', 'step', and
! 'finialize' calls ... but this requires some rewriting of the C++ code

module CFE_module

  use, intrinsic :: ISO_C_Binding, only: C_int, C_ptr, C_NULL_ptr
  implicit none

  private
  
! A new type representing the FiniteElement class
  type CFE_type
    private
    type(C_ptr) :: object = C_NULL_ptr
  end type CFE_type

  interface
! Interface to instantiate a new class
    function C_CFE__new() result(this) bind(C,name="CFE__new")
      import
      type(C_ptr) :: this
      !integer(C_int), value :: a, b
    end function C_CFE__new

! Interface to call the function 'run'
    subroutine C_CFE__delete(this) bind(C,name="CFE__delete")
      import
      type(C_ptr), value :: this
    end subroutine C_CFE__delete

! Interface to delete an instance of the FiniteElement class
    subroutine C_CFE__run(this) bind(C,name="CFE__run")
      import
      type(C_ptr), value :: this
    end subroutine C_CFE__run

  end interface

! Module procedures
  interface new
    module procedure CFE__new
  end interface new

  interface delete
    module procedure CFE__delete
  end interface delete

  interface run
    module procedure CFE__run
  end interface run

  public :: new, delete, run, CFE_type

contains

! Fortran wrapper routines to interface C wrappers
  subroutine CFE__new(this)
    type(CFE_type), intent(out) :: this
    !integer :: a,b
    this%object = C_CFE__new()
  end subroutine CFE__new

  subroutine CFE__delete(this)
    type(CFE_type), intent(inout) :: this
    call C_CFE__delete(this%object)
    this%object = C_NULL_ptr
  end subroutine CFE__delete

  subroutine CFE__run(this)
    type(CFE_type), intent(in) :: this
    call C_CFE__run(this%object)
  end subroutine CFE__run

end module CFE_module
