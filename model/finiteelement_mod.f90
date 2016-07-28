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
! Interface to instantiate a new instance of the FiniteElement class
    function C_CFE__new() result(this) bind(C,name="CFE__new")
      import
      type(C_ptr) :: this
      !integer(C_int), value :: a, b
    end function C_CFE__new

! Interface to instantiate a new instance of the Environment class
    function C_Cenv__new() result(this) bind(C,name="Cenv__new")
      import
      type(C_ptr) :: this
      !integer(C_int), value :: a, b
    end function C_Cenv__new

! Interface to delete an instance of the FiniteElement class
    subroutine C_CFE__delete(this) bind(C,name="CFE__delete")
      import
      type(C_ptr), value :: this
    end subroutine C_CFE__delete

! Interface to delete an instance of the Environment class
    subroutine C_Cenv__delete(this) bind(C,name="Cenv__delete")
      import
      type(C_ptr), value :: this
    end subroutine C_Cenv__delete

! Interface to call the function 'run'
    subroutine C_CFE__run(this) bind(C,name="CFE__run")
      import
      type(C_ptr), value :: this
    end subroutine C_CFE__run

! Interface to call the function 'init'
    subroutine C_CFE__init(this, pcpt) bind(C,name="CFE__init")
      import
      type(C_ptr), value :: this
      integer(C_int), value :: pcpt
    end subroutine C_CFE__init

! Interface to call the function 'step'
    subroutine C_CFE__step(this, pcpt) bind(C,name="CFE__step")
      import
      type(C_ptr), value :: this
      integer(C_int), value :: pcpt
    end subroutine C_CFE__step

! Interface to call the function 'finalise'
    subroutine C_CFE__finalise(this) bind(C,name="CFE__finalise")
      import
      type(C_ptr), value :: this
    end subroutine C_CFE__finalise

  end interface

! Module procedures
  interface new
    module procedure CFE__new
  end interface new

  interface new_env
    module procedure Cenv__new
  end interface new_env

  interface delete
    module procedure CFE__delete
  end interface delete

  interface delete_env
    module procedure Cenv__delete
  end interface delete_env

  interface run
    module procedure CFE__run
  end interface run

  interface init
    module procedure CFE__init
  end interface init

  interface step
    module procedure CFE__step
  end interface step

  interface finalise
    module procedure CFE__finalise
  end interface finalise

  public :: new, new_env, delete, delete_env, finalise, step, init, run, CFE_type

contains

! Fortran wrapper routines to interface C wrappers
  function CFE__new() result(this)
    type(CFE_type) :: this
    this%object = C_CFE__new()
  end function CFE__new

  function Cenv__new() result(this)
    type(CFE_type) :: this
    this%object = C_Cenv__new()
  end function Cenv__new

  subroutine CFE__delete(this)
    type(CFE_type), intent(inout) :: this
    call C_CFE__delete(this%object)
    this%object = C_NULL_ptr
  end subroutine CFE__delete
  
  subroutine Cenv__delete(this)
    type(CFE_type), intent(inout) :: this
    call C_Cenv__delete(this%object)
    this%object = C_NULL_ptr
  end subroutine Cenv__delete

  subroutine CFE__run(this)
    type(CFE_type), intent(in) :: this
    call C_CFE__run(this%object)
  end subroutine CFE__run

  subroutine CFE__init(this, pcpt)
    type(CFE_type), intent(inout) :: this
    integer, intent(out) :: pcpt
    call C_CFE__init(this%object, pcpt)
  end subroutine CFE__init

  subroutine CFE__step(this, pcpt)
    type(CFE_type), intent(inout) :: this
    integer, intent(in) :: pcpt
    call C_CFE__step(this%object, pcpt)
  end subroutine CFE__step

  subroutine CFE__finalise(this)
    type(CFE_type), intent(inout) :: this
    call C_CFE__finalise(this%object)
  end subroutine CFE__finalise

end module CFE_module
