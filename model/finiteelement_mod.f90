! A module file that represents an interface (Fortran 2003, standard)
! between the C++ code and a Fortran (2003) program.

module neXtSIM

  use, intrinsic :: ISO_C_Binding, only: C_int, C_ptr, C_NULL_ptr
  implicit none

  private
  
!========================================================================
! A new type representing the Environment and FiniteElement classes
!========================================================================

  type CXXClass
    private
    type(C_ptr) :: object = C_NULL_ptr
  end type CXXClass

!========================================================================
  interface
!========================================================================

! Interface to instantiate a new instance of the FiniteElement class
    function FiniteElementNew() result(this) bind(C,name="FiniteElementNew")
      import
      type(C_ptr) :: this
    end function FiniteElementNew

! Interface to instantiate a new instance of the Environment class
    function EnvironmentNew() result(this) bind(C,name="EnvironmentNew")
      import
      type(C_ptr) :: this
      !integer(C_int), value :: a, b
    end function EnvironmentNew

! Interface to delete an instance of the FiniteElement class
    subroutine FiniteElementDelete(this) bind(C,name="FiniteElementDelete")
      import
      type(C_ptr), value :: this
    end subroutine FiniteElementDelete

! Interface to delete an instance of the Environment class
    subroutine EnvironmentDelete(this) bind(C,name="EnvironmentDelete")
      import
      type(C_ptr), value :: this
    end subroutine EnvironmentDelete

! Interface to call the function 'run' - to be depricated
    subroutine FiniteElementRun(this) bind(C,name="FiniteElementRun")
      import
      type(C_ptr), value :: this
    end subroutine FiniteElementRun

! Interface to call the function 'init'
    subroutine FiniteElementInit(this, pcpt) bind(C,name="FiniteElementInit")
      import
      type(C_ptr), value :: this
      integer(C_int), value :: pcpt
    end subroutine FiniteElementInit

! Interface to call the function 'step'
    subroutine FiniteElementStep(this, pcpt) bind(C,name="FiniteElementStep")
      import
      type(C_ptr), value :: this
      integer(C_int), value :: pcpt
    end subroutine FiniteElementStep

! Interface to call the function 'finalise'
    subroutine FiniteElementFinalise(this) bind(C,name="FiniteElementFinalise")
      import
      type(C_ptr), value :: this
    end subroutine FiniteElementFinalise

  end interface

!========================================================================
! Module procedures
!========================================================================

  interface new
    module procedure Env_FE_new
  end interface new

  interface delete
    module procedure Env_FE_del
  end interface delete

  interface run
    module procedure FE_run
  end interface run

  interface init
    module procedure FE_init
  end interface init

  interface step
    module procedure FE_step
  end interface step

  interface finalise
    module procedure FE_finalise
  end interface finalise

  public :: new, delete, finalise, step, init, run, CXXClass

!========================================================================
contains
! Fortran wrapper routines to interface C wrappers
!========================================================================

  ! Instantiate a new instance of the Environment class first, then the
  ! FiniteElement class afterwards. This must always be done like this
  ! so I'm combining the two operations into one Fortran subroutine
  subroutine Env_FE_new(env, FE)
    type(CXXclass), intent(out) :: env, FE
    env%object = EnvironmentNew()
    FE%ovject  = FiniteElementNew()
  end subroutine Env_FE_new

  ! One delete subroutine for both the Environment and FiniteElement
  ! class instances to be consistent with Env_FE_new
  subroutine Env_FE_del(env, FE)
    type(CXXClass), intent(inout) :: env, FE
    call FiniteElementDelete(FE%object)
    call EnvironmentDelete(env%object)
  end subroutine Env_FE_del

  ! To be depricated
  subroutine FE_run(this)
    type(CXXClass), intent(in) :: this
    call FiniteElementRun(this%object)
  end subroutine FE_run

  subroutine FE_init(this, pcpt)
    type(CXXClass), intent(inout) :: this
    integer, intent(out) :: pcpt
    call FiniteElementInit(this%object, pcpt)
  end subroutine FE_init

  subroutine FE_step(this, pcpt)
    type(CXXClass), intent(inout) :: this
    integer, intent(in) :: pcpt
    call FiniteElementStep(this%object, pcpt)
  end subroutine FE_step

  subroutine FE_finalise(this)
    type(CXXClass), intent(inout) :: this
    call FiniteElementFinalise(this%object)
  end subroutine FE_finalise

end module neXtSIM
