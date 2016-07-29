! A module file that represents an interface (Fortran 2003, standard)
! between the C++ code and a Fortran (2003) program.

module neXtSIM

  use, intrinsic :: ISO_C_Binding, only: C_int, C_ptr, C_NULL_ptr, C_char, C_loc, C_null_char
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
    function EnvironmentNew(argc, argv) result(this) bind(C,name="EnvironmentNew")
      import
      type(C_ptr) :: this
      integer(C_int), value :: argc
      type(C_ptr) :: argv(*)
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

  public :: new, delete, finalise, step, init, run, CXXClass

!========================================================================
contains
! Fortran wrapper routines to interface C wrappers
!========================================================================

  ! Instantiate a new instance of the Environment class first, then the
  ! FiniteElement class afterwards. This must always be done like this
  ! so I'm combining the two operations into one Fortran subroutine
  subroutine new(env, FE, argc, argv)
    type(CXXclass), intent(out) :: env, FE
    integer(C_int), intent(in) :: argc
    character(len=*), intent(inout), target :: argv(argc+1)
    type(C_ptr) :: argvp(argc+2)
    integer :: i
    character(len=128), target :: buffy

    do i=1,argc+1
      argv(i) = trim(argv(i))//C_null_char
      argvp(i) = C_loc(argv(i))
    end do

    ! A strange work-around for a strange bug!
	! I simply pass a larger array than the one that should be needed
	! ... it appears that boost::mpi is addressing something beyond argc
	! :^/
    buffy = ""//C_null_char
    argvp(argc+2) = C_loc(buffy(1:1))

    env%object = EnvironmentNew(argc+1, argvp)
    FE%object  = FiniteElementNew()
  end subroutine new

  ! One delete subroutine for both the Environment and FiniteElement
  ! class instances to be consistent with new
  subroutine delete(env, FE)
    type(CXXClass), intent(inout) :: env, FE
    call FiniteElementDelete(FE%object)
    call EnvironmentDelete(env%object)
  end subroutine delete

  ! To be depricated
  subroutine run(this)
    type(CXXClass), intent(in) :: this
    call FiniteElementRun(this%object)
  end subroutine run

  subroutine init(this, pcpt)
    type(CXXClass), intent(inout) :: this
    integer, intent(out) :: pcpt
    call FiniteElementInit(this%object, pcpt)
  end subroutine init

  subroutine step(this, pcpt)
    type(CXXClass), intent(inout) :: this
    integer, intent(in) :: pcpt
    call FiniteElementStep(this%object, pcpt)
  end subroutine step

  subroutine finalise(this)
    type(CXXClass), intent(inout) :: this
    call FiniteElementFinalise(this%object)
  end subroutine finalise

end module neXtSIM
