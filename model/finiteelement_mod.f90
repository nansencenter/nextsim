! A module file that represents an interface (Fortran 2003, standard)
! between the C++ code and a Fortran (2003) program.

module neXtSIM

  use, intrinsic :: ISO_C_Binding, only: C_int, C_ptr, C_NULL_ptr, C_char, C_loc, C_null_char, c_f_pointer, C_double
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
      type(C_ptr), value :: argv
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
    function FiniteElementInit(this) result(pcpt) bind(C,name="FiniteElementInit")
      import
      type(C_ptr), value :: this
      integer(C_int) :: pcpt
    end function FiniteElementInit

! Interface to call the function 'step'
    subroutine FiniteElementStep(this, pcpt) bind(C,name="FiniteElementStep")
      import
      type(C_ptr), value :: this
      integer(C_int) :: pcpt
    end subroutine FiniteElementStep

! Interface to call the function 'finalise'
    subroutine FiniteElementFinalise(this) bind(C,name="FiniteElementFinalise")
      import
      type(C_ptr), value :: this
    end subroutine FiniteElementFinalise

! Interface to access variables on the grid
    function FiniteElementGetNCols(this) result(ncols) bind(C,name="FiniteElementGetNCols")
      import
      type(C_ptr), value :: this
      integer(C_int) :: ncols
    end function FiniteElementGetNCols

    function FiniteElementGetNRows(this) result(nrows) bind(C,name="FiniteElementGetNRows")
      import
      type(C_ptr), value :: this
      integer(C_int) :: nrows
    end function FiniteElementGetNRows

!    subroutine FiniteElementUpdateMoorings(this) bind(C,name="FiniteElementUpdateMoorings")
!      import
!      type(C_ptr), value :: this
!    end subroutine FiniteElementUpdateMoorings
!
!    subroutine FiniteElementResetMoorings(this) bind(C,name="FiniteElementResetMoorings")
!      import
!      type(C_ptr), value :: this
!    end subroutine FiniteElementResetMoorings
!
!    function FiniteElementGetConc(this) result(conc) bind(C,name="FiniteElementGetConc")
!      import
!      type(C_ptr), value :: this
!      type(C_ptr) :: conc
!    end function FiniteElementGetConc

  end interface

!========================================================================
! Module procedures
!========================================================================

  public :: new, delete, finalise, step, init, run, CXXClass
  !public :: getConc

!========================================================================
contains
! Fortran wrapper routines to interface C wrappers
!========================================================================

  ! Instantiate a new instance of the Environment class first, then the
  ! FiniteElement class afterwards. This must always be done like this
  ! so I'm combining the two operations into one Fortran subroutine
  ! NB!  ARGV MUST BE CONTIGUOUS ... or bad things may happen
  subroutine new(env, FE, argc, argv)
    type(CXXclass), intent(out) :: env, FE
    integer(C_int), intent(in) :: argc
    character(len=*), intent(inout), target :: argv(argc+1)
    type(C_ptr), target :: argvp(argc+1)
    type(C_ptr) :: argvdp
    integer :: i
    character(len=128), target :: buffy

    do i=1,argc+1
      argv(i) = trim(argv(i))//C_null_char
      argvp(i) = C_loc(argv(i))
    end do

    ! Must give the C-routine a double pointer!
    argvdp = C_loc(argvp(1))
    env%object = EnvironmentNew(argc+1, argvdp)
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
    pcpt = FiniteElementInit(this%object)
  end subroutine init

  subroutine step(this, pcpt)
    type(CXXClass), intent(inout) :: this
    integer, intent(inout) :: pcpt
    call FiniteElementStep(this%object, pcpt)
  end subroutine step

  subroutine finalise(this)
    type(CXXClass), intent(inout) :: this
    call FiniteElementFinalise(this%object)
  end subroutine finalise

!  subroutine getConc(this, conc)
!    type(CXXClass), intent(in) :: this
!    real(C_double), pointer, intent(out) :: conc(:,:)
!
!    integer :: ncols, nrows
!    type(C_ptr) :: conc_cptr
!
!    call FiniteElementUpdateMoorings(this%object)
!    ncols = FiniteElementGetNCols(this%object)
!    nrows = FiniteElementGetNRows(this%object)
!
!
!    conc_cptr = FiniteElementGetConc(this%object)
!    call c_f_pointer(conc_cptr, conc, [ncols, nrows])
!
!    call FiniteElementResetMoorings(this%object)
!  end subroutine getConc

end module neXtSIM
