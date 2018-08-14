! ------------------------------------------------------------------------------
! *** int oasis3_init_comp_(int *, char *);
! ------------------------------------------------------------------------------
integer function oasis3_init_comp(compid,model_name,lname) bind(c)

    use mod_oasis
    use iso_c_binding

    implicit none
  
    ! Arguments
    integer(c_int),   intent(out) :: compid
    integer(c_int),   intent(in)  :: lname
    character(c_char),intent(in)  :: model_name(lname+1)
    
    ! Locals
    integer :: ierror, i
    character(len=lname) :: name_str

    ! Loop over the input
    do i=1,lname
      name_str(i:i) = model_name(i)
    enddo

    call oasis_init_comp(compid,name_str,ierror)
    
    oasis3_init_comp = ierror    
    return
end function oasis3_init_comp

! ------------------------------------------------------------------------------
! *** int oasis3_terminate_(void);
! ------------------------------------------------------------------------------
integer function oasis3_terminate() bind(c)

    use mod_oasis
    use iso_c_binding

    implicit none

    ! Locals
    integer :: ierror

    call oasis_terminate(ierror)

    oasis3_terminate = ierror
    return
end function oasis3_terminate

! ------------------------------------------------------------------------------
! *** int oasis3_def_partition_(int *, int *, int *);
! ------------------------------------------------------------------------------
integer function oasis3_def_partition(il_part_id,ig_paral,ig_paral_len) bind(c)

    use mod_oasis
    use iso_c_binding

    implicit none

    ! Arguments
    integer(c_int),intent(out) :: il_part_id
    integer(c_int),intent(in)  :: ig_paral_len
    integer(c_int),intent(in)  :: ig_paral(ig_paral_len)
    
    ! Locals
    integer :: ierror

    call oasis_def_partition(il_part_id,ig_paral,ierror)
    
    oasis3_def_partition = ierror
    return
end function  oasis3_def_partition

! ------------------------------------------------------------------------------
! *** int oasis3_def_var_(int *, char *, int *, int *, int *, int *, int *);
! ------------------------------------------------------------------------------
integer function oasis3_def_var(var_id,name,il_part_id,var_nodims,kinout,var_actual_shape,var_type,lname) bind(c)

    use mod_oasis
    use iso_c_binding

    implicit none
    ! Arguments
    integer(c_int),   intent(out) :: var_id
    integer(c_int),   intent(in)  :: lname
    character(c_char),intent(in)  :: name(lname+1)
    integer(c_int),   intent(in)  :: il_part_id
    integer(c_int),   intent(in)  :: var_nodims(2)
    integer(c_int),   intent(in)  :: kinout
    integer(c_int),   intent(in)  :: var_actual_shape(2*var_nodims(1))
    integer(c_int),   intent(in)  :: var_type

    ! Locals
    integer :: ierror, i
    character(len=lname) :: name_str

    ! Loop over the input
    do i=1,lname
      name_str(i:i) = name(i)
    enddo

    call oasis_def_var(var_id,name_str,il_part_id,var_nodims,kinout,var_actual_shape,var_type,ierror)

    oasis3_def_var = ierror
    return
end function oasis3_def_var

! ------------------------------------------------------------------------------
! *** int oasis3_enddef_(void);
! ------------------------------------------------------------------------------
integer function oasis3_enddef() bind(c)

    use mod_oasis
    use iso_c_binding

    implicit none

    ! Locals
    integer :: ierror

    call oasis_enddef(ierror)

    oasis3_enddef = ierror
    return
end function oasis3_enddef

! ------------------------------------------------------------------------------
! *** int oasis3_put_1d_(int *, int *, double *, int *);
! ------------------------------------------------------------------------------
integer function oasis3_put_1d(var_id,date,field_array,field_array_len) bind(c)

    use mod_oasis
    use iso_c_binding

    implicit none

    ! Arguments
    integer(c_int),            intent(in) :: var_id
    integer(c_int),            intent(in) :: date
    integer(c_int),            intent(in) :: field_array_len
    real(c_double),intent(in) :: field_array(field_array_len)

    ! Locals
    integer :: info

    call oasis_put(var_id,date,field_array,info)

    oasis3_put_1d = info
    return
end function oasis3_put_1d

! ------------------------------------------------------------------------------
! *** int oasis3_put_2d_(int *, int *, double *, int *, int *);
! ------------------------------------------------------------------------------
integer function oasis3_put_2d(var_id,date,field_array,field_array_len_x,field_array_len_y) bind(c)

    use mod_oasis
    use iso_c_binding

    implicit none

    ! Arguments
    integer(c_int),            intent(in) :: var_id
    integer(c_int),            intent(in) :: date
    integer(c_int),            intent(in) :: field_array_len_x
    integer(c_int),            intent(in) :: field_array_len_y
    real(c_double),intent(in) :: field_array(field_array_len_x,field_array_len_y)

    ! Locals
    integer :: info

    call oasis_put(var_id,date,transpose(field_array),info)

    oasis3_put_2d = info
    return
end function oasis3_put_2d

! ------------------------------------------------------------------------------
! *** int oasis3_get_1d_(int *, int *, double *, int *);
! ------------------------------------------------------------------------------
integer function oasis3_get_1d(var_id,date,field_array,field_array_len) bind(c)

    use mod_oasis
    use iso_c_binding

    implicit none

    ! Arguments
    integer(c_int),            intent(in)  :: var_id
    integer(c_int),            intent(in)  :: date
    integer(c_int),            intent(in)  :: field_array_len
    real(c_double),intent(out) :: field_array(field_array_len)

    ! Locals
    integer :: info

    call oasis_get(var_id,date,field_array,info)

    oasis3_get_1d = info
    return
end function oasis3_get_1d

! ------------------------------------------------------------------------------
! *** int oasis3_get_2d_(int *, int *, double *, int *);
! ------------------------------------------------------------------------------
integer function oasis3_get_2d(var_id,date,field_array,field_array_len_x,field_array_len_y) bind(c)

    use mod_oasis
    use iso_c_binding

    implicit none

    ! Arguments
    integer(c_int),            intent(in)  :: var_id
    integer(c_int),            intent(in)  :: date
    integer(c_int),            intent(in)  :: field_array_len_x,field_array_len_y
    real(c_double),intent(out) :: field_array(field_array_len_y,field_array_len_x)

    ! Locals
    integer :: info
    real    :: tmp_array(field_array_len_x,field_array_len_y)

    call oasis_get(var_id,date,tmp_array,info)
    field_array = transpose(tmp_array)

    oasis3_get_2d = info
    return
end function oasis3_get_2d

! ------------------------------------------------------------------------------
! *** oasis3_abort(int *, const char *, const char *, int *, int *);
! ------------------------------------------------------------------------------
integer function oasis3_abort(comp_id, routine_name, abort_message, lname, lmess) bind(c)

  use mod_oasis
  use iso_c_binding

  implicit none
  
  ! Arguments
  integer(c_int),   intent(in)  :: comp_id, lname, lmess
  character(c_char),intent(in)  :: routine_name(lname+1)
  character(c_char),intent(in)  :: abort_message(lmess+1)
  
  ! Locals
  integer :: info, i
  character(len=lname) routine_name_str
  character(len=lmess) abort_message_str

  ! Loop over the input
  do i=1,lname
    routine_name_str(i:i) = routine_name(i)
  enddo
  do i=1,lmess
    abort_message_str(i:i) = abort_message(i)
  enddo

  call oasis_abort(comp_id, routine_name_str, abort_message_str)
  
  oasis3_abort = OASIS_Ok
  return
end function oasis3_abort

! ------------------------------------------------------------------------------
! *** oasis3_get_localcomm(MPI_Comm *);
! ------------------------------------------------------------------------------
integer function oasis3_get_localcomm(localcomm) bind(c)

  use mod_oasis
  use iso_c_binding

  implicit none
  
  ! Arguments
  integer(c_int),         intent(out)  :: localcomm
  
  ! Locals
  integer :: info

  call oasis_get_localcomm(localcomm, info)
  oasis3_get_localcomm = info

  return 
end function oasis3_get_localcomm

! ------------------------------------------------------------------------------
! *** oasis3_terminate_grids_writing
! ------------------------------------------------------------------------------
integer function oasis3_start_grids_writing(flag) bind(c)

  use mod_oasis
  use iso_c_binding

  implicit none

  ! Arguments
  integer(c_int),         intent(out)  :: flag
  
  call oasis_start_grids_writing(flag)
  
  oasis3_start_grids_writing = OASIS_Ok

  return 
end function oasis3_start_grids_writing

! ------------------------------------------------------------------------------
! *** oasis3_write_grid
! ------------------------------------------------------------------------------
integer function oasis3_write_grid(gridname, nx, ny, lon,lat, lgridname) bind(c)

  use mod_oasis
  use iso_c_binding

  implicit none
  
  ! Arguments
  integer(c_int),             intent(in) :: nx, ny, lgridname
  character(c_char),          intent(in) :: gridname(lgridname+1)
  real(c_double), intent(in) :: lon(ny,nx), lat(ny,nx)
  
  ! Locals
  integer :: i
  character(len=lgridname) :: gridname_str

  ! Loop over the input
  do i=1,lgridname
    gridname_str(i:i) = gridname(i)
  enddo

  call oasis_write_grid(gridname_str, nx, ny, transpose(lon), transpose(lat))
  
  oasis3_write_grid = OASIS_Ok

  return 
end function oasis3_write_grid

! ------------------------------------------------------------------------------
! *** oasis3_write_corner
! ------------------------------------------------------------------------------
integer function oasis3_write_corner(gridname, nx, ny, nc, clon,clat, lgridname) bind(c)

  use mod_oasis
  use iso_c_binding

  implicit none
  
  ! Arguments
  integer(c_int),             intent(in) :: nx, ny, nc, lgridname
  character(c_char),          intent(in) :: gridname(lgridname+1)
  real(c_double), intent(in) :: clon(nc,ny,nx), clat(nc,ny,nx)
  
  ! Locals
  integer :: i
  character(len=lgridname) :: gridname_str

  ! Loop over the input
  do i=1,lgridname
    gridname_str(i:i) = gridname(i)
  enddo

  call oasis_write_corner(gridname_str, nx, ny, nc, reshape(clon, (/3, 2, 1/)), reshape(clat, (/3, 2, 1/)))
  
  oasis3_write_corner = OASIS_Ok

  return 
end function oasis3_write_corner



! ------------------------------------------------------------------------------
! *** oasis3_write_area
! ------------------------------------------------------------------------------
integer function oasis3_write_area(gridname, nx, ny, area, lgridname) bind(c)

  use mod_oasis
  use iso_c_binding

  implicit none
  
  ! Arguments
  integer(c_int),             intent(in) :: nx, ny, lgridname
  character(c_char),          intent(in) :: gridname(lgridname+1)
  real(c_double), intent(in) :: area(ny,nx)
  
  ! Locals
  integer :: i
  character(len=lgridname) :: gridname_str

  ! Loop over the input
  do i=1,lgridname
    gridname_str(i:i) = gridname(i)
  enddo

  call oasis_write_area(gridname_str, nx, ny, transpose(area))
  
  oasis3_write_area = OASIS_Ok

  return 
end function oasis3_write_area



! ------------------------------------------------------------------------------
! *** oasis3_write_mask
! ------------------------------------------------------------------------------
integer function oasis3_write_mask(gridname, nx, ny, mask, lgridname) bind(c)

  use mod_oasis
  use iso_c_binding

  implicit none
  
  ! Arguments
  integer(c_int),         intent(in) :: nx, ny, lgridname
  character(c_char),      intent(in) :: gridname(lgridname+1)
  integer(c_int),         intent(in) :: mask(ny,nx)
  
  ! Locals
  integer :: i
  character(len=lgridname) :: gridname_str

  ! Loop over the input
  do i=1,lgridname
    gridname_str(i:i) = gridname(i)
  enddo

  call oasis_write_mask(gridname_str, nx, ny, transpose(mask))
  
  oasis3_write_mask = OASIS_Ok

  return 
end function oasis3_write_mask


! ------------------------------------------------------------------------------
! *** oasis3_terminate_grids_writing
! ------------------------------------------------------------------------------
integer function oasis3_terminate_grids_writing() bind(c)

  use mod_oasis
  use iso_c_binding

  implicit none
  
  call oasis_terminate_grids_writing()
  
  oasis3_terminate_grids_writing = OASIS_Ok

  return 
end function oasis3_terminate_grids_writing



