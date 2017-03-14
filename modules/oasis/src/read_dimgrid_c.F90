SUBROUTINE read_dimgrid_c (nlon,nlat,data_filename,dfsize,w_unit)

! Simple wrapper for read_dimgrid.F90 so that we can pass a string (or array of
! characters) from C to Fortran

implicit none

integer, intent(out) :: nlon,nlat
integer, intent(in)  :: dfsize
character(len=1), intent(in) :: data_filename(dfsize)
integer :: w_unit

integer :: i, iend
character(len=30)  data_filename_str

! Check the size of string from C
if ( dfsize > 30 ) then
  write(*,*) "Warning: data_filename too long in read_dimgrid-interface.F90. Truncating to 30"
  iend=30
else
  iend=dfsize
endif

! Loop over the input
do i=1,iend
  data_filename_str(i:i) = data_filename(i)
enddo
! Padd with spaces
do i=iend+1,30
  data_filename_str(i:i) = " "
enddo

CALL read_dimgrid(nlon,nlat,data_filename_str,w_unit)

END SUBROUTINE read_dimgrid_c
