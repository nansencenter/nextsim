SUBROUTINE read_grid_c (nlon,nlat,corners_ij_lus, &
                                       data_filename,dfsize, w_unit, &
                                       gridlon,gridlat, &
                                       gridclo,gridcla, &
                                       gridsrf, &
                                       indice_mask)

! Simple wrapper for read_grid.F90 so that we can pass a string (or array of
! characters) from C to Fortran

implicit none

INTEGER, INTENT(in)     :: nlon,nlat,corners_ij_lus
integer, intent(in)     :: dfsize
character(len=1), intent(in) :: data_filename(dfsize)
integer :: w_unit

DOUBLE PRECISION, DIMENSION(nlon,nlat)                :: gridlon,gridlat,gridsrf
DOUBLE PRECISION, DIMENSION(nlon,nlat,corners_ij_lus) :: gridclo,gridcla
INTEGER, DIMENSION(nlon,nlat)                      :: indice_mask

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

CALL read_grid (nlon,nlat,corners_ij_lus, &
                               data_filename_str, w_unit, &
                               gridlon,gridlat, &
                               gridclo,gridcla, &
                               gridsrf, &
                               indice_mask)

END SUBROUTINE read_grid_c
