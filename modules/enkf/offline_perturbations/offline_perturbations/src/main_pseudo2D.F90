program main_pseudo2D
  use netcdf
  use mod_pseudo
  use m_set_random_seed
  use mod_random_forcing
  implicit none
  integer, parameter:: xdim = 1024, ydim = 1024,xy_full = xdim*ydim
  real(8) :: synforc(xy_full, 6), randfld(xy_full, 6)
  integer :: i_step, i, id, file_exist
  integer :: ncid, stat
  integer :: xy_dimid
  integer :: varid(6)
  character(150) :: FILE_NAME, fileID
  logical :: start_from_restart=.false.
  !randfld_for_restart.dat saves randfld=(/ ran%slp,ran%wndspd,ran%snowfall,ran%dwlongw,ran%sss,ran%sst /), which can be used as a restart file to continue generating the perturbation series
  FILE_NAME='randfld_for_restart.dat'    ! 
  if(start_from_restart == .true.)then
      print*, start_from_restart, 'start from ',FILE_NAME
      open(13,file=trim(FILE_NAME))
      do id = 1,xy_full
        read(13,'(6f8.4)')  randfld(id, :)
      end do
      close(13)
  else
      print*, 'start perturbation from random values'
      randfld=0. ! note: randfld is not saved but variables in synforc are saved to file.
  endif
  !
  do i_step = 0,10*4 ! create a series of perturbations for one member
      synforc=0.
      !-------------------------------
      call limits_randf(xdim, ydim)  ! read in setting from pseudo2D.nml
      call init_fvars                ! init field variables
      call init_rand_update(synforc, randfld, 0)  ! call from mod_random_forcing.F90
      ! synforc is pure output overwrite everytime
      ! randfld saves previous result,used in load_ranffld()
      ! i_step=1 indicates to create a perturbation for the previous time.
      !-------------------------------

      ! ! save synforc and randfld in .netcdf 
      write(fileID,'(I4)') i_step
      FILE_NAME='synforc_'//trim(adjustl(fileID))//'.nc'      
      ! Create the netCDF file. The nf90_clobber parameter tells netCDF to overwrite this file, if it already exists.
      stat = nf90_create(FILE_NAME, NF90_CLOBBER, ncid) 
      ! Define the coordinate dimensions
      stat = nf90_def_dim(ncid, 'xy', xy_full, xy_dimid)
      ! Define the variables. The order of the variable list below is consistent with the order in function save_randfld_synforc()
      stat = nf90_def_var(ncid, "uwind",    NF90_FLOAT, xy_dimid, varid(1)) 
      stat = nf90_def_var(ncid, "vwind",    NF90_FLOAT, xy_dimid, varid(2)) 
      stat = nf90_def_var(ncid, "snowfall", NF90_FLOAT, xy_dimid, varid(3)) 
      stat = nf90_def_var(ncid, "Qlw_in",   NF90_FLOAT, xy_dimid, varid(4)) 
      stat = nf90_def_var(ncid, "sss",      NF90_FLOAT, xy_dimid, varid(5)) 
      stat = nf90_def_var(ncid, "sst",      NF90_FLOAT, xy_dimid, varid(6)) 
      ! End define mode. This tells netCDF we are done defining metadata.
      stat = nf90_enddef(ncid) 
      ! Write the pretend data to the file. Although netCDF supports
      ! reading and writing subsets of data, in this case we write all the
      ! data in one operation.
      do i = 1, 6 
         stat = nf90_put_var(ncid, varid(i), synforc(:,i)) 
      enddo
      ! Close the file. This frees up any internal netCDF resources associated with the file, and flushes any buffers.
      stat =  nf90_close(ncid) 
   enddo

! save final randfld in .dat
    write(fileID,'(I4)') i_step
    FILE_NAME='randfld_for_restart.dat'
    open(13,file=trim(FILE_NAME))
    do id = 1,xy_full
        write(13,'(6f8.4)')  randfld(id, :)
    end do
    close(13)
end program
