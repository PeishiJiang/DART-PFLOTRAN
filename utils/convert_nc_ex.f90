! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! Revised from convert_madis_profiler.f90 written by Nancy Colin
! $Id: convert_nc.f90 2019-09-09 15:48:00Z peishi.jiang@pnnl.gov $

program convert_nc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_nc - program that reads a netCDF file containing observation
!                          file and writes a DART obs_seq file using
!                          the DART library routines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8
use     utilities_mod, only : initialize_utilities, finalize_utilities
use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              increment_time, get_time, operator(-), GREGORIAN
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data
use          sort_mod, only : index_sort
use obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
                              create_3d_obs, getvar_int, getdimlen, getvar_real_2d, &
                              getvar_int_2d, query_varname

use           netcdf

! TODO
! How to make it in a more generic way to obtain the observation quantities in obs_kind_mod?
! Or should I generate a separate converter for each case?
use      obs_kind_mod, only : TEMPERATURE

implicit none

character(len=256),  parameter :: netcdf_file = '../data_files/obs_pflotran.nc'
character(len=256), parameter :: out_file    = '../data_files/obs_seq.pflotran'

logical, parameter :: use_input_qc           = .false.

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer  :: ncid, ntime, nloc
integer  :: n, i, oday, osec, nused, k, index
logical  :: file_exist, first_obs

real(r8) :: oerr, qc

real(r8), allocatable :: xloc(:), yloc(:), zloc(:), &
                         xlocu(:), ylocu(:), zlocu(:), &
                         tobs(:), tobsu(:)

! TODO
! Same here, the definition of observation variables should be more generic.
real(r8) :: temp_miss

! TODO
! Same here, the definition of observation variables should be more generic.
real(r8), allocatable :: temp(:,:)
integer,  allocatable :: qc_temp(:,:)

! TODO nvar should be changed according to the number of observation variables in .nc file
integer,  parameter   :: nvar = 1

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time

character(len=NF90_MAX_NAME) :: namelist(5)

!------------
! start of executable code
!------------

call initialize_utilities('convert_nc')


! TODO
! Should I change it to the start time in observation file?
! put the reference date into DART format
call set_calendar_type(GREGORIAN)
comp_day0 = set_date(2017, 4, 1, 0, 0, 0)

first_obs = .true.

ncid = nc_open_file_readonly(netcdf_file, 'convert_nc')

! Get the dimension
call getdimlen(ncid, "time", ntime)
call getdimlen(ncid, "location" , nloc)

! Allocate memory for variables in nc file
allocate(xloc(nloc)) ; allocate(xlocu(nloc*ntime))
allocate(yloc(nloc)) ; allocate(ylocu(nloc*ntime))
allocate(zloc(nloc)) ; allocate(zlocu(nloc*ntime))
allocate(tobs(ntime)); allocate(tobsu(nloc*ntime))

allocate(temp(ntime,nloc))    ;  allocate(qc_temp(ntime,nloc))

! read in the data arrays
call    getvar_real(ncid, "time",  tobs      ) ! time index
call    getvar_real(ncid, "x_location",  xloc) ! x location or easting
call    getvar_real(ncid, "y_location",  yloc) ! y location or northing
call    getvar_real(ncid, "z_location",  zloc) ! z location or latitude

call getvar_real_2d(ncid, "TEMPERATURE",  temp, temp_miss) ! temperature

! TODO
! if user says to use them, read in QCs if present
if (use_input_qc) then
   call getvar_int_2d(ncid, "TEMPERATUREQCR",   qc_temp) ! wind direction qc
else
   qc_temp = 0
endif

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=out_file, exist=file_exist)

! How to make it in a more generic way to obtain the observation quantities in obs_kind_mod?
! Or should I generate a separate converter for each case?
if ( file_exist ) then

  ! existing file found, append to it
  call read_obs_seq(out_file, 0, 0, nvar*ntime*nloc, obs_seq)

else

  ! create a new one
  call init_obs_sequence(obs_seq, num_copies, num_qc, nvar*ntime*nloc)
  do i = 1, num_copies
    call set_copy_meta_data(obs_seq, i, 'Observation')
  end do
  do i = 1, num_qc
    call set_qc_meta_data(obs_seq, i, 'Data QC')
  end do

endif

! TODO get the observation error from the nc file
oerr = 1.000_r8

! TODO Should it be provided by users, rather than a constant?
! Set the DART data quality control.  Be consistent with NCEP codes;
! 0 is 'must use', 1 is good, no reason not to use it.
qc = 1.0_r8

timeloop: do n = 1, ntime

  ! compute time of observation
  time_obs = increment_time(comp_day0, nint(tobs(n)))

  ! extract actual time of observation in file into oday, osec.
  call get_time(time_obs, osec, oday)

locloop: do k = 1, nloc
    do i = 1, nused
        ! Check for duplicate observations
      if ( xloc(k)  == xlocu(i) .and. &
           yloc(k)  == ylocu(i) .and. &
           zloc(k)  == zlocu(i) .and. &
           tobs(n)  == tobsu(i) ) cycle locloop
    end do

! TODO How to make it in a more generic way to obtain the observation quantities in obs_kind_mod?  Or should I generate a separate converter for each case?
  ! add wind component data to obs_seq
  if ( temp(n,k) /= temp_miss .and. qc_temp(n,k) == 0) then

   !if ( oerr == missing_r8)  cycle locloop

      call create_3d_obs(xloc(k), yloc(k), zloc(k), 0, temp(n,k), &
                         TEMPERATURE, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

  endif

  nused = nused + 1
  xlocu(nused) = xloc(k)
  ylocu(nused) = yloc(k)
  zlocu(nused) = zloc(k)
  tobsu(nused) = tobs(n)

end do locloop
end do timeloop

! need to wait to close file because in the loop it queries the
! report types.
call nc_close_file(ncid, 'convert_nc')

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, out_file)

! end of main program
call finalize_utilities()


end program

! <next few lines under version control, do not edit>
! $Id: convert_nc.f90 2019-09-09 15:48:00Z peishi.jiang@pnnl.gov $
! $Date: 2019-09-09 08:48:00 -0700 (Mon, 09 Sep 2019) $
