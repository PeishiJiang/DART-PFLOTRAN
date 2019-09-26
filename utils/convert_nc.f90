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
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR
use  netcdf_utilities_mod, only : nc_open_file_readonly, nc_close_file
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              increment_time, get_time, GREGORIAN, set_time, &
                              operator(-), operator(>), operator(<)
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data
use          sort_mod, only : index_sort
use obs_utilities_mod, only : getvar_real, get_or_fill_QC, add_obs_to_seq, &
                              create_3d_obs, getvar_int, getdimlen, getvar_real_2d, &
                              getvar_int_2d, query_varname

use           netcdf

use obs_kind_mod, only: TEMPERATURE

implicit none

!character(len=256),  parameter :: netcdf_file = '../data_files/obs_pflotran.nc'
!character(len=256), parameter :: out_file    = '../data_files/obs_seq.pflotran'
character(len=256) :: netcdf_file
character(len=256) :: out_file
integer :: obs_start_day,      &
           obs_start_second,   &
           obs_end_day,        &
           obs_end_second

logical, parameter :: use_input_qc           = .false.

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

integer  :: ncid, ios, ivar, ntime, nloc
integer  :: n, i, oday, osec, nused, k, index
logical  :: file_exist, first_obs

integer  :: iunit, io

real(r8) :: oerr, qc

real(r8), allocatable :: xloc(:), yloc(:), zloc(:), &
                         xlocu(:), ylocu(:), zlocu(:), &
                         tobs(:), tobsu(:)

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/model/pflotran/utils/convert_nc.f90 $"
character(len=32 ), parameter :: revision = "$Revision: Unknown $"
character(len=128), parameter :: revdate  = "$Date: 2019-09-17 08:48:00 -0700 (Tue, 17 Sept 2019) $"

character(len=512) :: string1, string2, string3

real(r8), allocatable :: temperature_val(:,:)

real(r8) :: temperature_miss

integer, allocatable :: qc_temperature(:,:)

integer, parameter :: nvar=1

namelist /convert_nc_nml/            &
    netcdf_file,                &
    out_file,                   &
    obs_start_day,              &
    obs_start_second,              &
    obs_end_day,              &
    obs_end_second

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time, start_time, end_time
character(len=256)      :: file_calendar
character(len=256)      :: unitstring
integer                 :: second, minute, hour, day, month, year

character(len=NF90_MAX_NAME) :: namelist(5)


!------------
! start of executable code
!------------

call initialize_utilities('convert_nc')

! Read in the convert_nc_nml
call find_namelist_in_file("input.nml", "convert_nc_nml", iunit)
read(iunit, nml = convert_nc_nml, iostat = io)
call check_namelist_read(iunit, io, "convert_nc_nml")

! Get the defined start and end of observation times
start_time = set_time(obs_start_second, obs_start_day)
end_time   = set_time(obs_end_second, obs_end_day)

first_obs = .true.

ncid = nc_open_file_readonly(netcdf_file, 'convert_nc')

! TODO
!call set_calendar_type(GREGORIAN)
!comp_day0 = set_date(2017, 4, 1, 0, 0, 0)
! Should I change it to the start time in observation file?
! put the reference date into DART format
ios = nf90_inq_varid(ncid, "time", ivar)
! Get the calendar type
ios = nf90_get_att(ncid, ivar, 'calendar', file_calendar)

!file_calendar = ''
if ( file_calendar == '' .or. file_calendar == 'None' ) then

   !> assumes time variable is real and fractional days.  if this isn't true,
   !> user has to supply their own read time routine.

   !comp_day0 = set_date(year, month, day, hour, minute, second)
   ios = nf90_get_att(ncid, ivar, 'units', unitstring)
   comp_day0 = set_time(0, 0)

else if ( file_calendar == 'GREGORIAN' .or. file_calendar == 'gregorian') then

   call set_calendar_type(GREGORIAN)
   ios = nf90_get_att(ncid, ivar, 'units', unitstring)

   if (unitstring(1:10) == 'days since') then

      read(unitstring,'(11x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
      if (ios /= 0) then
         write(string2,*)'expected "days since YYYY-MM-DD HH:MM:SS", got "'//trim(unitstring)//'"'
         write(string1,*)'Unable to interpret time unit attribute. Error status was ',ios
         call error_handler(E_ERR, 'convert_nc:', string1, &
                source, revision, revdate, text2=string2, text3=string3)
      endif

      ! This is the start of their calendar
      comp_day0 = set_date(year, month, day, hour, minute, second)

   else if (unitstring(1:13) == 'seconds since') then

      read(unitstring,'(14x,i4,5(1x,i2))',iostat=ios)year,month,day,hour,minute,second
      if (ios /= 0) then
         write(string1,*)'Unable to interpret time unit attribute. Error status was ',ios
         write(string2,*)'expected "seconds since YYYY-MM-DD HH:MM:SS", got "'//trim(unitstring)//'"'
         call error_handler(E_ERR, 'read_model_time:', string1, &
                source, revision, revdate, text2=string2)
      endif

      ! This is the start of their calendar
      comp_day0  = set_date(year, month, day, hour, minute, second)

   else

      write(string1, *) 'looking for "days since" or "seconds since" in the "units" attribute'
      call error_handler(E_ERR, 'read_model_time:', 'unable to set base time for gregorian calendar', &
                         source, revision, revdate, text2=string1)

   endif
else
   call error_handler(E_ERR, 'read_model_time:', &
    'calendar type "'//trim(file_calendar)//' unsupported by default read_model_time() routine', &
                      source, revision, revdate)
endif

! Get the dimension
call getdimlen(ncid, "time", ntime)
call getdimlen(ncid, "location" , nloc)

! Allocate memory for variables in nc file
allocate(xloc(nloc)) ; allocate(xlocu(nloc*ntime))
allocate(yloc(nloc)) ; allocate(ylocu(nloc*ntime))
allocate(zloc(nloc)) ; allocate(zlocu(nloc*ntime))
allocate(tobs(ntime)); allocate(tobsu(nloc*ntime))

allocate(temperature_val(ntime,nloc))
allocate(qc_temperature(ntime,nloc))

! read in the data arrays
call    getvar_real(ncid, "time",  tobs      ) ! time index
call    getvar_real(ncid, "x_location",  xloc) ! x location or easting
call    getvar_real(ncid, "y_location",  yloc) ! y location or northing
call    getvar_real(ncid, "z_location",  zloc) ! z location or latitude

call getvar_real_2d(ncid, 'TEMPERATURE',temperature_val,temperature_miss)

! Define or get the quality control value for each observation variable
if (use_input_qc) then
call getvar_int_2d(ncid, 'TEMPERATUREQCR', qc_temperature)
else
qc_temperature = 0
endif

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=out_file, exist=file_exist)

! How to make it in a more generic way to obtain the observation quantities in obs_kind_mod?
! Or should I generate a separate converter for each case?
if ( file_exist ) then

  write(string1,*)'Deleting the existing the DART observation file, with name', out_file
  open(unit=5, file=out_file, status="OLD")  ! root directory
  close(unit=5, status="DELETE")
  !call error_handler(E_ERR, 'convert_nc:', string1, source, revision, revdate)
  ! existing file found, append to it
  !call read_obs_seq(out_file, 0, 0, nvar*ntime*nloc, obs_seq)

end if

! create a new one
call init_obs_sequence(obs_seq, num_copies, num_qc, nvar*ntime*nloc)
do i = 1, num_copies
call set_copy_meta_data(obs_seq, i, 'observation')
end do
do i = 1, num_qc
call set_qc_meta_data(obs_seq, i, 'Data QC')
end do

! TODO get the observation error from the nc file
oerr = 1.000_r8

! TODO Should it be provided by users, rather than a constant?
! Set the DART data quality control.  Be consistent with NCEP codes;
! 0 is 'must use', 1 is good, no reason not to use it.
qc = 1.0_r8

timeloop: do n = 1, ntime

  ! compute time of observation
  if ( unitstring(1:3) == 'day') then
      time_obs = increment_time(comp_day0, nint(tobs(n)*86400))
  else if ( unitstring(1:6) == 'second') then
      time_obs = increment_time(comp_day0, nint(tobs(n)))
  end if

  ! extract actual time of observation in file into oday, osec.
  call get_time(time_obs, osec, oday)

  ! If time_obs is not within the obs_start and obs_end times, do not record the observation
  if ( (time_obs < start_time) .or. (time_obs > end_time)) then
      cycle
  end if
  !if ( (oday < obs_start_day) .or. &
       !(oday == obs_start_day .and. osec < obs_start_second) ) then
       !continue
  !else if ( (oday > obs_start_day) .or. &
            !(oday == obs_start_day .and. osec > obs_start_second) ) then
       !continue
  !end if

locloop: do k = 1, nloc
    do i = 1, nused
        ! Check for duplicate observations
      if ( xloc(k)  == xlocu(i) .and. &
           yloc(k)  == ylocu(i) .and. &
           zloc(k)  == zlocu(i) .and. &
           tobs(n)  == tobsu(i) ) cycle locloop
    end do

! Add each observation value here
if ( &
  temperature_val(n,k) /= temperature_miss .and. qc_temperature(n,k) == 0) then
   call create_3d_obs(xloc(k), yloc(k), zloc(k), 0, temperature_val(n,k), TEMPERATURE, oerr, oday, osec, qc, obs)
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
