! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod.f90 12591 2018-05-21 20:49:26Z nancy@ucar.edu $

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. The public interfaces listed
! must all be supported with the argument lists as indicated. Many of the interfaces
! are not required for minimal implementation (see the discussion of each
! interface and look for NULL INTERFACE).

! TODO (Something to think...)
! - Use % operator to put all the attributes under a variable, say pft
!   example: model_mod.f90 file in wrf


! Modules that are absolutely required for use are listed
use        types_mod, only : r8, i8, MISSING_R8, MISSING_I
use time_manager_mod, only : time_type, set_time, set_time_missing, set_date
use     location_mod, only : location_type, get_close_type, &
                             get_close_obs, get_close_state, &
                             convert_vertical_obs, convert_vertical_state, &
                             set_location, set_location_missing
use    utilities_mod, only : register_module, error_handler, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             E_ERR, E_MSG, &
                             find_namelist_in_file, check_namelist_read
use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode
use state_structure_mod, only : add_domain, get_index_start, get_index_end
use ensemble_manager_mod, only : ensemble_type
use dart_time_io_mod, only  : read_model_time, write_model_time
use default_model_mod, only : pert_model_copies, nc_write_model_vars
use obs_kind_mod, only: get_index_for_quantity

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          shortest_time_between_assimilations, &
          end_model,              &
          static_init_model,      &
          nc_write_model_atts,    &
          init_time,              &
          init_conditions

! public but in another module
public :: nc_write_model_vars,    &
          pert_model_copies,      &
          get_close_obs,          &
          get_close_state,        &
          convert_vertical_obs,   &
          read_model_time, &
          write_model_time, &
          convert_vertical_state


! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/models/template/model_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 12591 $"
character(len=128), parameter :: revdate  = "$Date: 2018-05-21 13:49:26 -0700 (Mon, 21 May 2018) $"

character(len=512) :: string1, string2, string3

type(location_type), allocatable :: state_loc(:)  ! state locations, compute once and store for speed

type(time_type) :: time_step

! EXAMPLE: perhaps a namelist here for anything you want to/can set at runtime.
! this is optional!  only add things which can be changed at runtime.
! TODO revise the time_step_days and time_step_seconds
integer, parameter :: MAX_STATE_VARIABLES = 40
integer, parameter :: MAX_STATE_NAME_LEN  = 256
integer  :: time_step_days      = 0
integer  :: time_step_seconds   = 3600

integer          :: model_size
! TODO revise var_names
character(len=MAX_STATE_NAME_LEN) :: var_names(MAX_STATE_VARIABLES)
character(len=MAX_STATE_NAME_LEN) :: var_qtynames(MAX_STATE_VARIABLES)

! TODO check what is NF90_MAX_NAME
! Everything needed to describe a variable
type progvartype
   private
   character(len=MAX_STATE_NAME_LEN) :: varname
   character(len=MAX_STATE_NAME_LEN) :: dartqtyname
   integer                           :: domain
   integer                           :: dartqtyind
end type progvartype

type(progvartype), dimension(MAX_STATE_VARIABLES) :: progvar

! uncomment this, the namelist related items in the 'use utilities' section above,
! and the namelist related items below in static_init_model() to enable the
! run-time namelist settings.
! TODO revise the model namelist
namelist /model_nml/            &
   model_size,                  &
   time_step_days,              &
   time_step_seconds,           &
   var_names
!namelist /model_nml/ model_size, time_step_days, time_step_seconds, pflotran_variables

! Define the grids/locations information
! For now, we assume structured cartesian coordinates
real(r8) :: x0, y0, z0  ! the lowest values
real(r8) :: dx, dy, dz  ! the grid sizes
integer  :: nx, ny, nz  ! the numbers of grids
namelist /grid_nml/  &
   x0,               &
   y0,               &
   z0,               &
   dx,               &
   dx,               &
   dy,               &
   dz,               &
   nx,               &
   ny,               &
   nz

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.
! Can be a NULL INTERFACE for the simplest models.

subroutine static_init_model()

 real(r8) :: x_loc, y_loc, z_loc
 integer  :: index_in, state_ind, i, j, k, dom_id
 integer  :: iunit, io, ivar

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the model information from model_nml in input.nml
! the model namelist includes model_size and the list of variable names
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! TODO
! Read a file about the spatial information of the model.
! For now, the information is read from input.nml file.
! Let's assume structured CARTESIAN coordinates for now.
! That is: x0, y0, z0, nx, ny, nz, dx, dy, dz
! There information can be read from the namelist file.
call find_namelist_in_file("input.nml", "grid_nml", iunit)
read(iunit, nml = grid_nml, iostat = io)
call check_namelist_read(iunit, io, "grid_nml")

! TODO
! Change model_size*nx*ny*nz to a more flexible count of
! the number of locations to be assimilated
! Create storage for locations
allocate(state_loc(model_size*nx*ny*nz))

! Define the locations of the model state variables
! naturally, this can be done VERY differently for more complicated models.
! set_location() is different for 1D vs. 3D models, not surprisingly.
index_in = 1
do state_ind = 1, model_size
    do k = 1, nz
        do j = 1, ny
            do i = 1,nx
                state_loc(index_in) = set_location(x0+dx*(i-1),y0+dy*(j-1),z0+dz*(k-1))
                index_in = index_in + 1
            end do
        end do
    end do
end do

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated. If this isn't settable at runtime
! feel free to hardcode it and not add it to a namelist.
time_step = set_time(time_step_seconds, &
                                  time_step_days)

! TODO
! Read the table of state vector/variables and DART variable quantity from obs_def_***_mod.f90 file
!call read_variable_info(iunit, model_size, var_names)

! TODO
! Add all the variable names to the domain by using add_domain()
! tell dart the size of the model
! Actually, we are using add_domain_from_spec() function
!dom_id = add_domain(int(model_size,i8))
dom_id = add_domain(model_size, var_names)

do ivar = 1, model_size
    progvar(ivar)%varname     = var_names(ivar)
    progvar(ivar)%domain      = dom_id
    progvar(ivar)%dartqtyname = var_qtynames(ivar)
    progvar(ivar)%dartqtyind  = get_index_for_quantity(var_qtynames(ivar))
end do

end subroutine static_init_model




!------------------------------------------------------------------
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

x = MISSING_R8

end subroutine init_conditions



!------------------------------------------------------------------
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to
! compute a timestep, for instance for radiation computations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a
! NULL INTERFACE.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

end subroutine adv_1step



!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer.
! This interface is required for all applications.

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size



!------------------------------------------------------------------
! Companion interface to init_conditions. Returns a time that is somehow
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no
! synthetic data experiments using perfect_model_obs are planned,
! this can be a NULL INTERFACE.

subroutine init_time(time)

type(time_type), intent(out) :: time

! for now, just set to 0
time = set_time(0,0)

end subroutine init_time




!------------------------------------------------------------------
! Given a state handle, a location, and a model state variable type,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the kind of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

! TODO
! Why do we need model_interpolate function in DA?
subroutine model_interpolate(state_handle, ens_size, location, obs_type, expected_obs, istatus)


type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

! This should be the result of the interpolation of a
! given kind (itype) of variable at the given location.
expected_obs(:) = MISSING_R8

! The return code for successful return should be 0.
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
istatus(:) = 1

end subroutine model_interpolate



!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.
! This interface is required for all applications.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

! TODO
! Revise it if the unit of time_step is not the same as desired
shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations



!------------------------------------------------------------------
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.
! TODO
subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

integer  :: n
character(len=32) :: varstring

! these should be set to the actual location and state quantity
location = state_loc(index_in)

! Get the variable type/variable name
if (present(var_type)) then

   var_type = MISSING_I

   FINDTYPE : do n = 1, model_size
      varstring = progvar(n)%varname
      if((index_in >= get_index_start(progvar(n)%domain, varstring)).and. &
         (index_in <= get_index_end(progvar(n)%domain, varstring)) ) then
         var_type = progvar(n)%dartqtyind
         !var_type = varstring
         !var_type = progvar(n)%dart_kind
         exit FINDTYPE
      endif
   enddo FINDTYPE

   if( var_type == MISSING_I ) then
      write(string1,*) 'Problem, cannot find base_offset, indx is: ', index_in
      call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
   endif

endif

end subroutine get_state_meta_data



!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()

! TODO
! More variables allocated in static_init_model() should be deallocated here
deallocate(state_loc)

end subroutine end_model


!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

! put file into define mode.

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate", revdate )

call nc_add_global_attribute(ncid, "model", "template")

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!!------------------------------------------------------------------
!!> read the time from the input file
!! TODO
!function read_model_time(filename)
!character(len=*),  intent(in) :: filename
!type(time_type)               :: read_model_time

!integer           :: year, month, day, hour, minute, second

!! Read the NetCDF file to get the time information

!! Assign the time information to read_model_time
!read_model_time = set_date(year, month, day, hour, minute, second)

!end function read_model_time

!!--------------------------------------------------------------------
!!> writes PFLOTRAN's model date and time of day into file.
!!>
!!> @param ncid         name of the file
!!> @param model_time   the current time of the model state
!!>
!! TODO
!subroutine write_model_time(ncid, dart_time)
!integer,         intent(in) :: ncid
!type(time_type), intent(in) :: dart_time

!end subroutine write_model_time

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

! <next few lines under version control, do not edit>
! $URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/models/template/model_mod.f90 $
! $Id: model_mod.f90 12591 2018-05-21 20:49:26Z nancy@ucar.edu $
! $Revision: 12591 $
! $Date: 2018-05-21 13:49:26 -0700 (Mon, 21 May 2018) $
