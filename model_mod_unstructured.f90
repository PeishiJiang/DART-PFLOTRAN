! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod.f90 12591 2019-09-05 20:49:26Z peishi.jiang@pnnl.gov $

module model_mod

! These are the interfaces required for PFLOTRAN to be compliant
! with the DART data assimilation infrastructure. The public interfaces listed
! must all be supported with the argument lists as indicated. Many of the interfaces
! are not required for minimal implementation (see the discussion of each
! interface and look for NULL INTERFACE).


! Modules that are absolutely required for use are listed
use        types_mod, only : r8, i8, MISSING_R8, MISSING_I
use time_manager_mod, only : time_type, set_time, set_time_missing, set_date
use     location_mod, only : location_type, get_close_type, &
                             get_close_obs, get_close_state, &
                             convert_vertical_obs, convert_vertical_state, &
                             set_location, get_location, set_location_missing, &
                             LocationDims
use    utilities_mod, only : register_module, error_handler, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             E_ERR, E_MSG, &
                             find_namelist_in_file, check_namelist_read
use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode, &
                                 nc_open_file_readonly, nc_close_file, nc_check
use distributed_state_mod, only : get_state
use state_structure_mod, only : add_domain, get_index_start, get_index_end, &
                                get_dart_vector_index, get_varid_from_kind, &
                                get_model_variable_indices, add_dimension_to_variable, &
                                finished_adding_domain, state_structure_info, &
                                get_io_num_unique_dims, get_dim_lengths, get_num_dims
use ensemble_manager_mod, only : ensemble_type
use dart_time_io_mod, only  : read_model_time, write_model_time
use default_model_mod, only : pert_model_copies, nc_write_model_vars
use obs_kind_mod, only: get_index_for_quantity
use xyz_location_mod, only : xyz_location_type, xyz_set_location, xyz_get_location,         &
                             xyz_get_close_type, xyz_get_close_init, xyz_get_close_destroy, &
                             xyz_find_nearest

use netcdf

!use obs_utilities_mod, only : getvar_real, getdimlen

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
   "$URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/models/pflotran/model_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 12591 $"
character(len=128), parameter :: revdate  = "$Date: 2018-05-21 13:49:26 -0700 (Mon, 21 May 2018) $"

integer :: dom_id

character(len=512) :: string1, string2, string3

type(location_type), allocatable :: state_loc(:)  ! state locations, compute once and store for speed

type(time_type) :: time_step

! EXAMPLE: perhaps a namelist here for anything you want to/can set at runtime.
! this is optional!  only add things which can be changed at runtime.
! TODO revise the time_step_days and time_step_seconds
integer, parameter :: MAX_STATE_VARIABLES = 40
integer, parameter :: MAX_STATE_NAME_LEN  = 256
integer  :: time_step_days
integer  :: time_step_seconds

integer          :: model_size
integer          :: nvar
! TODO revise var_names
character(len=MAX_STATE_NAME_LEN) :: var_names(MAX_STATE_VARIABLES)
character(len=MAX_STATE_NAME_LEN) :: var_qtynames(MAX_STATE_VARIABLES)
integer :: qty_list(MAX_STATE_VARIABLES)

! Everything needed to describe a variable
type progvartype
   private
   character(len=MAX_STATE_NAME_LEN) :: varname
   character(len=MAX_STATE_NAME_LEN) :: dartqtyname
   integer                           :: domain
   integer                           :: dartqtyind
end type progvartype

type(progvartype), dimension(MAX_STATE_VARIABLES) :: progvar
logical :: debug = .true.  ! turn up for more and more debug messages
character(len=256) :: template_file = 'null'   ! optional; sets sizes of arrays

! Structure for computing distances to cell centers, and assorted arrays
! needed for the get_close code.
type(xyz_get_close_type)             :: cc_gc
type(xyz_location_type), allocatable :: loc_set_xyz(:)
logical  :: search_initialized = .false.
real(r8) :: maxdist = 330000.0_r8

! uncomment this, the namelist related items in the 'use utilities' section above,
! and the namelist related items below in static_init_model() to enable the
! run-time namelist settings.
! TODO revise the model namelist
namelist /model_nml/            &
   nvar,                        &
   time_step_days,              &
   time_step_seconds,           &
   debug,                       &
   var_names,                   &
   template_file,               &
   var_qtynames

! Define the grids/locations information
integer  :: nloc
real(r8), allocatable :: x_loc_all(:), y_loc_all(:), z_loc_all(:)  ! the grid locations of each dimension at all points

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
 integer  :: index_in, state_ind, i, j, k
 integer  :: iunit, io, ivar
 integer  :: ncid, dimid, varid
 integer  :: num_dims

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the model information from model_nml in input.nml
! the model namelist includes and the list of variable names
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Get the variable quantity/kind indices
do ivar = 1, nvar
    qty_list(ivar) = get_index_for_quantity(var_qtynames(ivar))
end do

! TODO: add two domains for parameters and model states separately
! Add all the variable names to the domain by using add_domain()
! tell dart the size of the model
if (template_file /= 'null') then
    ! Use add_domain_from_file() function
    dom_id = add_domain(template_file, nvar, var_names, qty_list)
else
    write(string1,*) 'Problem, the template file cannot be null.'
    call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif

! Get the grid dimensions
! Read a file about the spatial information of the model.
ncid = nc_open_file_readonly(template_file, 'static_init_model')
! get the requested dimension size
call nc_check( nf90_inq_dimid(ncid, "location", dimid), &
               'static_init_model', 'inq dimid'//trim("location"))
call nc_check( nf90_inquire_dimension(ncid, dimid, len=nloc), &
               'static_init_model', 'inquire dimension'//trim("location"))

! Obtain the one-dimensional location in each dimension
allocate(x_loc_all(nloc))
allocate(y_loc_all(nloc))
allocate(z_loc_all(nloc))
call nc_check( nf90_inq_varid(ncid, "x_location", varid), &
               'static_init_model', 'inq varid'//trim("x_location"))
call nc_check( nf90_get_var(ncid, varid, x_loc_all), &
               'static_init_model', 'inquire variable'//trim("x_location"))
call nc_check( nf90_inq_varid(ncid, "y_location", varid), &
               'static_init_model', 'inq varid'//trim("y_location"))
call nc_check( nf90_get_var(ncid, varid, y_loc_all), &
               'static_init_model', 'inquire variable'//trim("y_location"))
call nc_check( nf90_inq_varid(ncid, "z_location", varid), &
               'static_init_model', 'inq varid'//trim("z_location"))
call nc_check( nf90_get_var(ncid, varid, z_loc_all), &
               'static_init_model', 'inquire variable'//trim("z_location"))

! Create storage for locations
model_size = nvar*nloc

! Define the locations of the model state variables
! naturally, this can be done VERY differently for more complicated models.
! set_location() is different for 1D vs. 3D models, not surprisingly.
allocate(state_loc(model_size))
index_in = 1
do ivar = 1, nvar
    do i = 1, nloc
        state_loc(index_in) = set_location(x_loc_all(i),y_loc_all(i),z_loc_all(i))
        index_in = index_in + 1
    end do
end do

! Assign the variable name information locally here
do ivar = 1, nvar
    progvar(ivar)%varname     = var_names(ivar)
    progvar(ivar)%domain      = dom_id
    progvar(ivar)%dartqtyname = var_qtynames(ivar)
    progvar(ivar)%dartqtyind  = qty_list(ivar)
end do

! Close the file
call nc_close_file(ncid, 'static_init_model')

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

subroutine model_interpolate(state_handle, ens_size, location, obs_qty, expected_obs, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_qty
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

! Local storage
real(r8), dimension(LocationDims) :: loc_array
logical  :: is_find=.false.
integer  :: loc_closest_ind
real(r8) :: loc_x, loc_y, loc_z

integer  :: e

!print *, 'ensemble', ens_size

! Let's assume failure.  Set return val to missing, then the code can
! just set istatus to something indicating why it failed, and return.
! If the interpolation is good, the expected_obs will be set to the
! good value, and the last line here sets istatus to 0.
! make any error codes set here be in the 10s
expected_obs(:) = MISSING_R8
istatus(:) = 0

! Get the individual location values
loc_array = get_location(location)
loc_x     = loc_array(1)
loc_y     = loc_array(2)
loc_z     = loc_array(3)

!if ((debug) .and. do_output()) print *, 'requesting interpolation at ', loc_x,loc_y,loc_z

! Get the nearest point by using xyz_location_mod first
call find_closest_loc(loc_array, is_find, loc_closest_ind)

if (is_find) then
    expected_obs = get_val(state_handle, ens_size, loc_closest_ind, obs_qty)

else
    write(string1,*) 'Problem, not able to get the nearest point for location: ', loc_x, loc_y, loc_z
    call error_handler(E_ERR,'model_interpolate',string1,source,revision,revdate)

end if

    ! if the forward operater failed set the value to missing_r8
do e = 1, ens_size
   if (istatus(e) /= 0) then
      expected_obs(e) = MISSING_R8
   endif
enddo

end subroutine model_interpolate


!------------------------------------------------------------
! Modified from models/mpas_atm/model_mod.f90
subroutine init_closest_center()

! initialize a GC structure
! to be used later in find_closest_cell_center().

! set up a GC in the locations mod

integer :: index_in, i

allocate(loc_set_xyz(nloc))

do i=1, nloc
   loc_set_xyz(i) = xyz_set_location(x_loc_all(i), y_loc_all(i), z_loc_all(i))
enddo

call xyz_get_close_init(cc_gc, maxdist, nloc, loc_set_xyz)

end subroutine init_closest_center


!------------------------------------------------------------------
! Modified from models/mpas_atm/model_mod.f90
subroutine find_closest_loc(loc_array, is_find, closest_loc_ind)

! Determine the index for the closest center to the given point

real(r8), dimension(LocationDims), intent(in) :: loc_array
logical,  intent(inout)                       :: is_find
integer,  intent(out)                         :: closest_loc_ind

real(r8)  :: xloc, yloc, zloc
type(xyz_location_type) :: pointloc
integer :: rc

! Get the individual location values
xloc = loc_array(1)
yloc = loc_array(2)
zloc = loc_array(3)

! do this exactly once.
if (.not. search_initialized) then
   call init_closest_center()
   search_initialized = .true.
endif

pointloc = xyz_set_location(xloc, yloc, zloc)

! Search the closet location
call xyz_find_nearest(cc_gc, pointloc, loc_set_xyz, closest_loc_ind, rc)

! decide what to do if we don't find anything.
if (rc /= 0 .or. closest_loc_ind < 0) then
    print *, 'cannot find nearest cell to the location (x,y,z): ', xloc, yloc, zloc
    is_find = .false.
! if we do find something, then get the nearest location
else
    is_find = .true.
    if (debug) then
        print *, "The nearest location is: ", x_loc_all(closest_loc_ind), y_loc_all(closest_loc_ind), &
                                              z_loc_all(closest_loc_ind)
        print *, "The indices of the nearest location are: ", closest_loc_ind
    end if
endif

end subroutine find_closest_loc


!------------------------------------------------------------------
function get_val(state_handle, ens_size, loc_ind, var_kind)

type(ensemble_type), intent(in) :: state_handle
integer, intent(in) :: var_kind
integer, intent(in) :: ens_size
integer  :: loc_ind
real(r8) :: get_val(ens_size)

integer :: i

character(len = 129) :: msg_string
integer :: var_id
integer(i8) :: state_index

var_id = get_varid_from_kind(dom_id, var_kind)

! TODO: probably need to revise this once state-space formulation is implemented.
! state_index = get_dart_vector_index(loc_x_ind, loc_y_ind, loc_z_ind, dom_id, var_id)
state_index = get_dart_vector_index(loc_ind, 1, 1, dom_id, var_id)
get_val     = get_state(state_index, state_handle)

!if (debug) then
    !print *, state_index
    !print *, get_val
    !!print *, get_index_start(dom_id, var_id)
    !!print *, get_num_dim(dom_id, var_id)
    !!print *, loc_x_ind, loc_y_ind, loc_z_ind
!end if

end function get_val

!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.
! This interface is required for all applications.

! Note that this function is unused in DART-PFLOTRAN, because we run
! PFLOTRAN as an external executable.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated. If this isn't settable at runtime
! feel free to hardcode it and not add it to a namelist.
! Note that time_step is unused in DART-PFLOTRAN, because we run
! PFLOTRAN as an external executable.
time_step = set_time(time_step_seconds, time_step_days)

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
subroutine get_state_meta_data(index_in, location, var_qty)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_qty

integer  :: n
integer  :: loc_x_ind, loc_y_ind, loc_z_ind
integer(i8) :: index_start, index_end
integer  :: var_id, dom_id
character(len=256) :: varstring
character(len=256) :: kind_string

! these should be set to the actual location and state quantity
location = state_loc(index_in)

!print *, get_location(location), index_in

! Get the indices of the location and the associated variable information
!call get_model_variable_indices(index_in, loc_x_ind, loc_y_ind, loc_z_ind, var_id, dom_id, var_qty, kind_string)

! Get the variable kind/variable name
if (present(var_qty)) then

   var_qty = MISSING_I

   FINDQTY : do n = 1, nvar
      dom_id      = progvar(n)%domain
      varstring   = progvar(n)%varname
      index_start = get_index_start(dom_id,varstring)
      index_end   = get_index_end(dom_id,varstring)
      !if (debug) then
          !print *, varstring, index_in, dom_id, index_start, index_end
      !endif
      if((index_in >= index_start).and. &
         (index_in <= index_end  ) ) then
         var_qty = progvar(n)%dartqtyind
         !var_qty = varstring
         !var_qty = progvar(n)%dart_kind
         exit FINDQTY
      endif
   enddo FINDQTY

   if( var_qty == MISSING_I ) then
      write(string1,*) 'Problem, cannot find base_offset, indx is: ', index_in
      call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
   endif

endif

end subroutine get_state_meta_data



!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()

! More variables allocated in static_init_model() should be deallocated here
!TODO: Modify it later once a full unstructured grid implementation is enabled.
deallocate(state_loc)

deallocate(x_loc_all)
deallocate(y_loc_all)
deallocate(z_loc_all)

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

call nc_add_global_attribute(ncid, "model", "pflotran")

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!!------------------------------------------------------------------
!!> read the time from the input file
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
