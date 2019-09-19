! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod.f90 12591 2019-09-05 20:49:26Z peishi.jiang@pnnl.gov $

module model_mod

! This is a template showing the interfaces required for a model to be compliant
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
                                 nc_begin_define_mode, nc_end_define_mode
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

integer :: dom_id

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

!namelist /model_nml/ model_size, time_step_days, time_step_seconds, pflotran_variables

! Define the grids/locations information
! TODO
! For now, we assume structured cartesian coordinates
real(r8) :: x0, y0, z0  ! the lowest values
real(r8) :: dx, dy, dz  ! the grid sizes
integer  :: nx, ny, nz  ! the numbers of grids
!real(r8) :: x_set(nx), y_set(ny), z_set(nz)  ! the grid locations in each dimension
real(r8), allocatable :: x_set(:), y_set(:), z_set(:)  ! the grid locations in each dimension
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
 integer  :: index_in, state_ind, i, j, k
 integer  :: iunit, io, ivar
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

! Get the variable quantity/kind indices
do ivar = 1, nvar
    qty_list(ivar) = get_index_for_quantity(var_qtynames(ivar))
end do

print *, 'Wait here'
! Add all the variable names to the domain by using add_domain()
! tell dart the size of the model
if (template_file /= 'null') then
    ! Use add_domain_from_file() function
    dom_id = add_domain(template_file, nvar, var_names, qty_list)

else
    ! Use add_domain_from_spec() function
    dom_id = add_domain(nvar, var_names, qty_list)

    ! TODO
    ! The size of the dimension time and ensemble member should be revised later on
    ! Add the dimension to each variable
    do ivar = 1, nvar
        !call add_dimension_to_variable(dom_id, ivar, "time", 6)
        !call add_dimension_to_variable(dom_id, ivar, "member", 1)
        call add_dimension_to_variable(dom_id, ivar, "z_location", nz)
        call add_dimension_to_variable(dom_id, ivar, "y_location", ny)
        call add_dimension_to_variable(dom_id, ivar, "x_location", nx)
    end do

    call finished_adding_domain(dom_id)

    print *, 'Test...'
    print *, get_index_start(dom_id,1), get_index_end(dom_id,1)
    print *, get_index_start(dom_id,2), get_index_end(dom_id,2)
    !num_dims = get_io_num_unique_dims(dom_id)
    !print *, num_dims
    !call state_structure_info(dom_id)
endif

!print *, 'Wait here b'
! TODO
! Read a file about the spatial information of the model.
! For now, the information is read from input.nml file.
! Let's assume structured CARTESIAN coordinates for now.
! That is: x0, y0, z0, nx, ny, nz, dx, dy, dz
! There information can be read from the namelist file.
call find_namelist_in_file("input.nml", "grid_nml", iunit)
read(iunit, nml = grid_nml, iostat = io)
call check_namelist_read(iunit, io, "grid_nml")

! Obtain the one-dimensional location in each dimension
allocate(x_set(nx))
allocate(y_set(ny))
allocate(z_set(nz))
x_set(1) = x0
do i = 2,nx
    x_set(i) = x_set(i-1)+dx
end do
y_set(1) = y0
do j = 2,ny
    y_set(j) = y_set(j-1)+dy
end do
z_set(1) = z0
do k = 2,nz
    z_set(k) = z_set(k-1)+dz
end do

! TODO
! Change model_size*nx*ny*nz to a more flexible count of
! the number of locations to be assimilated
! Create storage for locations
model_size = nvar*nx*ny*nz
allocate(state_loc(model_size))

! Define the locations of the model state variables
! naturally, this can be done VERY differently for more complicated models.
! set_location() is different for 1D vs. 3D models, not surprisingly.
index_in = 1
do ivar = 1, nvar
    do k = 1, nz
        do j = 1, ny
            do i = 1,nx
                state_loc(index_in) = set_location(x0+dx*(i-1),y0+dy*(j-1),z0+dz*(k-1))
                index_in = index_in + 1
            end do
        end do
    end do
end do


! Assign the variable name information locally here
do ivar = 1, nvar
    progvar(ivar)%varname     = var_names(ivar)
    progvar(ivar)%domain      = dom_id
    progvar(ivar)%dartqtyname = var_qtynames(ivar)
    progvar(ivar)%dartqtyind  = qty_list(ivar)

    !print *, varstring, index_in, get_index_start(progvar(n)%domain, varstring), get_index_end(progvar(n)%domain, varstring)
    !print *, progvar(ivar)
end do
print *, 'Wait here c'

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
subroutine model_interpolate(state_handle, ens_size, location, obs_qty, expected_obs, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_qty
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

! Local storage
real(r8), dimension(LocationDims) :: loc_array
real(r8) :: loc_x, loc_y, loc_z
integer  :: loc_x_ind, loc_y_ind, loc_z_ind
real(r8), dimension(LocationDims) :: loc_lrt, loc_llt, loc_urt, loc_ult
real(r8), dimension(LocationDims) :: loc_lrb, loc_llb, loc_urb, loc_ulb
integer, dimension(LocationDims)  :: loc_lrt_ind, loc_llt_ind, loc_urt_ind, loc_ult_ind
integer, dimension(LocationDims)  :: loc_lrb_ind, loc_llb_ind, loc_urb_ind, loc_ulb_ind
real(r8) :: w_lrt, w_llt, w_urt, w_ult
real(r8) :: w_lrb, w_llb, w_urb, w_ulb, w_sum
real(r8) :: val(2,2,2,ens_size)
integer  :: e
integer  :: i,j,k

!print *, 'ensemble', state_handle%num_copies
print *, 'ensemble', ens_size

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

if ((debug) .and. do_output()) print *, 'requesting interpolation at ', loc_x,loc_y,loc_z

! TODO
! For now, this is applied to structured cartesian grids. Later on, it should be modified to unstructured grids.
! Note that for interpolation in unstructured grids, functions and subroutines in location_mod (e.g., get_close_obs) can be used.
! Get the eight locations that surrounds the location to be interpolated
! Eight locations: lower right (top), lower left (top), upper right (top), upper left (top)
!                  lower right (bottom), lower left (bottom), upper right (bottom), upper left (bottom)
! And their indices.
call get_the_four_corners_locations(loc_array, loc_ult, loc_urt, loc_llt, loc_lrt, &
           loc_ulb, loc_urb, loc_llb, loc_lrb, &
           loc_ult_ind, loc_urt_ind, loc_llt_ind, loc_lrt_ind, &
           loc_ulb_ind, loc_urb_ind, loc_llb_ind, loc_lrb_ind)

! Get the weights of the four locations (based on the inverse distances of
! these locations to the location to be interpolated)
w_ult = sqrt(sum((loc_ult-loc_array)**2))
w_urt = sqrt(sum((loc_urt-loc_array)**2))
w_llt = sqrt(sum((loc_llt-loc_array)**2))
w_lrt = sqrt(sum((loc_lrt-loc_array)**2))
w_ulb = sqrt(sum((loc_ulb-loc_array)**2))
w_urb = sqrt(sum((loc_urb-loc_array)**2))
w_llb = sqrt(sum((loc_llb-loc_array)**2))
w_lrb = sqrt(sum((loc_lrb-loc_array)**2))
w_sum = w_ult+w_urt+w_llt+w_lrt+w_ulb+w_urb+w_llb+w_lrb
w_ult = w_ult / w_sum
w_urt = w_urt / w_sum
w_llt = w_llt / w_sum
w_lrt = w_lrt / w_sum
w_ulb = w_ulb / w_sum
w_urb = w_urb / w_sum
w_llb = w_llb / w_sum
w_lrb = w_lrb / w_sum

! Get the values of the four locations
! Four locations: lower right, lower left, upper right, upper left
val(1, 1, 1, :) =  get_val(state_handle, ens_size, loc_ult_ind, obs_qty)
val(1, 2, 1, :) =  get_val(state_handle, ens_size, loc_urt_ind, obs_qty)
val(2, 1, 1, :) =  get_val(state_handle, ens_size, loc_llt_ind, obs_qty)
val(2, 2, 1, :) =  get_val(state_handle, ens_size, loc_lrt_ind, obs_qty)
val(1, 1, 2, :) =  get_val(state_handle, ens_size, loc_ulb_ind, obs_qty)
val(1, 2, 2, :) =  get_val(state_handle, ens_size, loc_urb_ind, obs_qty)
val(2, 1, 2, :) =  get_val(state_handle, ens_size, loc_llb_ind, obs_qty)
val(2, 2, 2, :) =  get_val(state_handle, ens_size, loc_lrb_ind, obs_qty)

if (debug) then
    print *, 'The eight locations ....'
    print *, loc_ult
    print *, loc_urt
    print *, loc_llt
    print *, loc_lrt
    print *, loc_ulb
    print *, loc_urb
    print *, loc_llb
    print *, loc_lrb
    print *, "The values at the eight locations ...."
    print *, val(1,1,1,:),val(1,2,1,:),val(2,1,1,:),val(2,2,1,:),val(1,1,2,:),val(1,2,2,:),val(2,1,2,:),val(2,2,2,:)
end if

! Conduct the interpolation based on the weighted summation of the state values at the four locations
expected_obs = w_ult * val(1,1,1,:) + w_urt * val(1,2,1,:) + &
               w_llt * val(2,1,1,:) + w_lrt * val(2,2,1,:) + &
               w_ulb * val(1,1,2,:) + w_urb * val(1,2,2,:) + &
               w_llb * val(2,1,2,:) + w_lrb * val(2,2,2,:)

! if the forward operater failed set the value to missing_r8
do e = 1, ens_size
   if (istatus(e) /= 0) then
      expected_obs(e) = MISSING_R8
   endif
enddo

end subroutine model_interpolate

!------------------------------------------------------------------
! Get the four locations that surrounds the location to be interpolated
! Four locations: lower right, lower left, upper right, upper left
! And their indices.
subroutine get_the_four_corners_locations(loc_array, loc_ult, loc_urt, loc_llt, loc_lrt, &
           loc_ulb, loc_urb, loc_llb, loc_lrb, &
           loc_ult_ind, loc_urt_ind, loc_llt_ind, loc_lrt_ind, &
           loc_ulb_ind, loc_urb_ind, loc_llb_ind, loc_lrb_ind)

real(r8), dimension(LocationDims), intent(in)  :: loc_array
real(r8), dimension(LocationDims) :: loc_lrt, loc_llt, loc_urt, loc_ult
real(r8), dimension(LocationDims) :: loc_lrb, loc_llb, loc_urb, loc_ulb
integer, dimension(LocationDims)  :: loc_lrt_ind, loc_llt_ind, loc_urt_ind, loc_ult_ind
integer, dimension(LocationDims)  :: loc_lrb_ind, loc_llb_ind, loc_urb_ind, loc_ulb_ind
!real(r8), dimension(LocationDims), intent(out) :: loc_lr, loc_ll, loc_ur, loc_ul
!integer, dimension(LocationDims), intent(out)  :: loc_lr_ind, loc_ll_ind, loc_ur_ind, loc_ul_ind

real(r8) :: loc_x, loc_y, loc_z
integer  :: i,j,k

! Get the individual location values
loc_x     = loc_array(1)
loc_y     = loc_array(2)
loc_z     = loc_array(3)

! Get the location along x dimension
do i = 1,nx
    if (x_set(i) >= loc_x) then
        ! grid value
        loc_urt(1) = x_set(i)
        loc_lrt(1) = x_set(i)
        loc_urb(1) = x_set(i)
        loc_lrb(1) = x_set(i)
        ! grid index
        loc_urt_ind(1) = i
        loc_lrt_ind(1) = i
        loc_urb_ind(1) = i
        loc_lrb_ind(1) = i
        exit
    end if
end do
if (i == 1) then
    ! grid value
    loc_ult(1) = x_set(i)
    loc_llt(1) = x_set(i)
    loc_ulb(1) = x_set(i)
    loc_llb(1) = x_set(i)
    ! grid index
    loc_ult_ind(1) = i
    loc_llt_ind(1) = i
    loc_ulb_ind(1) = i
    loc_llb_ind(1) = i
else if (i == nx+1) then
    i = i-1
    ! grid value
    loc_ult(1) = x_set(i)
    loc_llt(1) = x_set(i)
    loc_urt(1) = x_set(i)
    loc_lrt(1) = x_set(i)
    loc_ulb(1) = x_set(i)
    loc_llb(1) = x_set(i)
    loc_urb(1) = x_set(i)
    loc_lrb(1) = x_set(i)
    ! grid index
    loc_ult_ind(1) = i
    loc_llt_ind(1) = i
    loc_urt_ind(1) = i
    loc_lrt_ind(1) = i
    loc_ulb_ind(1) = i
    loc_llb_ind(1) = i
    loc_urb_ind(1) = i
    loc_lrb_ind(1) = i
else
    ! grid value
    loc_ult(1) = x_set(i-1)
    loc_llt(1) = x_set(i-1)
    loc_ulb(1) = x_set(i-1)
    loc_llb(1) = x_set(i-1)
    ! grid index
    loc_ult_ind(1) = i-1
    loc_llt_ind(1) = i-1
    loc_ulb_ind(1) = i-1
    loc_llb_ind(1) = i-1
end if
!print *, 'check', i, loc_ul
!print *, 'check', i, loc_ll
!print *, 'check', i, loc_ur
!print *, 'check', i, loc_lr

! Get the location along y dimension
do j = 1,ny
    if (y_set(j) >= loc_y) then
        ! grid value
        loc_urt(2) = y_set(j)
        loc_ult(2) = y_set(j)
        loc_urb(2) = y_set(j)
        loc_ulb(2) = y_set(j)
        ! grid index
        loc_urt_ind(2) = j
        loc_ult_ind(2) = j
        loc_urb_ind(2) = j
        loc_ulb_ind(2) = j
        exit
    end if
end do
if (j == 1) then
    ! grid value
    loc_lrt(2) = y_set(j)
    loc_llt(2) = y_set(j)
    loc_lrb(2) = y_set(j)
    loc_llb(2) = y_set(j)
    ! grid index
    loc_lrt_ind(2) = j
    loc_llt_ind(2) = j
    loc_lrb_ind(2) = j
    loc_llb_ind(2) = j
else if (j == ny+1) then
    j = j-1
    ! grid value
    loc_ult(2) = y_set(j)
    loc_llt(2) = y_set(j)
    loc_urt(2) = y_set(j)
    loc_lrt(2) = y_set(j)
    loc_ulb(2) = y_set(j)
    loc_llb(2) = y_set(j)
    loc_urb(2) = y_set(j)
    loc_lrb(2) = y_set(j)
    ! grid index
    loc_lrt_ind(2) = j
    loc_llt_ind(2) = j
    loc_urt_ind(2) = j
    loc_ult_ind(2) = j
    loc_lrb_ind(2) = j
    loc_llb_ind(2) = j
    loc_urb_ind(2) = j
    loc_ulb_ind(2) = j
else
    ! grid value
    loc_lrt(2) = y_set(j-1)
    loc_llt(2) = y_set(j-1)
    loc_lrb(2) = y_set(j-1)
    loc_llb(2) = y_set(j-1)
    ! grid index
    loc_lrt_ind(2) = j-1
    loc_llt_ind(2) = j-1
    loc_lrb_ind(2) = j-1
    loc_llb_ind(2) = j-1
end if
!print *, 'check', j, loc_ul
!print *, 'check', j, loc_ll
!print *, 'check', j, loc_ur
!print *, 'check', j, loc_lr
!print *, 'check', j, y_set(j), y_set(j-1)
!print *, 'check', y_set
!print *, 'check', x_set

! Get the location along z dimension
do k = 1,nz
    if (z_set(k) >= loc_z) then
        ! grid value
        loc_ult(3) = z_set(k)
        loc_llt(3) = z_set(k)
        loc_urt(3) = z_set(k)
        loc_lrt(3) = z_set(k)
        ! grid index
        loc_lrt_ind(3) = k
        loc_llt_ind(3) = k
        loc_urt_ind(3) = k
        loc_ult_ind(3) = k
        exit
    end if
end do
if (k == 1) then
    ! grid value
    loc_ulb(3) = z_set(k)
    loc_llb(3) = z_set(k)
    loc_urb(3) = z_set(k)
    loc_lrb(3) = z_set(k)
    ! grid index
    loc_lrb_ind(3) = k
    loc_llb_ind(3) = k
    loc_urb_ind(3) = k
    loc_ulb_ind(3) = k
else if (k == nz+1) then
    k = k-1
    ! grid value
    loc_ult(3) = z_set(k)
    loc_llt(3) = z_set(k)
    loc_urt(3) = z_set(k)
    loc_lrt(3) = z_set(k)
    loc_ulb(3) = z_set(k)
    loc_llb(3) = z_set(k)
    loc_urb(3) = z_set(k)
    loc_lrb(3) = z_set(k)
    ! grid index
    loc_lrt_ind(3) = k
    loc_llt_ind(3) = k
    loc_urt_ind(3) = k
    loc_ult_ind(3) = k
    loc_lrb_ind(3) = k
    loc_llb_ind(3) = k
    loc_urb_ind(3) = k
    loc_ulb_ind(3) = k
else
    ! grid value
    loc_ulb(3) = z_set(k-1)
    loc_llb(3) = z_set(k-1)
    loc_urb(3) = z_set(k-1)
    loc_lrb(3) = z_set(k-1)
    ! grid index
    loc_lrb_ind(3) = k-1
    loc_llb_ind(3) = k-1
    loc_urb_ind(3) = k-1
    loc_ulb_ind(3) = k-1
end if
end subroutine get_the_four_corners_locations


!------------------------------------------------------------------
function get_val(state_handle, ens_size, loc_array_ind, var_kind)

type(ensemble_type), intent(in) :: state_handle
integer, dimension(LocationDims) :: loc_array_ind
!integer, intent(in) :: lon_index, lat_index, level, var_kind
integer, intent(in) :: var_kind
integer, intent(in) :: ens_size
real(r8) :: get_val(ens_size)

integer :: loc_x_ind, loc_y_ind, loc_z_ind

character(len = 129) :: msg_string
integer :: var_id
integer(i8) :: state_index

! Get the individual location values
loc_x_ind = loc_array_ind(1)
loc_y_ind = loc_array_ind(2)
loc_z_ind = loc_array_ind(3)

var_id = get_varid_from_kind(dom_id, var_kind)

! Find the index into state array and return this value
!dom_id = progvar(var_id)%domain
state_index = get_dart_vector_index(loc_x_ind, loc_y_ind, loc_z_ind, dom_id, var_id)
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
      if (debug) then
          print *, varstring, index_in, dom_id, index_start, index_end
      endif
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
deallocate(state_loc)

deallocate(x_set)
deallocate(y_set)
deallocate(z_set)

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
