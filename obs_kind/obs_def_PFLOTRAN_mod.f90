! This is PFLOTRAN observation definition file for DART.
!
! $Id: obs_def_PFLOTRAN_mod.f90 2020-03-09 peishi.jiang@pnnl.gov $

! BEGIN DART PREPROCESS KIND LIST
!TEMPERATURE,  QTY_PFLOTRAN_TEMPERATURE
!WATER_LEVEL_SENSOR, QTY_PFLOTRAN_WATER_LEVEL
!NORMALIZED_SPC_SENSOR, QTY_PFLOTRAN_GROUNDWATER_TRACER
! END DART PREPROCESS KIND LIST


! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_PFLOTRAN_mod, only : get_expected_pflotran_ref
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE


! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!         case(TEMPERATURE)
!            call get_expected_pflotran_ref(state_handle, ens_size, location, obs_def%time, obs_def%kind, expected_obs, istatus)
!         case(NORMALIZED_SPC_SENSOR)
!            call get_expected_pflotran_ref(state_handle, ens_size, location, obs_def%time, obs_def%kind, expected_obs, istatus)
!         case(WATER_LEVEL_SENSOR)
!            call get_expected_pflotran_ref(state_handle, ens_size, location, obs_def%time, obs_def%kind, expected_obs, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF


! BEGIN DART PREPROCESS READ_OBS_DEF
!         case(TEMPERATURE)
!            continue
!         case(NORMALIZED_SPC_SENSOR)
!            continue
!         case(WATER_LEVEL_SENSOR)
!            continue
! END DART PREPROCESS READ_OBS_DEF


! BEGIN DART PREPROCESS WRITE_OBS_DEF
!         case(TEMPERATURE)
!            continue
!         case(NORMALIZED_SPC_SENSOR)
!            continue
!         case(WATER_LEVEL_SENSOR)
!            continue
! END DART PREPROCESS WRITE_OBS_DEF


! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!         case(TEMPERATURE)
!            continue
!         case(NORMALIZED_SPC_SENSOR)
!            continue
!         case(WATER_LEVEL_SENSOR)
!            continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE

module obs_def_PFLOTRAN_mod

use        types_mod, only : r8, missing_r8, RAD2DEG, DEG2RAD, PI
use    utilities_mod, only : register_module, error_handler, E_ERR, &
                             nmlfileunit, check_namelist_read,      &
                             find_namelist_in_file, do_nml_file, do_nml_term, &
                             ascii_file_format
use time_manager_mod, only : time_type, print_time
use     location_mod, only : location_type, set_location, get_location
use  assim_model_mod, only : interpolate

use     obs_kind_mod, only : QTY_PFLOTRAN_TEMPERATURE, get_quantity_for_type_of_obs

use  ensemble_manager_mod, only : ensemble_type
use obs_def_utilities_mod, only : track_status

implicit none
private

public :: get_expected_pflotran_ref

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$pflotran/obs_kind/obs_def_PFLOTRAN_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 13137 $"
character(len=128), parameter :: revdate  = "$Date: 2020-03-09 (Mon, 09 Mar 2020) $"

character(len=129) :: string1, string2
integer  :: ii
integer  :: keycount 

logical, save :: module_initialized = .false.

contains

! ---------------------------------------------------
subroutine initialize_module

    ! Handle any module initialization tasks
    
    if (module_initialized) return
    
    call register_module(source, revision, revdate)
    module_initialized = .true.
    
end subroutine initialize_module

!> Distributed version of get_expected_pflotran_def
subroutine get_expected_pflotran_ref(state_handle, ens_size, location, obs_time, obs_type, expected_obs, istatus)
!------------------------------------------------------------------------------
! Purpose: Get the ensemble state vector for an observation given a time.
!------------------------------------------------------------------------------

type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
type(location_type), intent(in)  :: location
type(time_type),     intent(in)  :: obs_time
integer,             intent(in)  :: obs_type
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

! local variables
integer :: obs_qty
logical :: return_now
integer :: this_istatus(ens_size)

if ( .not. module_initialized ) call initialize_module

istatus = 0   ! to use track_status, it must start out 0

! Obtain the variable kind/quantity given the type
obs_qty = get_quantity_for_type_of_obs(obs_type)

! print *, obs_qty, obs_type
! call print_time(obs_time)

! this calls the model_interpolate in model_mod code.
call interpolate(state_handle, ens_size, location, obs_time, obs_qty, expected_obs, this_istatus)
call track_status(ens_size, this_istatus, expected_obs, istatus, return_now)
if (return_now) return

end subroutine get_expected_pflotran_ref

end module obs_def_PFLOTRAN_mod

! END DART PREPROCESS MODULE CODE

! <next few lines under version control, do not edit>
! $Id: obs_def_PFLOTRAN_mod.f90 2020-03-09 peishi.jiang@pnnl.gov $
! $Date: 2020-03-09$
