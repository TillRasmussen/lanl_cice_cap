!--------------- CICE Ocean solo model -----------------
! This is the CICE ocean solo model component that's NUOPC compiant.
! The public ocn_register method sets up all the model services such as
! initialize, run and finalize.
!
! Author:  Fei.Liu@gmail.com
!
! 5/10/13
! This is now acting as a cap/connector between NUOPC driver and GFDL CICE code.
! Right now it's working in solo ocean mode where it does not export/import any data
!

module cice_cap_mod

  use ice_blocks, only: nx_block, ny_block
  use ice_domain_size, only: max_blocks
  use ice_flux
  use CICE_RunMod
  use CICE_InitMod
  use CICE_FinalMod 

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS      => routine_SetServices, &
    model_label_SetClock  => label_SetClock, &
    model_label_Advance   => label_Advance

  implicit none
  private
  public SetServices

  type cice_internalstate_type
  end type

  type cice_internalstate_wrapper
    type(cice_internalstate_type), pointer :: ptr
  end type

  integer   :: import_slice = 1
  integer   :: export_slice = 1

  type CICE_Field_Definition
    character(len=64)                             :: short_name
    character(len=128)                            :: long_name
    character(len=64)                             :: standard_name
    character(len=64)                             :: unit
    logical                                       :: connected
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: farrayPtr => null()
  end type CICE_Field_Definition

  type(CICE_Field_Definition) :: import_from_atmos(18)
  type(CICE_Field_Definition) :: import_from_ocean(12)
  type(CICE_Field_Definition) :: export_to_atmos(21)
  type(CICE_Field_Definition) :: export_to_ocean(8)
  type(CICE_Field_Definition) :: pass_thr_ocn_to_atm(12)

  contains
  !-----------------------------------------------------------------------
  !------------------- CICE code starts here -----------------------
  !-----------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call model_routine_SS(gcomp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP1, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP2, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_FINALIZE, &
      userRoutine=cice_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
    ! No need to change clock settings
    call ESMF_MethodAdd(gcomp, label=model_label_SetClock, &
      userRoutine=SetClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_MethodAdd(gcomp, label=model_label_Advance, &
      userRoutine=ModelAdvance_slow, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call CICE_FieldsSetup()

    call CICE_BuildFieldDictionary(import_from_atmos, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call CICE_BuildFieldDictionary(import_from_ocean, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call CICE_BuildFieldDictionary(export_to_atmos, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call CICE_BuildFieldDictionary(export_to_ocean, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  subroutine InitializeP1(gcomp, importState, exportState, clock, rc)

    type(ESMF_GridComp)                    :: gcomp
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer, intent(out)                   :: rc

    ! Local Variables
    type(ESMF_VM)                          :: vm
    integer                                :: mpi_comm

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, mpiCommunicator=mpi_comm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call CICE_Initialize(mpi_comm)

    !call CICE_AdvertiseFields(importState, import_from_atmos, rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !call CICE_AdvertiseFields(importState, import_from_ocean, rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !call CICE_AdvertiseFields(exportState, export_to_atmos, rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !call CICE_AdvertiseFields(exportState, export_to_ocean, rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out

    write(*,*) '----- CICE initialization phase 1 completed'

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP2(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local Variables
    type(ESMF_VM)                          :: vm
    type(ESMF_Grid)                        :: gridIn
    type(ESMF_Grid)                        :: gridOut
    integer                                :: npet
    
    rc = ESMF_SUCCESS

    ! We can check if npet is 4 or some other value to make sure
    ! CICE is configured to run on the correct number of processors.

    ! create a Grid object for Fields
    ! we are going to create a single tile displaced pole grid from a gridspec
    ! file. We also use the exact decomposition in CICE so that the Fields
    ! created can wrap on the data pointers in internal part of CICE
    gridIn = ESMF_GridCreate('global_gx3_gridspec.nc', ESMF_FILEFORMAT_GRIDSPEC, &
      (/2,2/), isSphere=.true., coordNames=(/'ulon', 'ulat'/), &
        indexflag=ESMF_INDEX_DELOCAL, addCornerStagger=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    gridOut = gridIn ! for now out same as in

    !call CICE_RealizeFields(importState, gridIn, import_from_atmos, "Atmos import", rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !call CICE_RealizeFields(importState, gridIn, import_from_ocean, "Ocean import", rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !call CICE_RealizeFields(exportState, gridOut, export_to_atmos, "Atmos export", rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !call CICE_RealizeFields(exportState, gridOut, export_to_ocean, "Ocean export", rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out

    write(*,*) '----- CICE initialization phase 2 completed'

  end subroutine
  
  !-----------------------------------------------------------------------------

  ! CICE model uses same clock as parent gridComp
  subroutine SetClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_TimeInterval)       :: stabilityTimeStep, timestep

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeIntervalSet(timestep, m=60, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockSet(clock, timestep=timestep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
    call ESMF_TimeIntervalSet(stabilityTimeStep, m=60, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_GridCompSetClock(gcomp, clock, stabilityTimeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance_slow(gcomp, rc)
    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc
    
    ! local variables
    type(ESMF_Clock)                       :: clock
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Time)                        :: currTime
    type(ESMF_TimeInterval)                :: timeStep

    rc = ESMF_SUCCESS
    write(*,*) 'CICE: --- run phase called --- 1'
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.
    
    call NUOPC_ClockPrintCurrTime(clock, &
      "------>Advancing CICE from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call NUOPC_TimePrint(currTime + timeStep, &
      "--------------------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


    !call NUOPC_StateWrite(importState, filePrefix='field_ice_import_', &
    !  timeslice=import_slice, rc=rc) 
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !import_slice = import_slice + 1

    write(*,*) 'CICE: --- run phase called --- 2'
    call CICE_Run
    write(*,*) 'CICE: --- run phase called --- 3'

    !call NUOPC_StateWrite(exportState, filePrefix='field_ice_export_', &
    !  timeslice=export_slice, rc=rc) 
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    !export_slice = export_slice + 1

    write(*,*) 'CICE: --- run phase called ---'

  end subroutine 

  subroutine cice_model_finalize(gcomp, importState, exportState, clock, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Time)                        :: currTime

    rc = ESMF_SUCCESS


    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call CICE_Finalize

    write(*,*) 'CICE: --- completed ---'

  end subroutine cice_model_finalize

  subroutine CICE_BuildFieldDictionary(field_defs, rc)

    type(CICE_Field_Definition), intent(inout) :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: nfields, i

    rc = ESMF_SUCCESS

    nfields = size(field_defs)

    do i = 1, nfields

      ! importable field: i-directed wind stress into ocean
      ! Available from GSM atmosphere model: YES
      ! Corresponding GSM atmosphere output field name: mean_zonal_moment_flx
      if(.not. NUOPC_FieldDictionaryHasEntry(field_defs(i)%standard_name, rc=rc)) then
        call NUOPC_FieldDictionaryAddEntry(standardName=field_defs(i)%standard_name, &
          canonicalUnits=field_defs(i)%unit, &
          defaultLongName=field_defs(i)%long_name, &
          defaultShortName=field_defs(i)%short_name, &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

    enddo

  end subroutine CICE_BuildFieldDictionary

  subroutine CICE_AdvertiseFields(state, field_defs, rc)

    type(ESMF_State), intent(inout)             :: state
    type(CICE_Field_Definition), intent(inout) :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: nfields, i

    rc = ESMF_SUCCESS

    nfields = size(field_defs)

    do i = 1, nfields

      call NUOPC_StateAdvertiseField(state, &
        standardName=field_defs(i)%standard_name, &
        longname=field_defs(i)%long_name, &
        shortname=field_defs(i)%short_name, &
        name=field_defs(i)%short_name, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    enddo

  end subroutine CICE_AdvertiseFields

  subroutine CICE_RealizeFields(state, grid, field_defs, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Grid), intent(in)                 :: grid
    type(CICE_Field_Definition), intent(inout)  :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc

    integer                                     :: nfields, i
    type(ESMF_Field)                            :: field
    integer                                     :: npet, nx, ny, pet, elb(2), eub(2), clb(2), cub(2), tlb(2), tub(2)
    type(ESMF_VM)                               :: vm
    character(len=2048)                         :: info
 
    rc = ESMF_SUCCESS

      !call ESMF_VMGetCurrent(vm, rc=rc)
      !if (rc /= ESMF_SUCCESS) call ESMF_Finalize()

      !call ESMF_VMGet(vm, petcount=npet, localPet=pet, rc=rc)
      !if (rc /= ESMF_SUCCESS) call ESMF_Finalize()

      !call ESMF_GridGet(grid, exclusiveLBound=elb, exclusiveUBound=eub, &
      !                        computationalLBound=clb, computationalUBound=cub, &
      !                        totalLBound=tlb, totalUBound=tub, rc=rc)
      !if (rc /= ESMF_SUCCESS) call ESMF_Finalize()

      !write(info, *) pet, 'exc', elb, eub, 'comp', clb, cub, 'total', tlb, tub
      !call ESMF_LogWrite(tag // " Grid "// info, &
      !  ESMF_LOGMSG_INFO, &
      !  line=__LINE__, &
      !  file=__FILE__, &
      !  rc=rc)

    nfields = size(field_defs)

    do i = 1, nfields

      write(info, *) '/', field_defs(i)%short_name, ':', &
        lbound(field_defs(i)%farrayPtr,1), ubound(field_defs(i)%farrayPtr,1), &
        lbound(field_defs(i)%farrayPtr,2), ubound(field_defs(i)%farrayPtr,2)
      call ESMF_LogWrite(tag // " Field "// field_defs(i)%standard_name // info, &
        ESMF_LOGMSG_INFO, &
        line=__LINE__, &
        file=__FILE__, &
        rc=rc)

      if(associated(field_defs(i)%farrayPtr)) then
        field = ESMF_FieldCreate(grid=grid, &
          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_DELOCAL, &
          totalLWidth=(/1,1/), totalUWidth=(/1,1/),&
          ungriddedLBound=(/1/), ungriddedUBound=(/max_blocks/), &
          name=field_defs(i)%short_name, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (NUOPC_StateIsFieldConnected(state, fieldName=field_defs(i)%short_name)) then
        call NUOPC_StateRealizeField(state, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite(tag // " Field "// field_defs(i)%standard_name // " is connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        call ESMF_LogWrite(tag // " Field "// field_defs(i)%standard_name // " is not connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
        !if(associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
        ! remove a not connected Field from State
        call ESMF_StateRemove(state, (/field_defs(i)%short_name/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

    enddo


  end subroutine CICE_RealizeFields

  subroutine CICE_FieldsSetup

    import_from_atmos(1)%short_name = 'strax'
    import_from_atmos(1)%long_name = 'wind stress components'
    import_from_atmos(1)%standard_name = 'mean_zonal_moment_flx'
    import_from_atmos(1)%unit = 'N/m^2'
    import_from_atmos(1)%connected = .false.
    import_from_atmos(1)%farrayPtr => strax


    import_from_atmos(2)%short_name = 'stray'
    import_from_atmos(2)%long_name = 'wind stress components'
    import_from_atmos(2)%standard_name = 'mean_merid_moment_flx'
    import_from_atmos(2)%unit = 'N/m^2'
    import_from_atmos(2)%connected = .false.
    import_from_atmos(2)%farrayPtr => stray


    import_from_atmos(3)%short_name = 'zlvl'
    import_from_atmos(3)%long_name = 'atm level height'
    import_from_atmos(3)%standard_name = 'atm_level_height'
    import_from_atmos(3)%unit = 'm'
    import_from_atmos(3)%connected = .false.
    import_from_atmos(3)%farrayPtr => zlvl


    import_from_atmos(4)%short_name = 'uatm'
    import_from_atmos(4)%long_name = 'wind velocity components'
    import_from_atmos(4)%standard_name = 'wind_x'
    import_from_atmos(4)%unit = 'm/s'
    import_from_atmos(4)%connected = .false.
    import_from_atmos(4)%farrayPtr => uatm


    import_from_atmos(5)%short_name = 'vatm'
    import_from_atmos(5)%long_name = 'wind velocity components'
    import_from_atmos(5)%standard_name = 'wind_y'
    import_from_atmos(5)%unit = 'm/s'
    import_from_atmos(5)%connected = .false.
    import_from_atmos(5)%farrayPtr => vatm


    import_from_atmos(6)%short_name = 'wind'
    import_from_atmos(6)%long_name = 'wind speed'
    import_from_atmos(6)%standard_name = 'wind_speed'
    import_from_atmos(6)%unit = 'm/s'
    import_from_atmos(6)%connected = .false.
    import_from_atmos(6)%farrayPtr => wind


    import_from_atmos(7)%short_name = 'potT'
    import_from_atmos(7)%long_name = 'air potential temperature'
    import_from_atmos(7)%standard_name = 'air_potT'
    import_from_atmos(7)%unit = 'K'
    import_from_atmos(7)%connected = .false.
    import_from_atmos(7)%farrayPtr => potT


    import_from_atmos(8)%short_name = 'Tair'
    import_from_atmos(8)%long_name = 'air temperature'
    import_from_atmos(8)%standard_name = 'air_T'
    import_from_atmos(8)%unit = 'K'
    import_from_atmos(8)%connected = .false.
    import_from_atmos(8)%farrayPtr => Tair


    import_from_atmos(9)%short_name = 'Qa'
    import_from_atmos(9)%long_name = 'air specific humidity'
    import_from_atmos(9)%standard_name = 'air_q'
    import_from_atmos(9)%unit = 'kg/kg'
    import_from_atmos(9)%connected = .false.
    import_from_atmos(9)%farrayPtr => Qa


    import_from_atmos(10)%short_name = 'rhoa'
    import_from_atmos(10)%long_name = 'air density'
    import_from_atmos(10)%standard_name = 'air_density'
    import_from_atmos(10)%unit = 'kg/m^3'
    import_from_atmos(10)%connected = .false.
    import_from_atmos(10)%farrayPtr => rhoa   


    import_from_atmos(11)%short_name = 'swvdr'
    import_from_atmos(11)%long_name = 'sw down, visible, direct'
    import_from_atmos(11)%standard_name = 'mean_sw_flux_vis_dir'
    import_from_atmos(11)%unit = 'W/m^2'
    import_from_atmos(11)%connected = .true.
    import_from_atmos(11)%farrayPtr => swvdr   


    import_from_atmos(12)%short_name = 'swvdf'
    import_from_atmos(12)%long_name = 'sw down, visible, diffuse'
    import_from_atmos(12)%standard_name = 'mean_sw_flux_vis_dif'
    import_from_atmos(12)%unit = 'W/m^2'
    import_from_atmos(12)%connected = .true.
    import_from_atmos(12)%farrayPtr => swvdf


    import_from_atmos(13)%short_name = 'swidr'
    import_from_atmos(13)%long_name = 'sw down, near IR, direct'
    import_from_atmos(13)%standard_name = 'mean_sw_flux_nir_dir'
    import_from_atmos(13)%unit = 'W/m^2'
    import_from_atmos(13)%connected = .true.
    import_from_atmos(13)%farrayPtr => swidr


    import_from_atmos(14)%short_name = 'swidf'
    import_from_atmos(14)%long_name = 'sw down, near IR, diffuse'
    import_from_atmos(14)%standard_name = 'mean_sw_flux_nir_dif'
    import_from_atmos(14)%unit = 'W/m^2'
    import_from_atmos(14)%connected = .true.
    import_from_atmos(14)%farrayPtr => swidf


    import_from_atmos(15)%short_name = 'flw'
    import_from_atmos(15)%long_name = 'incoming longwave radiation'
    import_from_atmos(15)%standard_name = 'mean_down_lw_flx'
    import_from_atmos(15)%unit = 'W/m^2'
    import_from_atmos(15)%connected = .true.
    import_from_atmos(15)%farrayPtr => flw


    import_from_atmos(16)%short_name = 'frain'
    import_from_atmos(16)%long_name = 'rainfall rate'
    import_from_atmos(16)%standard_name = 'mean_prec_rate'
    import_from_atmos(16)%unit = 'kg/m^2 s'
    import_from_atmos(16)%connected = .true.
    import_from_atmos(16)%farrayPtr => frain


    import_from_atmos(17)%short_name = 'fsnow'
    import_from_atmos(17)%long_name = 'snowfall rate'
    import_from_atmos(17)%standard_name = 'mean_prec_rate'
    import_from_atmos(17)%unit = 'kg/m^2 s'
    import_from_atmos(17)%connected = .true.
    import_from_atmos(17)%farrayPtr => frain


    import_from_atmos(18)%short_name = 'faero_atm'
    import_from_atmos(18)%long_name = 'aerosol deposition rate'
    import_from_atmos(18)%standard_name = 'faero_atm'
    import_from_atmos(18)%unit = 'kg/m^2 s'
    import_from_atmos(18)%connected = .false.
    !import_from_atmos(18)%farrayPtr => faero_atm


    import_from_ocean(1)%short_name = 'uocn'
    import_from_ocean(1)%long_name = 'ocean current, x-direction'
    import_from_ocean(1)%standard_name = 'ocn_current_x'
    import_from_ocean(1)%unit = 'm/s'
    import_from_ocean(1)%connected = .false.
    import_from_ocean(1)%farrayPtr => uocn


    import_from_ocean(2)%short_name = 'vocn'
    import_from_ocean(2)%long_name = 'ocean current, y-direction'
    import_from_ocean(2)%standard_name = 'ocn_current_y'
    import_from_ocean(2)%unit = 'm/s'
    import_from_ocean(2)%connected = .false.
    import_from_ocean(2)%farrayPtr => vocn


    import_from_ocean(3)%short_name = 'ss_tltx'
    import_from_ocean(3)%long_name = 'sea surface slope, x-direction'
    import_from_ocean(3)%standard_name = 'sss_x'
    import_from_ocean(3)%unit = 'm/m'
    import_from_ocean(3)%connected = .false.
    import_from_ocean(3)%farrayPtr => ss_tltx


    import_from_ocean(4)%short_name = 'ss_tlty'
    import_from_ocean(4)%long_name = 'sea surface slope, y-direction'
    import_from_ocean(4)%standard_name = 'sss_y'
    import_from_ocean(4)%unit = 'm/m'
    import_from_ocean(4)%connected = .false.
    import_from_ocean(4)%farrayPtr => ss_tlty


    import_from_ocean(5)%short_name = 'sss'
    import_from_ocean(5)%long_name = 'sea surface salinity'
    import_from_ocean(5)%standard_name = 's_surf_ppt'
    import_from_ocean(5)%unit = 'ppt'
    import_from_ocean(5)%connected = .true.
    import_from_ocean(5)%farrayPtr => sss


    import_from_ocean(6)%short_name = 'sst'
    import_from_ocean(6)%long_name = 'sea surface temperature'
    import_from_ocean(6)%standard_name = 'sea_surface_temperature'
    import_from_ocean(6)%unit = 'C'
    import_from_ocean(6)%connected = .false.
    import_from_ocean(6)%farrayPtr => sst


    import_from_ocean(7)%short_name = 'frzmlt'
    import_from_ocean(7)%long_name = 'freezing/melting potential'
    import_from_ocean(7)%standard_name = 'frz_mlt_potential'
    import_from_ocean(7)%unit = 'W/m^2'
    import_from_ocean(7)%connected = .false.
    import_from_ocean(7)%farrayPtr => frzmlt


    import_from_ocean(8)%short_name = 'frzmlt_init'
    import_from_ocean(8)%long_name = 'frzmlt used in current time step'
    import_from_ocean(8)%standard_name = 'inst_frz_mlt_potential'
    import_from_ocean(8)%unit = 'W/m^2'
    import_from_ocean(8)%connected = .false.
    import_from_ocean(8)%farrayPtr => frzmlt_init 


    import_from_ocean(9)%short_name = 'Tf'
    import_from_ocean(9)%long_name = 'freezing temperature'
    import_from_ocean(9)%standard_name = 'freezing_temperature'
    import_from_ocean(9)%unit = 'C'
    import_from_ocean(9)%connected = .false.
    import_from_ocean(9)%farrayPtr => Tf


    import_from_ocean(10)%short_name = 'qdp'
    import_from_ocean(10)%long_name = 'deep ocean heat flux, negative upward'
    import_from_ocean(10)%standard_name = 'deep_ocean_heat_flux'
    import_from_ocean(10)%unit = 'W/m^2'
    import_from_ocean(10)%connected = .false.
    import_from_ocean(10)%farrayPtr => qdp


    import_from_ocean(11)%short_name = 'hmix'
    import_from_ocean(11)%long_name = 'mixed layer depth'
    import_from_ocean(11)%standard_name = 'mixed_layer_depth'
    import_from_ocean(11)%unit = 'm'
    import_from_ocean(11)%connected = .false.
    import_from_ocean(11)%farrayPtr => hmix


    import_from_ocean(12)%short_name = 'daice_da'
    import_from_ocean(12)%long_name = 'data assimilation concentration increment rate'
    import_from_ocean(12)%standard_name = 'daice_da'
    import_from_ocean(12)%unit = 'M/s'
    import_from_ocean(12)%connected = .false.
    import_from_ocean(12)%farrayPtr => daice_da


    export_to_atmos(1)%short_name = 'strairxT'
    export_to_atmos(1)%long_name = 'stress on ice by air, x-direction'
    export_to_atmos(1)%standard_name = 'stress_on_air_ice_x'
    export_to_atmos(1)%unit = 'N/m^2'
    export_to_atmos(1)%connected = .false.
    export_to_atmos(1)%farrayPtr => strairxT


    export_to_atmos(2)%short_name = 'strairyT'
    export_to_atmos(2)%long_name = 'stress on ice by air, y-direction'
    export_to_atmos(2)%standard_name = 'stress_on_air_ice_y'
    export_to_atmos(2)%unit = 'N/m^2'
    export_to_atmos(2)%connected = .false.
    export_to_atmos(2)%farrayPtr => strairyT


    export_to_atmos(3)%short_name = 'fsens'
    export_to_atmos(3)%long_name = 'sensible heat flux'
    export_to_atmos(3)%standard_name = 'mean_sensi_heat_flx'
    export_to_atmos(3)%unit = 'W/w^2'
    export_to_atmos(3)%connected = .false.
    export_to_atmos(3)%farrayPtr => fsens


    export_to_atmos(4)%short_name = 'flat'
    export_to_atmos(4)%long_name = 'latent heat flux'
    export_to_atmos(4)%standard_name = 'mean_laten_heat_flx'
    export_to_atmos(4)%unit = 'W/w^2'
    export_to_atmos(4)%connected = .false.
    export_to_atmos(4)%farrayPtr => flat


    export_to_atmos(5)%short_name = 'fswabs'
    export_to_atmos(5)%long_name = 'shortwave flux absorbed in ice and ocean'
    export_to_atmos(5)%standard_name = 'fswabs'
    export_to_atmos(5)%unit = 'W/w^2'
    export_to_atmos(5)%connected = .false.
    export_to_atmos(5)%farrayPtr => fswabs


    export_to_atmos(6)%short_name = 'flwout'
    export_to_atmos(6)%long_name = 'outgoing longwave radiation'
    export_to_atmos(6)%standard_name = 'flwout'
    export_to_atmos(6)%unit = 'W/w^2'
    export_to_atmos(6)%connected = .false.
    export_to_atmos(6)%farrayPtr => flwout


    export_to_atmos(7)%short_name = 'Tref'
    export_to_atmos(7)%long_name = '2m atm reference temperature'
    export_to_atmos(7)%standard_name = '2m_atm_ref_temperature'
    export_to_atmos(7)%unit = 'K'
    export_to_atmos(7)%connected = .false.
    export_to_atmos(7)%farrayPtr => Tref


    export_to_atmos(8)%short_name = 'Qref'
    export_to_atmos(8)%long_name = '2m atm reference spec humidity'
    export_to_atmos(8)%standard_name = '2m_atm_ref_spec_humidity'
    export_to_atmos(8)%unit = 'K'
    export_to_atmos(8)%connected = .false.
    export_to_atmos(8)%farrayPtr => Qref


    export_to_atmos(9)%short_name = 'evap'
    export_to_atmos(9)%long_name = 'evaporative water flux'
    export_to_atmos(9)%standard_name = 'evap_water_flux'
    export_to_atmos(9)%unit = 'kg/m^2 s'
    export_to_atmos(9)%connected = .false.
    export_to_atmos(9)%farrayPtr => evap


    export_to_atmos(10)%short_name = 'alvdr'
    export_to_atmos(10)%long_name = 'visible, direct'
    export_to_atmos(10)%standard_name = 'albedo_vis_dir'
    export_to_atmos(10)%unit = ''
    export_to_atmos(10)%connected = .false.
    export_to_atmos(10)%farrayPtr => alvdr


    export_to_atmos(11)%short_name = 'alidr'
    export_to_atmos(11)%long_name = 'near-ir, direct'
    export_to_atmos(11)%standard_name = 'albedo_nir_dir'
    export_to_atmos(11)%unit = ''
    export_to_atmos(11)%connected = .false.
    export_to_atmos(11)%farrayPtr => alidr


    export_to_atmos(12)%short_name = 'alvdf'
    export_to_atmos(12)%long_name = 'visible, diffuse'
    export_to_atmos(12)%standard_name = 'albedo_vis_dif'
    export_to_atmos(12)%unit = ''
    export_to_atmos(12)%connected = .false.
    export_to_atmos(12)%farrayPtr => alvdf


    export_to_atmos(13)%short_name = 'alidf'
    export_to_atmos(13)%long_name = 'near-ir, diffuse'
    export_to_atmos(13)%standard_name = 'albedo_nir_dif'
    export_to_atmos(13)%unit = ''
    export_to_atmos(13)%connected = .false.
    export_to_atmos(13)%farrayPtr => alidf


    export_to_atmos(14)%short_name = 'alvdr_ai'
    export_to_atmos(14)%long_name = 'grid-box-mean visible, direct'
    export_to_atmos(14)%standard_name = 'mean_albedo_vis_dir'
    export_to_atmos(14)%unit = ''
    export_to_atmos(14)%connected = .false.
    export_to_atmos(14)%farrayPtr => alvdr_ai 


    export_to_atmos(15)%short_name = 'alidr_ai'
    export_to_atmos(15)%long_name = 'grid-box-mean near-ir, direct'
    export_to_atmos(15)%standard_name = 'mean_albedo_nir_dir'
    export_to_atmos(15)%unit = ''
    export_to_atmos(15)%connected = .false.
    export_to_atmos(15)%farrayPtr => alidr_ai


    export_to_atmos(16)%short_name = 'alvdf_ai'
    export_to_atmos(16)%long_name = 'grid-box-mean visible, diffuse'
    export_to_atmos(16)%standard_name = 'mean_albedo_vis_dif'
    export_to_atmos(16)%unit = ''
    export_to_atmos(16)%connected = .false.
    export_to_atmos(16)%farrayPtr => alvdf_ai


    export_to_atmos(17)%short_name = 'alidf_ai'
    export_to_atmos(17)%long_name = 'grid-box-mean near-ir, diffuse'
    export_to_atmos(17)%standard_name = 'mean_albedo_nir_dif'
    export_to_atmos(17)%unit = ''
    export_to_atmos(17)%connected = .false.
    export_to_atmos(17)%farrayPtr => alidf_ai


    export_to_atmos(18)%short_name = 'albice'
    export_to_atmos(18)%long_name = 'bare ice albedo'
    export_to_atmos(18)%standard_name = 'bare_ice_albedo'
    export_to_atmos(18)%unit = ''
    export_to_atmos(18)%connected = .false.
    export_to_atmos(18)%farrayPtr => albice


    export_to_atmos(19)%short_name = 'albsno'
    export_to_atmos(19)%long_name = 'snow albedo'
    export_to_atmos(19)%standard_name = 'snow_albedo'
    export_to_atmos(19)%unit = ''
    export_to_atmos(19)%connected = .false.
    export_to_atmos(19)%farrayPtr => albsno


    export_to_atmos(20)%short_name = 'albpnd'
    export_to_atmos(20)%long_name = 'melt pond albedo'
    export_to_atmos(20)%standard_name = 'melt_pond_albedo'
    export_to_atmos(20)%unit = ''
    export_to_atmos(20)%connected = .false.
    export_to_atmos(20)%farrayPtr => albpnd


    export_to_atmos(21)%short_name = 'apeff_ai'
    export_to_atmos(21)%long_name = 'effective pond area used for radiation calculation'
    export_to_atmos(21)%standard_name = 'apeff_ai'
    export_to_atmos(21)%unit = 'm^2'
    export_to_atmos(21)%connected = .false.
    export_to_atmos(21)%farrayPtr => apeff_ai


    export_to_ocean(1)%short_name = 'strocnxT'
    export_to_ocean(1)%long_name = 'stress on ice by ocn, x-direction'
    export_to_ocean(1)%standard_name = 'stress_on_air_ocn_x'
    export_to_ocean(1)%unit = 'N/m^2'
    export_to_ocean(1)%connected = .false.
    export_to_ocean(1)%farrayPtr => strocnxT


    export_to_ocean(2)%short_name = 'strocnyT'
    export_to_ocean(2)%long_name = 'stress on ice by ocn, y-direction'
    export_to_ocean(2)%standard_name = 'stress_on_air_ocn_y'
    export_to_ocean(2)%unit = 'N/m^2'
    export_to_ocean(2)%connected = .false.
    export_to_ocean(2)%farrayPtr => strocnyT


    export_to_ocean(3)%short_name = 'fpond'
    export_to_ocean(3)%long_name = 'fresh water flux to ponds'
    export_to_ocean(3)%standard_name = 'fresh_water_flx_to_ponds'
    export_to_ocean(3)%unit = 'kg/m^2/s'
    export_to_ocean(3)%connected = .false.
    export_to_ocean(3)%farrayPtr => fpond


    export_to_ocean(4)%short_name = 'fresh'
    export_to_ocean(4)%long_name = 'fresh water flux to ocean'
    export_to_ocean(4)%standard_name = 'fresh_water_flx_ to_ocn'
    export_to_ocean(4)%unit = 'kg/m^2/s'
    export_to_ocean(4)%connected = .false.
    export_to_ocean(4)%farrayPtr => fresh


    export_to_ocean(5)%short_name = 'fsalt'
    export_to_ocean(5)%long_name = 'salt flux to ocean'
    export_to_ocean(5)%standard_name = 'salt_flx_to_ocn'
    export_to_ocean(5)%unit = 'kg/m^2/s'
    export_to_ocean(5)%connected = .false.
    export_to_ocean(5)%farrayPtr => fsalt


    export_to_ocean(6)%short_name = 'fhocn'
    export_to_ocean(6)%long_name = 'net heat flux to ocean'
    export_to_ocean(6)%standard_name = 'net_heat_flx_to_ocn'
    export_to_ocean(6)%unit = 'W/m^2'
    export_to_ocean(6)%connected = .false.
    export_to_ocean(6)%farrayPtr => fhocn


    export_to_ocean(7)%short_name = 'fswthru'
    export_to_ocean(7)%long_name = 'shortwave penetrating to ocean'
    export_to_ocean(7)%standard_name = 'sw_pen_to_ocean'
    export_to_ocean(7)%unit = 'W/m^2'
    export_to_ocean(7)%connected = .false.
    export_to_ocean(7)%farrayPtr => fswthru


    export_to_ocean(8)%short_name = 'faero_ocn'
    export_to_ocean(8)%long_name = 'aerosol flux to ocean'
    export_to_ocean(8)%standard_name = 'faero_ocn'
    export_to_ocean(8)%unit = 'kg/m^2/s'
    export_to_ocean(8)%connected = .false.
    !export_to_ocean(8)%farrayPtr => faero_ocn


    pass_thr_ocn_to_atm(1)%short_name = 'strairx_ocn'
    pass_thr_ocn_to_atm(1)%long_name = 'stress on ocean by air, x-direction'
    pass_thr_ocn_to_atm(1)%standard_name = 'strairx_ocn'
    pass_thr_ocn_to_atm(1)%unit = 'N/m^2'
    pass_thr_ocn_to_atm(1)%connected = .false.
    pass_thr_ocn_to_atm(1)%farrayPtr => strairx_ocn 


    pass_thr_ocn_to_atm(2)%short_name = 'strairy_ocn'
    pass_thr_ocn_to_atm(2)%long_name = 'stress on ocean by air, y-direction'
    pass_thr_ocn_to_atm(2)%standard_name = 'strairy_ocn'
    pass_thr_ocn_to_atm(2)%unit = 'N/m^2'
    pass_thr_ocn_to_atm(2)%connected = .false.
    pass_thr_ocn_to_atm(2)%farrayPtr => strairy_ocn


    pass_thr_ocn_to_atm(3)%short_name = 'fsens_ocn'
    pass_thr_ocn_to_atm(3)%long_name = 'sensible heat flux'
    pass_thr_ocn_to_atm(3)%standard_name = 'mean_sensi_heat_flx'
    pass_thr_ocn_to_atm(3)%unit = 'W/m^2'
    pass_thr_ocn_to_atm(3)%connected = .false.
    pass_thr_ocn_to_atm(3)%farrayPtr => fsens_ocn


    pass_thr_ocn_to_atm(4)%short_name = 'flat_ocn'
    pass_thr_ocn_to_atm(4)%long_name = 'latent heat flux'
    pass_thr_ocn_to_atm(4)%standard_name = 'mean_laten_heat_flx'
    pass_thr_ocn_to_atm(4)%unit = 'W/m^2'
    pass_thr_ocn_to_atm(4)%connected = .false.
    pass_thr_ocn_to_atm(4)%farrayPtr => flat_ocn


    pass_thr_ocn_to_atm(5)%short_name = 'flwout_ocn'
    pass_thr_ocn_to_atm(5)%long_name = 'outgoing longwave radiation'
    pass_thr_ocn_to_atm(5)%standard_name = 'flwout_ocn'
    pass_thr_ocn_to_atm(5)%unit = 'W/m^2'
    pass_thr_ocn_to_atm(5)%connected = .false.
    pass_thr_ocn_to_atm(5)%farrayPtr => flwout_ocn


    pass_thr_ocn_to_atm(6)%short_name = 'evap_ocn'
    pass_thr_ocn_to_atm(6)%long_name = 'evaporative water flux'
    pass_thr_ocn_to_atm(6)%standard_name = 'evap_ocn'
    pass_thr_ocn_to_atm(6)%unit = 'kg/m^2/s'
    pass_thr_ocn_to_atm(6)%connected = .false.
    pass_thr_ocn_to_atm(6)%farrayPtr => evap_ocn


    pass_thr_ocn_to_atm(7)%short_name = 'alvdr_ocn'
    pass_thr_ocn_to_atm(7)%long_name = 'visible, direct'
    pass_thr_ocn_to_atm(7)%standard_name = 'albedo_vis_dir'
    pass_thr_ocn_to_atm(7)%unit = ''
    pass_thr_ocn_to_atm(7)%connected = .false.
    pass_thr_ocn_to_atm(7)%farrayPtr => alvdr_ocn


    pass_thr_ocn_to_atm(8)%short_name = 'alidr_ocn'
    pass_thr_ocn_to_atm(8)%long_name = 'near-ir, direct'
    pass_thr_ocn_to_atm(8)%standard_name = 'albedo_nir_dir'
    pass_thr_ocn_to_atm(8)%unit = ''
    pass_thr_ocn_to_atm(8)%connected = .false.
    pass_thr_ocn_to_atm(8)%farrayPtr => alidr_ocn


    pass_thr_ocn_to_atm(9)%short_name = 'alvdf_ocn'
    pass_thr_ocn_to_atm(9)%long_name = 'visible, diffuse'
    pass_thr_ocn_to_atm(9)%standard_name = 'albedo_vis_dif'
    pass_thr_ocn_to_atm(9)%unit = ''
    pass_thr_ocn_to_atm(9)%connected = .false.
    pass_thr_ocn_to_atm(9)%farrayPtr => alvdf_ocn


    pass_thr_ocn_to_atm(10)%short_name = 'alidf_ocn'
    pass_thr_ocn_to_atm(10)%long_name = 'near-ir, diffuse'
    pass_thr_ocn_to_atm(10)%standard_name = 'albedo_nir_dif'
    pass_thr_ocn_to_atm(10)%unit = ''
    pass_thr_ocn_to_atm(10)%connected = .false.
    pass_thr_ocn_to_atm(10)%farrayPtr => alidf_ocn


    pass_thr_ocn_to_atm(11)%short_name = 'Tref_ocn'
    pass_thr_ocn_to_atm(11)%long_name = '2m atm reference temperature'
    pass_thr_ocn_to_atm(11)%standard_name = '2m_atm_ref_temperature'
    pass_thr_ocn_to_atm(11)%unit = 'K'
    pass_thr_ocn_to_atm(11)%connected = .false.
    pass_thr_ocn_to_atm(11)%farrayPtr => Tref_ocn


    pass_thr_ocn_to_atm(12)%short_name = 'Qref_ocn'
    pass_thr_ocn_to_atm(12)%long_name = '2m atm reference spec humidity'
    pass_thr_ocn_to_atm(12)%standard_name = '2m_atm_ref_spec_humidity'
    pass_thr_ocn_to_atm(12)%unit = 'kg/kg'
    pass_thr_ocn_to_atm(12)%connected = .false.
    pass_thr_ocn_to_atm(12)%farrayPtr => Qref_ocn

  end subroutine CICE_FieldsSetup

end module cice_cap_mod
