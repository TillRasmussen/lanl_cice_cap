!--------------- LANL CICE NUOPC CAP -----------------
! This is the LANL CICE model cap component that's NUOPC compiant.
!
! Author:  Fei.Liu@gmail.com
!
! 5/10/13
! This is now acting as a cap/connector between NUOPC driver and LANL CICE code.
!

module cice_cap_mod

  use ice_blocks, only: nx_block, ny_block, block, get_block
  use ice_domain_size, only: max_blocks, nx_global, ny_global
  use ice_domain, only: nblocks, blocks_ice
  use ice_flux
  use ice_state
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

  integer   :: import_slice = 0
  integer   :: export_slice = 0

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: transferOffer
    logical           :: assoc    ! is the farrayPtr associated with internal data
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: farrayPtr
  end type fld_list_type

  integer,parameter :: fldsMax = 100
  integer :: fldsToIce_num = 0
  type (fld_list_type) :: fldsToIce(fldsMax)
  integer :: fldsFrIce_num = 0
  type (fld_list_type) :: fldsFrIce(fldsMax)

  integer :: lsize    ! local number of gridcells for coupling
  real(ESMF_KIND_R8), pointer :: aice_cpl(:,:,:)
  character(len=256) :: tmpstr

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

    call CICE_AdvertiseFields(importState, fldsToIce_num, fldsToIce, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call CICE_AdvertiseFields(exportState, fldsFrIce_num, fldsFrIce, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

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
    type(ESMF_DistGrid)                    :: distgrid
    integer                                :: npet
    integer                                :: i,j,iblk, n
    integer                                :: ilo,ihi,jlo,jhi
    integer                                :: ig,jg,cnt
    type(block)                            :: this_block
    integer, pointer                       :: indexList(:)
    real(ESMF_KIND_R8), pointer            :: tarray(:,:)     


    rc = ESMF_SUCCESS

    ! We can check if npet is 4 or some other value to make sure
    ! CICE is configured to run on the correct number of processors.

    ! create a Grid object for Fields
    ! we are going to create a single tile displaced pole grid from a gridspec
    ! file. We also use the exact decomposition in CICE so that the Fields
    ! created can wrap on the data pointers in internal part of CICE

    write(tmpstr,'(a,2i8)') 'InitializeP2 ice nx,ny = ',nx_global,ny_global
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/nx_global,ny_global/), &
       regDecomp=(/2,2/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_DistGridPrint(distgrid=distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, elementCount=cnt, rc=rc)
    allocate(indexList(cnt))
    write(tmpstr,'(a,i8)') 'InitializeP2 distgrid cnt= ',cnt
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, seqIndexList=indexList, rc=rc)
    write(tmpstr,'(a,4i8)') 'InitializeP2 distgrid list= ',indexList(1),indexList(cnt),minval(indexList), maxval(indexList)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
    deallocate(IndexList)

    gridIn = ESMF_GridCreate('global_gx3_gridspec.nc', ESMF_FILEFORMAT_GRIDSPEC, &
!      (/2,2/), isSphere=.true., coordNames=(/'ulon', 'ulat'/), &
      distgrid=distgrid, isSphere=.true., coordNames=(/'ulon', 'ulat'/), &
      indexflag=ESMF_INDEX_DELOCAL, addCornerStagger=.true., rc=rc)
!      indexflag=ESMF_INDEX_GLOBAL, addCornerStagger=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridGetCoord(gridIn, coordDim=1, localDE=0,  &
       staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=tarray, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(tmpstr,'(a,2g15.7)') 'initializeP2 gridIn center1 = ',minval(tarray),maxval(tarray)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_GridGetCoord(gridIn, coordDim=2, localDE=0,  &
       staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=tarray, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(tmpstr,'(a,2g15.7)') 'initializeP2 gridIn center2 = ',minval(tarray),maxval(tarray)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    gridOut = gridIn ! for now out same as in

    call CICE_RealizeFields(importState, gridIn , fldsToIce_num, fldsToIce, "Ice import", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call CICE_RealizeFields(exportState, gridOut, fldsFrIce_num, fldsFrIce, "Ice export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

! Have to be careful with reset since states are pointing directly into cice arrays
!    call state_reset(ImportState, value=-99._ESMF_KIND_R8, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!    call state_reset(ExportState, value=-99._ESMF_KIND_R8, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

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
    type(block)                            :: this_block
    integer                                :: i,j,iblk,n,i2,j2
    integer                                :: ilo,ihi,jlo,jhi

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

    call state_diagnose(importState, 'cice_import', rc)

    import_slice = import_slice + 1
    call NUOPC_StateWrite(importState, filePrefix='field_ice_import_', &
      timeslice=import_slice, relaxedFlag=.true., rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(*,*) 'CICE: --- run phase called --- 2'
    call CICE_Run
    write(*,*) 'CICE: --- run phase called --- 3'

    !---- local modifications to coupling fields -----
    !---- good place to rotate vectors           -----

    aice_cpl = 0.
    do iblk = 1,nblocks
       do j = 1,ny_block
       do i = 1,nx_block
          aice_cpl(i,j,iblk) = aice(i,j,iblk)
!          write(tmpstr,'(a,3i6,2x,g17.7)') 'tcx aice = ',i,j,iblk,aice_cpl(i,j,iblk)
!          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
       enddo
       enddo
    enddo

    !-------------------------------------------------

    call state_diagnose(exportState, 'cice_export', rc)

    export_slice = export_slice + 1
    call NUOPC_StateWrite(exportState, filePrefix='field_ice_export_', &
      timeslice=export_slice, relaxedFlag=.true., rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(*,*) 'CICE: --- run phase called --- 4'

  end subroutine 

  subroutine cice_model_finalize(gcomp, importState, exportState, clock, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Time)                        :: currTime

    write(*,*) 'CICE: --- finalize called ---'
    rc = ESMF_SUCCESS

    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call CICE_Finalize

    write(*,*) 'CICE: --- completed ---'

  end subroutine cice_model_finalize

  subroutine CICE_AdvertiseFields(state, nfields, field_defs, rc)

    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i

    rc = ESMF_SUCCESS

    do i = 1, nfields

      call NUOPC_StateAdvertiseField(state, &
        standardName=field_defs(i)%stdname, &
        name=field_defs(i)%shortname, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    enddo

  end subroutine CICE_AdvertiseFields

  subroutine CICE_RealizeFields(state, grid, nfields, field_defs, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Grid), intent(in)                 :: grid
    integer, intent(in)                         :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc

    integer                                     :: i
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

    do i = 1, nfields

      if (field_defs(i)%assoc) then
        write(info, *) '/', field_defs(i)%shortname, ':', &
          lbound(field_defs(i)%farrayPtr,1), ubound(field_defs(i)%farrayPtr,1), &
          lbound(field_defs(i)%farrayPtr,2), ubound(field_defs(i)%farrayPtr,2), &
          lbound(field_defs(i)%farrayPtr,3), ubound(field_defs(i)%farrayPtr,3)
        call ESMF_LogWrite(tag // " Field "// field_defs(i)%stdname // info, &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=rc)
        field = ESMF_FieldCreate(grid=grid, &
          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_DELOCAL, &
!          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_GLOBAL, &
          totalLWidth=(/1,1/), totalUWidth=(/1,1/),&
          ungriddedLBound=(/1/), ungriddedUBound=(/max_blocks/), &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
          totalLWidth=(/1,1/), totalUWidth=(/1,1/),&
          ungriddedLBound=(/1/), ungriddedUBound=(/max_blocks/), &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (NUOPC_StateIsFieldConnected(state, fieldName=field_defs(i)%shortname)) then
        call NUOPC_StateRealizeField(state, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite(tag // " Field "// field_defs(i)%stdname // " is connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_FieldPrint(field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        call ESMF_LogWrite(tag // " Field "// field_defs(i)%stdname // " is not connected.", &
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
        call ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

    enddo


  end subroutine CICE_RealizeFields

  !-----------------------------------------------------------------------------

  subroutine state_diagnose(State, string, rc)
    ! ----------------------------------------------
    ! Diagnose status of fieldBundle
    ! ----------------------------------------------
    type(ESMF_State), intent(inout) :: State
    character(len=*), intent(in), optional :: string
    integer, intent(out), optional  :: rc

    ! local variables
    integer                     :: i,j,n
    integer                     :: fieldCount
    character(len=64) ,pointer  :: fieldNameList(:)
    character(len=64)           :: lstring
    real(ESMF_KIND_R8), pointer :: dataPtr(:,:,:)
    character(len=256)          :: tmpstr
    integer                     :: lrc, dbrc
    character(len=*),parameter  :: subname='(cice_cap:state_diagnose)'

    lstring = ''
    if (present(string)) then
       lstring = trim(string)
    endif

    call ESMF_StateGet(State, itemCount=fieldCount, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(State, itemNameList=fieldNameList, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    do n = 1, fieldCount
      call State_GetFldPtr(State, fieldNameList(n), dataPtr, rc=lrc)
      if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      write(tmpstr,'(A,3g14.7)') trim(subname)//' '//trim(lstring)//':'//trim(fieldNameList(n)), &
        minval(dataPtr),maxval(dataPtr),sum(dataPtr)
!      write(tmpstr,'(A)') trim(subname)//' '//trim(lstring)//':'//trim(fieldNameList(n))
      call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo
    deallocate(fieldNameList)

    if (present(rc)) rc = lrc

  end subroutine state_diagnose

  !-----------------------------------------------------------------------------

  subroutine state_reset(State, value, rc)
    ! ----------------------------------------------
    ! Set all fields to value in State
    ! If value is not provided, reset to 0.0
    ! ----------------------------------------------
    type(ESMF_State), intent(inout) :: State
    real(ESMF_KIND_R8), intent(in), optional :: value
    integer, intent(out), optional  :: rc

    ! local variables
    integer                     :: i,j,k,n
    integer                     :: fieldCount
    character(len=64) ,pointer  :: fieldNameList(:)
    real(ESMF_KIND_R8)          :: lvalue
    real(ESMF_KIND_R8), pointer :: dataPtr(:,:,:)
    character(len=*),parameter :: subname='(cice_cap:state_reset)'

    if (present(rc)) rc = ESMF_SUCCESS

    lvalue = 0._ESMF_KIND_R8
    if (present(value)) then
      lvalue = value
    endif

    call ESMF_StateGet(State, itemCount=fieldCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    allocate(fieldNameList(fieldCount))
    call ESMF_StateGet(State, itemNameList=fieldNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    do n = 1, fieldCount
      call State_GetFldPtr(State, fieldNameList(n), dataPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

      do k=lbound(dataPtr,3),ubound(dataPtr,3)
      do j=lbound(dataPtr,2),ubound(dataPtr,2)
      do i=lbound(dataPtr,1),ubound(dataPtr,1)
         dataPtr(i,j,k) = lvalue
      enddo
      enddo
      enddo

    enddo
    deallocate(fieldNameList)

  end subroutine state_reset

  !-----------------------------------------------------------------------------

  subroutine State_GetFldPtr(ST, fldname, fldptr, rc)
    type(ESMF_State), intent(in) :: ST
    character(len=*), intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:,:,:)
    integer, intent(out), optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(cice_cap:State_GetFldPtr)'

    call ESMF_StateGet(ST, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (present(rc)) rc = lrc

  end subroutine State_GetFldPtr

  !-----------------------------------------------------------------------------

  subroutine CICE_FieldsSetup

!   call fld_list_add(fldsToIce_num, fldsToIce, "inst_zonal_wind_height10m", "will provide", strax)
!   call fld_list_add(fldsToIce_num, fldsToIce, "inst_merid_wind_height10m", "will provide", stray)
!   call fld_list_add(fldsToIce_num, fldsToIce, "inst_pressure_height_surface", "will provide", zlvl)
!   call fld_list_add(fldsToIce_num, fldsToIce, "xx_pot_air_temp", "will provide", potT)
!   call fld_list_add(fldsToIce_num, fldsToIce, "inst_temp_height2m", "will provide", Tair)
!   call fld_list_add(fldsToIce_num, fldsToIce, "inst_spec_humid_height2m", "will provide", Qa)
!   call fld_list_add(fldsToIce_num, fldsToIce, "xx_inst_air_density", "will provide", rhoa)
!   call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_vis_dir_flx", "will provide", swvdr)
!   call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_vis_dif_flx", "will provide", swvdf)
!   call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_ir_dir_flx", "will provide", swidr)
!   call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_ir_dif_flx", "will provide", swidf)
!   call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_lw_flx", "will provide", flw)
!   call fld_list_add(fldsToIce_num, fldsToIce, "mean_prec_rate", "will provide", frain)
!   call fld_list_add(fldsToIce_num, fldsToIce, "xx_mean_fprec_rate", "will provide", frain)
!   call fld_list_add(fldsToIce_num, fldsToIce, "xx_faero_atm", "will provide", faero_atm)
!   call fld_list_add(fldsToIce_num, fldsToIce, "ocn_current_zonal", "will provide", uocn)
!   call fld_list_add(fldsToIce_num, fldsToIce, "ocn_current_merid", "will provide", vocn)
!   call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_slope_zonal", "will provide", ss_tltx)
!   call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_slope_merid", "will provide", ss_tlty)
!   call fld_list_add(fldsToIce_num, fldsToIce, "s_surf", "will provide", sss)
!   call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_temperature", "will provide", sst)
!   call fld_list_add(fldsToIce_num, fldsToIce, "freezing_melting_potential", "will provide", frzmlt)
!   call fld_list_add(fldsToIce_num, fldsToIce, "xx_inst_frz_mlt_potential", "will provide", frzmlt_init)
!   call fld_list_add(fldsToIce_num, fldsToIce, "freezing_temp", "will provide", Tf)
!   call fld_list_add(fldsToIce_num, fldsToIce, "mean_deep_ocean_down_heat_flx", "will provide", qdp)
!   call fld_list_add(fldsToIce_num, fldsToIce, "mixed_layer_depth", "will provide", hmix)
!   call fld_list_add(fldsToIce_num, fldsToIce, "xx_daice_da", "will provide", daice_da)

! tcraig, don't point directly into cice data YET (last field is optional in interface)
! instead, create space for the field when it's "realized".
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_zonal_wind_height10m", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_merid_wind_height10m", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_pressure_height_surface", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_temp_height2m", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_spec_humid_height2m", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_vis_dir_flx", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_vis_dif_flx", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_ir_dir_flx", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_ir_dif_flx", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_lw_flx", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_prec_rate", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "ocn_current_zonal", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "ocn_current_merid", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_slope_zonal", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_slope_merid", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "s_surf", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_temperature", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "freezing_melting_potential", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "freezing_temp", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_deep_ocean_down_heat_flx", "will provide")
    call fld_list_add(fldsToIce_num, fldsToIce, "mixed_layer_depth", "will provide")


    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_air_ice_zonal", "will provide", strairxT)
    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_air_ice_merid", "will provide", strairyT)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_sensi_heat_flx", "will provide", fsens)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_laten_heat_flx", "will provide", flat)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_fswabs", "will provide", fswabs)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_up_lw_flx", "will provide", flwout)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_inst_temp_height2m", "will provide", Tref)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_inst_spec_humid_height2m", "will provide", Qref)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_evap_rate", "will provide", evap)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_inst_vis_dir_albedo", "will provide", alvdr)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_inst_ir_dir_albedo", "will provide", alidr)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_inst_vis_dif_albedo", "will provide", alvdf)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_inst_ir_dif_albedo", "will provide", alidf)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_albedo_vis_dir", "will provide", alvdr_ai)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_albedo_nir_dir", "will provide", alidr_ai)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_albedo_vis_dif", "will provide", alvdf_ai)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_albedo_nir_dif", "will provide", alidf_ai)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_bare_ice_albedo", "will provide", albice)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_snow_albedo", "will provide", albsno)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_melt_pond_albedo", "will provide", albpnd)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_apeff_ai", "will provide", apeff_ai)

    allocate(aice_cpl(nx_block,ny_block,max_blocks))
    call fld_list_add(fldsFrIce_num, fldsFrIce, "ice_fraction", "will provide", aice_cpl)

    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_ocn_ice_zonal", "will provide", strocnxT)
    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_ocn_ice_merid", "will provide", strocnyT)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_fresh_water_flx_to_ponds", "will provide", fpond)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_fresh_water_to_ocean_rate", "will provide", fresh)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_salt_rate", "will provide", fsalt)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_net_heat_flx_to_ocn", "will provide", fhocn)
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_sw_pen_to_ocean", "will provide", fswthru)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_faero_ocn", "will provide", faero_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_strairx_ocn", "will provide", strairx_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_strairy_ocn", "will provide", strairy_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_sensi_heat_flx", "will provide", fsens_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_laten_heat_flx", "will provide", flat_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_flwout_ocn", "will provide", flwout_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_evap_ocn", "will provide", evap_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_albedo_vis_dir", "will provide", alvdr_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_albedo_nir_dir", "will provide", alidr_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_albedo_vis_dif", "will provide", alvdf_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_albedo_nir_dif", "will provide", alidf_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_2m_atm_ref_temperature", "will provide", Tref_ocn)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_2m_atm_ref_spec_humidity", "will provide", Qref_ocn)

  end subroutine CICE_FieldsSetup

  !-----------------------------------------------------------------------------

  subroutine fld_list_add(num, fldlist, stdname, transferOffer, data)
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    integer,             intent(inout)  :: num
    type(fld_list_type), intent(inout)  :: fldlist(:)
    character(len=*),    intent(in)     :: stdname
    character(len=*),    intent(in)     :: transferOffer
    real(ESMF_KIND_R8), dimension(:,:,:), optional, target :: data

    ! local variables
    integer :: rc
    character(len=256)          :: shortname
    character(len=*), parameter :: subname='(cice_cap:fld_list_add)'

    ! make sure that stdname is in the NUOPC Field Dictionary 
    call NUOPC_FieldDictionaryGetEntry(stdname, defaultShortName=shortname, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=trim(subname)//&
      ": invalid stdname: "//trim(stdname), &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! fill in the new entry

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//&
      ": ERROR num gt fldsMax "//trim(stdname), ESMF_LOGMSG_ERROR)
      return
    endif

    fldlist(num)%stdname        = trim(stdname)
    fldlist(num)%shortname      = trim(shortname)
    fldlist(num)%transferOffer  = trim(transferOffer)
    if (present(data)) then
      fldlist(num)%assoc        = .true.
      fldlist(num)%farrayPtr    => data
    else
      fldlist(num)%assoc        = .false.
    endif

  end subroutine fld_list_add

  !-----------------------------------------------------------------------------
end module cice_cap_mod
