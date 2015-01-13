!--------------- LANL CICE NUOPC CAP -----------------
! This is the LANL CICE model cap component that's NUOPC compiant.
!
! Author:  Fei.Liu@gmail.com
!
! 5/10/13
! This is now acting as a cap/connector between NUOPC driver and LANL CICE code.
!

module cice_cap_mod

  use ice_blocks, only: nx_block, ny_block, nblocks_tot, block, get_block, &
                        get_block_parameter
  use ice_domain_size, only: max_blocks, nx_global, ny_global
  use ice_domain, only: nblocks, blocks_ice, distrb_info
  use ice_distribution, only: ice_distributiongetblockloc
  use ice_constants, only: Tffresh, rad_to_deg
  use ice_flux
  use ice_grid, only: TLAT, TLON, ULAT, ULON, hm, tarea
  use ice_state
  use CICE_RunMod
  use CICE_InitMod
  use CICE_FinalMod 

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS      => SetServices, &
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
  character(len=256) :: tmpstr
  character(len=2048):: info
  logical :: isPresent
  integer :: dbrc     ! temporary debug rc value

  contains
  !-----------------------------------------------------------------------
  !------------------- CICE code starts here -----------------------
  !-----------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname='(cice_cap:SetServices)'

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
    character(len=*),parameter  :: subname='(cice_cap:InitializeP1)'

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

    write(info,*) subname,' --- initialization phase 1 completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

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
    integer                                :: i,j,iblk, n, i1,j1, DE
    integer                                :: ilo,ihi,jlo,jhi
    integer                                :: ig,jg,cnt
    integer                                :: peID,locID
    integer, pointer                       :: indexList(:)
    integer, pointer                       :: deLabelList(:)
    integer, pointer                       :: deBlockList(:,:,:)
    integer, pointer                       :: petMap(:)
    integer, pointer                       :: i_glob(:),j_glob(:)
    integer                                :: lbnd(2),ubnd(2)
    type(block)                            :: this_block
    type(ESMF_DELayout)                    :: delayout
    real(ESMF_KIND_R8), pointer            :: tarray(:,:)     
    real(ESMF_KIND_R8), pointer :: coordXcenter(:,:)
    real(ESMF_KIND_R8), pointer :: coordYcenter(:,:)
    real(ESMF_KIND_R8), pointer :: coordXcorner(:,:)
    real(ESMF_KIND_R8), pointer :: coordYcorner(:,:)
    integer(ESMF_KIND_I4), pointer :: gridmask(:,:)
    real(ESMF_KIND_R8), pointer :: gridarea(:,:)
    character(len=*),parameter  :: subname='(cice_cap:InitializeP2)'

    rc = ESMF_SUCCESS

    ! We can check if npet is 4 or some other value to make sure
    ! CICE is configured to run on the correct number of processors.

    ! create a Grid object for Fields
    ! we are going to create a single tile displaced pole grid from a gridspec
    ! file. We also use the exact decomposition in CICE so that the Fields
    ! created can wrap on the data pointers in internal part of CICE

    write(tmpstr,'(a,2i8)') subname//' ice nx,ny = ',nx_global,ny_global
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

!    distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/nx_global,ny_global/), &
!       regDecomp=(/2,2/), rc=rc)

    allocate(deBlockList(2,2,nblocks_tot))
    allocate(petMap(nblocks_tot))
    allocate(deLabelList(nblocks_tot))

    write(tmpstr,'(a,1i8)') subname//' nblocks = ',nblocks_tot
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    do n = 1, nblocks_tot
       deLabelList(n) = n
       call get_block_parameter(n,ilo=ilo,ihi=ihi,jlo=jlo,jhi=jhi, &
          i_glob=i_glob,j_glob=j_glob)
       deBlockList(1,1,n) = i_glob(ilo)
       deBlockList(1,2,n) = i_glob(ihi)
       deBlockList(2,1,n) = j_glob(jlo)
       deBlockList(2,2,n) = j_glob(jhi)
       call ice_distributionGetBlockLoc(distrb_info,n,peID,locID)
       petMap(n) = peID - 1
       write(tmpstr,'(a,2i8)') subname//' IDs  = ',n,peID
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(tmpstr,'(a,3i8)') subname//' iglo = ',n,deBlockList(1,1,n),deBlockList(1,2,n)
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(tmpstr,'(a,3i8)') subname//' jglo = ',n,deBlockList(2,1,n),deBlockList(2,2,n)
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    enddo

    delayout = ESMF_DELayoutCreate(petMap, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/nx_global,ny_global/), &
!        indexflag = ESMF_INDEX_DELOCAL, &
        deBlockList=deBlockList, &
!        deLabelList=deLabelList, &
        delayout=delayout, &
        rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    deallocate(deLabelList)
    deallocate(deBlockList)
    deallocate(petMap)

!    call ESMF_DistGridPrint(distgrid=distgrid, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, elementCount=cnt, rc=rc)
    allocate(indexList(cnt))
    write(tmpstr,'(a,i8)') subname//' distgrid cnt= ',cnt
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    call ESMF_DistGridGet(distgrid=distgrid, localDE=0, seqIndexList=indexList, rc=rc)
    write(tmpstr,'(a,4i8)') subname//' distgrid list= ',indexList(1),indexList(cnt),minval(indexList), maxval(indexList)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    deallocate(IndexList)

!    gridIn = ESMF_GridCreate('global_gx3_gridspec.nc', ESMF_FILEFORMAT_GRIDSPEC, &
!!      (/2,2/), isSphere=.true., coordNames=(/'ulon', 'ulat'/), &
!      distgrid=distgrid, isSphere=.true., coordNames=(/'ulon', 'ulat'/), &
!      indexflag=ESMF_INDEX_DELOCAL, addCornerStagger=.true., rc=rc)
!!      indexflag=ESMF_INDEX_GLOBAL, addCornerStagger=.true., rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    gridIn = ESMF_GridCreate(distgrid=distgrid, &
       coordSys = ESMF_COORDSYS_SPH_DEG, &
       rc = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CORNER, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_MASK, itemTypeKind=ESMF_TYPEKIND_I4, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_GridAddItem(gridIn, itemFlag=ESMF_GRIDITEM_AREA, itemTypeKind=ESMF_TYPEKIND_R8, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    do iblk = 1,nblocks
       DE = iblk-1
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi

       call ESMF_GridGetCoord(gridIn, coordDim=1, localDE=DE, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           computationalLBound=lbnd, computationalUBound=ubnd, &
           farrayPtr=coordXcenter, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       call ESMF_GridGetCoord(gridIn, coordDim=2, localDE=DE, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           farrayPtr=coordYcenter, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       write(tmpstr,'(a,5i8)') subname//' iblk center bnds ',iblk,lbnd,ubnd
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       if (lbnd(1) /= 1 .or. lbnd(2) /= 1 .or. ubnd(1) /= ihi-ilo+1 .or. ubnd(2) /= jhi-jlo+1) then
          write(tmpstr,'(a,5i8)') subname//' iblk bnds ERROR '
          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
          rc = ESMF_FAILURE
          return
       endif

       call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_MASK, localDE=DE, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           farrayPtr=gridmask, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       call ESMF_GridGetItem(gridIn, itemflag=ESMF_GRIDITEM_AREA, localDE=DE, &
           staggerloc=ESMF_STAGGERLOC_CENTER, &
           farrayPtr=gridarea, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       do j1 = lbnd(2),ubnd(2)
       do i1 = lbnd(1),ubnd(1)
          i = i1 + ilo - lbnd(1)
          j = j1 + jlo - lbnd(2)
          coordXcenter(i1,j1) = TLON(i,j,iblk) * rad_to_deg
          coordYcenter(i1,j1) = TLAT(i,j,iblk) * rad_to_deg
          gridmask(i1,j1) = nint(hm(i,j,iblk))
          gridarea(i1,j1) = tarea(i,j,iblk)
       enddo
       enddo

       call ESMF_GridGetCoord(gridIn, coordDim=1, localDE=DE, &
           staggerloc=ESMF_STAGGERLOC_CORNER, &
           computationalLBound=lbnd, computationalUBound=ubnd, &
           farrayPtr=coordXcorner, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       call ESMF_GridGetCoord(gridIn, coordDim=2, localDE=DE, &
           staggerloc=ESMF_STAGGERLOC_CORNER, &
           farrayPtr=coordYcorner, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

       write(tmpstr,'(a,5i8)') subname//' iblk corner bnds ',iblk,lbnd,ubnd
       call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

       do j1 = lbnd(2),ubnd(2)
       do i1 = lbnd(1),ubnd(1)
          i = i1 + ilo - lbnd(1)
          j = j1 + jlo - lbnd(2)
          coordXcorner(i1,j1) = ULON(i,j,iblk) * rad_to_deg
          coordYcorner(i1,j1) = ULAT(i,j,iblk) * rad_to_deg
       enddo
       enddo

    enddo

    call ESMF_GridGetCoord(gridIn, coordDim=1, localDE=0,  &
       staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=tarray, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(tmpstr,'(a,2g15.7)') subname//' gridIn center1 = ',minval(tarray),maxval(tarray)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_GridGetCoord(gridIn, coordDim=2, localDE=0,  &
       staggerLoc=ESMF_STAGGERLOC_CENTER, farrayPtr=tarray, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(tmpstr,'(a,2g15.7)') subname//' gridIn center2 = ',minval(tarray),maxval(tarray)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_GridGetCoord(gridIn, coordDim=1, localDE=0,  &
       staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=tarray, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(tmpstr,'(a,2g15.7)') subname//' gridIn corner1 = ',minval(tarray),maxval(tarray)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_GridGetCoord(gridIn, coordDim=2, localDE=0,  &
       staggerLoc=ESMF_STAGGERLOC_CORNER, farrayPtr=tarray, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(tmpstr,'(a,2g15.7)') subname//' gridIn corner2 = ',minval(tarray),maxval(tarray)
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

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
    call state_reset(ExportState, value=-99._ESMF_KIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!    call State_getFldPtr(exportState,'ifrac'    ,dataPtr_ifrac,rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
!    call State_getFldPtr(exportState,'sit'      ,dataPtr_itemp,rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
!    dataPtr_ifrac = -99._ESMF_KIND_R8
!    dataPtr_itemp = -99._ESMF_KIND_R8

    write(info,*) subname,' --- initialization phase 2 completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)

  end subroutine
  
  !-----------------------------------------------------------------------------

  ! CICE model uses same clock as parent gridComp
  subroutine SetClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_TimeInterval)       :: stabilityTimeStep, timestep
    character(len=*),parameter  :: subname='(cice_cap:SetClock)'

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
    call NUOPC_CompSetClock(gcomp, clock, stabilityTimeStep, rc=rc)
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
    type(ESMF_Field)                       :: lfield
    type(block)                            :: this_block
    character(len=64)                      :: fldname
    integer                                :: i,j,iblk,n,i1,i2,j1,j2
    integer                                :: ilo,ihi,jlo,jhi
    type(ESMF_StateItem_Flag)              :: itemType
    ! imports
    real(ESMF_KIND_R8), pointer :: dataPtr_ith2m(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ishh2m(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_izwh10m(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_imwh10m(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ips(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_mdlwfx(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swvr(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swvf(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swir(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_swif(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_lprec(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fprec(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sst(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sss(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sssz(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_sssm(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ocncz(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_ocncm(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fmpot(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_mld(:,:,:)
    ! exports
    real(ESMF_KIND_R8), pointer :: dataPtr_ifrac(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_itemp(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alvdr(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alidr(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alvdf(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_alidf(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strairxT(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strairyT(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strocnxT(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_strocnyT(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fswthru(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_flwout(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_fsens(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_flat(:,:,:)
    real(ESMF_KIND_R8), pointer :: dataPtr_evap(:,:,:)
    character(len=*),parameter  :: subname='(cice_cap:ModelAdvance_slow)'

    rc = ESMF_SUCCESS
    write(info,*) subname,' --- run phase 1 called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
    
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
!    call NUOPC_StateWrite(importState, filePrefix='field_ice_import_', &
!      timeslice=import_slice, relaxedFlag=.true., rc=rc) 
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    ! tcraig, skip first field (dummyfield) due to StateWrite errors
    do i = 2,fldsToice_num
      fldname = fldsToice(i)%shortname
      call ESMF_StateGet(importState, itemName=trim(fldname), itemType=itemType, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (itemType /= ESMF_STATEITEM_NOTFOUND) then
        call ESMF_StateGet(importState, itemName=trim(fldname), field=lfield, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        call NUOPC_FieldWrite(lfield, file='field_ice_import_'//trim(fldname)//'.nc', &
          timeslice=import_slice, relaxedFlag=.true., rc=rc) 
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
    enddo

    call State_getFldPtr(importState,'ith2m'  ,dataPtr_ith2m,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'ishh2m' ,dataPtr_ishh2m,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'izwh10m',dataPtr_izwh10m,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'imwh10m',dataPtr_imwh10m,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'ips'    ,dataPtr_ips,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mdlwfx' ,dataPtr_mdlwfx,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'sw_flux_vis_dir',dataPtr_swvr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'sw_flux_vis_dif',dataPtr_swvf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'sw_flux_nir_dir',dataPtr_swir,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'sw_flux_nir_dif',dataPtr_swif,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'lprec'  ,dataPtr_lprec,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'fprec'  ,dataPtr_fprec,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'sst'    ,dataPtr_sst,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'sss'    ,dataPtr_sss,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'sssz'   ,dataPtr_sssz,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'sssm'   ,dataPtr_sssm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'ocncz'  ,dataPtr_ocncz,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'ocncm'  ,dataPtr_ocncm,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'fmpot'  ,dataPtr_fmpot,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(importState,'mld'    ,dataPtr_mld,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return

    do iblk = 1,nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo,jhi
       do i = ilo,ihi
          i1 = i - ilo + 1
          j1 = j - jlo + 1
!          Tair   (i,j,iblk) = dataPtr_ith2m  (i1,j1,iblk)
!          Qa     (i,j,iblk) = dataPtr_ishh2m (i1,j1,iblk)
!          strax  (i,j,iblk) = dataPtr_izwh10m(i1,j1,iblk)
!          stray  (i,j,iblk) = dataPtr_imwh10m(i1,j1,iblk)
!          zlvl   (i,j,iblk) = dataPtr_ips    (i1,j1,iblk)
!          flw    (i,j,iblk) = dataPtr_mdlwfx (i1,j1,iblk)
!          swvdr  (i,j,iblk) = dataPtr_swvr   (i1,j1,iblk)
!          swvdf  (i,j,iblk) = dataPtr_swvf   (i1,j1,iblk)
!          swidr  (i,j,iblk) = dataPtr_swir   (i1,j1,iblk)
!          swidf  (i,j,iblk) = dataPtr_swif   (i1,j1,iblk)
!          frain  (i,j,iblk) = dataPtr_lprec  (i1,j1,iblk)
!          ??     (i,j,iblk) = dataPtr_fprec  (i1,j1,iblk)
!          sst    (i,j,iblk) = dataPtr_sst    (i1,j1,iblk)
!          sss    (i,j,iblk) = dataPtr_sss    (i1,j1,iblk)
!          ss_tltx(i,j,iblk) = dataPtr_sssz   (i1,j1,iblk)
!          ss_tlty(i,j,iblk) = dataPtr_sssm   (i1,j1,iblk)
!          uocn   (i,j,iblk) = dataPtr_ocncz  (i1,j1,iblk)
!          vocn   (i,j,iblk) = dataPtr_ocncm  (i1,j1,iblk)
!          frzmlt (i,j,iblk) = dataPtr_fmpot  (i1,j1,iblk)
!          hmix   (i,j,iblk) = dataPtr_mld    (i1,j1,iblk)
!!          write(tmpstr,'(a,3i6,2x,g17.7)') subname//' sst = ',i,j,iblk,dataPtr_sst(i,j,iblk)
!!          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       enddo
       enddo
    enddo

    write(info,*) subname,' --- run phase 2 called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
    call CICE_Run
    write(info,*) subname,' --- run phase 3 called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    !---- local modifications to coupling fields -----
    !---- good place to rotate vectors           -----

    call State_getFldPtr(exportState,'ifrac'    ,dataPtr_ifrac,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'sit'      ,dataPtr_itemp,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'iivisdira',dataPtr_alvdr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'iivisdifa',dataPtr_alvdf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'iiirdira' ,dataPtr_alidr,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'iiirdifa' ,dataPtr_alidf,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'strairxT' ,dataPtr_strairxT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'strairyT' ,dataPtr_strairyT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'strocnxT' ,dataPtr_strocnxT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'strocnyT' ,dataPtr_strocnyT,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mswpenocn',dataPtr_fswthru,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mulwfxi'  ,dataPtr_flwout,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mshfxai'  ,dataPtr_fsens,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mlhfxai'  ,dataPtr_flat,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return
    call State_getFldPtr(exportState,'mevapai'  ,dataPtr_evap,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) return

    write(info, *) subname//' ifrac size :', &
      lbound(dataPtr_ifrac,1), ubound(dataPtr_ifrac,1), &
      lbound(dataPtr_ifrac,2), ubound(dataPtr_ifrac,2), &
      lbound(dataPtr_ifrac,3), ubound(dataPtr_ifrac,3)
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    dataPtr_ifrac = 0.
    dataPtr_itemp = 0.
    do iblk = 1,nblocks
       this_block = get_block(blocks_ice(iblk),iblk)
       ilo = this_block%ilo
       ihi = this_block%ihi
       jlo = this_block%jlo
       jhi = this_block%jhi
       do j = jlo,jhi
       do i = ilo,ihi
          i1 = i - ilo + 1
          j1 = j - jlo + 1
          dataPtr_ifrac   (i1,j1,iblk) = aice(i,j,iblk)
          dataPtr_itemp   (i1,j1,iblk) = Tffresh + trcr(i,j,1,iblk)
          dataPtr_alvdr   (i1,j1,iblk) = alvdr(i,j,iblk)
          dataPtr_alidr   (i1,j1,iblk) = alidr(i,j,iblk)
          dataPtr_alvdf   (i1,j1,iblk) = alvdf(i,j,iblk)
          dataPtr_alidf   (i1,j1,iblk) = alidf(i,j,iblk)
          dataPtr_strairxT(i1,j1,iblk) = strairxT(i,j,iblk)
          dataPtr_strairyT(i1,j1,iblk) = strairyT(i,j,iblk)
          dataPtr_strocnxT(i1,j1,iblk) = strocnxT(i,j,iblk)
          dataPtr_strocnyT(i1,j1,iblk) = strocnyT(i,j,iblk)
          dataPtr_fswthru (i1,j1,iblk) = fswthru(i,j,iblk)
          dataPtr_flwout  (i1,j1,iblk) = flwout(i,j,iblk)
          dataPtr_fsens   (i1,j1,iblk) = fsens(i,j,iblk)
          dataPtr_flat    (i1,j1,iblk) = flat(i,j,iblk)
          dataPtr_evap    (i1,j1,iblk) = evap(i,j,iblk)
!!          write(tmpstr,'(a,3i6,2x,g17.7)') subname//' aice = ',i,j,iblk,dataPtr_ifrac(i,j,iblk)
!!          call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       enddo
       enddo
    enddo

    !-------------------------------------------------

    call state_diagnose(exportState, 'cice_export', rc)

    export_slice = export_slice + 1
!    call NUOPC_StateWrite(exportState, filePrefix='field_ice_export_', &
!      timeslice=export_slice, relaxedFlag=.true., rc=rc) 
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
    ! tcraig, skip first field (dummyfield) due to StateWrite errors
    do i = 2,fldsFrIce_num
      fldname = fldsFrIce(i)%shortname
      call ESMF_StateGet(exportState, itemName=trim(fldname), itemType=itemType, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (itemType /= ESMF_STATEITEM_NOTFOUND .and. (fldname=='ifrac' .or. fldname=='sit')) then
!      if (itemType /= ESMF_STATEITEM_NOTFOUND) then
        call ESMF_StateGet(exportState, itemName=trim(fldname), field=lfield, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
!        call NUOPC_FieldWrite(lfield, file='field_ice_export_'//trim(fldname)//'.nc', &
!          timeslice=export_slice, relaxedFlag=.true., rc=rc) 
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
      endif
    enddo

    write(info,*) subname,' --- run phase 4 called --- ',rc
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine 

  subroutine cice_model_finalize(gcomp, importState, exportState, clock, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Time)                        :: currTime
    character(len=*),parameter  :: subname='(cice_cap:cice_model_finalize)'

    rc = ESMF_SUCCESS

    write(info,*) subname,' --- finalize called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call CICE_Finalize

    write(info,*) subname,' --- finalize completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine cice_model_finalize

  subroutine CICE_AdvertiseFields(state, nfields, field_defs, rc)

    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i
    character(len=*),parameter  :: subname='(cice_cap:CICE_AdvertiseFields)'

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
    character(len=*),parameter  :: subname='(cice_cap:CICE_RealizeFields)'
 
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
      !call ESMF_LogWrite(subname // tag // " Grid "// info, &
      !  ESMF_LOGMSG_INFO, &
      !  line=__LINE__, &
      !  file=__FILE__, &
      !  rc=dbrc)

    do i = 1, nfields

      if (field_defs(i)%assoc) then
        write(info, *) subname, tag, ' Field ', field_defs(i)%shortname, ':', &
          lbound(field_defs(i)%farrayPtr,1), ubound(field_defs(i)%farrayPtr,1), &
          lbound(field_defs(i)%farrayPtr,2), ubound(field_defs(i)%farrayPtr,2), &
          lbound(field_defs(i)%farrayPtr,3), ubound(field_defs(i)%farrayPtr,3)
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
        field = ESMF_FieldCreate(grid=grid, &
          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_DELOCAL, &
!          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_GLOBAL, &
!          totalLWidth=(/1,1/), totalUWidth=(/1,1/),&
          ungriddedLBound=(/1/), ungriddedUBound=(/max_blocks/), &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
!          totalLWidth=(/1,1/), totalUWidth=(/1,1/),&
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
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=dbrc)
!        call ESMF_FieldPrint(field=field, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
      else
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is not connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=dbrc)
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
    ! Diagnose status of state
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
    integer                     :: lrc
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
  logical function FieldBundle_FldChk(FB, fldname, rc)
    type(ESMF_FieldBundle), intent(in) :: FB
    character(len=*)      ,intent(in) :: fldname
    integer, intent(out), optional :: rc

    ! local variables
    integer :: lrc
    character(len=*),parameter :: subname='(module_MEDIATOR:FieldBundle_FldChk)'

    if (present(rc)) rc = ESMF_SUCCESS

    FieldBundle_FldChk = .false.

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), isPresent=isPresent, rc=lrc)
    if (present(rc)) rc = lrc
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (isPresent) then
       FieldBundle_FldChk = .true.
    endif

  end function FieldBundle_FldChk

  !-----------------------------------------------------------------------------

  subroutine FieldBundle_GetFldPtr(FB, fldname, fldptr, rc)
    type(ESMF_FieldBundle), intent(in) :: FB
    character(len=*)      , intent(in) :: fldname
    real(ESMF_KIND_R8), pointer, intent(in) :: fldptr(:,:)
    integer, intent(out), optional :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(module_MEDIATOR:FieldBundle_GetFldPtr)'

    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_FieldBundleGet(FB, fieldName=trim(fldname), field=lfield, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    if (present(rc)) rc = lrc

  end subroutine FieldBundle_GetFldPtr

  !-----------------------------------------------------------------------------

  subroutine CICE_FieldsSetup
    character(len=*),parameter  :: subname='(cice_cap:CICE_FieldsSetup)'

!--------- import fields to Sea Ice -------------

! tcraig, don't point directly into cice data YET (last field is optional in interface)
! instead, create space for the field when it's "realized".
    call fld_list_add(fldsToIce_num, fldsToIce, "dummyfield"               , "will provide", shortname="dummyfield")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_temp_height2m"       , "will provide", shortname="ith2m")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_spec_humid_height2m" , "will provide", shortname="ishh2m")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_zonal_wind_height10m", "will provide", shortname="izwh10m")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_merid_wind_height10m", "will provide", shortname="imwh10m")
    call fld_list_add(fldsToIce_num, fldsToIce, "inst_pres_height_surface" , "will provide", shortname="ips")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_lw_flx"         , "will provide", shortname="mdlwfx")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_vis_dir_flx" , "will provide", shortname="sw_flux_vis_dir")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_vis_dif_flx" , "will provide", shortname="sw_flux_vis_dif")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_ir_dir_flx"  , "will provide", shortname="sw_flux_nir_dir")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_ir_dif_flx"  , "will provide", shortname="sw_flux_nir_dif")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_prec_rate"           , "will provide", shortname="lprec")
    call fld_list_add(fldsToIce_num, fldsToIce, "mean_fprec_rate"          , "will provide", shortname="fprec")
    call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_temperature"  , "will provide", shortname="sst")
    call fld_list_add(fldsToIce_num, fldsToIce, "s_surf"                   , "will provide", shortname="sss")
    call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_slope_zonal"  , "will provide", shortname="sssz")
    call fld_list_add(fldsToIce_num, fldsToIce, "sea_surface_slope_merid"  , "will provide", shortname="sssm")
    call fld_list_add(fldsToIce_num, fldsToIce, "ocn_current_zonal"        , "will provide", shortname="ocncz")
    call fld_list_add(fldsToIce_num, fldsToIce, "ocn_current_merid"        , "will provide", shortname="ocncm")
    call fld_list_add(fldsToIce_num, fldsToIce, "freezing_melting_potential", "will provide", shortname="fmpot")
    call fld_list_add(fldsToIce_num, fldsToIce, "mixed_layer_depth"        , "will provide", shortname="mld")

!   call fld_list_add(fldsToIce_num, fldsToIce, "inst_zonal_wind_height10m", "will provide", strax)
!   call fld_list_add(fldsToIce_num, fldsToIce, "inst_merid_wind_height10m", "will provide", stray)
!   call fld_list_add(fldsToIce_num, fldsToIce, "inst_pres_height_surface" , "will provide", zlvl)
!   call fld_list_add(fldsToIce_num, fldsToIce, "xx_pot_air_temp"          , "will provide", potT)
!   call fld_list_add(fldsToIce_num, fldsToIce, "inst_temp_height2m"       , "will provide", Tair)
!   call fld_list_add(fldsToIce_num, fldsToIce, "inst_spec_humid_height2m" , "will provide", Qa)
!   call fld_list_add(fldsToIce_num, fldsToIce, "xx_inst_air_density"      , "will provide", rhoa)
!   call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_vis_dir_flx" , "will provide", swvdr)
!   call fld_list_add(fldsToIce_num, fldsToIce, "mean_down_sw_vis_dif_flx" , "will provide", swvdf)
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

!--------- export fields from Sea Ice -------------

!tcx    call fld_list_add(fldsFrIce_num, fldsFrIce, "dummyfield"                      , "will provide", dummyfield)
!tcx    call fld_list_add(fldsFrIce_num, fldsFrIce, "sea_ice_temperature"             , "will provide", icetemp_cpl)
!tcx    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_vis_dir_albedo"         , "will provide", alvdr)
!tcx    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_ir_dir_albedo"          , "will provide", alidr)
!tcx    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_vis_dif_albedo"         , "will provide", alvdf)
!tcx    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_ir_dif_albedo"          , "will provide", alidf)
!tcx    call fld_list_add(fldsFrIce_num, fldsFrIce, "ice_fraction"                    , "will provide", aice_cpl)
!    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_air_ice_zonal"         , "will provide", strairxT)
!    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_air_ice_merid"         , "will provide", strairyT)
!    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_ocn_ice_zonal"         , "will provide", strocnxT)
!    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_ocn_ice_merid"         , "will provide", strocnyT)
!    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_sw_pen_to_ocn"              , "will provide", fswthru)
!    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_up_lw_flx_ice"              , "will provide", flwout)
!    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_sensi_heat_flx_atm_into_ice", "will provide", fsens)
!    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_laten_heat_flx_atm_into_ice", "will provide", flat)
!tcx    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_evap_rate_atm_into_ice"     , "will provide", evap)
    call fld_list_add(fldsFrIce_num, fldsFrIce, "dummyfield"                      , "will provide", shortname="dummyfield")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "sea_ice_temperature"             , "will provide", shortname="sit")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_vis_dir_albedo"         , "will provide", shortname="iivisdira")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_ir_dir_albedo"          , "will provide", shortname="iiirdira")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_vis_dif_albedo"         , "will provide", shortname="iivisdifa")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "inst_ice_ir_dif_albedo"          , "will provide", shortname="iiirdifa")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "ice_fraction"                    , "will provide", shortname="ifrac")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_air_ice_zonal"         , "will provide", shortname="strairxT")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_air_ice_merid"         , "will provide", shortname="strairyT")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_ocn_ice_zonal"         , "will provide", shortname="strocnxT")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "stress_on_ocn_ice_merid"         , "will provide", shortname="strocnyT")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_sw_pen_to_ocn"              , "will provide", shortname="mswpenocn")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_up_lw_flx_ice"              , "will provide", shortname="mulwfxi")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_sensi_heat_flx_atm_into_ice", "will provide", shortname="mshfxai")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_laten_heat_flx_atm_into_ice", "will provide", shortname="mlhfxai")
    call fld_list_add(fldsFrIce_num, fldsFrIce, "mean_evap_rate_atm_into_ice"     , "will provide", shortname="mevapai")

!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_inst_temp_height2m", "will provide", Tref)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_inst_spec_humid_height2m", "will provide", Qref)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_albedo_vis_dir", "will provide", alvdr_ai)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_albedo_nir_dir", "will provide", alidr_ai)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_albedo_vis_dif", "will provide", alvdf_ai)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_albedo_nir_dif", "will provide", alidf_ai)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_bare_ice_albedo", "will provide", albice)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_snow_albedo", "will provide", albsno)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_melt_pond_albedo", "will provide", albpnd)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_apeff_ai", "will provide", apeff_ai)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_fresh_water_flx_to_ponds", "will provide", fpond)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_fresh_water_to_ocean_rate", "will provide", fresh)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_mean_salt_rate", "will provide", fsalt)
!   call fld_list_add(fldsFrIce_num, fldsFrIce, "xx_net_heat_flx_to_ocn", "will provide", fhocn)
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

  subroutine fld_list_add(num, fldlist, stdname, transferOffer, data, shortname)
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    integer,             intent(inout)  :: num
    type(fld_list_type), intent(inout)  :: fldlist(:)
    character(len=*),    intent(in)     :: stdname
    character(len=*),    intent(in)     :: transferOffer
    real(ESMF_KIND_R8), dimension(:,:,:), optional, target :: data
    character(len=*),    intent(in),optional :: shortname

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(cice_cap:fld_list_add)'

    ! fill in the new entry

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num gt fldsMax "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
      return
    endif

    fldlist(num)%stdname        = trim(stdname)
    if (present(shortname)) then
       fldlist(num)%shortname   = trim(shortname)
    else
       fldlist(num)%shortname   = trim(stdname)
    endif
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
