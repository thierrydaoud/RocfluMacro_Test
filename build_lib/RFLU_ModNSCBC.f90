










!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
! ******************************************************************************
!
! Purpose: Collection of routines for NCSBC boundary treatment
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModNSCBC.F90,v 1.3 2016/02/08 20:03:05 fred Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModNSCBC

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModError
  USE ModMPI

  USE RFLU_ModAxisymmetry, ONLY: RFLU_AXI_SourceTermsNSCBC
  USE RFLU_ModBoundConvertCv, ONLY: RFLU_BXV_ConvertCvCons2Prim, &
                                    RFLU_BXV_ConvertCvPrim2Cons
  USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                               RFLU_ConvertCvPrim2Cons, &
                               RFLU_ScalarConvertCvCons2Prim, &
                               RFLU_ScalarConvertCvPrim2Cons

  USE RFLU_ModDifferentiationBFaces
  USE RFLU_ModEntropyFixes  
  USE RFLU_ModRindStates, ONLY: RFLU_SetRindStateFarfieldPerf
  
  USE ModInterfaces, ONLY: BcondFarfieldPerf, &
                           BcondInflowPerf, &
                           BcondOutflowPerf, &
                           MixtPerf_C_DGP, & 
                           MixtPerf_C_GHoVm2, &
                           MixtPerf_D_PRT, &                                                      
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_Ho_CpTUVW, &
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_P_DRT, &
                           MixtPerf_R_CpG, &                           
                           MixtPerf_R_M, &
                           MixtPerf_T_DPR, & 
                           RFLU_DecideNeedBGradFace, &
                           RFLU_DecidePrint, &
                           RFLU_PrintLocInfo
  
  USE RFLU_ModJWL

   USE ModInterfacesSpecies, ONLY: SPEC_GetSpeciesIndex
   USE ModSpecies,           ONLY: t_spec_input, &
                                   t_spec_type     


  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModNSCBC.F90,v $ $Revision: 1.3 $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_NSCBC_CompFirstPatchFlux, &
            RFLU_NSCBC_CompPatchFlux, &
            RFLU_NSCBC_CompSecondPatchFlux, &
            RFLU_NSCBC_CompRhs, &
            RFLU_NSCBC_DecideHaveNSCBC, &
            RFLU_NSCBC_InitFF, &            
            RFLU_NSCBC_InitIFTotAng, &
            RFLU_NSCBC_InitIFVelTemp, &
            RFLU_NSCBC_InitIJ, &             
            RFLU_NSCBC_InitNSWHeat, &
            RFLU_NSCBC_InitNSWTemp, &    
            RFLU_NSCBC_InitOF, &            
            RFLU_NSCBC_InitSW                    

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: RFLU_NSCBC_CompRhsFF, &
             RFLU_NSCBC_CompRhsIFTotAng, &
             RFLU_NSCBC_CompRhsIFVelTemp, &             
             RFLU_NSCBC_CompRhsIJ, & 
             RFLU_NSCBC_CompRhsNSWHeat, &
             RFLU_NSCBC_CompRhsNSWTemp, &             
             RFLU_NSCBC_CompRhsOF, &             
             RFLU_NSCBC_CompRhsSW

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
 
 
  
  
  
! ******************************************************************************
!
! Purpose: Compute patch flux with first-order accuracy.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompFirstPatchFlux(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: decidePrintFlag
  INTEGER, PARAMETER :: MAX_INOUTFLOW_LOCS = 10
  INTEGER :: c1,ifl,indGs,nLocs,loc(MAX_INOUTFLOW_LOCS,MIN_VAL:MAX_VAL)
  REAL(RFREAL) :: ah,cpgas,de,dissFact,dp,dq,dr,du,dv1,dv2,dv5,dw,efc,epsentr, &
                  fs,ggas,Hh,Hl,Hr,iCmassRef,iCpRef,irl,irr,l1,l2,l5,nm,nx,ny, &
                  nz,pl,pr,pRef,qh,ql,qr,rgas,rh,rl,rr,rRef,sh,tl,tr,t1,t2,t3, &
                  t5,uh,ul,ur,vh,vl,vr,vRef,wh,wl,wr,wt,xc,yc,zc,Yproducts,al,el,&
                 ar,er,irh,eh,b,c,d,f

  REAL(RFREAL) :: flxConv(5),flxDiss(5),flx(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: mfMixt
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv,rhs,vals,sd
  TYPE(t_global), POINTER :: global
 
  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecProducts
  TYPE(t_spec_input), POINTER :: pSpecInput
!FRED - Added JWL EOS capability - 9/11/15


! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_CompFirstPatchFlux',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  indGs    = pRegion%grid%indGs
  epsentr  = pRegion%mixtInput%epsentr
  dissFact = pRegion%mixtInput%dissFact

  cv   => pRegion%mixt%cv
  dv   => pRegion%mixt%dv
  gv   => pRegion%mixt%gv
  rhs  => pRegion%mixt%rhs
  sd   => pRegion%mixt%sd

  IF ( global%specUsed .EQV. .TRUE. ) THEN
  pSpecInput => pRegion%specInput
  END IF

  nLocs = 0

  decidePrintFlag = RFLU_DecidePrint(global)

  pRef = global%refPressure
  rRef = global%refDensity
  vRef = global%refVelocity

  iCpRef    = 2.0_RFREAL/(rRef*vRef*vRef)
  iCmassRef = 1.0_RFREAL/(rRef*vRef)

! ******************************************************************************
! Define constants
! ******************************************************************************

  cpgas = global%refCp
  ggas  = global%refGamma
  rgas  = MixtPerf_R_CpG(cpgas,ggas)

! ******************************************************************************
! Main Body
! ******************************************************************************

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    nx = pPatch%fn(XCOORD,ifl)
    ny = pPatch%fn(YCOORD,ifl)
    nz = pPatch%fn(ZCOORD,ifl)
    nm = pPatch%fn(XYZMAG,ifl)

    xc = pPatch%fc(XCOORD,ifl)
    yc = pPatch%fc(YCOORD,ifl)
    zc = pPatch%fc(ZCOORD,ifl)

    fs = pPatch%gs(indGs*ifl)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    rl  = cv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl

    ul  = cv(CV_MIXT_XMOM,c1)*irl
    vl  = cv(CV_MIXT_YMOM,c1)*irl
    wl  = cv(CV_MIXT_ZMOM,c1)*irl
    pl  = dv(DV_MIXT_PRES,c1)

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
     IF (global%specUsed .EQV. .TRUE.) THEN
       IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
          scalarConvFlag = .FALSE.
       ELSE
          scalarConvFlag = .TRUE.
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
       END IF

    iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
    Yproducts = pRegion%spec%cv(iCvSpecProducts,c1)
    CALL RFLU_JWL_ComputeEnergyMixt(pRegion,c1,ggas,rgas,pl,rl,Yproducts,al,el,tl)
    Hl = el + 0.5_RFREAL * (ul*ul + vl*vl + wl*wl) + (pl/rl)

    IF (scalarConvFlag .EQV. .TRUE.) THEN
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
    END IF
   END IF !specUsed
    ELSE
       tl = MixtPerf_T_DPR(rl,pl,rgas)
       Hl = MixtPerf_Ho_CpTUVW(cpgas,tl,ul,vl,wl)
    END IF

    ql = ul*nx + vl*ny + wl*nz - fs

! ------------------------------------------------------------------------------
!   Right state : get right values from boundary arrays cv and dv
! ------------------------------------------------------------------------------

    rr   = pPatch%mixt%cv(CV_MIXT_DENS,ifl)
    irr  = 1.0_RFREAL/rr

    ur   = pPatch%mixt%cv(CV_MIXT_XMOM,ifl)*irr
    vr   = pPatch%mixt%cv(CV_MIXT_YMOM,ifl)*irr
    wr   = pPatch%mixt%cv(CV_MIXT_ZMOM,ifl)*irr
    pr   = pPatch%mixt%dv(DV_MIXT_PRES,ifl)

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
     IF (global%specUsed .EQV. .TRUE.) THEN
       IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
          scalarConvFlag = .FALSE.
       ELSE
          scalarConvFlag = .TRUE.
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
       END IF
    Yproducts = pRegion%spec%cv(iCvSpecProducts,c1)!Fred - see above

    CALL RFLU_JWL_ComputeEnergyMixt(pRegion,c1,ggas,rgas,pr,rr,Yproducts,ar,er,tr)
    Hr = er + 0.5_RFREAL * (ur*ur + vr*vr + wr*wr) + (pr/rr)

    IF (scalarConvFlag .EQV. .TRUE.) THEN
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
    END IF
   END IF !specUsed
    ELSE
    tr = MixtPerf_T_DPR(rr,pr,rgas)
    Hr = MixtPerf_Ho_CpTUVW(cpgas,tr,ur,vr,wr)
    END IF

    qr = ur*nx + vr*ny + wr*nz - fs

! ==============================================================================
!   Compute Roe-averaged face variables (h for hat)
! ==============================================================================

    rh = SQRT(rl*rr)
    wt = rl/(rl + rh)

    uh = wt*ul + (1.0_RFREAL-wt)*ur
    vh = wt*vl + (1.0_RFREAL-wt)*vr
    wh = wt*wl + (1.0_RFREAL-wt)*wr
    Hh = wt*Hl + (1.0_RFREAL-wt)*Hr

    qh = uh*nx + vh*ny + wh*nz - fs
    sh = 0.5_RFREAL*(uh*uh + vh*vh + wh*wh)

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
     IF (global%specUsed .EQV. .TRUE.) THEN
       IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
          scalarConvFlag = .FALSE.
       ELSE
          scalarConvFlag = .TRUE.
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
       END IF

    irh = 1.0_RFREAL/rh
    Yproducts = pRegion%spec%cv(iCvSpecProducts,c1) !Fred - see above
    eh = wt*el + (1.0_RFREAL-wt)*er
    CALL RFLU_JWL_ComputePressureMixt(pRegion,c1,ggas,rgas,eh,rh,Yproducts,ah,b,d,c,f)

    IF (scalarConvFlag .EQV. .TRUE.) THEN
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
    END IF
   END IF !specUsed
    ELSE
    ah = MixtPerf_C_GHoVm2(ggas,Hh,2.0_RFREAL*sh) ! NOTE factor of 2
    END IF

! ==============================================================================
!   Compute eigenvalues
! ==============================================================================

    l1 = ABS(qh - ah)
    l2 = ABS(qh     )
    l5 = ABS(qh + ah)

! ==============================================================================
!   Compute entropy fix
! ==============================================================================

    efc = epsentr*l5
    l1  = EntropyFixHartenHyman(l1,efc)
    l2  = EntropyFixHartenHyman(l2,efc)
    l5  = EntropyFixHartenHyman(l5,efc)

! ==============================================================================
!   Compute wavestrengths
! ==============================================================================

    dr = rr - rl
    du = ur - ul
    de = vr - vl ! use de to denote dv because dv is already used
    dw = wr - wl
    dp = pr - pl
    dq = du*nx + de*ny + dw*nz ! Note: no gs contribution

    dv1 = (dp - rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv5 = (dp + rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv2 = dr - dp/(ah*ah)

    t1 = l1*dv1
    t2 = l2*dv2
    t3 = l2*rh
    t5 = l5*dv5

! =============================================================================
!   Compute fluxes
! =============================================================================

! ------------------------------------------------------------------------------
!   Central part
! ------------------------------------------------------------------------------

    flxConv(1) = 0.5_RFREAL*(ql*rl            + qr*rr           )*nm
    flxConv(2) = 0.5_RFREAL*(ql*rl*ul + pl*nx + qr*rr*ur + pr*nx)*nm
    flxConv(3) = 0.5_RFREAL*(ql*rl*vl + pl*ny + qr*rr*vr + pr*ny)*nm
    flxConv(4) = 0.5_RFREAL*(ql*rl*wl + pl*nz + qr*rr*wr + pr*nz)*nm
    flxConv(5) = 0.5_RFREAL*(ql*rl*Hl + pl*fs + qr*rr*Hr + pr*fs)*nm

! ------------------------------------------------------------------------------
!   Dissipative part
! ------------------------------------------------------------------------------

    flxDiss(1) = 0.5_RFREAL*dissFact*(t1            + t2                                   + t5           )*nm
    flxDiss(2) = 0.5_RFREAL*dissFact*(t1*(uh-nx*ah) + t2*uh + t3*(du-nx*dq)                + t5*(uh+nx*ah))*nm
    flxDiss(3) = 0.5_RFREAL*dissFact*(t1*(vh-ny*ah) + t2*vh + t3*(de-ny*dq)                + t5*(vh+ny*ah))*nm
    flxDiss(4) = 0.5_RFREAL*dissFact*(t1*(wh-nz*ah) + t2*wh + t3*(dw-nz*dq)                + t5*(wh+nz*ah))*nm
    flxDiss(5) = 0.5_RFREAL*dissFact*(t1*(Hh-qh*ah) + t2*sh + t3*(uh*du+vh*de+wh*dw-qh*dq) + t5*(Hh+qh*ah))*nm

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    rhs(CV_MIXT_DENS,c1) = rhs(CV_MIXT_DENS,c1) + (flxConv(1) - flxDiss(1))
    rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + (flxConv(2) - flxDiss(2))
    rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + (flxConv(3) - flxDiss(3))
    rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + (flxConv(4) - flxDiss(4))
    rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + (flxConv(5) - flxDiss(5))

    flx(1) = flxConv(1) - flxDiss(1)
    flx(2) = flxConv(2) - flxDiss(2)
    flx(3) = flxConv(3) - flxDiss(3)
    flx(4) = flxConv(4) - flxDiss(4)
    flx(5) = flxConv(5) - flxDiss(5)

    pPatch%cp(ifl)          = iCpRef*(0.5_RFREAL*(pl + pr) - pRef)
    pPatch%cmass(ifl)       = iCmassRef*flx(1)/nm
    pPatch%cmom(XCOORD,ifl) = iCpRef*0.5_RFREAL*(ul+ur)*flx(1)/nm
    pPatch%cmom(YCOORD,ifl) = iCpRef*0.5_RFREAL*(vl+vr)*flx(1)/nm
    pPatch%cmom(ZCOORD,ifl) = iCpRef*0.5_RFREAL*(wl+wr)*flx(1)/nm

    IF ( pRegion%irkStep == 1 ) THEN
      global%massIn  = global%massIn  - MIN(flx(1),0.0_RFREAL)
      global%massOut = global%massOut + MAX(flx(1),0.0_RFREAL)
    END IF ! pRegion

    IF ( (global%checkLevel == CHECK_HIGH) .AND. &
         (global%verbLevel == VERBOSE_HIGH) .AND. &
         (global%myProcid == MASTERPROC) .AND. &
         (decidePrintFlag .EQV. .TRUE.) ) THEN
      IF ( pPatch%bcType == BC_OUTFLOW ) THEN
        IF ( flx(1) < 0.0_RFREAL ) THEN
          nLocs = nLocs + 1

          IF ( nLocs == 1 ) THEN
            global%warnCounter = global%warnCounter + 1
  
            WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, &
                  '*** WARNING *** Inflow detected at outflow boundary!'
            WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                           pRegion%iRegionGlobal
            WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME, &
                    'RFLU_NSCBC_CompFirstPatchFlux'
            IF ( global%flowType == FLOW_UNSTEADY ) THEN
              WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME, &
                                                  'Current time:', &
                                                  global%currentTime
            ELSE
              WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME, &
                                               'Current iteration number:', &
                                               global%currentIter
            END IF ! global%flowType

            WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME, &
                                           'Runge-Kutta stage:', &
                                           pRegion%irkStep
          END IF ! nLocs

          IF ( nLocs <= MAX_INOUTFLOW_LOCS ) THEN
            loc(nLocs,MIN_VAL:MAX_VAL) = c1
          END IF ! nLocs
        END IF ! flx(1)
      END IF ! pPatch%bcType
    END IF ! global%checkLevel
  END DO ! ifl

! ------------------------------------------------------------------------------
! Write info on inflow at outflow boundary
! ------------------------------------------------------------------------------

  IF ( (global%checkLevel == CHECK_HIGH) .AND. &
       (global%verbLevel == VERBOSE_HIGH) .AND. &
       (global%myProcid == MASTERPROC) .AND. &
       (decidePrintFlag .EQV. .TRUE.) .AND. &
       (nLocs > 0) ) THEN
    IF ( pPatch%bcType == BC_OUTFLOW ) THEN
      IF ( nLocs > MAX_INOUTFLOW_LOCS ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, &
               'Only wrote the first',MAX_INOUTFLOW_LOCS,'of',nLocs, &
               'outflow faces with inflow.'
        CALL RFLU_PrintLocInfo(pRegion,loc,MAX_INOUTFLOW_LOCS, &
                             LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
      ELSE
        CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, &
                               LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
      END IF ! nLocs
    END IF ! pPatch%bcType
  END IF ! nLocs

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompFirstPatchFlux
  









 
! ******************************************************************************
!
! Purpose: Compute patch flux.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompPatchFlux(pRegion,pPatch)

  IMPLICIT NONE
 
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: decidePrintFlag
  INTEGER, PARAMETER :: MAX_INOUTFLOW_LOCS = 10
  INTEGER :: c1,ifl,indGs,nLocs,loc(MAX_INOUTFLOW_LOCS,MIN_VAL:MAX_VAL)
  REAL(RFREAL) :: ah,cpgas,de,dissFact,dp,dq,dr,du,dv1,dv2,dv5,dw,efc,epsentr, &
                  fs,ggas,Hh,Hl,Hr,iCmassRef,iCpRef,irl,irr,l1,l2,l5,nm,nx,ny, &
                  nz,pl,pr,pRef,qh,ql,qr,rgas,rh,rl,rr,rRef,sh,tl,tr,t1,t2,t3, &
                  t5,uh,ul,ur,vh,vl,vr,vRef,wh,wl,wr,wt,xc,yc,zc,Yproducts,eh, &
                  b,c,d,f,irh,ar,er

  REAL(RFREAL) :: flxConv(5),flxDiss(5),flx(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: mfMixt
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv,rhs,vals,sd
  TYPE(t_global), POINTER :: global

  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecProducts
  TYPE(t_spec_input), POINTER :: pSpecInput
  !Fred - Added JWL EOS Capabilities - 9/11/15

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_CompPatchFlux',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  indGs    = pRegion%grid%indGs
  epsentr  = pRegion%mixtInput%epsentr
  dissFact = pRegion%mixtInput%dissFact

  cv   => pRegion%mixt%cv
  dv   => pRegion%mixt%dv
  gv   => pRegion%mixt%gv
  rhs  => pRegion%mixt%rhs
  sd   => pRegion%mixt%sd

  IF ( global%specUsed .EQV. .TRUE. ) THEN
  pSpecInput => pRegion%specInput
  END IF


  nLocs = 0

  decidePrintFlag = RFLU_DecidePrint(global)

  pRef = global%refPressure
  rRef = global%refDensity
  vRef = global%refVelocity

  iCpRef    = 2.0_RFREAL/(rRef*vRef*vRef)
  iCmassRef = 1.0_RFREAL/(rRef*vRef)

! ******************************************************************************
! Define constants
! ******************************************************************************

  cpgas = global%refCp
  ggas  = global%refGamma
  rgas  = MixtPerf_R_CpG(cpgas,ggas)

! ******************************************************************************
! Main Body
! ******************************************************************************

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    nx = pPatch%fn(XCOORD,ifl)
    ny = pPatch%fn(YCOORD,ifl)
    nz = pPatch%fn(ZCOORD,ifl)
    nm = pPatch%fn(XYZMAG,ifl)

    xc = pPatch%fc(XCOORD,ifl)
    yc = pPatch%fc(YCOORD,ifl)
    zc = pPatch%fc(ZCOORD,ifl)

    fs = pPatch%gs(indGs*ifl)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!  get right values from boundary arrays cv and dv
! ------------------------------------------------------------------------------

    rr   = pPatch%mixt%cv(CV_MIXT_DENS,ifl)
    irr  = 1.0_RFREAL/rr

    ur   = pPatch%mixt%cv(CV_MIXT_XMOM,ifl)*irr
    vr   = pPatch%mixt%cv(CV_MIXT_YMOM,ifl)*irr
    wr   = pPatch%mixt%cv(CV_MIXT_ZMOM,ifl)*irr
    pr   = pPatch%mixt%dv(DV_MIXT_PRES,ifl)

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
     IF (global%specUsed .EQV. .TRUE.) THEN
       IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
          scalarConvFlag = .FALSE.
       ELSE
          scalarConvFlag = .TRUE.
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
       END IF

    iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
    Yproducts = pRegion%spec%cv(iCvSpecProducts,c1) !Fred - removed *irr

    Yproducts = pRegion%spec%cv(iCvSpecProducts,c1) !Fred - removed *irr

    CALL RFLU_JWL_ComputeEnergyMixt(pRegion,c1,ggas,rgas,pr,rr,Yproducts,ar,er,tr)
    Hr = er + 0.5_RFREAL * (ur*ur + vr*vr + wr*wr) + (pr/rr)

    IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
    END IF
   END IF !specUsed
    ELSE
    tr = MixtPerf_T_DPR(rr,pr,rgas)
    Hr = MixtPerf_Ho_CpTUVW(cpgas,tr,ur,vr,wr)
   END IF

    qr = ur*nx + vr*ny + wr*nz - fs

! ------------------------------------------------------------------------------
!  get Left state also from boundary arrays cv and dv
! ------------------------------------------------------------------------------

    rl  = rr 

    ul  = ur
    vl  = vr
    wl  = wr
    pl  = pr

    tl = tr
    Hl = Hr

    ql = qr 

! ==============================================================================
!   Compute Roe-averaged face variables (h for hat)
! ==============================================================================

    rh = SQRT(rl*rr)
    wt = rl/(rl + rh)

    uh = wt*ul + (1.0_RFREAL-wt)*ur
    vh = wt*vl + (1.0_RFREAL-wt)*vr
    wh = wt*wl + (1.0_RFREAL-wt)*wr
    Hh = wt*Hl + (1.0_RFREAL-wt)*Hr

    qh = uh*nx + vh*ny + wh*nz - fs
    sh = 0.5_RFREAL*(uh*uh + vh*vh + wh*wh)


    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
      IF (global%specUsed .EQV. .TRUE.) THEN
       IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
          scalarConvFlag = .FALSE.
       ELSE
          scalarConvFlag = .TRUE.
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
       END IF

    irh = 1.0_RFREAL/rh
    Yproducts = pRegion%spec%cv(iCvSpecProducts,c1) !Fred - see above
    CALL RFLU_JWL_ComputePressureMixt(pRegion,c1,ggas,rgas,er,rh,Yproducts,ah,b,d,c,f)

    IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
    END IF
   END IF !specUsed
    ELSE
    ah = MixtPerf_C_GHoVm2(ggas,Hh,2.0_RFREAL*sh) ! NOTE factor of 2
    END IF
! ==============================================================================
!   Compute eigenvalues
! ==============================================================================

    l1 = ABS(qh - ah)
    l2 = ABS(qh     )
    l5 = ABS(qh + ah)

! ==============================================================================
!   Compute entropy fix
! ==============================================================================

    efc = epsentr*l5
    l1  = EntropyFixHartenHyman(l1,efc)
    l2  = EntropyFixHartenHyman(l2,efc)
    l5  = EntropyFixHartenHyman(l5,efc)

! ==============================================================================
!   Compute wavestrengths
! ==============================================================================

    dr = rr - rl
    du = ur - ul
    de = vr - vl ! use de to denote dv because dv is already used
    dw = wr - wl
    dp = pr - pl
    dq = du*nx + de*ny + dw*nz ! Note: no gs contribution

    dv1 = (dp - rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv5 = (dp + rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv2 = dr - dp/(ah*ah)

    t1 = l1*dv1
    t2 = l2*dv2
    t3 = l2*rh
    t5 = l5*dv5

! =============================================================================
!   Compute fluxes
! =============================================================================

! ------------------------------------------------------------------------------
!   Central part
! ------------------------------------------------------------------------------

    flxConv(1) = 0.5_RFREAL*(ql*rl            + qr*rr           )*nm
    flxConv(2) = 0.5_RFREAL*(ql*rl*ul + pl*nx + qr*rr*ur + pr*nx)*nm
    flxConv(3) = 0.5_RFREAL*(ql*rl*vl + pl*ny + qr*rr*vr + pr*ny)*nm
    flxConv(4) = 0.5_RFREAL*(ql*rl*wl + pl*nz + qr*rr*wr + pr*nz)*nm
    flxConv(5) = 0.5_RFREAL*(ql*rl*Hl + pl*fs + qr*rr*Hr + pr*fs)*nm

! ------------------------------------------------------------------------------
!   Dissipative part
! ------------------------------------------------------------------------------

    flxDiss(1) = 0.5_RFREAL*dissFact*(t1            + t2                                   + t5           )*nm
    flxDiss(2) = 0.5_RFREAL*dissFact*(t1*(uh-nx*ah) + t2*uh + t3*(du-nx*dq)                + t5*(uh+nx*ah))*nm
    flxDiss(3) = 0.5_RFREAL*dissFact*(t1*(vh-ny*ah) + t2*vh + t3*(de-ny*dq)                + t5*(vh+ny*ah))*nm
    flxDiss(4) = 0.5_RFREAL*dissFact*(t1*(wh-nz*ah) + t2*wh + t3*(dw-nz*dq)                + t5*(wh+nz*ah))*nm
    flxDiss(5) = 0.5_RFREAL*dissFact*(t1*(Hh-qh*ah) + t2*sh + t3*(uh*du+vh*de+wh*dw-qh*dq) + t5*(Hh+qh*ah))*nm

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    rhs(CV_MIXT_DENS,c1) = rhs(CV_MIXT_DENS,c1) + (flxConv(1) - flxDiss(1))
    rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + (flxConv(2) - flxDiss(2))
    rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + (flxConv(3) - flxDiss(3))
    rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + (flxConv(4) - flxDiss(4))
    rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + (flxConv(5) - flxDiss(5))

    flx(1) = flxConv(1) - flxDiss(1)
    flx(2) = flxConv(2) - flxDiss(2)
    flx(3) = flxConv(3) - flxDiss(3)
    flx(4) = flxConv(4) - flxDiss(4)
    flx(5) = flxConv(5) - flxDiss(5)

    pPatch%cp(ifl)          = iCpRef*(0.5_RFREAL*(pl + pr) - pRef)
    pPatch%cmass(ifl)       = iCmassRef*flx(1)/nm
    pPatch%cmom(XCOORD,ifl) = iCpRef*0.5_RFREAL*(ul+ur)*flx(1)/nm
    pPatch%cmom(YCOORD,ifl) = iCpRef*0.5_RFREAL*(vl+vr)*flx(1)/nm
    pPatch%cmom(ZCOORD,ifl) = iCpRef*0.5_RFREAL*(wl+wr)*flx(1)/nm

    IF ( pRegion%irkStep == 1 ) THEN
      global%massIn  = global%massIn  - MIN(flx(1),0.0_RFREAL)
      global%massOut = global%massOut + MAX(flx(1),0.0_RFREAL)
    END IF ! pRegion

    IF ( (global%checkLevel == CHECK_HIGH) .AND. &
         (global%verbLevel == VERBOSE_HIGH) .AND. &
         (global%myProcid == MASTERPROC) .AND. &
         (decidePrintFlag .EQV. .TRUE.) ) THEN
      IF ( pPatch%bcType == BC_OUTFLOW ) THEN
        IF ( flx(1) < 0.0_RFREAL ) THEN
          nLocs = nLocs + 1

          IF ( nLocs == 1 ) THEN
            global%warnCounter = global%warnCounter + 1

            WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, &
                  '*** WARNING *** Inflow detected at outflow boundary!'
            WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                             pRegion%iRegionGlobal
            WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME, &
                    'RFLU_NSCBC_CompPatchFlux'

            IF ( global%flowType == FLOW_UNSTEADY ) THEN
              WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME, &
                                                'Current time:', &
                                                  global%currentTime
            ELSE
              WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME, &
                                             'Current iteration number:', &
                                               global%currentIter
            END IF ! global%flowType

            WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME, &
                                           'Runge-Kutta stage:', &
                                           pRegion%irkStep
          END IF ! nLocs

          IF ( nLocs <= MAX_INOUTFLOW_LOCS ) THEN
            loc(nLocs,MIN_VAL:MAX_VAL) = c1
          END IF ! nLocs
        END IF ! flx(1)
      END IF ! pPatch%bcType
    END IF ! global%checkLevel
  END DO ! ifl

! ------------------------------------------------------------------------------
! Write info on inflow at outflow boundary
! ------------------------------------------------------------------------------

  IF ( (global%checkLevel == CHECK_HIGH) .AND. &
       (global%verbLevel == VERBOSE_HIGH) .AND. &
       (global%myProcid == MASTERPROC) .AND. &
       (decidePrintFlag .EQV. .TRUE.) .AND. &
       (nLocs > 0) ) THEN
    IF ( pPatch%bcType == BC_OUTFLOW ) THEN
      IF ( nLocs > MAX_INOUTFLOW_LOCS ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, &
               'Only wrote the first',MAX_INOUTFLOW_LOCS,'of',nLocs, &
               'outflow faces with inflow.'
        CALL RFLU_PrintLocInfo(pRegion,loc,MAX_INOUTFLOW_LOCS, &
                               LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
      ELSE
        CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, &
                               LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
      END IF ! nLocs
    END IF ! pPatch%bcType
  END IF ! nLocs

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompPatchFlux
  







! ******************************************************************************
!
! Purpose: Compute patch fluxes with second-order accuracy.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompSecondPatchFlux(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: decidePrintFlag
  INTEGER, PARAMETER :: MAX_INOUTFLOW_LOCS = 10
  INTEGER :: c1,ifl,indGs,nLocs,loc(MAX_INOUTFLOW_LOCS,MIN_VAL:MAX_VAL)
  REAL(RFREAL) :: ah,cpgas,de,dissFact,dp,dq,dr,dx,dy,dz,du,dv1,dv2,dv5,dw,efc, &
                  epsentr,fs,ggas,Hh,Hl,Hr,iCmassRef,iCpRef,irl,irr,l1,l2,l5,  &
                  nm,nx,ny,nz,pl,pr,pRef,qh,ql,qr,rgas,rh,rl,rr,rRef,sh,tl,tr, &
                  t1,t2,t3,t5,uh,ul,ur,vh,vl,vr,vRef,wh,wl,wr,wt,xc,yc,zc,b,c,d,f,&
                 irh,eh,Yproducts,al,el,ar,er
  REAL(RFREAL) :: flxConv(5),flxDiss(5),flx(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: mfMixt
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv,rhs,vals,sd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
  TYPE(t_global), POINTER :: global

  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecProducts
  TYPE(t_spec_input), POINTER :: pSpecInput
!Fred - Added JWL EOS Capabilities - 9/11/15

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_CompSecondPatchFlux',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  indGs    = pRegion%grid%indGs
  epsentr  = pRegion%mixtInput%epsentr
  dissFact = pRegion%mixtInput%dissFact

  cv   => pRegion%mixt%cv
  dv   => pRegion%mixt%dv
  gv   => pRegion%mixt%gv
  grad => pRegion%mixt%gradCell
  rhs  => pRegion%mixt%rhs
  sd   => pRegion%mixt%sd

  IF ( global%specUsed .EQV. .TRUE. ) THEN
  pSpecInput => pRegion%specInput
  END IF

  nLocs = 0

  decidePrintFlag = RFLU_DecidePrint(global)

  pRef = global%refPressure
  rRef = global%refDensity
  vRef = global%refVelocity

  iCpRef    = 2.0_RFREAL/(rRef*vRef*vRef)
  iCmassRef = 1.0_RFREAL/(rRef*vRef)

! ******************************************************************************
! Define constants
! ******************************************************************************

  cpgas = global%refCp
  ggas  = global%refGamma
  rgas  = MixtPerf_R_CpG(cpgas,ggas)

! ******************************************************************************
! Main Body
! ******************************************************************************

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    nx = pPatch%fn(XCOORD,ifl)
    ny = pPatch%fn(YCOORD,ifl)
    nz = pPatch%fn(ZCOORD,ifl)
    nm = pPatch%fn(XYZMAG,ifl)

    xc = pPatch%fc(XCOORD,ifl)
    yc = pPatch%fc(YCOORD,ifl)
    zc = pPatch%fc(ZCOORD,ifl)

    fs = pPatch%gs(indGs*ifl)

! ==============================================================================
!   Compute left and right states
! ==============================================================================

! ------------------------------------------------------------------------------
!   Left state
! ------------------------------------------------------------------------------

    rl  = cv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl

    ul  = cv(CV_MIXT_XMOM,c1)*irl
    vl  = cv(CV_MIXT_YMOM,c1)*irl
    wl  = cv(CV_MIXT_ZMOM,c1)*irl
    pl  = dv(DV_MIXT_PRES,c1)

    dx  = xc - pRegion%grid%cofg(XCOORD,c1)
    dy  = yc - pRegion%grid%cofg(YCOORD,c1)
    dz  = zc - pRegion%grid%cofg(ZCOORD,c1)

    rl = rl + grad(XCOORD,GRC_MIXT_DENS,c1)*dx &
            + grad(YCOORD,GRC_MIXT_DENS,c1)*dy &
            + grad(ZCOORD,GRC_MIXT_DENS,c1)*dz
    ul = ul + grad(XCOORD,GRC_MIXT_XVEL,c1)*dx &
            + grad(YCOORD,GRC_MIXT_XVEL,c1)*dy &
            + grad(ZCOORD,GRC_MIXT_XVEL,c1)*dz
    vl = vl + grad(XCOORD,GRC_MIXT_YVEL,c1)*dx &
            + grad(YCOORD,GRC_MIXT_YVEL,c1)*dy &
            + grad(ZCOORD,GRC_MIXT_YVEL,c1)*dz
    wl = wl + grad(XCOORD,GRC_MIXT_ZVEL,c1)*dx &
            + grad(YCOORD,GRC_MIXT_ZVEL,c1)*dy &
            + grad(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
    pl = pl + grad(XCOORD,GRC_MIXT_PRES,c1)*dx &
            + grad(YCOORD,GRC_MIXT_PRES,c1)*dy &
            + grad(ZCOORD,GRC_MIXT_PRES,c1)*dz

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
     IF (global%specUsed .EQV. .TRUE.) THEN
       IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
          scalarConvFlag = .FALSE.
       ELSE
          scalarConvFlag = .TRUE.
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
       END IF

    iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
    Yproducts = pRegion%spec%cv(iCvSpecProducts,c1) !Fred - removed *irl

    CALL RFLU_JWL_ComputeEnergyMixt(pRegion,c1,ggas,rgas,pl,rl,Yproducts,al,el,tl)
    Hl = el + 0.5_RFREAL * (ul*ul + vl*vl + wl*wl) + (pl/rl)

    IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
    END IF
   END IF !specUsed
    ELSE
    tl = MixtPerf_T_DPR(rl,pl,rgas)
    Hl = MixtPerf_Ho_CpTUVW(cpgas,tl,ul,vl,wl)
    END IF

    ql = ul*nx + vl*ny + wl*nz - fs

! ------------------------------------------------------------------------------
!   Get right values from boundary arrays cv and dv
! ------------------------------------------------------------------------------

    rr   = pPatch%mixt%cv(CV_MIXT_DENS,ifl)
    irr  = 1.0_RFREAL/rr

    ur   = pPatch%mixt%cv(CV_MIXT_XMOM,ifl)*irr
    vr   = pPatch%mixt%cv(CV_MIXT_YMOM,ifl)*irr
    wr   = pPatch%mixt%cv(CV_MIXT_ZMOM,ifl)*irr
    pr   = pPatch%mixt%dv(DV_MIXT_PRES,ifl)

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
      IF (global%specUsed .EQV. .TRUE.) THEN
       IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
          scalarConvFlag = .FALSE.
       ELSE
          scalarConvFlag = .TRUE.
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
       END IF

    Yproducts = pRegion%spec%cv(iCvSpecProducts,c1) !Fred - see above

    CALL RFLU_JWL_ComputeEnergyMixt(pRegion,c1,ggas,rgas,pr,rr,Yproducts,ar,er,tr)
    Hr = er + 0.5_RFREAL * (ur*ur + vr*vr + wr*wr) + (pr/rr)

    IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
    END IF
   END IF !specUsed
    ELSE
    tr = MixtPerf_T_DPR(rr,pr,rgas)
    Hr = MixtPerf_Ho_CpTUVW(cpgas,tr,ur,vr,wr)
    END IF

    qr = ur*nx + vr*ny + wr*nz - fs

! ==============================================================================
!   Compute Roe-averaged face variables (h for hat)
! ==============================================================================

    rh = SQRT(rl*rr)
    wt = rl/(rl + rh)

    uh = wt*ul + (1.0_RFREAL-wt)*ur
    vh = wt*vl + (1.0_RFREAL-wt)*vr
    wh = wt*wl + (1.0_RFREAL-wt)*wr
    Hh = wt*Hl + (1.0_RFREAL-wt)*Hr

    qh = uh*nx + vh*ny + wh*nz - fs
    sh = 0.5_RFREAL*(uh*uh + vh*vh + wh*wh)

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
      IF (global%specUsed .EQV. .TRUE.) THEN
       IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
          scalarConvFlag = .FALSE.
       ELSE
          scalarConvFlag = .TRUE.
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
       END IF

    irh = 1.0_RFREAL/rh
    Yproducts = pRegion%spec%cv(iCvSpecProducts,c1) !Fred - see above
    eh = wt*el + (1.0_RFREAL-wt)*er
    CALL RFLU_JWL_ComputePressureMixt(pRegion,c1,ggas,rgas,eh,rh,Yproducts,ah,b,d,c,f)

    IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
    END IF
   END IF !specUsed
    ELSE
    ah = MixtPerf_C_GHoVm2(ggas,Hh,2.0_RFREAL*sh) ! NOTE factor of 2
    END IF

! ==============================================================================
!   Compute eigenvalues
! ==============================================================================

    l1 = ABS(qh - ah)
    l2 = ABS(qh     )
    l5 = ABS(qh + ah)

! ==============================================================================
!   Compute entropy fix
! ==============================================================================

    efc = epsentr*l5
    l1  = EntropyFixHartenHyman(l1,efc)
    l2  = EntropyFixHartenHyman(l2,efc)
    l5  = EntropyFixHartenHyman(l5,efc)

! ==============================================================================
!   Compute wavestrengths
! ==============================================================================

    dr = rr - rl
    du = ur - ul
    de = vr - vl ! use de to denote dv because dv is already used
    dw = wr - wl
    dp = pr - pl
    dq = du*nx + de*ny + dw*nz ! Note: no gs contribution

    dv1 = (dp - rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv5 = (dp + rh*ah*dq)/(2.0_RFREAL*ah*ah)
    dv2 = dr - dp/(ah*ah)

    t1 = l1*dv1
    t2 = l2*dv2
    t3 = l2*rh
    t5 = l5*dv5

! =============================================================================
!   Compute fluxes
! =============================================================================

! ------------------------------------------------------------------------------
!   Central part
! ------------------------------------------------------------------------------

    flxConv(1) = 0.5_RFREAL*(ql*rl            + qr*rr           )*nm
    flxConv(2) = 0.5_RFREAL*(ql*rl*ul + pl*nx + qr*rr*ur + pr*nx)*nm
    flxConv(3) = 0.5_RFREAL*(ql*rl*vl + pl*ny + qr*rr*vr + pr*ny)*nm
    flxConv(4) = 0.5_RFREAL*(ql*rl*wl + pl*nz + qr*rr*wr + pr*nz)*nm
    flxConv(5) = 0.5_RFREAL*(ql*rl*Hl + pl*fs + qr*rr*Hr + pr*fs)*nm

! ------------------------------------------------------------------------------
!   Dissipative part
! ------------------------------------------------------------------------------

    flxDiss(1) = 0.5_RFREAL*dissFact*(t1            + t2                                   + t5           )*nm
    flxDiss(2) = 0.5_RFREAL*dissFact*(t1*(uh-nx*ah) + t2*uh + t3*(du-nx*dq)                + t5*(uh+nx*ah))*nm
    flxDiss(3) = 0.5_RFREAL*dissFact*(t1*(vh-ny*ah) + t2*vh + t3*(de-ny*dq)                + t5*(vh+ny*ah))*nm
    flxDiss(4) = 0.5_RFREAL*dissFact*(t1*(wh-nz*ah) + t2*wh + t3*(dw-nz*dq)                + t5*(wh+nz*ah))*nm
    flxDiss(5) = 0.5_RFREAL*dissFact*(t1*(Hh-qh*ah) + t2*sh + t3*(uh*du+vh*de+wh*dw-qh*dq) + t5*(Hh+qh*ah))*nm

! ==============================================================================
!   Accumulate into residual
! ==============================================================================

    rhs(CV_MIXT_DENS,c1) = rhs(CV_MIXT_DENS,c1) + (flxConv(1) - flxDiss(1))
    rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + (flxConv(2) - flxDiss(2))
    rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + (flxConv(3) - flxDiss(3))
    rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + (flxConv(4) - flxDiss(4))
    rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + (flxConv(5) - flxDiss(5))

    flx(1) = flxConv(1) - flxDiss(1)
    flx(2) = flxConv(2) - flxDiss(2)
    flx(3) = flxConv(3) - flxDiss(3)
    flx(4) = flxConv(4) - flxDiss(4)
    flx(5) = flxConv(5) - flxDiss(5)

    pPatch%cp(ifl)          = iCpRef*(0.5_RFREAL*(pl + pr) - pRef)
    pPatch%cmass(ifl)       = iCmassRef*flx(1)/nm
    pPatch%cmom(XCOORD,ifl) = iCpRef*0.5_RFREAL*(ul+ur)*flx(1)/nm
    pPatch%cmom(YCOORD,ifl) = iCpRef*0.5_RFREAL*(vl+vr)*flx(1)/nm
    pPatch%cmom(ZCOORD,ifl) = iCpRef*0.5_RFREAL*(wl+wr)*flx(1)/nm

    IF ( pRegion%irkStep == 1 ) THEN
      global%massIn  = global%massIn  - MIN(flx(1),0.0_RFREAL)
      global%massOut = global%massOut + MAX(flx(1),0.0_RFREAL)
    END IF ! pRegion

    IF ( (global%checkLevel == CHECK_HIGH) .AND. &
         (global%verbLevel == VERBOSE_HIGH) .AND. &
         (global%myProcid == MASTERPROC) .AND. &
         (decidePrintFlag .EQV. .TRUE.) ) THEN
      IF ( pPatch%bcType == BC_OUTFLOW ) THEN
        IF ( flx(1) < 0.0_RFREAL ) THEN
          nLocs = nLocs + 1

          IF ( nLocs == 1 ) THEN
            global%warnCounter = global%warnCounter + 1

            WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, &
                  '*** WARNING *** Reverse flow detected at global patch', &
                  pPatch%iPatchGlobal
            WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                             pRegion%iRegionGlobal
            IF ( global%flowType == FLOW_UNSTEADY ) THEN
              WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME, &
                                                  'Current time:', &
                                                  global%currentTime
            ELSE
              WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME, &
                                               'Current iteration number:', &
                                               global%currentIter
            END IF ! global%flowType

            WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME, &
                                           'Runge-Kutta stage:', &
                                           pRegion%irkStep
          END IF ! nLocs

          IF ( nLocs <= MAX_INOUTFLOW_LOCS ) THEN
            loc(nLocs,MIN_VAL:MAX_VAL) = c1
          END IF ! nLocs
        END IF ! flx(1)
      END IF ! pPatch%bcType
    END IF ! global%checkLevel
  END DO ! ifl

! ------------------------------------------------------------------------------
! Write info on inflow at outflow boundary
! ------------------------------------------------------------------------------

  IF ( (global%checkLevel == CHECK_HIGH) .AND. &
       (global%verbLevel == VERBOSE_HIGH) .AND. &
       (global%myProcid == MASTERPROC) .AND. &
       (decidePrintFlag .EQV. .TRUE.) .AND. &
       (nLocs > 0) ) THEN
    IF ( pPatch%bcType == BC_OUTFLOW ) THEN
      IF ( nLocs > MAX_INOUTFLOW_LOCS ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, &
               'Only wrote the first',MAX_INOUTFLOW_LOCS,'of',nLocs, &
               'outflow faces with inflow.'
        CALL RFLU_PrintLocInfo(pRegion,loc,MAX_INOUTFLOW_LOCS, &
                               LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
      ELSE
        CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, &
                               LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
      END IF ! nLocs
    END IF ! pPatch%bcType
  END IF ! nLocs

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompSecondPatchFlux




 

  

! ******************************************************************************
!
! Purpose: Wrapper for computing RHS.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompRhs(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,iPatch
  INTEGER :: varInfo(CV_MIXT_XVEL:CV_MIXT_TEMP)
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_CompRhs',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Compute gradients at boundary faces
! ******************************************************************************

  CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
      CALL RFLU_BXV_ConvertCvCons2Prim(pRegion,pPatch,CV_MIXT_STATE_DUVWP)

      CALL RFLU_ComputeGradBFacesWrapper(pRegion,pPatch,CV_MIXT_DENS,CV_MIXT_PRES, &
                                         GRBF_MIXT_DENS,GRBF_MIXT_PRES, &
                                         pRegion%mixt%cv,pPatch%mixt%gradFace)
  
      CALL RFLU_BXV_ConvertCvPrim2Cons(pRegion,pPatch,CV_MIXT_STATE_CONS)
    END IF ! RFLU_DecideNeedBGradFace
  END DO ! iPatch

  CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)

! ******************************************************************************
! Loop over patches and call Rhs computation routine based on type of BC
! ******************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      SELECT CASE( pPatch%bcType )
        CASE (BC_SLIPWALL)
          CALL RFLU_NSCBC_CompRhsSW(pRegion,pPatch)
        CASE (BC_NOSLIPWALL_HFLUX)
          CALL RFLU_NSCBC_CompRhsNSWHeat(pRegion,pPatch)
        CASE (BC_NOSLIPWALL_TEMP)
          CALL RFLU_NSCBC_CompRhsNSWTemp(pRegion,pPatch)
        CASE (BC_INFLOW_TOTANG)
          CALL RFLU_NSCBC_CompRhsIFTotAng(pRegion,pPatch)
        CASE (BC_INFLOW_VELTEMP)
          CALL RFLU_NSCBC_CompRhsIFVelTemp(pRegion,pPatch)
        CASE (BC_OUTFLOW)
          CALL RFLU_NSCBC_CompRhsOF(pRegion,pPatch)
        CASE (BC_FARFIELD)
          CALL RFLU_NSCBC_CompRhsFF(pRegion,pPatch)
        CASE (BC_INJECTION)
          CALL RFLU_NSCBC_CompRhsIJ(pRegion,pPatch)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,1498)
      END SELECT ! pPatch%bcType

      CALL RFLU_AXI_SourceTermsNSCBC(pRegion,pPatch)
    END IF ! pPatch%bcKind
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompRhs  
  
    
  
  
  
  
  
   
  
  
  
! ******************************************************************************
!
! Purpose: Compute RHS for farfield.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompRhsFF(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: farFieldInflow,farFieldSupersonic
  INTEGER :: c1,distrib,errorFlag,ifl,indCp,indMol
  REAL(RFREAL) :: a,aoa,aos,cp,d1,d2,d3,d4,d5,dotProduct,dndx,dndy,dndz,dpdn, &
                  dpds,dpdt,dpdx,dpdy,dpdz,drdn,drds,drdt,drdx,drdy,drdz,dsdx, &
                  dsdy,dsdz,dtdx,dtdy,dtdz,dudn,duds,dudt,dudx,dudy,dudz,dundn, &
                  dunds,dundt,dusdn,dusds,dusdt,dutdn,dutds,dutdt,dvdn,dvds, &
                  dvdt,dvdx,dvdy,dvdz,dwdn,dwds,dwdt,dwdx,dwdy,dwdz,dxdn,dxds, &
                  dxdt,dydn,dyds,dydt,dzdn,dzds,dzdt,g,L1,L2,L3,L4,L5,lambda1, &
                  lambda2,lambda3,lambda4,lambda5,mm,mf,mfNormal,nx,ny,nz,p,r, &
                  rgas,rhs1,rhs2,rhs3,rhs4,rhs5,ru,rv,rw,sx,sy,sz,tx,ty,tz,u, &
                  un,us,ut,v,w
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv,pDv,pRhs,pGv,vals
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGradFace
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_CompRhsFF',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main Body - loop over all boundary faces
! ******************************************************************************

  vals => pPatch%mixt%vals

  distrib   = pPatch%mixt%distrib

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  pCv  => pPatch%mixt%cv
  pDv  => pPatch%mixt%dv
  pGv  => pRegion%mixt%gv ! NOTE get gv from volume
  pRhs => pPatch%mixt%rhs

  pGradFace => pPatch%mixt%gradFace

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    cp   = pGv(GV_MIXT_CP ,indCp *c1)
    mm   = pGv(GV_MIXT_MOL,indMol*c1)
    rgas = MixtPerf_R_M(mm)
    g    = MixtPerf_G_CpR(cp,rgas)

    r  = pCv(CV_MIXT_DENS,ifl)
    ru = pCv(CV_MIXT_XMOM,ifl)
    rv = pCv(CV_MIXT_YMOM,ifl)
    rw = pCv(CV_MIXT_ZMOM,ifl)
    p  = pDv(DV_MIXT_PRES,ifl)
    a  = pDv(DV_MIXT_SOUN,ifl)

    u = ru/r
    v = rv/r
    w = rw/r

    drdx = pGradFace(XCOORD,GRBF_MIXT_DENS,ifl) 
    drdy = pGradFace(YCOORD,GRBF_MIXT_DENS,ifl) 
    drdz = pGradFace(ZCOORD,GRBF_MIXT_DENS,ifl) 

    dudx = pGradFace(XCOORD,GRBF_MIXT_XVEL,ifl) 
    dudy = pGradFace(YCOORD,GRBF_MIXT_XVEL,ifl) 
    dudz = pGradFace(ZCOORD,GRBF_MIXT_XVEL,ifl) 

    dvdx = pGradFace(XCOORD,GRBF_MIXT_YVEL,ifl) 
    dvdy = pGradFace(YCOORD,GRBF_MIXT_YVEL,ifl) 
    dvdz = pGradFace(ZCOORD,GRBF_MIXT_YVEL,ifl) 

    dwdx = pGradFace(XCOORD,GRBF_MIXT_ZVEL,ifl) 
    dwdy = pGradFace(YCOORD,GRBF_MIXT_ZVEL,ifl) 
    dwdz = pGradFace(ZCOORD,GRBF_MIXT_ZVEL,ifl) 

    dpdx = pGradFace(XCOORD,GRBF_MIXT_PRES,ifl) 
    dpdy = pGradFace(YCOORD,GRBF_MIXT_PRES,ifl) 
    dpdz = pGradFace(ZCOORD,GRBF_MIXT_PRES,ifl) 

! ==============================================================================
!   Computation of normal and tangential vectors
! ==============================================================================

    nx = pPatch%fn(XCOORD,ifl)
    ny = pPatch%fn(YCOORD,ifl)
    nz = pPatch%fn(ZCOORD,ifl)    

    SELECT CASE ( pRegion%mixtInput%dimens )
      CASE ( 2 )
        sx = -ny
        sy =  nx
        sz =  nz    

        tx = 0.0_RFREAL
        ty = 0.0_RFREAL
        tz = 1.0_RFREAL    
      CASE ( 3 )
        sx = -ny
        sy =  nx
        sz =  nz    

        tx = 0.0_RFREAL
        ty = 0.0_RFREAL
        tz = 1.0_RFREAL    
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,1664)
    END SELECT ! pRegion%mixtInput%dimens    

    dxdn = nx
    dydn = ny
    dzdn = nz

    dxds = sx
    dyds = sy
    dzds = sz

    dxdt = tx
    dydt = ty
    dzdt = tz

    dndx = dxdn
    dndy = dydn
    dndz = dzdn

    dsdx = dxds
    dsdy = dyds
    dsdz = dzds

    dtdx = dxdt
    dtdy = dydt
    dtdz = dzdt

! ==============================================================================
!   Determine if current face is inflow or outflow  AND if supersonic in normal 
!   face direction
! ==============================================================================

    aoa = vals(BCDAT_FARF_ATTACK,distrib*ifl)
    aos = vals(BCDAT_FARF_SLIP  ,distrib*ifl)

    dotProduct = nx*(COS(aoa)*COS(aos)) + ny*(SIN(aoa)*COS(aos)) + nz*SIN(aos)

    IF ( dotProduct < 0.0_RFREAL ) THEN
      farFieldInflow = .TRUE.
    ELSE
      farFieldInflow = .FALSE.
    END IF ! dotProduct

    mf = vals(BCDAT_FARF_MACH  ,distrib*ifl)
    mfNormal = mf*ABS(dotProduct)

    IF ( mfNormal > 1.0_RFREAL ) THEN
      farFieldSupersonic = .TRUE.
    ELSE
      farFieldSupersonic = .FALSE.
    END IF ! dotProduct

! ==============================================================================
!   Computation of normal and tangential derivatives
! ==============================================================================

    drdn = drdx*dxdn + drdy*dydn + drdz*dzdn
    drds = drdx*dxds + drdy*dyds + drdz*dzds
    drdt = drdx*dxdt + drdy*dydt + drdz*dzdt

    dudn = dudx*dxdn + dudy*dydn + dudz*dzdn
    duds = dudx*dxds + dudy*dyds + dudz*dzds
    dudt = dudx*dxdt + dudy*dydt + dudz*dzdt

    dvdn = dvdx*dxdn + dvdy*dydn + dvdz*dzdn
    dvds = dvdx*dxds + dvdy*dyds + dvdz*dzds
    dvdt = dvdx*dxdt + dvdy*dydt + dvdz*dzdt

    dwdn = dwdx*dxdn + dwdy*dydn + dwdz*dzdn
    dwds = dwdx*dxds + dwdy*dyds + dwdz*dzds
    dwdt = dwdx*dxdt + dwdy*dydt + dwdz*dzdt

    dpdn = dpdx*dxdn + dpdy*dydn + dpdz*dzdn
    dpds = dpdx*dxds + dpdy*dyds + dpdz*dzds
    dpdt = dpdx*dxdt + dpdy*dydt + dpdz*dzdt

! ==============================================================================
!   Computation of normal and tangential vector components
! ==============================================================================
 
    un = u*dndx + v*dndy + w*dndz
    us = u*dsdx + v*dsdy + w*dsdz 
    ut = u*dtdx + v*dtdy + w*dtdz 

! ==============================================================================
!   Computation of normal and tangential derivatives
! ==============================================================================

    dundn = dudn*dndx + dvdn*dndy + dwdn*dndz
    dunds = duds*dndx + dvds*dndy + dwds*dndz
    dundt = dudt*dndx + dvdt*dndy + dwdt*dndz

    dusdn = dudn*dsdx + dvdn*dsdy + dwdn*dsdz
    dusds = duds*dsdx + dvds*dsdy + dwds*dsdz
    dusdt = dudt*dsdx + dvdt*dsdy + dwdt*dsdz

    dutdn = dudn*dtdx + dvdn*dtdy + dwdn*dtdz
    dutds = duds*dtdx + dvds*dtdy + dwds*dtdz
    dutdt = dudt*dtdx + dvdt*dtdy + dwdt*dtdz

! ==============================================================================
!   Computation of eigenvalues 
! ==============================================================================

    lambda1 = un - a
    lambda2 = un
    lambda3 = un
    lambda4 = un
    lambda5 = un + a

! ==============================================================================
!   Computations of L
! ==============================================================================

    IF ( farFieldInflow .EQV. .TRUE. ) THEN 
      IF ( farFieldSupersonic .EQV. .TRUE. ) THEN
        L1 = 0.0_RFREAL
        L2 = 0.0_RFREAL
        L3 = 0.0_RFREAL
        L4 = 0.0_RFREAL
        L5 = 0.0_RFREAL
      ELSE
        L5 = lambda5*(dpdn + r*a*dundn)
        L1 = 0.0_RFREAL
        L2 = 0.0_RFREAL 
        L3 = 0.0_RFREAL 
        L4 = 0.0_RFREAL
      END IF ! farFieldSupersonic 
    ELSE
      IF ( farFieldSupersonic .EQV. .TRUE. ) THEN
        L2 = lambda2*(a*a*drdn - dpdn)
        L3 = lambda3*dusdn
        L4 = lambda4*dutdn
        L5 = lambda5*(dpdn + r*a*dundn)
        L1 = lambda1*(dpdn - r*a*dundn)
      ELSE
        L2 = lambda2*(a*a*drdn - dpdn)
        L3 = lambda3*dusdn
        L4 = lambda4*dutdn
        L5 = lambda5*(dpdn + r*a*dundn)
        L1 = 0.0_RFREAL 
      END IF ! farFieldSupersonic 
    END IF ! farFieldInflow

! ==============================================================================
!   Computations of d
! ==============================================================================

    d1 = (L2 + 0.5_RFREAL*(L1+L5))/(a*a)
    d2 = 0.5_RFREAL*(L1+L5)
    d3 = (L5-L1)/(2.0_RFREAL*r*a)
    d4 = L3
    d5 = L4 

! ==============================================================================
!   Computations of Rhs
! ==============================================================================

    rhs1 = d1 + r*dusds + us*drds + r*dutdt + ut*drdt
    rhs2 = un*d1 + r*d3 + r*un*dusds + r*us*dunds + un*us*drds &
                        + r*un*dutdt + r*ut*dundt + un*ut*drdt
    rhs3 = us*d1 + r*d4 + dpds &
         + 2.0_RFREAL*r*us*dusds + us*us*drds &
         + r*us*dutdt + r*ut*dusdt + us*ut*drdt
    rhs4 = ut*d1 + r*d5 + dpdt &
         + r*ut*dusds + r*us*dutds + ut*us*drds &
         + 2.0_RFREAL*r*ut*dutdt + ut*ut*drdt
    rhs5 = 0.5_RFREAL*(un*un+us*us+ut*ut)*d1 + d2/(g-1.0_RFREAL) &
         + r*un*d3 + r*us*d4 + r*ut*d5 &
         + us*(  0.5_RFREAL*(un*un+us*us+ut*ut)*drds &
              + r*(un*dunds+us*dusds+ut*dutds) &
              + g*dpds/(g-1.0_RFREAL) ) &
         + (0.5_RFREAL*r*(un*un+us*us+ut*ut) &
              + g*p/(g-1.0_RFREAL))*dusds &
         + ut*(  0.5_RFREAL*(un*un+us*us+ut*ut)*drdt &
              + r*(un*dundt+us*dusdt+ut*dutdt) &
              + g*dpdt/(g-1.0_RFREAL) ) &
         + (0.5_RFREAL*r*(un*un+us*us+ut*ut) &
              + g*p/(g-1.0_RFREAL))*dutdt

    pRhs(CV_MIXT_DENS,ifl) = rhs1 
    pRhs(CV_MIXT_XMOM,ifl) = rhs2*dxdn + rhs3*dxds + rhs4*dxdt 
    pRhs(CV_MIXT_YMOM,ifl) = rhs2*dydn + rhs3*dyds + rhs4*dydt 
    pRhs(CV_MIXT_ZMOM,ifl) = rhs2*dzdn + rhs3*dzds + rhs4*dzdt
    pRhs(CV_MIXT_ENER,ifl) = rhs5 

! ------------------------------------------------------------------------------
!   Adding source terms for moving reference frame formulation
! ------------------------------------------------------------------------------

    IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN

      pRhs(CV_MIXT_XMOM,ifl) = pRhs(CV_MIXT_XMOM,ifl) + r*pRegion%mvfAcc(XCOORD)
      pRhs(CV_MIXT_YMOM,ifl) = pRhs(CV_MIXT_YMOM,ifl) + r*pRegion%mvfAcc(YCOORD)
      pRhs(CV_MIXT_ZMOM,ifl) = pRhs(CV_MIXT_ZMOM,ifl) + r*pRegion%mvfAcc(ZCOORD)
      pRhs(CV_MIXT_ENER,ifl) = pRhs(CV_MIXT_ENER,ifl) &
                             + r*u*pRegion%mvfAcc(XCOORD) &
                             + r*v*pRegion%mvfAcc(YCOORD) &
                             + r*w*pRegion%mvfAcc(ZCOORD)
    END IF ! global%mvFrameFlag
  END DO ! ifl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompRhsFF


  
  
  
  

! ******************************************************************************
!
! Purpose: Compute RHS for inflow with imposed total quantities.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompRhsIFTotAng(pRegion,pPatch) 

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_CompRhsIFTotAng',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main body
! ******************************************************************************

  pPatch%mixt%rhs = 0.0_RFREAL

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompRhsIFTotAng 


  
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Compute RHS for inflow with imposed velocity and temperature.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompRhsIFVelTemp(pRegion,pPatch) 

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,distrib,errorFlag,ifl,indCp,indMol
  REAL(RFREAL) :: a,aoa,aos,BC_dudt,BC_dvdt,BC_dwdt,BC_dTdt,BC_dundt,BC_dusdt, &
                  BC_dutdt,cp,d1,d2,d3,d4,d5,dotProduct,dndx,dndy,dndz,dpdn, &
                  dpds,dpdt,dpdx,dpdy,dpdz,drdn,drds,drdt,drdx,drdy,drdz,dsdx, &
                  dsdy,dsdz,dtdx,dtdy,dtdz,dudn,duds,dudt,dudx,dudy,dudz,dundn, &
                  dunds,dundt,dusdn,dusds,dusdt,dutdn,dutds,dutdt,dvdn,dvds, &
                  dvdt,dvdx,dvdy,dvdz,dwdn,dwds,dwdt,dwdx,dwdy,dwdz,dxdn,dxds, &
                  dxdt,dydn,dyds,dydt,dzdn,dzds,dzdt,g,L1,L2,L3,L4,L5,lambda1, &
                  lambda2,lambda3,lambda4,lambda5,mm,mf,mfNormal,nx,ny,nz,p,r, &
                  rgas,rhs1,rhs2,rhs3,rhs4,rhs5,ru,rv,rw,sx,sy,sz,T,tx,ty,tz, &
                  u,un,us,ut,v,w
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv,pDv,pRhs,pGv
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGradFace
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_CompRhsIFVelTemp',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main body - loop over all boundary faces
! ******************************************************************************

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  pCv  => pPatch%mixt%cv
  pDv  => pPatch%mixt%dv
  pGv  => pRegion%mixt%gv ! NOTE get gv from volume
  pRhs => pPatch%mixt%rhs

  pGradFace => pPatch%mixt%gradFace

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    cp   = pGv(GV_MIXT_CP ,indCp *c1)
    mm   = pGv(GV_MIXT_MOL,indMol*c1)
    rgas = MixtPerf_R_M(mm)
    g    = MixtPerf_G_CpR(cp,rgas)

    r  = pCv(CV_MIXT_DENS,ifl)
    ru = pCv(CV_MIXT_XMOM,ifl)
    rv = pCv(CV_MIXT_YMOM,ifl)
    rw = pCv(CV_MIXT_ZMOM,ifl)
    p  = pDv(DV_MIXT_PRES,ifl)
    T  = pDv(DV_MIXT_TEMP,ifl)
    a  = pDv(DV_MIXT_SOUN,ifl)

    u = ru/r
    v = rv/r
    w = rw/r

    drdx = pGradFace(XCOORD,GRBF_MIXT_DENS,ifl) 
    drdy = pGradFace(YCOORD,GRBF_MIXT_DENS,ifl) 
    drdz = pGradFace(ZCOORD,GRBF_MIXT_DENS,ifl) 

    dudx = pGradFace(XCOORD,GRBF_MIXT_XVEL,ifl) 
    dudy = pGradFace(YCOORD,GRBF_MIXT_XVEL,ifl) 
    dudz = pGradFace(ZCOORD,GRBF_MIXT_XVEL,ifl) 

    dvdx = pGradFace(XCOORD,GRBF_MIXT_YVEL,ifl) 
    dvdy = pGradFace(YCOORD,GRBF_MIXT_YVEL,ifl) 
    dvdz = pGradFace(ZCOORD,GRBF_MIXT_YVEL,ifl) 

    dwdx = pGradFace(XCOORD,GRBF_MIXT_ZVEL,ifl) 
    dwdy = pGradFace(YCOORD,GRBF_MIXT_ZVEL,ifl) 
    dwdz = pGradFace(ZCOORD,GRBF_MIXT_ZVEL,ifl) 

    dpdx = pGradFace(XCOORD,GRBF_MIXT_PRES,ifl) 
    dpdy = pGradFace(YCOORD,GRBF_MIXT_PRES,ifl) 
    dpdz = pGradFace(ZCOORD,GRBF_MIXT_PRES,ifl) 

! ==============================================================================
!   Computation of boundary condition at inflow
! ==============================================================================

    BC_dudt = 0.0_RFREAL
    BC_dvdt = 0.0_RFREAL
    BC_dwdt = 0.0_RFREAL
    BC_dTdt = 0.0_RFREAL

! ==============================================================================
!   Computation of normal and tangential vectors
! ==============================================================================

    nx = pPatch%fn(XCOORD,ifl)
    ny = pPatch%fn(YCOORD,ifl)
    nz = pPatch%fn(ZCOORD,ifl)    

    SELECT CASE ( pRegion%mixtInput%dimens )
      CASE ( 2 )
        sx = -ny
        sy =  nx
        sz =  nz    

        tx = 0.0_RFREAL
        ty = 0.0_RFREAL
        tz = 1.0_RFREAL    
      CASE ( 3 )
        sx = -ny
        sy =  nx
        sz =  nz    

        tx = 0.0_RFREAL
        ty = 0.0_RFREAL
        tz = 1.0_RFREAL    
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,2094)
    END SELECT ! pRegion%mixtInput%dimens    

    dxdn = nx
    dydn = ny
    dzdn = nz

    dxds = sx
    dyds = sy
    dzds = sz

    dxdt = tx
    dydt = ty
    dzdt = tz

    dndx = dxdn
    dndy = dydn
    dndz = dzdn

    dsdx = dxds
    dsdy = dyds
    dsdz = dzds

    dtdx = dxdt
    dtdy = dydt
    dtdz = dzdt

! ==============================================================================
!   Computation of normal and tangential derivatives
! ==============================================================================

    drdn = drdx*dxdn + drdy*dydn + drdz*dzdn
    drds = drdx*dxds + drdy*dyds + drdz*dzds
    drdt = drdx*dxdt + drdy*dydt + drdz*dzdt

    dudn = dudx*dxdn + dudy*dydn + dudz*dzdn
    duds = dudx*dxds + dudy*dyds + dudz*dzds
    dudt = dudx*dxdt + dudy*dydt + dudz*dzdt

    dvdn = dvdx*dxdn + dvdy*dydn + dvdz*dzdn
    dvds = dvdx*dxds + dvdy*dyds + dvdz*dzds
    dvdt = dvdx*dxdt + dvdy*dydt + dvdz*dzdt

    dwdn = dwdx*dxdn + dwdy*dydn + dwdz*dzdn
    dwds = dwdx*dxds + dwdy*dyds + dwdz*dzds
    dwdt = dwdx*dxdt + dwdy*dydt + dwdz*dzdt

    dpdn = dpdx*dxdn + dpdy*dydn + dpdz*dzdn
    dpds = dpdx*dxds + dpdy*dyds + dpdz*dzds
    dpdt = dpdx*dxdt + dpdy*dydt + dpdz*dzdt

! ==============================================================================
!   Computation of normal and tangential vector components
! ==============================================================================
 
    un = u*dndx + v*dndy + w*dndz
    us = u*dsdx + v*dsdy + w*dsdz 
    ut = u*dtdx + v*dtdy + w*dtdz 

! ==============================================================================
!   Computation of normal and tangential derivatives
! ==============================================================================

    dundn = dudn*dndx + dvdn*dndy + dwdn*dndz
    dunds = duds*dndx + dvds*dndy + dwds*dndz
    dundt = dudt*dndx + dvdt*dndy + dwdt*dndz

    dusdn = dudn*dsdx + dvdn*dsdy + dwdn*dsdz
    dusds = duds*dsdx + dvds*dsdy + dwds*dsdz
    dusdt = dudt*dsdx + dvdt*dsdy + dwdt*dsdz

    dutdn = dudn*dtdx + dvdn*dtdy + dwdn*dtdz
    dutds = duds*dtdx + dvds*dtdy + dwds*dtdz
    dutdt = dudt*dtdx + dvdt*dtdy + dwdt*dtdz

! ==============================================================================
!   Computation of boundary condition at inflow in normal and tangential direc.
! ==============================================================================

    BC_dundt = BC_dudt*dndx + BC_dvdt*dndy + BC_dwdt*dndz
    BC_dusdt = BC_dudt*dsdx + BC_dvdt*dsdy + BC_dwdt*dsdz
    BC_dutdt = BC_dudt*dtdx + BC_dvdt*dtdy + BC_dwdt*dtdz 

! ==============================================================================
!   Computation of eigenvalues 
! ==============================================================================

    lambda1 = un - a
    lambda2 = un
    lambda3 = un
    lambda4 = un
    lambda5 = un + a

! ==============================================================================
!   Computations of L
! ==============================================================================

    IF ( pPatch%mixt%switches(BCSWI_INFLOW_TYPE) /= BCOPT_SUPERSONIC ) THEN
      IF ( pPatch%reflect == BC_REFLECTING ) THEN
        L5 = lambda5*(dpdn + r*a*dundn)
        L1 = L5 + 2.0_RFREAL*r*a*BC_dundt
        L2 = 0.5_RFREAL*(g-1.0_RFREAL)*(L1+L5) + (r*a*a/T)*BC_dTdt
        L3 = -BC_dusdt 
        L4 = -BC_dutdt
      ELSE
        L5 = lambda5*(dpdn + r*a*dundn)
        L1 = 0.0_RFREAL
        L2 = 0.0_RFREAL 
        L3 = 0.0_RFREAL 
        L4 = 0.0_RFREAL
      END IF ! pPatch%reflect
    ELSE
      L1 = 0.0_RFREAL
      L2 = 0.0_RFREAL
      L3 = 0.0_RFREAL
      L4 = 0.0_RFREAL
      L5 = 0.0_RFREAL
    END IF ! pPatch%mixt%switches 

! ==============================================================================
!   Computations of d
! ==============================================================================

    d1 = (L2 + 0.5_RFREAL*(L1+L5))/(a*a)
    d2 = 0.5_RFREAL*(L1+L5)
    d3 = (L5-L1)/(2.0_RFREAL*r*a)
    d4 = L3
    d5 = L4 

! ==============================================================================
!   Computations of Rhs
! ==============================================================================

    rhs1 = d1 + r*dusds + us*drds + r*dutdt + ut*drdt
    rhs2 = un*d1 + r*d3 + r*un*dusds + r*us*dunds + un*us*drds &
                        + r*un*dutdt + r*ut*dundt + un*ut*drdt
    rhs3 = us*d1 + r*d4 + dpds &
         + 2.0_RFREAL*r*us*dusds + us*us*drds &
         + r*us*dutdt + r*ut*dusdt + us*ut*drdt
    rhs4 = ut*d1 + r*d5 + dpdt &
         + r*ut*dusds + r*us*dutds + ut*us*drds &
         + 2.0_RFREAL*r*ut*dutdt + ut*ut*drdt
    rhs5 = 0.5_RFREAL*(un*un+us*us+ut*ut)*d1 + d2/(g-1.0_RFREAL) &
         + r*un*d3 + r*us*d4 + r*ut*d5 &
         + us*(  0.5_RFREAL*(un*un+us*us+ut*ut)*drds &
              + r*(un*dunds+us*dusds+ut*dutds) &
              + g*dpds/(g-1.0_RFREAL) ) &
         + (0.5_RFREAL*r*(un*un+us*us+ut*ut) &
              + g*p/(g-1.0_RFREAL))*dusds &
         + ut*(  0.5_RFREAL*(un*un+us*us+ut*ut)*drdt &
              + r*(un*dundt+us*dusdt+ut*dutdt) &
              + g*dpdt/(g-1.0_RFREAL) ) &
         + (0.5_RFREAL*r*(un*un+us*us+ut*ut) &
              + g*p/(g-1.0_RFREAL))*dutdt

    pRhs(CV_MIXT_DENS,ifl) = rhs1 
    pRhs(CV_MIXT_XMOM,ifl) = rhs2*dxdn + rhs3*dxds + rhs4*dxdt 
    pRhs(CV_MIXT_YMOM,ifl) = rhs2*dydn + rhs3*dyds + rhs4*dydt 
    pRhs(CV_MIXT_ZMOM,ifl) = rhs2*dzdn + rhs3*dzds + rhs4*dzdt
    pRhs(CV_MIXT_ENER,ifl) = rhs5 

    IF ( pPatch%mixt%switches(BCSWI_INFLOW_TYPE) /= BCOPT_SUPERSONIC ) THEN
      IF ( pPatch%reflect == BC_REFLECTING ) THEN
        pRhs(CV_MIXT_XMOM,ifl) = 0.0_RFREAL 
        pRhs(CV_MIXT_YMOM,ifl) = 0.0_RFREAL 
        pRhs(CV_MIXT_ZMOM,ifl) = 0.0_RFREAL 
        pRhs(CV_MIXT_ENER,ifl) = 0.0_RFREAL 
      END IF ! pPatch%reflect
    END IF ! pPatch%mixt%switches
  END DO ! ifl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompRhsIFVelTemp 
    
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Compute RHS for injection boundary.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompRhsIJ(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_CompRhsIJ',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main body
! ******************************************************************************

  pPatch%mixt%rhs = 0.0_RFREAL

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompRhsIJ


  
  
  


! ******************************************************************************
!
! Purpose: Compute RHS for no-slip wall with imposed heat flux.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompRhsNSWHeat(pRegion,pPatch) 

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_CompRhsNSWHeat',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main body
! ******************************************************************************

  pPatch%mixt%rhs = 0.0_RFREAL

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompRhsNSWHeat 


  
  
  
  
  
! ******************************************************************************
!
! Purpose: Compute RHS for no-slip wall with imposed temperature.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompRhsNSWTemp(pRegion,pPatch) 

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_CompRhsNSWTemp',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main body
! ******************************************************************************

  pPatch%mixt%rhs = 0.0_RFREAL


! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompRhsNSWTemp   
  
  
  
  
  
  

! ******************************************************************************
!
! Purpose: Compute RHS for outflow.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompRhsOF(pRegion,pPatch) 

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: bcOptType,c1,distrib,errorFlag,ifl,indCp,indMol
  REAL(RFREAL) :: a,aoa,aos,cp,d1,d2,d3,d4,d5,dotProduct,dndx,dndy,dndz,dpdn, &
                  dpds,dpdt,dpdx,dpdy,dpdz,drdn,drds,drdt,drdx,drdy,drdz,dsdx, &
                  dsdy,dsdz,dtdx,dtdy,dtdz,dudn,duds,dudt,dudx,dudy,dudz,dundn, &
                  dunds,dundt,dusdn,dusds,dusdt,dutdn,dutds,dutdt,dvdn,dvds, &
                  dvdt,dvdx,dvdy,dvdz,dwdn,dwds,dwdt,dwdx,dwdy,dwdz,dxdn,dxds, &
                  dxdt,dydn,dyds,dydt,dzdn,dzds,dzdt,g,L1,L2,L3,L4,L5,lambda1, &
                  lambda2,lambda3,lambda4,lambda5,mm,mf,mfNormal,nscbcK,nx,ny, &
                  nz,p,Pinf,r,rgas,rhs1,rhs2,rhs3,rhs4,rhs5,ru,rv,rw,sx,sy,sz, &
                  T,tx,ty,tz,u,un,us,ut,v,w,omega,H,ie,E,    &
                  mp,ra,wTNT,rho0, &
                  xc,yc,zc,smag,eps             !3D tangent vector calculations
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv,pDv,pRhs,pGv,vals
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGradFace
  TYPE(t_global), POINTER :: global
  ! Josh
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
  TYPE(t_spec_input), POINTER :: pSpecInput
  TYPE(t_spec_type), POINTER :: pSpecType
  REAL(RFREAL) :: YProducts
  INTEGER :: iCvSpecAir,iCvSpecExplosive,iCvSpecProducts,iSpec

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

! Josh
  pCvSpec => pRegion%spec%cv
  pSpecInput => pRegion%specInput


  CALL RegisterFunction(global,'RFLU_NSCBC_CompRhsOF',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main Body - loop over all boundary faces
! ******************************************************************************

  distrib   = pPatch%mixt%distrib
  bcOptType = pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE)

  IF ( bcOptType /= BCOPT_SUPERSONIC ) THEN
    vals => pPatch%mixt%vals
  END IF ! bcOptType

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  pCv  => pPatch%mixt%cv
  pDv  => pPatch%mixt%dv
  pGv  => pRegion%mixt%gv ! NOTE get gv from volume
  pRhs => pPatch%mixt%rhs

  pGradFace => pPatch%mixt%gradFace

  eps  = 1.0E-15_RFREAL !Tiny epsilon used in 3D vector calcs to 
                        !avoid divide by zero

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    cp   = pGv(GV_MIXT_CP ,indCp *c1)
    mm   = pGv(GV_MIXT_MOL,indMol*c1)
    rgas = MixtPerf_R_M(mm)
    g    = MixtPerf_G_CpR(cp,rgas)

    r  = pCv(CV_MIXT_DENS,ifl)
    ru = pCv(CV_MIXT_XMOM,ifl)
    rv = pCv(CV_MIXT_YMOM,ifl)
    rw = pCv(CV_MIXT_ZMOM,ifl)
    p  = pDv(DV_MIXT_PRES,ifl)
    a  = pDv(DV_MIXT_SOUN,ifl)

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
      ra   = pRegion%mixtInput%prepRealVal3   !Ambient farfield density
      rho0 = pRegion%mixtInput%prepRealVal24  !Reference density for HE material
      wTNT = pRegion%mixtInput%prepRealVal17  !omega from the used JWL equation
      IF (r .LE. ra) THEN
       mp = 0.0_RFREAL
      ELSE
       mp = (wTNT-(g-1.0_RFREAL))/(rho0-ra) 
      END IF

    omega = max(mp*(r-ra)+(g-1.0_RFREAL),wTNT) 

    ! TLJ added from ModJWL to be consistent; TLJ
    IF (r > rho0) THEN
       mp = 0.0_RFREAL
       omega = wTNT
    ENDIF
    END IF !Modified omega for one-equation model

    u = ru/r
    v = rv/r
    w = rw/r

    drdx = pGradFace(XCOORD,GRBF_MIXT_DENS,ifl) 
    drdy = pGradFace(YCOORD,GRBF_MIXT_DENS,ifl) 
    drdz = pGradFace(ZCOORD,GRBF_MIXT_DENS,ifl) 

    dudx = pGradFace(XCOORD,GRBF_MIXT_XVEL,ifl) 
    dudy = pGradFace(YCOORD,GRBF_MIXT_XVEL,ifl) 
    dudz = pGradFace(ZCOORD,GRBF_MIXT_XVEL,ifl) 

    dvdx = pGradFace(XCOORD,GRBF_MIXT_YVEL,ifl) 
    dvdy = pGradFace(YCOORD,GRBF_MIXT_YVEL,ifl) 
    dvdz = pGradFace(ZCOORD,GRBF_MIXT_YVEL,ifl) 

    dwdx = pGradFace(XCOORD,GRBF_MIXT_ZVEL,ifl) 
    dwdy = pGradFace(YCOORD,GRBF_MIXT_ZVEL,ifl) 
    dwdz = pGradFace(ZCOORD,GRBF_MIXT_ZVEL,ifl) 

    dpdx = pGradFace(XCOORD,GRBF_MIXT_PRES,ifl) 
    dpdy = pGradFace(YCOORD,GRBF_MIXT_PRES,ifl) 
    dpdz = pGradFace(ZCOORD,GRBF_MIXT_PRES,ifl) 

! ==============================================================================
!   Computation of normal and tangential vectors
! ==============================================================================

    nx = pPatch%fn(XCOORD,ifl)
    ny = pPatch%fn(YCOORD,ifl)
    nz = pPatch%fn(ZCOORD,ifl)    

    SELECT CASE ( pRegion%mixtInput%dimens )
      CASE ( 2 )
        sx = -ny
        sy =  nx
        sz =  nz    

        tx = 0.0_RFREAL
        ty = 0.0_RFREAL
        tz = 1.0_RFREAL    
      CASE ( 3 ) !Fred - fix this tangential vector computation...it doesn't
                 !work for 3D... 

        xc = pPatch%fc(XCOORD,ifl)
        yc = pPatch%fc(YCOORD,ifl)
        zc = pPatch%fc(ZCOORD,ifl) 
        !Face centroid coordinates for use in calculating tangent vectors
 
        smag = (xc*xc + yc*yc + (nx*xc+ny*yc)*(nx*xc+ny*yc)/(nz+eps)/(nz+eps))**0.5_RFREAL
        !To normalize s 

        IF (ABS(smag) <= eps) THEN  ! smag = 0 so the face centroid is the 
                                    ! origin and general solution is singular
         smag = (1.0_RFREAL + nx*nx/(nz+eps)/(nz+eps))**0.5_RFREAL
         sx = -1.0_RFREAL/smag
         sy = 0.0_RFREAL/smag
         sz = nx/(nz + eps)/smag
         
        ELSE  
         sx = -xc/smag 
         sy = -yc/smag
         sz = (nx*xc+ny*yc)/(nz+eps)/smag       
         !Note, if nz=0, then this will cause sz to basically be 1, 
         !which reverts back to the originally coded, non-generic method...in
         !some sense
        END IF

        tx = ny*sz-nz*sy
        ty = nz*sx-nx*sz
        tz = nx*sy-ny*sx
        !t = cross(n,s)   
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,2680)
    END SELECT ! pRegion%mixtInput%dimens    

    dxdn = nx
    dydn = ny
    dzdn = nz

    dxds = sx
    dyds = sy
    dzds = sz

    dxdt = tx
    dydt = ty
    dzdt = tz

    dndx = dxdn
    dndy = dydn
    dndz = dzdn

    dsdx = dxds
    dsdy = dyds
    dsdz = dzds

    dtdx = dxdt
    dtdy = dydt
    dtdz = dzdt

! ==============================================================================
!   Computation of normal and tangential derivatives
! ==============================================================================

    drdn = drdx*dxdn + drdy*dydn + drdz*dzdn
    drds = drdx*dxds + drdy*dyds + drdz*dzds
    drdt = drdx*dxdt + drdy*dydt + drdz*dzdt

    dudn = dudx*dxdn + dudy*dydn + dudz*dzdn
    duds = dudx*dxds + dudy*dyds + dudz*dzds
    dudt = dudx*dxdt + dudy*dydt + dudz*dzdt

    dvdn = dvdx*dxdn + dvdy*dydn + dvdz*dzdn
    dvds = dvdx*dxds + dvdy*dyds + dvdz*dzds
    dvdt = dvdx*dxdt + dvdy*dydt + dvdz*dzdt

    dwdn = dwdx*dxdn + dwdy*dydn + dwdz*dzdn
    dwds = dwdx*dxds + dwdy*dyds + dwdz*dzds
    dwdt = dwdx*dxdt + dwdy*dydt + dwdz*dzdt

    dpdn = dpdx*dxdn + dpdy*dydn + dpdz*dzdn
    dpds = dpdx*dxds + dpdy*dyds + dpdz*dzds
    dpdt = dpdx*dxdt + dpdy*dydt + dpdz*dzdt

! ==============================================================================
!   Computation of normal and tangential vector components
! ==============================================================================
 
    un = u*dndx + v*dndy + w*dndz
    us = u*dsdx + v*dsdy + w*dsdz 
    ut = u*dtdx + v*dtdy + w*dtdz 

! ==============================================================================
!   Computation of normal and tangential derivatives
! ==============================================================================

    dundn = dudn*dndx + dvdn*dndy + dwdn*dndz
    dunds = duds*dndx + dvds*dndy + dwds*dndz
    dundt = dudt*dndx + dvdt*dndy + dwdt*dndz

    dusdn = dudn*dsdx + dvdn*dsdy + dwdn*dsdz
    dusds = duds*dsdx + dvds*dsdy + dwds*dsdz
    dusdt = dudt*dsdx + dvdt*dsdy + dwdt*dsdz

    dutdn = dudn*dtdx + dvdn*dtdy + dwdn*dtdz
    dutds = duds*dtdx + dvds*dtdy + dwds*dtdz
    dutdt = dudt*dtdx + dvdt*dtdy + dwdt*dtdz

! ==============================================================================
!   Computation of eigenvalues 
! ==============================================================================

    nscbcK = pPatch%nscbcK

    IF ( bcOptType /= BCOPT_SUPERSONIC ) THEN
      Pinf = vals(BCDAT_OUTFLOW_PRESS,distrib*ifl)
    ELSE
      Pinf = p
    END IF ! bcOptType

    ! TLJ added from dupe_meso2
    ! Speed of sound for ideal gas computed above; here it
    !   is modified for JWL
    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
       ! Compute Y here
       iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput, &
                                                   'PRODUCTS')
       YProducts = pCvSpec(iCvSpecProducts,ifl)
       ie = RFLU_JWL_E_PR(pRegion,p,r,YProducts)
       a = RFLU_JWL_C_ER(pRegion,ie,r,YProducts)
    END IF

    lambda1 = un - a
    lambda2 = un
    lambda3 = un
    lambda4 = un
    lambda5 = un + a

! ==============================================================================
!   Computations of L
! ==============================================================================

    IF ( pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE) == BCOPT_SUPERSONIC ) THEN
      L2 = lambda2*(a*a*drdn - dpdn)
      L3 = lambda3*dusdn
      L4 = lambda4*dutdn
      L5 = lambda5*(dpdn + r*a*dundn)
      L1 = lambda1*(dpdn - r*a*dundn)
    ELSE
      IF ( pPatch%reflect == BC_REFLECTING ) THEN
        L2 = lambda2*(a*a*drdn - dpdn)
        L3 = lambda3*dusdn
        L4 = lambda4*dutdn
        L5 = lambda5*(dpdn + r*a*dundn)
        L1 = -L5
      ELSE ! non reflecting BC
        L2 = lambda2*(a*a*drdn - dpdn)
        L3 = lambda3*dusdn
        L4 = lambda4*dutdn
        L5 = lambda5*(dpdn + r*a*dundn)
        L1 = nscbcK*(p-Pinf)
      END IF ! pPatch%reflect
    END IF ! pPatch%mixt%switches

! ==============================================================================
!   Computations of d
! ==============================================================================

    d1 = (L2 + 0.5_RFREAL*(L1+L5))/(a*a)
    d2 = 0.5_RFREAL*(L1+L5)
    d3 = (L5-L1)/(2.0_RFREAL*r*a)
    d4 = L3
    d5 = L4 

! ==============================================================================
!   Computations of Rhs
! ==============================================================================

    rhs1 = d1 + r*dusds + us*drds + r*dutdt + ut*drdt
    rhs2 = un*d1 + r*d3 + r*un*dusds + r*us*dunds + un*us*drds &
                        + r*un*dutdt + r*ut*dundt + un*ut*drdt
    rhs3 = us*d1 + r*d4 + dpds &
         + 2.0_RFREAL*r*us*dusds + us*us*drds &
         + r*us*dutdt + r*ut*dusdt + us*ut*drdt
    rhs4 = ut*d1 + r*d5 + dpdt &
         + r*ut*dusds + r*us*dutds + ut*us*drds &
         + 2.0_RFREAL*r*ut*dutdt + ut*ut*drdt

   IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
    
    E  = pCv(CV_MIXT_ENER,ifl)/r !ie + 0.5_RFREAL*(un*un+us*us+ut*ut)
    H  = E+(p/r)

    rhs5 = (H-(a*a)/omega)*d1 + d2/omega + r*un*d3 &
           + r*us*d4 + r*ut*d5 + &
           us*(H-(a*a)/omega)*drds + r*(H+us*us)*dusds + &
           r*us*un*dunds + r*us*ut*dutds + ((omega+1.0_RFREAL)*us/omega)*dpds + &
           ut*(H-(a*a)/omega)*drdt + r*(H+ut*ut)*dutdt + &
           r*ut*un*dundt + r*ut*us*dusdt + ((omega+1.0_RFREAL)*ut/omega)*dpdt

   ELSE
    rhs5 = 0.5_RFREAL*(un*un+us*us+ut*ut)*d1 + d2/(g-1.0_RFREAL) &
         + r*un*d3 + r*us*d4 + r*ut*d5 &
         + us*(  0.5_RFREAL*(un*un+us*us+ut*ut)*drds &
              + r*(un*dunds+us*dusds+ut*dutds) &
              + g*dpds/(g-1.0_RFREAL) ) &
         + (0.5_RFREAL*r*(un*un+us*us+ut*ut) &
              + g*p/(g-1.0_RFREAL))*dusds &
         + ut*(  0.5_RFREAL*(un*un+us*us+ut*ut)*drdt &
              + r*(un*dundt+us*dusdt+ut*dutdt) &
              + g*dpdt/(g-1.0_RFREAL) ) &
         + (0.5_RFREAL*r*(un*un+us*us+ut*ut) &
              + g*p/(g-1.0_RFREAL))*dutdt
    END IF

    pRhs(CV_MIXT_DENS,ifl) = rhs1 
    pRhs(CV_MIXT_XMOM,ifl) = rhs2*dxdn + rhs3*dxds + rhs4*dxdt 
    pRhs(CV_MIXT_YMOM,ifl) = rhs2*dydn + rhs3*dyds + rhs4*dydt 
    pRhs(CV_MIXT_ZMOM,ifl) = rhs2*dzdn + rhs3*dzds + rhs4*dzdt

    IF ( pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE) == BCOPT_SUPERSONIC ) THEN
      pRhs(CV_MIXT_ENER,ifl) = rhs5 
    ELSE IF ( pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE) == BCOPT_SUBSONIC ) THEN
      IF ( pPatch%reflect == BC_REFLECTING ) THEN
        pRhs(CV_MIXT_ENER,ifl) = 0.0_RFREAL
      ELSE ! non reflecting BC
        pRhs(CV_MIXT_ENER,ifl) = rhs5 
      END IF ! pPatch%reflect
    END IF ! pPatch%mixt%switches
  END DO ! ifl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompRhsOF
 

  
  

  
  
  



! ******************************************************************************
!
! Purpose: Compute RHS for slip wall.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_CompRhsSW(pRegion,pPatch) 

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: c1,errorFlag,ifl,indCp,indMol
  REAL(RFREAL) :: a,aoa,aos,cp,d1,d2,d3,d4,d5,dotProduct,dndx,dndy,dndz,dpdn, &
                  dpds,dpdt,dpdx,dpdy,dpdz,drdn,drds,drdt,drdx,drdy,drdz,dsdx, &
                  dsdy,dsdz,dtdx,dtdy,dtdz,dudn,duds,dudt,dudx,dudy,dudz, &
                  dundn,dunds,dundt,dusdn,dusds,dusdt,dutdn,dutds,dutdt,dvdn, &
                  dvds,dvdt,dvdx,dvdy,dvdz,dwdn,dwds,dwdt,dwdx,dwdy,dwdz,dxdn, &
                  dxds,dxdt,dydn,dyds,dydt,dzdn,dzds,dzdt,g,L1,L2,L3,L4,L5, &
                  lambda1,lambda2,lambda3,lambda4,lambda5,mm,mf,mfNormal,nx, &
                  ny,nz,p,Pinf,r,rgas,rhs1,rhs2,rhs3,rhs4,rhs5,ru,rv,rw,sx, &
                  sy,sz,T,tx,ty,tz,u,un,us,ut,v,w
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv,pDv,pRhs,pGv
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGradFace
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_CompRhsSW',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main Body - loop over all boundary faces
! ******************************************************************************

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  pCv  => pPatch%mixt%cv
  pDv  => pPatch%mixt%dv
  pGv  => pRegion%mixt%gv ! NOTE get gv from volume
  pRhs => pPatch%mixt%rhs

  pGradFace => pPatch%mixt%gradFace

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    cp   = pGv(GV_MIXT_CP ,indCp *c1)
    mm   = pGv(GV_MIXT_MOL,indMol*c1)
    rgas = MixtPerf_R_M(mm)
    g    = MixtPerf_G_CpR(cp,rgas)

    r  = pCv(CV_MIXT_DENS,ifl)
    ru = pCv(CV_MIXT_XMOM,ifl)
    rv = pCv(CV_MIXT_YMOM,ifl)
    rw = pCv(CV_MIXT_ZMOM,ifl)
    p  = pDv(DV_MIXT_PRES,ifl)
    a  = pDv(DV_MIXT_SOUN,ifl)

    u = ru/r
    v = rv/r
    w = rw/r

    drdx = pGradFace(XCOORD,GRBF_MIXT_DENS,ifl)
    drdy = pGradFace(YCOORD,GRBF_MIXT_DENS,ifl)
    drdz = pGradFace(ZCOORD,GRBF_MIXT_DENS,ifl)

    dudx = pGradFace(XCOORD,GRBF_MIXT_XVEL,ifl)
    dudy = pGradFace(YCOORD,GRBF_MIXT_XVEL,ifl)
    dudz = pGradFace(ZCOORD,GRBF_MIXT_XVEL,ifl)

    dvdx = pGradFace(XCOORD,GRBF_MIXT_YVEL,ifl)
    dvdy = pGradFace(YCOORD,GRBF_MIXT_YVEL,ifl)
    dvdz = pGradFace(ZCOORD,GRBF_MIXT_YVEL,ifl)

    dwdx = pGradFace(XCOORD,GRBF_MIXT_ZVEL,ifl)
    dwdy = pGradFace(YCOORD,GRBF_MIXT_ZVEL,ifl)
    dwdz = pGradFace(ZCOORD,GRBF_MIXT_ZVEL,ifl)

    dpdx = pGradFace(XCOORD,GRBF_MIXT_PRES,ifl)
    dpdy = pGradFace(YCOORD,GRBF_MIXT_PRES,ifl)
    dpdz = pGradFace(ZCOORD,GRBF_MIXT_PRES,ifl)

! ==============================================================================
!   Computation of normal and tangential vectors
! ==============================================================================

    nx = pPatch%fn(XCOORD,ifl)
    ny = pPatch%fn(YCOORD,ifl)
    nz = pPatch%fn(ZCOORD,ifl)

    SELECT CASE ( pRegion%mixtInput%dimens )
      CASE ( 2 )
        sx = -ny
        sy =  nx
        sz =  nz

        tx = 0.0_RFREAL
        ty = 0.0_RFREAL
        tz = 1.0_RFREAL
      CASE ( 3 )
        sx = -ny
        sy =  nx
        sz =  nz

        tx = 0.0_RFREAL
        ty = 0.0_RFREAL
        tz = 1.0_RFREAL
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,3032)
    END SELECT ! pRegion%mixtInput%dimens

    dxdn = nx
    dydn = ny
    dzdn = nz

    dxds = sx
    dyds = sy
    dzds = sz

    dxdt = tx
    dydt = ty
    dzdt = tz

    dndx = dxdn
    dndy = dydn
    dndz = dzdn

    dsdx = dxds
    dsdy = dyds
    dsdz = dzds

    dtdx = dxdt
    dtdy = dydt
    dtdz = dzdt

! ==============================================================================
!   Computation of normal and tangential derivatives
! ==============================================================================

    drdn = drdx*dxdn + drdy*dydn + drdz*dzdn
    drds = drdx*dxds + drdy*dyds + drdz*dzds
    drdt = drdx*dxdt + drdy*dydt + drdz*dzdt

    dudn = dudx*dxdn + dudy*dydn + dudz*dzdn
    duds = dudx*dxds + dudy*dyds + dudz*dzds
    dudt = dudx*dxdt + dudy*dydt + dudz*dzdt

    dvdn = dvdx*dxdn + dvdy*dydn + dvdz*dzdn
    dvds = dvdx*dxds + dvdy*dyds + dvdz*dzds
    dvdt = dvdx*dxdt + dvdy*dydt + dvdz*dzdt

    dwdn = dwdx*dxdn + dwdy*dydn + dwdz*dzdn
    dwds = dwdx*dxds + dwdy*dyds + dwdz*dzds
    dwdt = dwdx*dxdt + dwdy*dydt + dwdz*dzdt

    dpdn = dpdx*dxdn + dpdy*dydn + dpdz*dzdn
    dpds = dpdx*dxds + dpdy*dyds + dpdz*dzds
    dpdt = dpdx*dxdt + dpdy*dydt + dpdz*dzdt

! ==============================================================================
!   Computation of normal and tangential vector components
! ==============================================================================

    un = u*dndx + v*dndy + w*dndz
    us = u*dsdx + v*dsdy + w*dsdz
    ut = u*dtdx + v*dtdy + w*dtdz

! ==============================================================================
!   Normal component of velocity should be zero.
! ==============================================================================

    un = 0.0_RFREAL

! ==============================================================================
!   Computation of normal and tangential derivatives
! ==============================================================================

    dundn = dudn*dndx + dvdn*dndy + dwdn*dndz
    dunds = duds*dndx + dvds*dndy + dwds*dndz
    dundt = dudt*dndx + dvdt*dndy + dwdt*dndz

    dusdn = dudn*dsdx + dvdn*dsdy + dwdn*dsdz
    dusds = duds*dsdx + dvds*dsdy + dwds*dsdz
    dusdt = dudt*dsdx + dvdt*dsdy + dwdt*dsdz

    dutdn = dudn*dtdx + dvdn*dtdy + dwdn*dtdz
    dutds = duds*dtdx + dvds*dtdy + dwds*dtdz
    dutdt = dudt*dtdx + dvdt*dtdy + dwdt*dtdz

! ==============================================================================
!   Computation of eigenvalues
! ==============================================================================

    lambda1 = un - a
    lambda2 = un
    lambda3 = un
    lambda4 = un
    lambda5 = un + a

! ==============================================================================
!   Computations of L
! ==============================================================================

    L2 = lambda2*(a*a*drdn - dpdn)
    L3 = lambda3*dusdn
    L4 = lambda4*dutdn
    L5 = lambda5*(dpdn + r*a*dundn)
    L1 = L5

! ==============================================================================
!   Computations of d
! ==============================================================================

    d1 = (L2 + 0.5_RFREAL*(L1+L5))/(a*a)
    d2 = 0.5_RFREAL*(L1+L5)
    d3 = (L5-L1)/(2.0_RFREAL*r*a)
    d4 = L3
    d5 = L4

! ==============================================================================
!   Computations of Rhs
! ==============================================================================

    rhs1 = d1 + r*dusds + us*drds + r*dutdt + ut*drdt
! Normal velocity is zero, no need to solve euler equaiton ---------------------
    rhs2 = 0.0_RFREAL
    rhs3 = us*d1 + r*d4 + dpds &
         + 2.0_RFREAL*r*us*dusds + us*us*drds &
         + r*us*dutdt + r*ut*dusdt + us*ut*drdt
    rhs4 = ut*d1 + r*d5 + dpdt &
         + r*ut*dusds + r*us*dutds + ut*us*drds &
         + 2.0_RFREAL*r*ut*dutdt + ut*ut*drdt
    rhs5 = 0.5_RFREAL*(un*un+us*us+ut*ut)*d1 + d2/(g-1.0_RFREAL) &
         + r*un*d3 + r*us*d4 + r*ut*d5 &
         + us*(  0.5_RFREAL*(un*un+us*us+ut*ut)*drds &
              + r*(un*dunds+us*dusds+ut*dutds) &
              + g*dpds/(g-1.0_RFREAL) ) &
         + (0.5_RFREAL*r*(un*un+us*us+ut*ut) &
              + g*p/(g-1.0_RFREAL))*dusds &
         + ut*(  0.5_RFREAL*(un*un+us*us+ut*ut)*drdt &
              + r*(un*dundt+us*dusdt+ut*dutdt) &
              + g*dpdt/(g-1.0_RFREAL) ) &
         + (0.5_RFREAL*r*(un*un+us*us+ut*ut) &
              + g*p/(g-1.0_RFREAL))*dutdt

    pRhs(CV_MIXT_DENS,ifl) = rhs1
    pRhs(CV_MIXT_XMOM,ifl) = rhs2*dxdn + rhs3*dxds + rhs4*dxdt
    pRhs(CV_MIXT_YMOM,ifl) = rhs2*dydn + rhs3*dyds + rhs4*dydt
    pRhs(CV_MIXT_ZMOM,ifl) = rhs2*dzdn + rhs3*dzds + rhs4*dzdt
    pRhs(CV_MIXT_ENER,ifl) = rhs5
  END DO ! ifl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_CompRhsSW 
  
  
  
  
  
  
 


! ******************************************************************************
!
! Purpose: Determine whether need boundary face gradients and related functions.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

LOGICAL FUNCTION RFLU_NSCBC_DecideHaveNSCBC(pRegion)

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  INTEGER :: iPatch
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_DecideHaveNSCBC',"../modflu/RFLU_ModNSCBC.F90")

! *****************************************************************************
! Initialize
! *****************************************************************************

  RFLU_NSCBC_DecideHaveNSCBC = .FALSE.

! *****************************************************************************
! Determine whether need boundary face gradients
! *****************************************************************************

  patchLoop: DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      RFLU_NSCBC_DecideHaveNSCBC = .TRUE.
  
      EXIT patchLoop
    END IF ! pPatch%bcKind
  END DO patchLoop

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_NSCBC_DecideHaveNSCBC







! ******************************************************************************
!
! Purpose: Initialize farfield boundary.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_InitFF(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: corrFlag
  INTEGER :: errorFlag,ifl,c1,bcOptType,distrib,gasModel,indCp,indGs,indMol
  REAL(RFREAL) :: aoa,aos,corr,cp,irl,Eo,g,gc,liftCoef,mf,mw,nm,nx,ny,nz,pf, &
                  pl,pr,rel,rer,rl,rr,rul,rur,rvl,rvr,rwl,rwr,tf,ul,vl,wl,xc, &
                  yc,zc
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv,vals
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_InitFF',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main Body 
! ******************************************************************************

  corr = pPatch%mixt%switches(BCSWI_FARF_CORR)

  vals => pPatch%mixt%vals
      
  IF ( corr == BCOPT_CORR_YES ) THEN 
    corrFlag = .TRUE. 
  ELSE 
    corrFlag = .FALSE.
    liftCoef = 0.0_RFREAL
  END IF ! corr

  indCp    = pRegion%mixtInput%indCp
  indMol   = pRegion%mixtInput%indMol
  gasModel = pRegion%mixtInput%gasModel

  distrib  =  pPatch%mixt%distrib

  pPatch%mixt%cvState = CV_MIXT_STATE_CONS

  IF (pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,3346, &
                   'region mixt cvState is incompatible...')
  END IF ! pRegion%mixt%cvState

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    nx = pPatch%fn(XCOORD,ifl)
    ny = pPatch%fn(YCOORD,ifl)
    nz = pPatch%fn(ZCOORD,ifl)
    nm = pPatch%fn(XYZMAG,ifl)

    xc = pPatch%fc(XCOORD,ifl)
    yc = pPatch%fc(YCOORD,ifl)
    zc = pPatch%fc(ZCOORD,ifl)

    rl  = pRegion%mixt%cv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl
    rul = pRegion%mixt%cv(CV_MIXT_XMOM,c1)*irl
    rvl = pRegion%mixt%cv(CV_MIXT_YMOM,c1)*irl
    rwl = pRegion%mixt%cv(CV_MIXT_ZMOM,c1)*irl
    rel = pRegion%mixt%cv(CV_MIXT_ENER,c1)*irl

    pl  = pRegion%mixt%dv(DV_MIXT_PRES,c1)

    mf  = vals(BCDAT_FARF_MACH  ,distrib*ifl)
    aoa = vals(BCDAT_FARF_ATTACK,distrib*ifl)
    aos = vals(BCDAT_FARF_SLIP  ,distrib*ifl)
    pf  = vals(BCDAT_FARF_PRESS ,distrib*ifl)
    tf  = vals(BCDAT_FARF_TEMP  ,distrib*ifl)

    IF ( gasModel == GAS_MODEL_TCPERF ) THEN

      mw = pRegion%mixt%gv(GV_MIXT_MOL,indMol*c1)
      cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *c1)
      gc = MixtPerf_R_M(mw)
      g  = MixtPerf_G_CpR(cp,gc)

      rel = rl*MixtPerf_Eo_DGPUVW(rl,g,pl,ul,vl,wl)
      
      CALL RFLU_SetRindStateFarfieldPerf(global,cp,mw,nx,ny,nz,mf,pf,tf, &
                                         aoa,aos,corrFlag,liftCoef,xc,yc, &
                                         zc,rl,rul,rvl,rwl,rel,rr,rur, & 
                                         rvr,rwr,rer,pr)     
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,3391)
    END IF ! gasModel

    pPatch%mixt%cv(CV_MIXT_DENS,ifl) = rr 
    pPatch%mixt%cv(CV_MIXT_XMOM,ifl) = rur
    pPatch%mixt%cv(CV_MIXT_YMOM,ifl) = rvr
    pPatch%mixt%cv(CV_MIXT_ZMOM,ifl) = rwr
    pPatch%mixt%cv(CV_MIXT_ENER,ifl) = rer
    
! TEMPORARY : overwriting above computation of farfield BC data by cell data
!             this is done to avoid waves at farfield just in begining of run
    pPatch%mixt%cv(CV_MIXT_DENS,ifl) = pRegion%mixt%cv(CV_MIXT_DENS,c1)
    pPatch%mixt%cv(CV_MIXT_XMOM,ifl) = pRegion%mixt%cv(CV_MIXT_XMOM,c1)
    pPatch%mixt%cv(CV_MIXT_YMOM,ifl) = pRegion%mixt%cv(CV_MIXT_YMOM,c1)
    pPatch%mixt%cv(CV_MIXT_ZMOM,ifl) = pRegion%mixt%cv(CV_MIXT_ZMOM,c1)
    pPatch%mixt%cv(CV_MIXT_ENER,ifl) = pRegion%mixt%cv(CV_MIXT_ENER,c1) 
! END TEMPORARY
  END DO ! ifl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_InitFF







! ******************************************************************************
!
! Purpose: Initialize inflow boundary with imposed total quantities.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_InitIFTotAng(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,ifl,c1
  INTEGER :: bcOptFixed,bcOptType,distrib,gasModel,indCp,indGs,indMol
  REAL(RFREAL) :: cp,betah,betav,g,gc,mach,mw,nx,ny,nz,nm,pr,ptot,rer,rl,rr, &
                  rul,rvl,rwl,rur,rvr,rwr,ttot
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv,vals
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_InitIFTotAng',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main Body 
! ******************************************************************************

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  distrib    = pPatch%mixt%distrib
  bcOptType  = pPatch%mixt%switches(BCSWI_INFLOW_TYPE)
  bcOptFixed = pPatch%mixt%switches(BCSWI_INFLOW_FIXED)

  vals => pPatch%mixt%vals

  pPatch%mixt%cvState = CV_MIXT_STATE_CONS

  IF (pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,3490, &
                   'region mixt cvState is incompatible...')
  END IF ! pRegion%mixt%cvState

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    mw = pRegion%mixt%gv(GV_MIXT_MOL,indMol*c1)
    cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *c1)

    gc = MixtPerf_R_M(mw)
    g  = MixtPerf_G_CpR(cp,gc)

    nx = pPatch%fn(XCOORD,ifl)
    ny = pPatch%fn(YCOORD,ifl)
    nz = pPatch%fn(ZCOORD,ifl)
    nm = pPatch%fn(XYZMAG,ifl)

    rl  = pRegion%mixt%cv(CV_MIXT_DENS,c1)
    rul = pRegion%mixt%cv(CV_MIXT_XMOM,c1)
    rvl = pRegion%mixt%cv(CV_MIXT_YMOM,c1)
    rwl = pRegion%mixt%cv(CV_MIXT_ZMOM,c1)

    ptot  = vals(BCDAT_INFLOW_PTOT, distrib*ifl)
    ttot  = vals(BCDAT_INFLOW_TTOT, distrib*ifl)
    betah = vals(BCDAT_INFLOW_BETAH,distrib*ifl)
    betav = vals(BCDAT_INFLOW_BETAV,distrib*ifl)

    IF ( bcOptType == BCOPT_SUPERSONIC ) THEN
      mach = vals(BCDAT_INFLOW_MACH,distrib*ifl)
    ELSE
      mach = 0.0_RFREAL
    END IF ! bcOptType

    CALL BcondInflowPerf(bcOptType,bcOptFixed,ptot,ttot,betah,betav, &
                         mach,nx,ny,nz,cp,mw,rl,rul,rvl,rwl,rr,rur, &
                         rvr,rwr,rer,pr)

    pPatch%mixt%cv(CV_MIXT_DENS,ifl) = rr 
    pPatch%mixt%cv(CV_MIXT_XMOM,ifl) = rur
    pPatch%mixt%cv(CV_MIXT_YMOM,ifl) = rvr
    pPatch%mixt%cv(CV_MIXT_ZMOM,ifl) = rwr
    pPatch%mixt%cv(CV_MIXT_ENER,ifl) = rer
  END DO ! ifl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_InitIFTotAng





! ******************************************************************************
!
! Purpose: Initialize inflow with imposed velocity and temperature.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_InitIFVelTemp(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,ifl,c1,bcOptType,distrib,gasModel,indCp,indGs,indMol
  REAL(RFREAL) :: cp,Eo,g,gc,mw,pr,rl,rr,rur,rvr,rwr,tr,ur,vr,wr
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv,vals
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_InitIFVelTemp',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main body 
! ******************************************************************************

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  distrib   = pPatch%mixt%distrib
  bcOptType = pPatch%mixt%switches(BCSWI_INFLOW_TYPE)

  vals => pPatch%mixt%vals

  pPatch%mixt%cvState = CV_MIXT_STATE_CONS

  IF (pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,3610, &
                   'region mixt cvState is incompatible...')
  END IF ! pRegion%mixt%cvState

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    mw = pRegion%mixt%gv(GV_MIXT_MOL,indMol*c1)
    cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *c1)

    gc = MixtPerf_R_M(mw)
    g  = MixtPerf_G_CpR(cp,gc)

    rl  = pRegion%mixt%cv(CV_MIXT_DENS,c1)

    ur = vals(BCDAT_INFLOW_U,distrib*ifl)
    vr = vals(BCDAT_INFLOW_V,distrib*ifl)
    wr = vals(BCDAT_INFLOW_W,distrib*ifl)
    tr = vals(BCDAT_INFLOW_T,distrib*ifl)

    IF ( bcOptType == BCOPT_SUPERSONIC ) THEN
      pr = vals(BCDAT_INFLOW_P,distrib*ifl)
      rr = MixtPerf_D_PRT(pr,gc,tr)
    ELSE
      rr = rl
      pr = MixtPerf_P_DRT(rr,gc,tr)
    END IF ! bcOptType

    rur = rr*ur
    rvr = rr*vr
    rwr = rr*wr

    pPatch%mixt%cv(CV_MIXT_DENS,ifl) = rr 
    pPatch%mixt%cv(CV_MIXT_XMOM,ifl) = rur
    pPatch%mixt%cv(CV_MIXT_YMOM,ifl) = rvr
    pPatch%mixt%cv(CV_MIXT_ZMOM,ifl) = rwr

    Eo = MixtPerf_Eo_DGPUVW(rr,g,Pr,ur,vr,wr)

    pPatch%mixt%cv(CV_MIXT_ENER,ifl) = rr*Eo
  END DO ! ifl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_InitIFVelTemp







! ******************************************************************************
!
! Purpose: Initialize injection boundary.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_InitIJ(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_InitIJ',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main body 
! ******************************************************************************

  pPatch%mixt%cv = 0.0_RFREAL

  pPatch%mixt%cvState = CV_MIXT_STATE_CONS

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_InitIJ







! ******************************************************************************
!
! Purpose: Initialize no-slip wall with imposed heat flux.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_InitNSWHeat(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_InitNSWHeat',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main body 
! ******************************************************************************

  pPatch%mixt%cv = 0.0_RFREAL

  pPatch%mixt%cvState = CV_MIXT_STATE_CONS

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_InitNSWHeat






! ******************************************************************************
!
! Purpose: Initialize no-slip wall with imposed temperature.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_InitNSWTemp(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_InitNSWTemp',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main body 
! ******************************************************************************

  pPatch%mixt%cv = 0.0_RFREAL

  pPatch%mixt%cvState = CV_MIXT_STATE_CONS

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_InitNSWTemp






! ******************************************************************************
!
! Purpose: Initialize outflow.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_InitOF(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,ifl,c1,bcOptType,distrib,gasModel,indCp,indGs,indMol
  REAL(RFREAL) :: cp,Eo,g,gc,mw,pl,pr,rl,rul,rvl,rwl,rel,ul,vl,wl
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv,vals
  TYPE(t_global), POINTER :: global


  INTEGER :: iCvSpecProducts
  REAL(RFREAL) :: YProducts
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
  TYPE(t_spec_input), POINTER :: pSpecInput

!Fred - 12/29/20 - Adding in species module to NSCBC init

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_NSCBC_InitOF',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main Body 
! ******************************************************************************
  
  IF ( global%specUsed .EQV. .TRUE. ) THEN
  pCvSpec => pRegion%spec%cv
  pSpecInput => pRegion%specInput
  END IF
  
  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  distrib =  pPatch%mixt%distrib
  bcOptType = pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE)

  IF ( bcOptType /= BCOPT_SUPERSONIC ) THEN
    vals => pPatch%mixt%vals
  END IF ! bcOptType

  pPatch%mixt%cvState = CV_MIXT_STATE_CONS

  IF (pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,3947, &
                   'region mixt cvState is incompatible...')
  END IF ! pRegion%mixt%cvState

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    mw = pRegion%mixt%gv(GV_MIXT_MOL,indMol*c1)
    cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *c1)

    gc = MixtPerf_R_M(mw)
    g  = MixtPerf_G_CpR(cp,gc)

    rl  = pRegion%mixt%cv(CV_MIXT_DENS,c1)
    rul = pRegion%mixt%cv(CV_MIXT_XMOM,c1)
    rvl = pRegion%mixt%cv(CV_MIXT_YMOM,c1)
    rwl = pRegion%mixt%cv(CV_MIXT_ZMOM,c1)
    rel = pRegion%mixt%cv(CV_MIXT_ENER,c1)

    
    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
    IF ( global%specUsed .EQV. .TRUE. ) THEN
       iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
       YProducts       = pCvSpec(iCvSpecProducts,c1) !Get local species fraction

    ELSE
       WRITE(*,*) 'Error - Must invoke Species module when using JWL NSCBC!!!'
       CALL ErrorStop(global,ERR_REACHED_DEFAULT,3975) ! Defensive coding    

    END IF !SpecUsed
    END IF !JWL

    ul = rul/rl
    vl = rvl/rl
    wl = rwl/rl

    pl  = pRegion%mixt%dv(DV_MIXT_PRES,c1)

    IF ( bcOptType /= BCOPT_SUPERSONIC ) THEN
      pr = vals(BCDAT_OUTFLOW_PRESS,distrib*ifl)
    ELSE
      pr = pl
    END IF ! bcOptType

    pPatch%mixt%cv(CV_MIXT_DENS,ifl) = rl 
    pPatch%mixt%cv(CV_MIXT_XMOM,ifl) = rul
    pPatch%mixt%cv(CV_MIXT_YMOM,ifl) = rvl
    pPatch%mixt%cv(CV_MIXT_ZMOM,ifl) = rwl

   IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
    Eo = 0.0_RFREAL !Fred - forcing this to zero if using JWL w/o species...will
                    !cause code to crash at first time step 
    IF ( global%specUsed .EQV. .TRUE. ) THEN
        Eo = RFLU_JWL_E_PR(pRegion,pl,rl,YProducts)
        Eo = Eo + 0.5_RFREAL*(ul*ul+vl*vl+wl*wl)  
    ELSE
        WRITE(*,*) 'Error - Must invoke Species module when using JWL NSCBC!!!'
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,4007) ! Defensive coding
    END IF !SpecUsed 
   ELSE 
    Eo = MixtPerf_Eo_DGPUVW(rl,g,Pr,ul,vl,wl)
   END IF !MixtJWL

    pPatch%mixt%cv(CV_MIXT_ENER,ifl) = rl*Eo
  END DO ! ifl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_InitOF







! ******************************************************************************
!
! Purpose: Initialize slip wall.
!
! Description: None.
!
! Input:
!   pRegion	Pointer to region
!   pPatch 	Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_NSCBC_InitSW(pRegion,pPatch)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,ifl,c1,bcOptType,distrib,gasModel,indCp,indGs,indMol
  REAL(RFREAL) :: cp,dndx,dndy,dndz,dsdx,dsdy,dsdz,dtdx,dtdy,dtdz,dxdn,dxds, &
                  dxdt,dydn,dyds,dydt,dzdn,dzds,dzdt,Eo,g,gc,mw,nx,ny,nz,rl, &
                  rr,rul,rvl,rwl,rel,sx,sy,sz,tx,ty,tz,pl,pr,ul,un,ur,us,ut, &
                  vl,vr,wl,wr,irl,al,tl,YProducts
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv,vals
  TYPE(t_global), POINTER :: global

  LOGICAL :: scalarConvFlag
  INTEGER :: iCvSpecProducts
  TYPE(t_spec_input), POINTER :: pSpecInput
  !Fred - Added JWL EOS Capabilities - 9/11/15

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  IF ( global%specUsed .EQV. .TRUE. ) THEN
  pSpecInput => pRegion%specInput
  END IF

  CALL RegisterFunction(global,'RFLU_NSCBC_InitSW',"../modflu/RFLU_ModNSCBC.F90")

! ******************************************************************************
! Main Body 
! ******************************************************************************

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

  distrib =  pPatch%mixt%distrib

  pPatch%mixt%cvState = CV_MIXT_STATE_CONS

  IF (pRegion%mixt%cvState /= CV_MIXT_STATE_CONS ) THEN
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,4107, &
                   'region mixt cvState is incompatible...')
  END IF ! pRegion%mixt%cvState

  DO ifl = 1,pPatch%nBFaces
    c1 = pPatch%bf2c(ifl)

    mw = pRegion%mixt%gv(GV_MIXT_MOL,indMol*c1)
    cp = pRegion%mixt%gv(GV_MIXT_CP ,indCp *c1)

    gc = MixtPerf_R_M(mw)
    g  = MixtPerf_G_CpR(cp,gc)

    rl  = pRegion%mixt%cv(CV_MIXT_DENS,c1)
    rul = pRegion%mixt%cv(CV_MIXT_XMOM,c1)
    rvl = pRegion%mixt%cv(CV_MIXT_YMOM,c1)
    rwl = pRegion%mixt%cv(CV_MIXT_ZMOM,c1)
    rel = pRegion%mixt%cv(CV_MIXT_ENER,c1)

    ul = rul/rl
    vl = rvl/rl
    wl = rwl/rl

    pl  = pRegion%mixt%dv(DV_MIXT_PRES,c1)

! ------------------------------------------------------------------------------
!   Computation of normal and tangential vectors
! ------------------------------------------------------------------------------

    nx = pPatch%fn(XCOORD,ifl)
    ny = pPatch%fn(YCOORD,ifl)
    nz = pPatch%fn(ZCOORD,ifl)

    SELECT CASE ( pRegion%mixtInput%dimens )
      CASE ( 2 )
        sx = -ny
        sy =  nx
        sz =  nz

        tx = 0.0_RFREAL
        ty = 0.0_RFREAL
        tz = 1.0_RFREAL
      CASE ( 3 )
        sx = -ny
        sy =  nx
        sz =  nz

        tx = 0.0_RFREAL
        ty = 0.0_RFREAL
        tz = 1.0_RFREAL
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,4158)
    END SELECT ! pRegion%mixtInput%dimens

    dxdn = nx
    dydn = ny
    dzdn = nz

    dxds = sx
    dyds = sy
    dzds = sz

    dxdt = tx
    dydt = ty
    dzdt = tz

    dndx = dxdn
    dndy = dydn
    dndz = dzdn

    dsdx = dxds
    dsdy = dyds
    dsdz = dzds

    dtdx = dxdt
    dtdy = dydt
    dtdz = dzdt

! ------------------------------------------------------------------------------
!   Computation of normal and tangential velocity components
! ------------------------------------------------------------------------------

    un = ul*dndx + vl*dndy + wl*dndz
    us = ul*dsdx + vl*dsdy + wl*dsdz
    ut = ul*dtdx + vl*dtdy + wl*dtdz

! ------------------------------------------------------------------------------
!   Normal component of velocity should be zero.
! ------------------------------------------------------------------------------

    un = 0.0_RFREAL

! ------------------------------------------------------------------------------
!   Computation of velocity components from normal and tangential components
! ------------------------------------------------------------------------------

    rr = rl
    ur = un*dxdn + us*dxds + ut*dxdt
    vr = un*dydn + us*dyds + ut*dydt
    wr = un*dzdn + us*dzds + ut*dzdt
    pr = pl

! ------------------------------------------------------------------------------
!   Computation of conservative variables
! ------------------------------------------------------------------------------

    pPatch%mixt%cv(CV_MIXT_DENS,ifl) = rr
    pPatch%mixt%cv(CV_MIXT_XMOM,ifl) = rr*ur
    pPatch%mixt%cv(CV_MIXT_YMOM,ifl) = rr*vr
    pPatch%mixt%cv(CV_MIXT_ZMOM,ifl) = rr*wr

    IF (pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL) THEN
     IF (global%specUsed .EQV. .TRUE.) THEN
       IF (pRegion%spec%cvState == CV_MIXT_STATE_PRIM) THEN
          scalarConvFlag = .FALSE.
       ELSE
          scalarConvFlag = .TRUE.
 CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
       END IF

    irl = 1.0_RFREAL/rl
    iCvSpecProducts = SPEC_GetSpeciesIndex(global,pSpecInput,'PRODUCTS')
    Yproducts = pRegion%spec%cv(iCvSpecProducts,c1)

    CALL RFLU_JWL_ComputeEnergyMixt(pRegion,c1,g,gc,pl,rl,Yproducts,al,Eo,tl)

     IF (scalarConvFlag .EQV. .TRUE.) THEN
  CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv,pRegion%spec%cvState)
     END IF
    END IF !specUsed
    ELSE
    Eo = MixtPerf_Eo_DGPUVW(rr,g,Pr,ur,vr,wr)
    END IF

    pPatch%mixt%cv(CV_MIXT_ENER,ifl) = rr*Eo
  END DO ! ifl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_NSCBC_InitSW







END MODULE RFLU_ModNSCBC

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModNSCBC.F90,v $
! Revision 1.3  2016/02/08 20:03:05  fred
! Fixing non-1 flag compiling error
!
! Revision 1.2  2016/02/04 19:58:07  fred
! Adding iterative JWL EOS capabilities for the cylindrical detonation case
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.6  2009/07/08 19:11:58  mparmar
! Adapted for axisymmetric computations
!
! Revision 1.5  2008/12/06 08:43:41  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:16:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/12/05 13:23:11  haselbac
! Rm unnecessary RCSIdentString decls, ifort compiler on vonkarman complained
!
! Revision 1.2  2007/06/18 18:00:59  mparmar
! Initial attempt to incorporate moving reference frame with NSCBC
!
! Revision 1.1  2007/04/09 18:49:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:00:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.5  2007/02/27 13:05:23  haselbac
! Cosmetics only
!
! Revision 1.4  2006/10/20 21:18:45  mparmar
! Added mass/momentum coeffs in flux routines and cosmetic clean-up
!
! Revision 1.3  2006/08/21 16:11:19  haselbac
! Clean-up: comments, order of routines, missing headers, indentation, etc
!
! Revision 1.2  2006/08/19 19:44:08  haselbac
! Significant clean-up and cosmetic changes
!
! Revision 1.1  2006/08/19 15:37:46  mparmar
! Initial revision
!
! ******************************************************************************

