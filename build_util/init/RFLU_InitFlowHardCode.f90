










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
! Purpose: Initialize flow field in a region using hard-coded values.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine assumes a perfect gas.
!
! ******************************************************************************
!
!!  Current cases:
!     barrelExp
!     cylds, rectshktb
!     wedge
!!  Legacy  cases:
!     acoustic
!     channelacoust
!     cyldet
!     ctint
!     shktb
!     cylpotential
!     cylwall
!     sphdet, detwav3D
!     ds (series)
!     vortexNSCBC
!     gcgc
!     gaussianpulse
!     gmpjet
!     gtlin
!     gttri
!     jet
!     kjet (series)
!     MShock_H2O_Air001
!     MShock_H2O_Air002
!     MShock_Air_Air_He001
!     mvsh
!     farf
!     nscbc1
!     nscbc2
!     nscbc3
!     nscbc4
!     nscbc5
!     nscbc6
!     nscbc7
!     nscbc8
!     ncavity
!     onera_c0, onera_c0_2d_100x50
!     onera_c0_3d
!     pipeacoust
!     radialpulse
!     rayleighproblem
!     ringleb
!     ShockBubble
!     skews_ms2p0, skews_ms3p0, skews_ms4p0
!     somm_spi
!     sphds
!     sphpotential
!     stg1d, stg2d
!     st_sod1, st_sod1_mp2
!     st_sod2, st_sod2_mp2
!     ssvorth (series)
!     taylorvortex
!     Two_Rarefaction
!     volcmod2dv3
!     vort
!     vortex_part
!     wcramp, wcrampsc, wcrampsm
!     2DShock001
!
!
! ******************************************************************************
!
! $Id: RFLU_InitFlowHardCode.F90,v 1.7 2016/02/05 20:01:55 fred Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InitFlowHardCode(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input, &
                        t_mixt
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters
    
  USE RFLU_ModBessel  
  USE RFLU_ModExactFlow  
  USE RFLU_ModFlowHardCode        
    
  USE ModInterfaces, ONLY: MixtPerf_Cv_CpR, &
                           MixtPerf_D_CGP, &
                           MixtPerf_D_PRT, &
                           MixtGasLiq_Eo_CvmTVm2, &
                           MixtPerf_Eo_DGPUVW, &
                           Mixt_Eo_eUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_P_DRT, &
                           MixtPerf_R_CpG, &  
                           MixtPerf_R_M, &
                           MixtPerf_C_DGP

  USE RFLU_ModJWL ! Subbu
         
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: errorString,iFileName1,RCSIdentString
  INTEGER :: errorFlag,flag,iBc,icg,im,in,indCp,indMol,iq,m,n,n1,n2,n3
  REAL(RFREAL) :: A,A_c,A_cl,A1,A2,Amp,aTot,c,const,cp,cpRef,Cvm,cvg,cvl,cvv,d, &
                  dInc,dMin,dOffs,dTot,dummyReal,etaqm,ExplEner,g,gaussAmp,gaussCenter, & 
                  gaussWidth,gc,gcRef,gRef,gx,gy,gz,H,height,L,LD,L_l,Lx,Ly,Lz,Mi, &
                  mInj,Mo,muo,mw,nx,ny,nz,omega,p,perturb,pg,pi,pl,pMin,po,pOffs, &
                  psi,pTot,pv,r,r1,r2,Rc,Rc_l,radius,refD,refL,refNu,refP,refU,rg, &
                  rGas,ri,rl,ro,rv,rVap,t,temp,tTot,the,theta,u,uo,uo_c,um,ur,ut,uz, &
                  v,vInj,Vm2,w,x,xMin,xo,xx,y,yp,yMin,yo,z,zMin,zz,eAir,ux,vy,wz
  REAL(RFREAL),ALLOCATABLE,DIMENSION(:) :: randNum
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld,pGv,pDv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput

  ! Saptarshi - 01/10/2015 - begin  
  REAL(RFREAL) :: dx,dy,dz,coordSinInt,coordSinIntEnd
  ! Saptarshi - 01/10/2015 - end
  ! TLJ - 05/14/2024 - begin  
  REAL(RFREAL) :: xs,ys,zs,rs
  REAL(RFREAL) :: emin, emax
  REAL(RFREAL) :: ea, pa, ra
  ! TLJ - 05/14/2024 - end  

! Rahul - Read RBA Input file
  CHARACTER(CHRLEN) :: endString
  INTEGER :: b,i,iLow,iHigh,inb
  REAL(RFREAL) :: wLow,wHigh,wInv,rLow,rHigh,e,e2,nearx
  REAL(RFREAL) :: tempYx, xQ, Rmax, Pmax, Umax
  REAL(RFREAL) :: etotal,TLJ_Cp_air,TLJ_Cv_air,TLJ_T0_air,TLJ_eamb,YR
  REAL(RFREAL),ALLOCATABLE,DIMENSION(:) :: xData,radData,rData,uData,eData
  REAL(RFREAL),ALLOCATABLE,DIMENSION(:) :: YData,eData2,YRdata
  CHARACTER(LEN=20) :: tmpstr
  INTEGER :: LIMIT_LOOP
  INTEGER :: ixend
  INTEGER :: inearx

! Rahul - End

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitFlowHardCode.F90,v $ $Revision: 1.7 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InitFlowHardCode', &
                        "../../utilities/init/RFLU_InitFlowHardCode.F90")
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Initializing flow field from hard code...'
                             
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME,'Case: ',TRIM(global%casename)
    END IF ! global%verbLevel                         
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid      => pRegion%grid
  pCv        => pRegion%mixt%cv
  pCvOld     => pRegion%mixt%cvOld
  pGv        => pRegion%mixt%gv    
  pDv        => pRegion%mixt%dv
  pMixtInput => pRegion%mixtInput

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol  

! ******************************************************************************
! Set constant gas properties. NOTE used only in some cases
! ******************************************************************************

  cpRef = global%refCp
  gRef  = global%refGamma  
  gcRef = MixtPerf_R_CpG(cpRef,gRef)

! ******************************************************************************
! Initialize flow field based on user input and fluid model
! ******************************************************************************

  SELECT CASE ( pMixtInput%fluidModel )
  
! ==============================================================================
!   Incompressible fluid model
! ==============================================================================  
  
    CASE ( FLUID_MODEL_INCOMP ) 
      pRegion%mixt%cvState = CV_MIXT_STATE_PRIM
      
! TEMPORARY, to be replaced by proper code
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,277) 
! END TEMPORARY

! ==============================================================================
!   Compressible fluid model
! ==============================================================================  
    
    CASE ( FLUID_MODEL_COMP )
      IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
        pRegion%mixt%cvState = CV_MIXT_STATE_DUVWP
      ELSE
        pRegion%mixt%cvState = CV_MIXT_STATE_CONS
      END IF ! global%solverType

        SELECT CASE ( global%casename )

!-------------------------------------------------------------------------------
!       Couture: Barrel with Explosive: Barrel Length = 40 mm, Diameter = 12mm 
!       Explosive region is located from:  z = 0mm:33.5mm
!       Particle bed is located from:      z = 34mm:40mm
!       TLJ - coded for PBXN-5 with JWL parameters for LX-10-0
!-------------------------------------------------------------------------------
        CASE ("barrelExp")

           print*,'TLJ RFLU_InitFlowHardCode : ',TRIM(global%casename)
           print*,'TLJ specUsed       : ',global%specUsed
           print*,'TLJ nSpecies       : ',global%nSpecies 
           print*,'TLJ Gas Model      : ',pRegion%mixtInput%gasModel
           print*,'TLJ RB value       : ',pMixtInput%prepRealVal26
           print*,'TLJ RB File        : ',TRIM(pMixtInput%rBurnFile)

           !***************
           IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_JWL ) THEN ! Normal Case

              DO icg = 1,pGrid%nCellsTot
                 x = pGrid%cofg(XCOORD,icg)
                 y = pGrid%cofg(YCOORD,icg)
                 z = pGrid%cofg(ZCOORD,icg)
                 radius = SQRT(x**2 + y**2)
                 IF ( radius < pMixtInput%prepRealVal1 .AND. &
                       (z < pMixtInput%prepRealVal2) ) THEN
                    ! High Pressure Region
                    d = pMixtInput%prepRealVal5
                    u = 0.0_RFREAL
                    v = 0.0_RFREAL
                    !w = 0.0_RFREAL
                    w = pMixtInput%prepRealVal7
                    p = pMixtInput%prepRealVal6
                 ELSE
                    ! Ambient Region
                    d = pMixtInput%prepRealVal3
                    u = 0.0_RFREAL
                    v = 0.0_RFREAL
                    w = 0.0_RFREAL
                    p = pMixtInput%prepRealVal4
                 END IF
                 mw = pGv(GV_MIXT_MOL,indMol*icg)
                 cp = pGv(GV_MIXT_CP ,indCp *icg)
                 gc = MixtPerf_R_M(mw)
                 g  = MixtPerf_G_CpR(cp,gc)
   
                 pCv(CV_MIXT_DENS,icg) = d
                 pCv(CV_MIXT_XMOM,icg) = d*u
                 pCv(CV_MIXT_YMOM,icg) = d*v
                 pCv(CV_MIXT_ZMOM,icg) = d*w
                 pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
    
              END DO ! icg
           ENDIF ! non-JWL case
           !***************

           !+++++++++++++++
           IF ( pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL ) THEN ! JWL Case

           IF ( pMixtInput%prepRealVal26 /= 1.0000 ) THEN ! JWL with Constant values
              DO icg = 1,pGrid%nCellsTot
                 x = pGrid%cofg(XCOORD,icg)  !Extract x coordinate of cellcentroid
                 y = pGrid%cofg(YCOORD,icg)  !Extract y coordinate of cellcentroid
                 z = pGrid%cofg(ZCOORD,icg)  !Extract z coordinate of cellcentroid

                 radius = SQRT(x**2 + y**2)

                 IF ( radius < pMixtInput%prepRealVal6 .AND. &
                    (z < pMixtInput%prepRealVal1) ) THEN
                    ! Explosive
                    d = pMixtInput%prepRealVal24
                    u = 0.0_RFREAL
                    v = 0.0_RFREAL
                    w = pMixtInput%prepRealVal5
                    e = pMixtInput%prepRealVal13/d
                 ELSE
                    ! Unshocked Air
                    d = pMixtInput%prepRealVal3
                    u = 0.0_RFREAL
                    v = 0.0_RFREAL
                    w = 0.0_RFREAL
                    p = pMixtInput%prepRealVal2
                 END IF

                 mw = pGv(GV_MIXT_MOL,indMol*icg)
                 cp = pGv(GV_MIXT_CP ,indCp *icg)
                 gc = MixtPerf_R_M(mw)
                 g  = MixtPerf_G_CpR(cp,gc)

                 pCv(CV_MIXT_DENS,icg) = d
                 pCv(CV_MIXT_XMOM,icg) = d*u
                 pCv(CV_MIXT_YMOM,icg) = d*v
                 pCv(CV_MIXT_ZMOM,icg) = d*w

                 ! Compute or read e - jgarno
                 IF (d == pMixtInput%prepRealVal3) THEN
                    ! Ideal gas law for Air
                    pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
                 ELSE
                    ! JWL for Products
                    ! Note e2 should be read in from onedx_hmx_mgeos.dat
                    ! READ(1990,*) (xData(icg),rData(icg),uData(icg), &
                    !   eData(icg),Ydata(icg),eData2(icg),YRdata(icg), icg=1,n-1)
                    ! wLow = (x-xData(inearx))/(xData(inearx+1)-xData(inearx))
                    ! d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                    ! e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                    ! e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                    ! u = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                    ! v = 0.0_RFREAL
                    ! w = 0.0_RFREAL
                    e2 = 189.86636969074604_RFREAL
                    p = RFLU_JWL_P_ER(pRegion,e2,d,1.00_RFREAL)
                    if ( e2 .lt. 2.00E+05_RFREAL ) then
                       e2 = 2.00E+05_RFREAL
                    endif
                    pCv(CV_MIXT_ENER,icg) = d*(e2 + 0.5_RFREAL*w*w)
                    !pCv(CV_MIXT_ENER,icg) = d*e
                 ENDIF
              END DO ! icg

           ENDIF ! Constant values case

  
           IF ( pMixtInput%prepRealVal26 == 1.0000 ) THEN ! JWL with Reactive Burn

           ! Read in RocSDT grid and dependent variables
           ! Read in onedx_hmx_mgeos.dat
           ! -0.02 < xData < 0.0399975  (meters)
           ! xData is the RocSDT grid
           ! For the barrelExp problem, x corresponds to z

           iFileName1 = pMixtInput%rBurnFile
           OPEN(unit=1990, FILE=iFileName1, FORM='FORMATTED', &
               IOSTAT=errorFlag)
           global%error = errorFlag
           IF ( global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_FILE_OPEN,428,iFileName1)
           END IF ! global%error

           n = 0
           DO
              n = n+1
              READ(1990,'(A)') endString
              IF ( endString == TRIM('# END') ) THEN
                 CLOSE(1990)
                 EXIT
              END IF ! name
              IF ( n >= LIMIT_INFINITE_LOOP ) THEN
                 ! Guard against infinite loop
                 CALL ErrorStop(global,ERR_INFINITE_LOOP,441)
              END IF ! n
           END DO ! Infinite DO
           n = n-1

           ! Allocate arrays to store column data from file
           ! xc  = xData  (data along the streamwise length of barrel)
           ! rho = rData
           ! u   = uData
           ! p   = eData
           ! phi = Ydata  (location of air/explosive material interface)
           ! ein = eData2 (specific internal energy)
           ! Y_R = YRdata (mass fraction of reactant)
           ALLOCATE (xData(n),rData(n),uData(n),eData(n),eData2(n),Ydata(n), &
                     YRdata(n),STAT=errorFlag)
           global%error = errorFlag
           IF ( global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_ALLOCATE,458,'radData')
           END IF ! global%error
   
           OPEN(unit=1990,FILE=iFileName1,FORM='FORMATTED', &
                    IOSTAT=errorFlag)
           global%error = errorFlag
           IF ( global%error /= ERR_NONE ) THEN
                  CALL ErrorStop(global,ERR_FILE_OPEN,465,iFileName1)
           END IF ! global%error
           READ(1990,*) (xData(icg),rData(icg),uData(icg), &
              eData(icg),Ydata(icg),eData2(icg),YRdata(icg), icg=1,n)
           CLOSE(1990)

           !open (unit=1990,file='onedx_hmx_mgeos.dat',FORM='FORMATTED')
           !READ(1990,*) (xData(icg),rData(icg),uData(icg), &
           !   eData(icg),Ydata(icg),eData2(icg),YRdata(icg), icg=1,n)
           !close(1990)


           ixend = maxloc(xData(:),1)
           Rmax  = maxval(rData(:))
           Pmax  = maxval(eData(:))
           Umax  = maxval(uData(:))
           emin  = minval(eData2(:))
           emax  = maxval(eData2(:))

           ! TLJ - this should be updated using input values
           TLJ_Cp_air = 1004.64_RFREAL      ! (J/kg-K)
           TLJ_Cv_air = 717.60_RFREAL       ! (J/kg-K)
           TLJ_T0_air = 288.15_RFREAL       ! (K)
           TLJ_eamb = TLJ_Cv_air*TLJ_T0_air ! (J/kg)
           ! TLJ - now updated using input values
           TLJ_eamb = pMixtInput%prepRealVal2/pMixtInput%prepRealVal3/0.4_RFREAL
           IF ( global%myProcid == MASTERPROC) THEN
              print*,'*****   Reactive Burn Profile Info   *****'
              print*,'PID ',global%myProcid
              print*,'Burn rate file: ',TRIM(iFileName1)
              print*,'Number of RocSDT steps: ',n
              print*,'ixend = ',ixend
              print*,'R (loc,max) = ',xData(maxloc(rData,1)),Rmax
              print*,'P (loc,max) = ',xData(maxloc(eData,1)),Pmax
              print*,'U (loc,max) = ',xData(maxloc(uData,1)),Umax
              print*,'e (loc,min) = ',xData(minloc(eData2,1)),emin
              print*,'e (loc,max) = ',xData(maxloc(eData2,1)),emax
              print*,'eamb_TLJ = ',TLJ_eamb
           ENDIF

           ! Josh Garno - For Reactive Burn Interpolation scheme
           ! TLJ - Do not change below unless you know what you are doing

           ! Step 1:
           !    Locate the x-point where the mass fraction of the
           !    reactant YR is one. This defines the end of the
           !    reaction zone.
           xQ = minloc(abs(YRdata(:)-1.00E-03_RFREAL),1)
           xQ = xData(nint(xQ))
           IF ( global%myProcid == MASTERPROC) THEN
              print*,'xQ = ',xQ
           ENDIF

           ! Step 2:
           !    Shift RocSDT domain so that the mass fraction of the
           !    reactant YR is at the end of the explosive, given
           !    by tempYx = pMixtInput%prepRealVal1.
           tempYx = pMixtInput%prepRealVal1
           IF ( global%myProcid == MASTERPROC) THEN
              print*,'xshift = ',xQ,tempYx,tempYx-xQ
           ENDIF
           xData(:) = xData(:) + (tempYx-xQ)

           ! Step 3:
           !    Now find the location of the material interface,
           !    and call this xQ (over-written from above defn.).
           !    This defines the beginning of the reaction zone.
           !    Thus, the explosive region will between xQ < x < prepRealVal(1)
           xQ = minloc(abs(Ydata(:)-0.5_RFREAL),1)
           xQ = xData(nint(xQ))
   
           IF ( global%myProcid == MASTERPROC) THEN
              WRITE(6,*) 'maxval(xData(:)) before shift', &
                      maxval(xData(:)),pRegion%iRegionGlobal
           ENDIF

           ! Step 4:
           !    Interpolate data onto Rocflu mesh and compute
           !    the initial conditions

           IF ( global%myProcid == MASTERPROC) THEN
              print*
              print*,xQ,' <= x <= ',pMixtInput%prepRealVal1
              print*
           ENDIF


           DO icg = 1,pGrid%nCellsTot
              x = pGrid%cofg(XCOORD,icg)  !Extract x coordinate of cellcentroid
              y = pGrid%cofg(YCOORD,icg)  !Extract y coordinate of cellcentroid
              z = pGrid%cofg(ZCOORD,icg)  !Extract z coordinate of cellcentroid
   
              radius = SQRT(x**2 + y**2)
   
              IF ( (radius < pMixtInput%prepRealVal6) .AND. &
                 ((xQ <= z) .AND. (z < pMixtInput%prepRealVal1) )) THEN
                 ! Explosive Region

                 ! Translate x->z, u->w
                 inearx = minloc(abs(xData(:)-z),1)
                 nearx  = z-xData(inearx)

                 IF ( (inearx .gt. 1) .AND. (inearx .lt. ixend) ) THEN
                    IF ( nearx .ge. 0.00_RFREAL ) THEN
                       wLow = (z-xData(inearx))/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = 0.0_RFREAL
                       v = 0.0_RFREAL
                       w = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                    ELSE
                       inearx = inearx-1
                       wLow = (z-xData(inearx))/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = 0.0_RFREAL
                       v = 0.0_RFREAL
                       w = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                    ENDIF
                 ELSEIF ( inearx .eq. 1 ) THEN
                    IF ( nearx .ge. 0.00_RFREAL ) THEN
                       wLow = (z-xData(inearx))/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = 0.0_RFREAL
                       v = 0.0_RFREAL
                       w = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                    ELSE
                       wLow = (xData(inearx)-z)/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (pMixtInput%prepRealVal3-rData(inearx))*wLow
                       e = eData(inearx) + (pMixtInput%prepRealVal2-eData(inearx))*wLow
                       e2 = eData2(inearx) + (2.00E+05_RFREAL-eData2(inearx))*wLow
                       u = 0.0_RFREAL
                       v = 0.0_RFREAL
                       w = uData(inearx) + (0.00_RFREAL-uData(inearx))*wLow
                       Y = 0.0_RFREAL
                       YR = 0.0_RFREAL
                    ENDIF
                 ELSEIF ( inearx .eq. ixend ) THEN
                    IF ( nearx .ge. 0.00_RFREAL ) THEN
                       wLow = (z-xData(inearx))/(xData(inearx)-xData(inearx-1))
                       d = rData(inearx) + (pMixtInput%prepRealVal3-rData(inearx))*wLow
                       e = eData(inearx) + (pMixtInput%prepRealVal2-eData(inearx))*wLow
                       e2 = eData2(inearx) + (2.00E+05_RFREAL-eData2(inearx))*wLow
                       u = 0.0_RFREAL
                       v = 0.0_RFREAL
                       w = uData(inearx) + (0.00_RFREAL-uData(inearx))*wLow
                       Y = 1.0_RFREAL
                       YR = 1.0_RFREAL
                    ELSE
                       inearx = inearx-1
                       wLow = (z-xData(inearx))/(xData(inearx)-xData(inearx-1))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = 0.0_RFREAL
                       v = 0.0_RFREAL
                       w = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                    ENDIF
                 ELSE
                    WRITE(*,*) 'Missed a spot'
                 ENDIF ! inearx
   
                 IF ( d > Rmax .OR. d < 0.00d0) THEN
                    WRITE(6,*) 'd > Rmax: d, icg, inearx, wLow, x, &
                      & xData(inearx-1:inearx+1), &
                      & max(xData(:)),r Data(inearx-1:inearx+1)',&
                      & d,icg,inearx,wLow,z,xData(inearx-1),xData(inearx),xData(inearx+1), &
                      & maxval(xData(:)),rData(inearx-1),rData(inearx),rData(inearx+1)
                    CALL ErrorStop(global,ERR_INVALID_VALUE,645,'Invalid den InitFlow')
                 ENDIF

                 IF (Y  .LE. 1.0E-20_RFREAL) Y  = 0.0_RFREAL
                 IF (YR .LE. 1.0E-20_RFREAL) YR = 0.0_RFREAL

                 ! We reset the minimum internal energy to the ambient
                 ! values because RocSDT assumes the strong shock approximation
                 if (e2 .lt. TLJ_eamb) e2 = TLJ_eamb
                 etotal = d*(e2+0.5_RFREAL*(u*u+v*v+w*w))
                 p = e

                 ! TLJ - treatment for open or closed barrel
                 !   Comment out below for Brad's geometry
                 !   Remove negative velocities in reactive burn for barrel case
                 !if (u <= 0.0_RFREAL) u = 0.0_RFREAL
                 !if (v <= 0.0_RFREAL) v = 0.0_RFREAL
                 !if (w <= 0.0_RFREAL) w = 0.0_RFREAL

                 pCv(CV_MIXT_DENS,icg) = d
                 pCv(CV_MIXT_XMOM,icg) = d*u
                 pCv(CV_MIXT_YMOM,icg) = d*v
                 pCv(CV_MIXT_ZMOM,icg) = d*w
                 pCv(CV_MIXT_ENER,icg) = d*(e2 + 0.5_RFREAL*(u*u+v*v+w*w))
   
                 ! We need these for plotting solutions
                 !pDv(DV_MIXT_PRES,icg) = p
                 !pDv(DV_MIXT_TEMP,icg) = RFLU_JWL_T_PR(pRegion,p,d,e2,1.0_RFREAL)
                 !pDv(DV_MIXT_SOUN,icg) = RFLU_JWL_C_ER(pRegion,e2,d,1.0_RFREAL)
                 T = RFLU_JWL_T_PR(pRegion,p,d,e2,1.0_RFREAL)
                 a = RFLU_JWL_C_ER(pRegion,e2,d,1.0_RFREAL)
   
              ELSE
                 ! Ambient Region
                 d = pMixtInput%prepRealVal3
                 u = 0.0_RFREAL
                 v = 0.0_RFREAL
                 w = 0.0_RFREAL
                 p = pMixtInput%prepRealVal2
   
                 mw = pGv(GV_MIXT_MOL,indMol*icg)
                 cp = pGv(GV_MIXT_CP ,indCp *icg)
                 gc = MixtPerf_R_M(mw)
                 g  = MixtPerf_G_CpR(cp,gc)
   
                 e2 = p/(d*0.4_RFREAL)  ! e = p/(\rho*(gamma-1))
   
                 pCv(CV_MIXT_DENS,icg) = d
                 pCv(CV_MIXT_XMOM,icg) = d*u
                 pCv(CV_MIXT_YMOM,icg) = d*v
                 pCv(CV_MIXT_ZMOM,icg) = d*w
                 pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
   
                 ! We need these for plotting solutions
                 !pDv(DV_MIXT_PRES,icg) = p
                 !pDv(DV_MIXT_TEMP,icg) = 300.0_RFREAL
                 !pDv(DV_MIXT_SOUN,icg) = MixtPerf_C_DGP(d,g,p)
                 T = p/(gc*d)           ! T = p/(\rho*Rgas)
                 a = MixtPerf_C_DGP(d,g,p)
                 Y = 0.0_RFREAL
                 YR = 0.0_RFREAL

              END IF

              IF (IsNaN(pCv(CV_MIXT_ENER,icg)) .EQV. .TRUE.) THEN
                 print*,'TLJ Input leads to NaN'
                 print*,'    e,d',pCv(CV_MIXT_ENER,icg),d,icg
                 CALL ErrorStop(global,ERR_INVALID_VALUE,712,'Invalid val InitFlow')
              ENDIF
              IF (IsNaN(p) .EQV. .TRUE.) THEN
                 print*,'TLJ Input leads to NaN'
                 print*,'    p,d',p,d,icg
                 CALL ErrorStop(global,ERR_INVALID_VALUE,717,'Invalid val InitFlow')
              ENDIF

              ! Josh Garno - WRITE INIT Check File
              IF ( radius < pMixtInput%prepRealVal6 .AND. &
                       (z < 0.04_RFREAL) ) THEN
                       !(z < pMixtInput%prepRealVal1) ) THEN
                 OPEN(1119,file='INIT_Centerline.dat',form='formatted',status='old' &
                          ,position='append')
                 WRITE(1119,'(11(E23.16,1X),I8,1X,E23.16)') &
                            x,y,z,d,w,e2,p,T,a,Y,YR,icg,maxval(xData(:))
                 CLOSE(1119)
              ENDIF ! y
   
           END DO ! icg

           DEALLOCATE(xData,rData,eData,uData,Ydata,eData2,YRdata)
   
           ENDIF ! RB case
           ENDIF ! JWL case
           !+++++++++++++++


! ------------------------------------------------------------------------------
!       Cylinder or rectangular shock tube
!       Thierry Daoud - Quaterplane case
!       Modified for quarterplane with JWL EOS 1-equation model
! ------------------------------------------------------------------------------

        CASE ( "cylds", "rectshktb" )

           print*,'TLJ RFLU_InitFlowHardCode : ',TRIM(global%casename)
           print*,'TLJ specUsed       : ',global%specUsed
           print*,'TLJ nSpecies       : ',global%nSpecies 
           print*,'TLJ Gas Model      : ',pRegion%mixtInput%gasModel
           print*,'TLJ RB value       : ',pMixtInput%prepRealVal26
           print*,'TLJ RB File        : ',TRIM(pMixtInput%rBurnFile)

           !***************
           IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_JWL ) THEN ! Normal Case

              DO icg = 1,pGrid%nCellsTot
                 x = pGrid%cofg(XCOORD,icg)
                 y = pGrid%cofg(YCOORD,icg)
                 z = pGrid%cofg(ZCOORD,icg)

                 IF (TRIM(global%casename)=="cylds") THEN
                    radius = SQRT(x**2 + y**2)
                    theta = ATAN2(y,x)
                    ur = pMixtInput%prepRealVal3 * cos(theta)
                    ut = pMixtInput%prepRealVal3 * sin(theta)
                 ELSEIF (TRIM(global%casename)=="rectshktb") THEN
                    radius = x
                    ur = pMixtInput%prepRealVal3
                    ut = 0.0_RFREAL
                 ENDIF

                 IF ( radius < pMixtInput%prepRealVal1 ) THEN
                    d = pMixtInput%prepRealVal2
                    u = ur
                    v = ut
                    w = 0.0_RFREAL
                    p = pMixtInput%prepRealVal4
                 ELSE
                    d = pMixtInput%prepRealVal5
                    u = 0.0_RFREAL
                    v = 0.0_RFREAL
                    w = 0.0_RFREAL
                    p = pMixtInput%prepRealVal6
                 END IF ! x
                 mw = pGv(GV_MIXT_MOL,indMol*icg)
                 cp = pGv(GV_MIXT_CP ,indCp *icg)
                 gc = MixtPerf_R_M(mw)
                 g  = MixtPerf_G_CpR(cp,gc)

                 pCv(CV_MIXT_DENS,icg) = d
                 pCv(CV_MIXT_XMOM,icg) = d*u
                 pCv(CV_MIXT_YMOM,icg) = d*v
                 pCv(CV_MIXT_ZMOM,icg) = d*w
                 pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)

              END DO ! icg

           ENDIF
           !***************

           !+++++++++++++++
           IF ( pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL ) THEN ! JWL Case


           IF ( pMixtInput%prepRealVal26 /= 1.0000 ) THEN ! JWL with Constant values
              DO icg = 1,pGrid%nCellsTot
                 x = pGrid%cofg(XCOORD,icg)  !Extract x coordinate of cellcentroid
                 y = pGrid%cofg(YCOORD,icg)  !Extract y coordinate of cellcentroid
                 z = pGrid%cofg(ZCOORD,icg)  !Extract z coordinate of cellcentroid

                 IF (TRIM(global%casename)=="cylds") THEN
                    radius = SQRT(x**2 + y**2)
                    theta = ATAN2(y,x)
                    ur = pMixtInput%prepRealVal5 * cos(theta)
                    ut = pMixtInput%prepRealVal5 * sin(theta)
                 ELSEIF (TRIM(global%casename)=="rectshktb") THEN
                    radius = x
                    ur = pMixtInput%prepRealVal5
                    ut = 0.0_RFREAL
                 ENDIF

                 IF ( radius < pMixtInput%prepRealVal1 ) THEN
                    ! Explosive
                    d = pMixtInput%prepRealVal24
                    u = ur
                    v = ut
                    w = 0.0_RFREAL
                    e = pMixtInput%prepRealVal13/d
                 ELSE
                    ! Unshocked Air
                    d = pMixtInput%prepRealVal3
                    u = 0.0_RFREAL
                    v = 0.0_RFREAL
                    w = 0.0_RFREAL
                    p = pMixtInput%prepRealVal2
                 END IF
                 mw = pGv(GV_MIXT_MOL,indMol*icg)
                 cp = pGv(GV_MIXT_CP ,indCp *icg)
                 gc = MixtPerf_R_M(mw)
                 g  = MixtPerf_G_CpR(cp,gc)

                 pCv(CV_MIXT_DENS,icg) = d
                 pCv(CV_MIXT_XMOM,icg) = d*u
                 pCv(CV_MIXT_YMOM,icg) = d*v
                 pCv(CV_MIXT_ZMOM,icg) = d*w

                 ! Compute or read e - jgarno
                 IF (d == pMixtInput%prepRealVal3) THEN
                    ! Ideal gas law for Air
                    pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
                 ELSE
                    ! JWL for Products
                    ! Note e2 should be read in from onedx_hmx_mgeos.dat
                    ! READ(1990,*) (xData(icg),rData(icg),uData(icg), &
                    !   eData(icg),Ydata(icg),eData2(icg),YRdata(icg), icg=1,n-1)
                    ! wLow = (x-xData(inearx))/(xData(inearx+1)-xData(inearx))
                    ! d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                    ! e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                    ! e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                    ! u = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                    ! v = 0.0_RFREAL
                    ! w = 0.0_RFREAL
                    e2 = e
                    p = RFLU_JWL_P_ER(pRegion,e2,d,1.00_RFREAL)

                    ra = pMixtInput%prepRealVal3
                    pa = pMixtInput%prepRealVal2
                    ea = pa/ra/0.4_RFREAL
                    if ( e2 .lt. ea ) then
                       e2 = ea
                       p  = pa
                    endif
                    pCv(CV_MIXT_ENER,icg) = d*(e2 + 0.5_RFREAL*(u*u+v*v+w*w))
                 ENDIF
              END DO ! icg

           ENDIF ! Constant values case


           IF ( pMixtInput%prepRealVal26 == 1.0000 ) THEN ! JWL with Reactive Burn

           ro = pMixtInput%prepRealVal1 ! Diaphragm loc

           ! Read in RocSDT grid and dependent variables
           ! Read in onedx_hmx_mgeos.dat
           ! -0.02 < xData < 0.0399975  (meters)
           ! xData is the RocSDT grid
           ! For any cylindrical problem, x corresponds to radius
           ! RocSDT currently assumes planar, so there are no
           !   curvature assumptions in the detonation front.
           !   Thus, RocSDT sould be modified to allow for curvature
           !   in the future.

           iFileName1 = pMixtInput%rBurnFile
           OPEN(unit=1990, FILE=iFileName1, FORM='FORMATTED', &
               IOSTAT=errorFlag)
           global%error = errorFlag
           IF ( global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_FILE_OPEN,901,iFileName1)
           END IF ! global%error
   
           n = 0
           DO
              n = n+1
              READ(1990,'(A)') endString
              IF ( endString == TRIM('# END') ) THEN
                 CLOSE(1990)
                 EXIT
              END IF ! name
              IF ( n >= LIMIT_INFINITE_LOOP ) THEN
                 ! Guard against infinite loop
                 CALL ErrorStop(global,ERR_INFINITE_LOOP,914)
              END IF ! n
           END DO ! Infinite DO
           n = n-1

           ! Allocate arrays to store column data from file
           ! xc  = xData  (radial coordinate)
           ! rho = rData  (density)
           ! u   = uData  (radial velocity)
           ! p   = eData  (pressure)
           ! phi = Ydata  (location of air/explosive material interface)
           ! ein = eData2 (specific internal energy)
           ! Y_R = YRdata (mass fraction of reactant)
           ALLOCATE (xData(n),rData(n),uData(n),eData(n),eData2(n),Ydata(n), &
                     YRdata(n),STAT=errorFlag)
           global%error = errorFlag
           OPEN(unit=1990,FILE=iFileName1,FORM='FORMATTED', &
                    IOSTAT=errorFlag)
           global%error = errorFlag
           IF ( global%error /= ERR_NONE ) THEN
                  CALL ErrorStop(global,ERR_FILE_OPEN,934,iFileName1)
           END IF ! global%error
           READ(1990,*) (xData(icg),rData(icg),uData(icg), &
              eData(icg),Ydata(icg),eData2(icg),YRdata(icg), icg=1,n)
           CLOSE(1990)
           
           ixend = maxloc(xData(:),1)
           Rmax  = maxval(rData(:))
           Pmax  = maxval(eData(:))
           Umax  = maxval(uData(:))
           emin  = minval(eData2(:))
           emax  = maxval(eData2(:))

           ! TLJ - this should be updated using input values
           TLJ_Cp_air = 1004.64_RFREAL      ! (J/kg-K)
           TLJ_Cv_air = 717.60_RFREAL       ! (J/kg-K)
           TLJ_T0_air = 288.15_RFREAL        ! (K)
           TLJ_eamb = TLJ_Cv_air*TLJ_T0_air ! (J/kg)
           ! TLJ - now updated using input values
           TLJ_eamb = pMixtInput%prepRealVal2/pMixtInput%prepRealVal3/0.4_RFREAL
           IF ( global%myProcid == MASTERPROC) THEN
              print*,'*****   Reactive Burn Profile Info   *****'
              print*,'Burn rate file: ',TRIM(iFileName1)
              print*,'Number of RocSDT steps: ',n
              print*,'ixend = ',ixend
              print*,'R (loc,max) = ',xData(maxloc(rData,1)),Rmax
              print*,'P (loc,max) = ',xData(maxloc(eData,1)),Pmax
              print*,'U (loc,max) = ',xData(maxloc(uData,1)),Umax
              print*,'e (loc,min) = ',xData(minloc(eData2,1)),emin
              print*,'e (loc,max) = ',xData(maxloc(eData2,1)),emax
              print*,'eamb_TLJ = ',TLJ_eamb
           ENDIF

           ! Josh Garno - For Reactive Burn Interpolation scheme
           ! TLJ - Do not change below unless you know what you are doing

           ! Step 1:
           !    Locate the x-point where the mass fraction of the
           !    reactant YR is one. This defines the end of the
           !    reaction zone.
           xQ = minloc(abs(YRdata(:)-1.00E-03_RFREAL),1)
           xQ = xData(nint(xQ))
           IF ( global%myProcid == MASTERPROC) THEN
              print*,'xQ = ',xQ
           ENDIF
   
           ! Step 2:
           !    Shift RocSDT domain so that the mass fraction of the
           !    reactant YR is at the end of the explosive, given
           !    by tempYx = pMixtInput%prepRealVal1.
           tempYx = pMixtInput%prepRealVal1
           IF ( global%myProcid == MASTERPROC) THEN
              print*,'xshift = ',xQ,tempYx,tempYx-xQ
           ENDIF
           xData(:) = xData(:) + (tempYx-xQ)
   
           ! Step 3:
           !    Now find the location of the material interface,
           !    and call this xQ (over-written from above defn.).
           !    This defines the beginning of the reaction zone.
           !    Thus, the explosive region will between xQ < x < prepRealVal(1)
           xQ = minloc(abs(Ydata(:)-0.5_RFREAL),1)
           xQ = xData(nint(xQ))
   
           IF ( global%myProcid == MASTERPROC) THEN
              WRITE(6,*) 'maxval(xData(:)) before shift', &
                      maxval(xData(:)),pRegion%iRegionGlobal
           ENDIF
   
           ! Step 4:
           !    Interpolate data onto Rocflu mesh and compute
           !    the initial conditions
   
           IF ( global%myProcid == MASTERPROC) THEN
              print*
              print*,xQ,' <= x <= ',pMixtInput%prepRealVal1
              print*
           ENDIF
   
           DO icg = 1,pGrid%nCellsTot
              x = pGrid%cofg(XCOORD,icg)  !Extract x coordinate of cellcentroid
              y = pGrid%cofg(YCOORD,icg)  !Extract y coordinate of cellcentroid
              z = pGrid%cofg(ZCOORD,icg)  !Extract z coordinate of cellcentroid
   
              radius = SQRT(x**2 + y**2)
              theta = DATAN2(y,x)
   
              IF ( (radius < ro) .AND. (xQ <= radius)) THEN
                 ! Explosive Region
                 inearx = minloc(abs(xData(:)-radius),1)
                 nearx  = radius-xData(inearx)
   
                 IF ( (inearx .gt. 1) .AND. (inearx .lt. ixend) ) THEN
                    IF ( nearx .ge. 0.00_RFREAL ) THEN
                       wLow = (radius-xData(inearx))/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                    ELSE
                       inearx = inearx-1
                       wLow = (radius-xData(inearx))/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                   ENDIF
                 ELSEIF ( inearx .eq. 1 ) THEN
                    IF ( nearx .ge. 0.00_RFREAL ) THEN
                       wLow = (radius-xData(inearx))/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                    ELSE
                       wLow = (xData(inearx)-radius)/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (pMixtInput%prepRealVal3-rData(inearx))*wLow
                       e = eData(inearx) + (pMixtInput%prepRealVal2-eData(inearx))*wLow
                       e2 = eData2(inearx) + (2.00E+05_RFREAL-eData2(inearx))*wLow
                       u = uData(inearx) + (0.00_RFREAL-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = 0.0_RFREAL
                       YR = 0.0_RFREAL
                    ENDIF
                 ELSEIF ( inearx .eq. ixend ) THEN
                    IF ( nearx .ge. 0.00_RFREAL ) THEN
                       wLow = (radius-xData(inearx))/(xData(inearx)-xData(inearx-1))
                       d = rData(inearx) + (pMixtInput%prepRealVal3-rData(inearx))*wLow
                       e = eData(inearx) + (pMixtInput%prepRealVal2-eData(inearx))*wLow
                       e2 = eData2(inearx) + (2.00E+05_RFREAL-eData2(inearx))*wLow
                       u = uData(inearx) + (0.00_RFREAL-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = 1.0_RFREAL
                       YR = 1.0_RFREAL
                    ELSE
                       inearx = inearx-1
                       wLow = (radius-xData(inearx))/(xData(inearx)-xData(inearx-1))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                    ENDIF
                 ELSE
                    WRITE(*,*) 'Missed a spot'
                 ENDIF ! inearx

                 IF (Y  .LE. 1.0E-20_RFREAL) Y  = 0.0_RFREAL
                 IF (YR .LE. 1.0E-20_RFREAL) YR = 0.0_RFREAL
   
                 ux = u*cos(theta)
                 vy = u*sin(theta)
                 wz = w

                 ! We reset the minimum internal energy to the ambient
                 ! values because RocSDT assumes the strong shock approximation
                 if (e2 .lt. TLJ_eamb) e2 = TLJ_eamb
                 etotal = d*(e2+0.5_RFREAL*(ux*ux + vy*vy + wz*wz))
                 p = e

                 pCv(CV_MIXT_DENS,icg) = d
                 pCv(CV_MIXT_XMOM,icg) = d*ux
                 pCv(CV_MIXT_YMOM,icg) = d*vy
                 pCv(CV_MIXT_ZMOM,icg) = d*wz
                 pCv(CV_MIXT_ENER,icg) = etotal
   
                 ! We need these for plotting solutions
                 !pDv(DV_MIXT_PRES,icg) = p
                 !pDv(DV_MIXT_TEMP,icg) = RFLU_JWL_T_PR(pRegion,p,d,e2,1.0_RFREAL)
                 !pDv(DV_MIXT_SOUN,icg) = RFLU_JWL_C_ER(pRegion,e2,d,1.0_RFREAL)
                 T = RFLU_JWL_T_PR(pRegion,p,d,e2,1.0_RFREAL)
                 a = RFLU_JWL_C_ER(pRegion,e2,d,1.0_RFREAL)
   
               ELSE
                 ! Unshocked Air
                 d = pMixtInput%prepRealVal3
                 u = 0.0_RFREAL
                 v = 0.0_RFREAL
                 w = 0.0_RFREAL
                 p = pMixtInput%prepRealVal2
   
                 mw = pGv(GV_MIXT_MOL,indMol*icg)
                 cp = pGv(GV_MIXT_CP ,indCp *icg)
                 gc = MixtPerf_R_M(mw)
                 g  = MixtPerf_G_CpR(cp,gc)
   
                 e2 = p/(d*0.4_RFREAL)
   
                 pCv(CV_MIXT_DENS,icg) = d
                 pCv(CV_MIXT_XMOM,icg) = d*u
                 pCv(CV_MIXT_YMOM,icg) = d*v
                 pCv(CV_MIXT_ZMOM,icg) = d*w
                 pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
   
                 ! We need these for plotting solutions
                 !pDv(DV_MIXT_PRES,icg) = p
                 !pDv(DV_MIXT_TEMP,icg) = 300.0_RFREAL
                 !pDv(DV_MIXT_SOUN,icg) = MixtPerf_C_DGP(d,g,p)
                 T = p/(gc*d)           ! T = p/(\rho*Rgas)
                 a = MixtPerf_C_DGP(d,g,p)
                 Y = 0.0_RFREAL
                 YR = 0.0_RFREAL
   
              ENDIF
   
              IF (IsNaN(pCv(CV_MIXT_ENER,icg)) .EQV. .TRUE.) THEN
                 print*,'TLJ Input leads to NaN'
                 print*,'    e,d',pCv(CV_MIXT_ENER,icg),d,icg
                 CALL ErrorStop(global,ERR_INVALID_VALUE,1159,'Invalid val InitFlow')
              ENDIF
              IF (IsNaN(p) .EQV. .TRUE.) THEN
                 print*,'TLJ Input leads to NaN'
                 print*,'    p,d',p,d,icg
                 CALL ErrorStop(global,ERR_INVALID_VALUE,1164,'Invalid val InitFlow')
              ENDIF
   
              ! Josh Garno - WRITE INIT Check File
              IF ( radius < ro ) THEN
                 OPEN(1119,file='INIT_Centerline.dat',form='formatted',status='old' &
                          ,position='append')
                 WRITE(1119,'(11(E23.16,1X),I8,1X,E23.16)') &
                            radius,theta,z,d,u,e2,p,T,a,Y,YR,icg,maxval(xData(:))
                 CLOSE(1119)
              ENDIF ! y

           END DO ! icg

           DEALLOCATE(xData,rData,eData,uData,Ydata,eData2,YRdata)

           ENDIF ! RB case
           ENDIF ! JWL case
           !+++++++++++++++


! ------------------------------------------------------------------------------
!       Wedge case - Josh Gillis
! ------------------------------------------------------------------------------

        CASE ( "wedge" )

          print*,'TLJ RFLU_InitFlowHardCode : ',TRIM(global%casename)
          print*,'TLJ specUsed       : ',global%specUsed
          print*,'TLJ nSpecies       : ',global%nSpecies 
          print*,'TLJ Gas Model      : ',pRegion%mixtInput%gasModel
          print*,'TLJ RB value       : ',pMixtInput%prepRealVal26
          print*,'TLJ RB File        : ',TRIM(pMixtInput%rBurnFile)

          !***************
          !***************
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_JWL ) THEN ! Normal Case

          ro = pMixtInput%prepRealVal1 ! Diaphragm loc

          DO icg = 1,pGrid%nCellsTot
             x = pGrid%cofg(XCOORD,icg)
             y = pGrid%cofg(YCOORD,icg)
             z = pGrid%cofg(ZCOORD,icg)

             radius = sqrt(x*x+y*y)

             IF ( radius < ro ) THEN
                d = pMixtInput%prepRealVal2
                u = pMixtInput%prepRealVal3
                v = 0.0_RFREAL
                w = 0.0_RFREAL
                p = pMixtInput%prepRealVal4
             ELSE
                d = pMixtInput%prepRealVal5
                u = 0.0_RFREAL
                v = 0.0_RFREAL
                w = 0.0_RFREAL
                p = pMixtInput%prepRealVal6
             END IF ! x
             mw = pGv(GV_MIXT_MOL,indMol*icg)
             cp = pGv(GV_MIXT_CP ,indCp *icg)
             gc = MixtPerf_R_M(mw)
             g  = MixtPerf_G_CpR(cp,gc)
 
             pCv(CV_MIXT_DENS,icg) = d
             pCv(CV_MIXT_XMOM,icg) = d*u
             pCv(CV_MIXT_YMOM,icg) = d*v
             pCv(CV_MIXT_ZMOM,icg) = d*w
             pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)

          END DO ! icg

          ENDIF ! End normal high-pressure/low-pressure AIR case
          !***************
  
        !+++++++++++++++
        IF ( pRegion%mixtInput%gasModel == GAS_MODEL_MIXT_JWL ) THEN ! JWL Case

        IF ( pMixtInput%prepRealVal26 /= 1.0000 ) THEN ! JWL with Constant values
           DO icg = 1,pGrid%nCellsTot
                 x = pGrid%cofg(XCOORD,icg)  !Extract x coordinate of cellcentroid
                 y = pGrid%cofg(YCOORD,icg)  !Extract y coordinate of cellcentroid
                 z = pGrid%cofg(ZCOORD,icg)  !Extract z coordinate of cellcentroid

                 radius = SQRT(x**2 + y**2)

                 IF ( radius < pMixtInput%prepRealVal6 .AND. &
                    (z < pMixtInput%prepRealVal1) ) THEN
                    ! Explosive
                    d = pMixtInput%prepRealVal24
                    u = 0.0_RFREAL
                    v = 0.0_RFREAL
                    w = pMixtInput%prepRealVal5
                    e = pMixtInput%prepRealVal13/d
                 ELSE
                    ! Unshocked Air
                    d = pMixtInput%prepRealVal3
                    u = 0.0_RFREAL
                    v = 0.0_RFREAL
                    w = 0.0_RFREAL
                    p = pMixtInput%prepRealVal2
                    END IF

                 mw = pGv(GV_MIXT_MOL,indMol*icg)
                 cp = pGv(GV_MIXT_CP ,indCp *icg)
                 gc = MixtPerf_R_M(mw)
                 g  = MixtPerf_G_CpR(cp,gc)

                 pCv(CV_MIXT_DENS,icg) = d
                 pCv(CV_MIXT_XMOM,icg) = d*u
                 pCv(CV_MIXT_YMOM,icg) = d*v
                 pCv(CV_MIXT_ZMOM,icg) = d*w

                 ! Compute or read e - jgarno
                 IF (d == pMixtInput%prepRealVal3) THEN
                    ! Ideal gas law for Air
                    pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
                 ELSE
                    ! JWL for Products
                    ! Note e2 should be read in from onedx_hmx_mgeos.dat
                    ! READ(1990,*) (xData(icg),rData(icg),uData(icg), &
                    !   eData(icg),Ydata(icg),eData2(icg),YRdata(icg), icg=1,n-1)
                    ! wLow = (x-xData(inearx))/(xData(inearx+1)-xData(inearx))
                    ! d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                    ! e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                    ! e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                    ! u = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                    ! v = 0.0_RFREAL
                    ! w = 0.0_RFREAL
                    !e2 = 189.86636969074604_RFREAL
                    e2 = e
                    p = RFLU_JWL_P_ER(pRegion,e2,d,1.00_RFREAL)
                    if ( e2 .lt. 2.00E+05_RFREAL ) then
                       e2 = 2.00E+05_RFREAL
                    endif
                    pCv(CV_MIXT_ENER,icg) = d*(e2 + 0.5_RFREAL*w*w)
                    !pCv(CV_MIXT_ENER,icg) = d*e
                 ENDIF
              END DO ! icg

           ENDIF ! Constant values case


           IF ( pMixtInput%prepRealVal26 == 1.0000 ) THEN ! JWL with Reactive Burn

           ro = pMixtInput%prepRealVal1 ! Diaphragm loc

           ! Read in RocSDT grid and dependent variables
           ! Read in onedx_hmx_mgeos.dat
           ! -0.02 < xData < 0.0399975  (meters)
           ! xData is the RocSDT grid
           ! For any cylindrical problem, x corresponds to radius
           ! RocSDT currently assumes planar, so there are no
           !   curvature assumptions in the detonation front.
           !   Thus, RocSDT sould be modified to allow for curvature
           !   in the future.

           iFileName1 = pMixtInput%rBurnFile
           OPEN(unit=1990, FILE=iFileName1, FORM='FORMATTED', &
               IOSTAT=errorFlag)
           global%error = errorFlag
           IF ( global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_FILE_OPEN,1327,iFileName1)
           END IF ! global%error
   
           n = 0
           DO
              n = n+1
              READ(1990,'(A)') endString
              IF ( endString == TRIM('# END') ) THEN
                 CLOSE(1990)
                 EXIT
              END IF ! name
              IF ( n >= LIMIT_INFINITE_LOOP ) THEN
                 ! Guard against infinite loop
                 CALL ErrorStop(global,ERR_INFINITE_LOOP,1340)
           END IF ! n
           END DO ! Infinite DO
           n = n-1

           ! Allocate arrays to store column data from file
           ! xc  = xData  (radial coordinate)
           ! rho = rData  (density)
           ! u   = uData  (radial velocity)
           ! p   = eData  (pressure)
           ! phi = Ydata  (location of air/explosive material interface)
           ! ein = eData2 (specific internal energy)
           ! Y_R = YRdata (mass fraction of reactant)
           ALLOCATE (xData(n),rData(n),uData(n),eData(n),eData2(n),Ydata(n), &
                     YRdata(n),STAT=errorFlag)
           global%error = errorFlag
           IF ( global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_ALLOCATE,1357,'radData')
           END IF ! global%error

           OPEN(unit=1990,FILE=iFileName1,FORM='FORMATTED', &
                    IOSTAT=errorFlag)
           global%error = errorFlag
           IF ( global%error /= ERR_NONE ) THEN
                  CALL ErrorStop(global,ERR_FILE_OPEN,1364,iFileName1)
           END IF ! global%error
           READ(1990,*) (xData(icg),rData(icg),uData(icg), &
              eData(icg),Ydata(icg),eData2(icg),YRdata(icg), icg=1,n)
           CLOSE(1990)

           ixend = maxloc(xData(:),1)
           Rmax  = maxval(rData(:))
           Pmax  = maxval(eData(:))
           Umax  = maxval(uData(:))
           emin  = minval(eData2(:))
           emax  = maxval(eData2(:))

           ! TLJ - this should be updated using input values
           TLJ_Cp_air = 1004.64_RFREAL      ! (J/kg-K)
           TLJ_Cv_air = 717.60_RFREAL       ! (J/kg-K)
           TLJ_T0_air = 288.15_RFREAL       ! (K)
           TLJ_eamb = TLJ_Cv_air*TLJ_T0_air ! (J/kg)
           ! TLJ - now updated using input values
           TLJ_eamb = pMixtInput%prepRealVal2/pMixtInput%prepRealVal3/0.4_RFREAL
           IF ( global%myProcid == MASTERPROC) THEN
              print*,'*****   Reactive Burn Profile Info   *****'
              print*,'Burn rate file: ',TRIM(iFileName1)
              print*,'Number of RocSDT steps: ',n
              print*,'ixend = ',ixend
              print*,'R (loc,max) = ',xData(maxloc(rData,1)),Rmax
              print*,'P (loc,max) = ',xData(maxloc(eData,1)),Pmax
              print*,'U (loc,max) = ',xData(maxloc(uData,1)),Umax
              print*,'e (loc,min) = ',xData(minloc(eData2,1)),emin
              print*,'e (loc,max) = ',xData(maxloc(eData2,1)),emax
              print*,'eamb_TLJ = ',TLJ_eamb
           ENDIF

           ! Josh Garno - For Reactive Burn Interpolation scheme
           ! TLJ - Do not change below unless you know what you are doing

           ! Step 1:
           !    Locate the x-point where the mass fraction of the
           !    reactant YR is one. This defines the end of the
           !    reaction zone.
           xQ = minloc(abs(YRdata(:)-1.00E-03_RFREAL),1)
           xQ = xData(nint(xQ))
           IF ( global%myProcid == MASTERPROC) THEN
              print*,'xQ = ',xQ
           ENDIF

           ! Step 2:
           !    Shift RocSDT domain so that the mass fraction of the
           !    reactant YR is at the end of the explosive, given
           !    by tempYx = pMixtInput%prepRealVal1.
           tempYx = pMixtInput%prepRealVal1
           IF ( global%myProcid == MASTERPROC) THEN
              print*,'xshift = ',xQ,tempYx,tempYx-xQ
           ENDIF
           xData(:) = xData(:) + (tempYx-xQ)

           ! Step 3:
           !    Now find the location of the material interface,
           !    and call this xQ (over-written from above defn.).
           !    This defines the beginning of the reaction zone.
           !    Thus, the explosive region will between xQ < x < prepRealVal(1)
           xQ = minloc(abs(Ydata(:)-0.5_RFREAL),1)
           xQ = xData(nint(xQ))

           IF ( global%myProcid == MASTERPROC) THEN
              WRITE(6,*) 'maxval(xData(:)) before shift', &
                      maxval(xData(:)),pRegion%iRegionGlobal
           ENDIF

           ! Step 4:
           !    Interpolate data onto Rocflu mesh and compute
           !    the initial conditions

           IF ( global%myProcid == MASTERPROC) THEN
              print*
              print*,xQ,' <= x <= ',pMixtInput%prepRealVal1
              print*
           ENDIF

           DO icg = 1,pGrid%nCellsTot
              x = pGrid%cofg(XCOORD,icg)  !Extract x coordinate of cellcentroid
              y = pGrid%cofg(YCOORD,icg)  !Extract y coordinate of cellcentroid
              z = pGrid%cofg(ZCOORD,icg)  !Extract z coordinate of cellcentroid

              radius = SQRT(x**2 + y**2)
              theta = DATAN2(y,x)

              IF ( (radius < ro) .AND. (xQ <= radius)) THEN
                 ! Explosive Region
                 inearx = minloc(abs(xData(:)-radius),1)
                 nearx  = radius-xData(inearx)

                 IF ( (inearx .gt. 1) .AND. (inearx .lt. ixend) ) THEN
                    IF ( nearx .ge. 0.00_RFREAL ) THEN
                       wLow = (radius-xData(inearx))/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                    ELSE
                       inearx = inearx-1
                       wLow = (radius-xData(inearx))/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                   ENDIF
                 ELSEIF ( inearx .eq. 1 ) THEN
                    IF ( nearx .ge. 0.00_RFREAL ) THEN
                       wLow = (radius-xData(inearx))/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                    ELSE
                       wLow = (xData(inearx)-radius)/(xData(inearx+1)-xData(inearx))
                       d = rData(inearx) + (pMixtInput%prepRealVal3-rData(inearx))*wLow
                       e = eData(inearx) + (pMixtInput%prepRealVal2-eData(inearx))*wLow
                       e2 = eData2(inearx) + (2.00E+05_RFREAL-eData2(inearx))*wLow
                       u = uData(inearx) + (0.00_RFREAL-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = 0.0_RFREAL
                       YR = 0.0_RFREAL
                    ENDIF
                 ELSEIF ( inearx .eq. ixend ) THEN
                    IF ( nearx .ge. 0.00_RFREAL ) THEN
                       wLow = (radius-xData(inearx))/(xData(inearx)-xData(inearx-1))
                       d = rData(inearx) + (pMixtInput%prepRealVal3-rData(inearx))*wLow
                       e = eData(inearx) + (pMixtInput%prepRealVal2-eData(inearx))*wLow
                       e2 = eData2(inearx) + (2.00E+05_RFREAL-eData2(inearx))*wLow
                       u = uData(inearx) + (0.00_RFREAL-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = 1.0_RFREAL
                       YR = 1.0_RFREAL
                    ELSE
                       inearx = inearx-1
                       wLow = (radius-xData(inearx))/(xData(inearx)-xData(inearx-1))
                       d = rData(inearx) + (rData(inearx+1)-rData(inearx))*wLow
                       e = eData(inearx) + (eData(inearx+1)-eData(inearx))*wLow
                       e2 = eData2(inearx) + (eData2(inearx+1)-eData2(inearx))*wLow
                       u = uData(inearx) + (uData(inearx+1)-uData(inearx))*wLow
                       v = 0.0_RFREAL
                       w = 0.0_RFREAL
                       Y = YData(inearx) + (YData(inearx+1)-YData(inearx))*wLow
                       YR = YRData(inearx) + (YRData(inearx+1)-YRData(inearx))*wLow
                    ENDIF
                 ELSE
                    WRITE(*,*) 'Missed a spot'
                 ENDIF ! inearx

                 IF (Y  .LE. 1.0E-20_RFREAL) Y  = 0.0_RFREAL
                 IF (YR .LE. 1.0E-20_RFREAL) YR = 0.0_RFREAL

                 ux = u*cos(theta)
                 vy = u*sin(theta)
                 wz = w

                 ! We reset the minimum internal energy to the ambient
                 ! values because RocSDT assumes the strong shock approximation
                 if (e2 .lt. TLJ_eamb) e2 = TLJ_eamb
                 etotal = d*(e2+0.5_RFREAL*(ux*ux + vy*vy + wz*wz))
                 p = e

                 pCv(CV_MIXT_DENS,icg) = d
                 pCv(CV_MIXT_XMOM,icg) = d*ux
                 pCv(CV_MIXT_YMOM,icg) = d*vy
                 pCv(CV_MIXT_ZMOM,icg) = d*wz
                 pCv(CV_MIXT_ENER,icg) = etotal

                 ! We need these for plotting solutions
                 !pDv(DV_MIXT_PRES,icg) = p
                 !pDv(DV_MIXT_TEMP,icg) = RFLU_JWL_T_PR(pRegion,p,d,e2,1.0_RFREAL)
                 !pDv(DV_MIXT_SOUN,icg) = RFLU_JWL_C_ER(pRegion,e2,d,1.0_RFREAL)
                 T = RFLU_JWL_T_PR(pRegion,p,d,e2,1.0_RFREAL)
                 a = RFLU_JWL_C_ER(pRegion,e2,d,1.0_RFREAL)

               ELSE
                 ! Unshocked Air
                 d = pMixtInput%prepRealVal3
                 u = 0.0_RFREAL
                 v = 0.0_RFREAL
                 w = 0.0_RFREAL
                 p = pMixtInput%prepRealVal2

                 mw = pGv(GV_MIXT_MOL,indMol*icg)
                 cp = pGv(GV_MIXT_CP ,indCp *icg)
                 gc = MixtPerf_R_M(mw)
                 g  = MixtPerf_G_CpR(cp,gc)

                 e2 = p/(d*0.4_RFREAL)

                 pCv(CV_MIXT_DENS,icg) = d
                 pCv(CV_MIXT_XMOM,icg) = d*u
                 pCv(CV_MIXT_YMOM,icg) = d*v
                 pCv(CV_MIXT_ZMOM,icg) = d*w
                 pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)

                 ! We need these for plotting solutions
                 !pDv(DV_MIXT_PRES,icg) = p
                 !pDv(DV_MIXT_TEMP,icg) = 300.0_RFREAL
                 !pDv(DV_MIXT_SOUN,icg) = MixtPerf_C_DGP(d,g,p)
                 T = p/(gc*d)           ! T = p/(\rho*Rgas)
                 a = MixtPerf_C_DGP(d,g,p)
                 Y = 0.0_RFREAL
                 YR = 0.0_RFREAL

              ENDIF

              IF (IsNaN(pCv(CV_MIXT_ENER,icg)) .EQV. .TRUE.) THEN
                 print*,'TLJ Input leads to NaN'
                 print*,'    e,d',pCv(CV_MIXT_ENER,icg),d,icg
                 CALL ErrorStop(global,ERR_INVALID_VALUE,1589,'Invalid val InitFlow')
              ENDIF
              IF (IsNaN(p) .EQV. .TRUE.) THEN
                 print*,'TLJ Input leads to NaN'
                 print*,'    p,d',p,d,icg
                 CALL ErrorStop(global,ERR_INVALID_VALUE,1594,'Invalid val InitFlow')
              ENDIF

              ! Josh Garno - WRITE INIT Check File
              IF ( radius < ro ) THEN
                 OPEN(1119,file='INIT_Centerline.dat',form='formatted',status='old' &
                          ,position='append')
                 WRITE(1119,'(11(E23.16,1X),I8,1X,E23.16)') &
                            radius,theta,z,d,u,e2,p,T,a,Y,YR,icg,maxval(xData(:))
                 CLOSE(1119)
              ENDIF ! y
   
           END DO ! icg

           DEALLOCATE(xData,rData,eData,uData,Ydata,eData2,YRdata)

           ENDIF ! RB case
           ENDIF ! JWL case
           !+++++++++++++++
   

! ------------------------------------------------------------------------------
!       Acoustic flow
! ------------------------------------------------------------------------------

        CASE ( "acoustic" ) 
          A1 = pMixtInput%prepRealVal1
          A2 = pMixtInput%prepRealVal2
          Mo = pMixtInput%prepRealVal3
          ro = pMixtInput%prepRealVal4
          po = pMixtInput%prepRealVal5

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE           
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              t = global%currentTime
              CALL RFLU_ComputeExactFlowAcoustic(global,x,t,ro,po,Mo,g, &
                                                 A1,A2,d,u,v,w,p)

              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w

              t = global%currentTime + global%dtImposed/2.0_RFREAL
              CALL RFLU_ComputeExactFlowAcoustic(global,x,t,ro,po,Mo,g, &
                                                 A1,A2,d,u,v,w,p)

              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_PRES,icg) = p

              t = global%currentTime - global%dtImposed/2.0_RFREAL
              CALL RFLU_ComputeExactFlowAcoustic(global,x,t,ro,po,Mo,g, &
                                                 A1,A2,d,u,v,w,p)

              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p
            ELSE
              t = global%currentTime
              CALL RFLU_ComputeExactFlowAcoustic(global,x,t,ro,po,Mo,g, &
                                                 A1,A2,d,u,v,w,p)

              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
            END IF ! solverType
          END DO ! icg


!
! ------------------------------------------------------------------------------
!       Channel acoustics. NOTE the channel is assumed to have the x-coordinate 
!       running down the axis. 
! ------------------------------------------------------------------------------

        CASE ( "channelacoust" )
          CALL RFLU_GetParamsHardCodePAcoust(pTot,aTot)
        
          dTot = MixtPerf_D_CGP(aTot,gRef,pTot)        
                                          
          Lx = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert)) 
          Ly = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert)) 
          Lz = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert)) 

! TEMPORARY : Hardcoded values for parallel
          Lx = 1.0_RFREAL 
          Ly = 1.0_RFREAL 
          Lz = 1.0_RFREAL 

          n1 = pMixtInput%prepIntVal1
          n2 = pMixtInput%prepIntVal2
          n3 = pMixtInput%prepIntVal3
                  
          const = MAX(pMixtInput%prepRealVal1,0.0_RFREAL)        
                                                                 
          omega = aTot*SQRT((n1*global%pi/Lx)**2 + (n2*global%pi/Ly)**2 + (n3*global%pi/Lz)**2)         
        
          IF ( global%verbLevel > VERBOSE_LOW ) THEN           
            WRITE(STDOUT,'(A,5X,A,3(1X,I2))') SOLVER_NAME, &
                  'Mode:',n1,n2,n3                   
            WRITE(STDOUT,'(A,5X,A,3(1X,E13.6))') SOLVER_NAME, &
                  'Size of domain (m):       ',Lx,Ly,Lz
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Total density (kg/m^3):   ',dTot
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Total pressure (N/m^2):   ',pTot            
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Angular frequency (rad/s):',omega
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Constant (-):             ',const                  
          END IF ! global%verbLevel                         
       
          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            DO icg = 1,pGrid%nCellsTot
              x = pGrid%cofg(XCOORD,icg)
              y = pGrid%cofg(YCOORD,icg)
              z = pGrid%cofg(ZCOORD,icg)

              t = global%currentTime

              CALL RFLU_ComputeExactFlowCAcoust(global,x,y,z,global%currentTime, &
                                                Lx,Ly,Lz,n1,n2,n3,omega, &
                                                dTot,pTot,aTot,const,d,u,v,w,p) 
                                    
              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w
            
              t = global%currentTime + global%dtImposed/2.0_RFREAL

              CALL RFLU_ComputeExactFlowCAcoust(global,x,y,z,global%currentTime, &
                                                Lx,Ly,Lz,n1,n2,n3,omega, &
                                                dTot,pTot,aTot,const,d,u,v,w,p) 
             
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_PRES,icg) = p

              t = global%currentTime - global%dtImposed/2.0_RFREAL

              CALL RFLU_ComputeExactFlowCAcoust(global,x,y,z,global%currentTime, &
                                                Lx,Ly,Lz,n1,n2,n3,omega, &
                                                dTot,pTot,aTot,const,d,u,v,w,p) 

              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p
            END DO ! icg    
          ELSE
            DO icg = 1,pGrid%nCellsTot
              x = pGrid%cofg(XCOORD,icg)
              y = pGrid%cofg(YCOORD,icg)
              z = pGrid%cofg(ZCOORD,icg)

              CALL RFLU_ComputeExactFlowCAcoust(global,x,y,z,global%currentTime, &
                                                Lx,Ly,Lz,n1,n2,n3,omega, &
                                                dTot,pTot,aTot,const,d,u,v,w,p) 
                                    
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,gRef,p,u,v,w)               
            END DO ! icg    
          END IF ! solverType
 
! Subbu - Init cyldet case
! ------------------------------------------------------------------------------
!       Detonation wave in 3D - Cylindrical domain
! ------------------------------------------------------------------------------
          CASE ( "cyldet" )

          IF (global%ifReadFromFile) THEN
            iFileName1 = pMixtInput%rBurnFile
            n = 0 
! Begin - Count no.of lines in the RBA data file
  
            OPEN(IF_READRB, FILE=iFileName1, FORM='FORMATTED', &
                 IOSTAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_FILE_OPEN,1785,iFileName1)
            END IF ! global%error
                
! Infinite DO
  
            DO 
              n = n+1
              READ(IF_READRB,'(A)') endString 
              IF ( endString == TRIM('# END') ) THEN
                CLOSE(IF_READRB)
                EXIT
              END IF ! name
! Guard against infinite loop 
              IF ( n >= LIMIT_INFINITE_LOOP ) THEN 
                CALL ErrorStop(global,ERR_INFINITE_LOOP,1799)
              END IF ! n
            END DO ! Infinite DO
  
! End - Counting lines in the file
           
! Allocate arrays to store column data from file
  
            ALLOCATE (radData(n-1),rData(n-1),uData(n-1),eData(n-1), &
                      STAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_ALLOCATE,1811,'radData')
            END IF ! global%error
  
! Begin - Read data file and store data in arrays
  
            OPEN(IF_READRB,FILE=iFileName1,FORM='FORMATTED', &
                 IOSTAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN
               CALL ErrorStop(global,ERR_FILE_OPEN,1820,iFileName1)
            END IF ! global%error
            READ(IF_READRB,*) (radData(icg),rData(icg),uData(icg), &
                                      eData(icg), icg=1,n-1)
            CLOSE(IF_READRB)
  
! End - Reading data from file
  
! Convert radius to m from micrometer 
! This is specific to RocSDT output
            radData = radData*1.0E-6_RFREAL
  
          END IF ! global%ifReadFromFile

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)

            IF (ABS(pMixtInput%prepIntVal1) == 0) THEN  !n=0 ; 2D(r-z)
             radius  = y
            ELSE                                        !2D(r-theta) or 3D(r-theta-z)
             radius  = SQRT(x**2+y**2)
             theta = DATAN2(y,x)
            END IF

            ro = pMixtInput%prepRealVal1 ! Diaphragm loc
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            IF (radius .GT. ro) THEN ! Air - Perfect gas
              p = pMixtInput%prepRealVal2
              d = pMixtInput%prepRealVal3
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)
              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
              e  = MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
            ELSE
              d = pMixtInput%prepRealVal4*pMixtInput%prepRealVal3
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              !ExplEner = 7.0e9_RFREAL
              pCv(CV_MIXT_ENER,icg) = pMixtInput%prepRealVal13

! Interpolate from the reactive burn data
              IF (global%ifReadFromFile) THEN
loop:           DO b=1,n-1
                 i = b
  
! Check if radius falls beyond the minimum bound of the input data 

                 IF (radData(1) .GT. radius) THEN
                  WRITE(*,*) 'Input value lies outside the range...rHigh'!,radius(1),i
                  rHigh = radData(1)
                  rLow  = radData(2)
                  iHigh = 2
                  iLow  = 1
                  EXIT loop
                 END IF
  
! End check
  
                 IF ( radData(i) .LT. radius ) THEN
                  rLow = radData(i)
                  iLow = i
                 ELSE IF ( ABS(radius-radData(i)) .LE. 1.0E-15_RFREAL ) THEN
                  rLow  = radius
                  rHigh = radData(i+1)
                  iLow  = i+1
                  iHigh = i
                  EXIT loop
                 END IF
  
                 i=i+1
  
! Check if radius falls beyond the maximum bound of the input data 
  
                 IF (i .GT. n-1) THEN
                  rHigh = radData(i-1)
                  rLow  = radData(i-2)
                  iHigh = i-1
                  iLow  = i-2
                  EXIT loop
                 END IF ! (i .GT. n-1)
  
! End check
                 IF (radData(i) .GT. radius)THEN
                  rHigh = radData(i)
                  iHigh = i
                  EXIT
                 ELSE IF ( ABS(radius-radData(i)) .LE. 1.0E-15_RFREAL ) THEN
                  rLow  = radius
                  rHigh = radData(i+1)
                  iHigh = i
                  EXIT loop
                 END IF
  
                END DO loop ! loop 
   
                wLow  = ABS(radius-rLow)
                wHigh = ABS(rHigh-radius) 
                wInv  = 1.0/(wLow+wHigh)
                    d = ( wLow *rData(iLow ) + &
                          wHigh*rData(iHigh) )*wInv
                   ur = ( wLow *uData(iLow ) + &
                          wHigh*uData(iHigh) )*wInv
                    u = ur*COS(theta) 
                    v = ur*SIN(theta)
                    w = 0.0_RFREAL
                    e = ( wLow *eData(iLow ) + &
                          wHigh*eData(iHigh) )*wInv
                pCv(CV_MIXT_ENER,icg) = d*Mixt_Eo_eUVW(e,u,v,w)
              END IF ! global%ifReadFromFile  

            END IF ! radius

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
!            pCv(CV_MIXT_ENER,icg) = d*Mixt_Eo_eUVW(e,u,v,w) 


          END DO ! icg

! Adding specific mode number/random mode numbers to base flow

          flag        = pMixtInput%prepIntVal3
          gaussAmp    = pMixtInput%prepRealVal6
          gaussCenter = pMixtInput%prepRealVal1
          gaussWidth  = pMixtInput%prepRealVal5*gaussCenter

          IF (flag == 1) THEN ! Random init
           ALLOCATE(randNum(pGrid%nCellsTot))
           DO icg = 1,pGrid%nCellsTot
             randNum(icg) = icg
           END DO
           CALL RANDOM_SEED()
           CALL RANDOM_NUMBER(randNum)

           DO icg = 1,pGrid%nCellsTot
             x = pGrid%cofg(XCOORD,icg)
             y = pGrid%cofg(YCOORD,icg)
             z = pGrid%cofg(ZCOORD,icg)
             radius  = SQRT(x**2+y**2)
             perturb = gaussAmp*EXP(-((radius-gaussCenter)**2) &
                      /(2.0_RFREAL*gaussWidth**2))*randNum(icg)

             pCv(CV_MIXT_DENS,icg) = pCv(CV_MIXT_DENS,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_XMOM,icg) = pCv(CV_MIXT_XMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_YMOM,icg) = pCv(CV_MIXT_YMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_ZMOM,icg) = pCv(CV_MIXT_ZMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_ENER,icg) = pCv(CV_MIXT_ENER,icg)*(1.0_RFREAL+perturb)
           END DO ! icg

          ELSE
           m = pMixtInput%prepIntVal1 ! Wave no in theta
           n = pMixtInput%prepIntVal2 ! Wave number in z
           DO icg = 1,pGrid%nCellsTot
             x = pGrid%cofg(XCOORD,icg)
             y = pGrid%cofg(YCOORD,icg)
             z = pGrid%cofg(ZCOORD,icg)
             IF (ABS(m) == 0) THEN ! m=0; 2D-(r-z) 
              radius  = y
              zz      = x
              the     = DATAN2(z,y)
             ELSE                  ! 2D- (r-theta) OR 3D-(r-theta-z)
              radius  = SQRT(x**2+y**2)
              zz      = z
              the     = DATAN2(y,x)
             END IF
             !perturb = gaussAmp*EXP(-((radius-gaussCenter)**2) &
             !        /(2.0_RFREAL*gaussWidth**2))*DCOS(m*the)*DCOS(n*zz/gaussCenter) ! gaussCenter = r0
             IF (radius .LE. pMixtInput%prepRealVal1) THEN
               perturb = gaussAmp*DCOS(m*the)
             ELSE
               perturb = 0
             END IF
             pCv(CV_MIXT_DENS,icg) = pCv(CV_MIXT_DENS,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_XMOM,icg) = pCv(CV_MIXT_XMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_YMOM,icg) = pCv(CV_MIXT_YMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_ZMOM,icg) = pCv(CV_MIXT_ZMOM,icg)*(1.0_RFREAL+perturb)
             pCv(CV_MIXT_ENER,icg) = pCv(CV_MIXT_ENER,icg)*(1.0_RFREAL+perturb)
           END DO ! icg
          END IF

! Rahul - Deallocate arrays storing the PBA file input
          IF (global%ifReadFromFile) THEN
            DEALLOCATE (radData,rData,uData,eData,STAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_DEALLOCATE,2019,'rBurnData')
            END IF ! global%error
          END IF ! global%ifReadFromFile

! Subbu - End init cyldet case


! ------------------------------------------------------------------------------
!       Contact Interface, Added by CRN
! ------------------------------------------------------------------------------

        CASE ( "ctint" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE
              d = pMixtInput%prepRealVal5
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            END IF ! x

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)

            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
! Saptarshi - 01/10/2015
!       Square ShockTube problems
!       A-B-A TYPE MODEL        
! ------------------------------------------------------------------------------

        CASE ( "shktb" )
          DO icg = 1,pGrid%nCellsTot           !Loop over all cells in domain
            x = pGrid%cofg(XCOORD,icg)       !Extract x coordinate of cellcentroid
            y = pGrid%cofg(YCOORD,icg)       !Extract y coordinate of cellcentroid
            z = pGrid%cofg(ZCOORD,icg)       !Extract z coordinate of cellcentroid

!**************************************************************************************
!IVAL = 0 EGG CRATE INTERFACE 
!     = 1 MULTIMODE COSINE SURFACE AT THE INTERFACE THOUGH SAME K VALUE 
!     = 2 SINGLE MODE COSINE SURFACE AT THE INTERFACE ALONG  Y
!     = 3 FLAT SURFACE
!     = 4 FRONT FLAT - INVERSED CHEVRON ENDS
!***************************************************************************************           

            IF( pMixtInput%prepIntVal1 == 0 ) THEN ! Egg Crate
               dx = pMixtInput%prepRealVal10 &
                    *abs(sin(pMixtInput%prepRealVal11*y) &
                    *sin(pMixtInput%prepRealVal11*z))
               coordSinInt    = pMixtInput%prepRealVal5 + dx
               coordSinIntEnd = pMixtInput%prepRealVal12 + dx
            ELSE IF ( pMixtInput%prepIntVal1 == 1 ) THEN  ! MultiMode Cosine surface
              dy = pMixtInput%prepRealVal10 &
                  *(cos(pMixtInput%prepRealVal11*y))
              dz = pMixtInput%prepRealVal10 &
                  *(cos(pMixtInput%prepRealVal11*z))
              dx = dy+dz
              coordSinInt    = pMixtInput%prepRealVal5 + dx
              coordSinIntEnd = pMixtInput%prepRealVal12 + dx
            ELSE IF ( pMixtInput%prepIntVal1 == 2 ) THEN  ! SingleMode Cosine Surface along Y
             dx = pMixtInput%prepRealVal10 &
                  *(cos(pMixtInput%prepRealVal11*y))
             coordSinInt    = pMixtInput%prepRealVal5 + dx
             coordSinIntEnd = pMixtInput%prepRealVal12 + dx
            ELSE IF ( pMixtInput%prepIntVal1 == 3) THEN ! For Flat Surface 
              dx = 0.000000000000000
              coordSinInt    = pMixtInput%prepRealVal5 + dx
              coordSinIntEnd = pMixtInput%prepRealVal12 + dx
            END IF

!*******************************************************************************************
!INTERFACE CONDITION ENDS

!*******************************************************************************************

            IF ( x .LT. pMixtInput%prepRealVal1 ) THEN !Shocked air portion
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE IF( x .LT. pMixtInput%prepRealVal5 ) THEN !Unshocked air portion
              d = pMixtInput%prepRealVal6
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal7
            ELSE IF( x .LT. coordSinInt ) THEN      !Unshocked air portion
              d = pMixtInput%prepRealVal6
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal7
            ELSE IF(x .LT. pMixtInput%prepRealVal12)THEN    !Unshocked SF6 portion
              d = pMixtInput%prepRealVal8
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal9
            ELSE IF (x .LT. coordSinIntEnd )THEN    !Unshocked SF6 portion
              d = pMixtInput%prepRealVal8
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal9
            ELSE                             !Unshocked Air
              d = pMixtInput%prepRealVal6
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal7
            END IF ! ( x < pMixtInput%prepRealVal1 

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)

            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)

          END DO ! icg

          ! Saptarshi - 01/10/2013 - end 

! ------------------------------------------------------------------------------
!       Potential flow around cylinder.
! ------------------------------------------------------------------------------

        CASE ( "cylpotential" ) 
          ro  = pMixtInput%prepRealVal1
          uo  = pMixtInput%prepRealVal2
          po  = pMixtInput%prepRealVal3
          R   = pMixtInput%prepRealVal4

!Subbu: Open data file

            ! ******************************************************************************
            ! Open file for data
            ! ******************************************************************************

            WRITE(iFileName1,'(A,1PE11.5,A)') TRIM(global%outDir)// &
                              TRIM(global%casename)// &
                              '.surf_',global%currentTime,'.dat'

            OPEN(IF_EXTR_DATA1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
            IOSTAT=errorFlag)
 !End

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

 !Subbu: BEWARE!!!!!!!!!! 
 !TEMPORARY ONLY :x and y ARE SWAPPED
            xo = x
            yo = y
            !temp = xo
            !xo = yo
            !yo = -temp
 !END TEMPORARY

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE           
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

            CALL RFLU_ComputeExactFlowCylPotential(global,xo,yo,ro,uo,po,R, &
                                                   d,u,v,w,p)

 !Subbu: BEWARE!!!!!!!!!! 
 !TEMPORARY ONLY :u and v ARE SWAPPED
            !temp = u
            !u = -v
            !v = temp
 !END TEMPORARY

 !Subbu: Write data to file
             
             WRITE(IF_EXTR_DATA1,'(6(1X,E23.16))') xo,yo,x,y,u,v
    
 !End Subbu
            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_PRES,icg) = p

              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p
            ELSE
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
            END IF ! solverType
          END DO ! icg

 ! Subbu: Close file

             !******************************************************************************
             ! Close file
             ! ******************************************************************************

              CLOSE(IF_EXTR_DATA1,IOSTAT=errorFlag)
 ! End Subbu

! ------------------------------------------------------------------------------
!       Incompressible flow around cylinder near a wall.
! ------------------------------------------------------------------------------

        CASE ( "cylwall" )
          ro  = pMixtInput%prepRealVal1 ! Uniform Density
          uo  = pMixtInput%prepRealVal2 ! Free-stream Velocity
          po  = pMixtInput%prepRealVal3 ! Free-stream Pressure
          R   = pMixtInput%prepRealVal4 ! Radius of cylinder
          LD  = pMixtInput%prepRealVal5 ! L/D ratio
          
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

         !   CALL RFLU_ComputeExactFlowCylWall(global,x,y,ro,uo,po,R,LD, &
         !                                          d,u,v,w,p)

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_PRES,icg) = p

              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p
            ELSE
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
            END IF ! solverType
          END DO ! icg


! Subbu - Sep 2011
! ------------------------------------------------------------------------------
!       Detonation wave in 3D
! ------------------------------------------------------------------------------
          CASE ( "sphdet", "detwav3D" )

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)

            !radius  = SQRT(x**2+y**2+z**2)
            radius = x

            IF ( radius > pMixtInput%prepRealVal1 ) THEN
              p = pMixtInput%prepRealVal2
              d = pMixtInput%prepRealVal3
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
            ELSE
              p = pMixtInput%prepRealVal4*pMixtInput%prepRealVal2
              d = pMixtInput%prepRealVal5*pMixtInput%prepRealVal3
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
            END IF ! rad
    

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)

            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! End Subbu - Sep 2011


! ------------------------------------------------------------------------------
!       Diffracting shock issuing from open-ended shock tube
! ------------------------------------------------------------------------------
! ----- Without backchamber ----------------------------------------------------

        CASE ( "ds_7p25"   , "ds_20p0"   , "ds_50p0"   , "ds_125p0"   , & 
               "ds_7p25_v2", "ds_20p0_v2", "ds_50p0_v2", "ds_125p0_v2", &
               "ds_7p25_v3", "ds_20p0_v3", "ds_50p0_v3", "ds_125p0_v3", &
               "ds_7p25_v4", "ds_20p0_v4", "ds_50p0_v4", "ds_125p0_v4", &
               "ds_7p25_v5", "ds_20p0_v5", "ds_50p0_v5", "ds_125p0_v5" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal6
            END IF ! x

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg
  
! ----- With backchamber -----------------------------------------------------	  

        CASE ( "ds_7p25_v6", "ds_20p0_v6", "ds_50p0_v6", "ds_125p0_v6" )
          DO icg = 1,pGrid%nCellsTot
            x =     pGrid%cofg(XCOORD,icg)
            y = ABS(pGrid%cofg(YCOORD,icg))

            IF ( (x < pMixtInput%prepRealVal1) .AND. & 
                 (y < pMixtInput%prepRealVal2) ) THEN
              d = pMixtInput%prepRealVal3
              u = pMixtInput%prepRealVal4
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal5
            ELSE
              d = pMixtInput%prepRealVal6
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal7
            END IF ! x

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg	  

! ------------------------------------------------------------------------------
!       hexahedral channel mesh
! ------------------------------------------------------------------------------

        CASE ( "vortexNSCBC" ) 

!         a = speed of sound
!         C = -0.0005_RFREAL*a
!         Rc = 0.15_RFREAL
!         u0 = 1.1_RFREAL*a
!         pinf = pMixtInput%prepRealVal3

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

!            radius = (x*x+y*y)**(0.5_RFREAL)
!         psi = C*exp(-(x*x+y*y)/(2.0_RFREAL*Rc*Rc))

!            IF ( radius < Rc ) THEN
!              d = pMixtInput%prepRealVal2
!              u = u0 - x*psi/(d*Rc*Rc)
!              v = y*psi/(d*Rc*Rc)
!              w = 0.0_RFREAL
!              p = pinf + d*C*psi/(Rc*Rc)
!            ELSE
!              d = pMixtInput%prepRealVal2
!              u = u0
!              v = 0.0_RFREAL 
!              w = 0.0_RFREAL
!              p = pinf
!            END IF ! radius

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Generic compressible gravity current
! ------------------------------------------------------------------------------

        CASE ( "gcgc" ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            IF ( (x < pMixtInput%prepRealVal1) .AND. & 
                 (y < pMixtInput%prepRealVal2) ) THEN 
              d = pMixtInput%prepRealVal3
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE 
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            END IF ! x
            
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)            
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg     

  
! ------------------------------------------------------------------------------
!       Gaussian pulse
! ------------------------------------------------------------------------------

        CASE ( "gaussianpulse" ) 
          A_c  = 0.00001_RFREAL
! TEMPORARY : Have to decide some standard means to initialize
!          uo_c = 0.1_RFREAL
          uo_c = 0.0_RFREAL
  
          L  = pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c
  
          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4
  
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE           
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

            u = uo
            p = po
            d = ro
            v = 0.0_RFREAL
            w = 0.0_RFREAL

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              IF ( x > -L/2.0_RFREAL .AND. x < L/2.0_RFREAL ) THEN
                x = x + L/2.0_RFREAL
                CALL RFLU_ComputeExactFlowGaussianPulse(global,x,ro,uo,po, &
                                                        c,g,L,A,d,u,v,w,p)
                x = x - L/2.0_RFREAL
              END IF ! x

              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w

              x = x - (uo+c)*global%dtImposed/2.0_RFREAL

              IF ( x > -L/2.0_RFREAL .AND. x < L/2.0_RFREAL ) THEN
                x = x + L/2.0_RFREAL
                CALL RFLU_ComputeExactFlowGaussianPulse(global,x,ro,uo,po, &
                                                        c,g,L,A,d,u,v,w,p)
                x = x - L/2.0_RFREAL
              END IF ! x

              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_PRES,icg) = p

              x = x + (uo+c)*global%dtImposed

              IF ( x > -L/2.0_RFREAL .AND. x < L/2.0_RFREAL ) THEN
                x = x + L/2.0_RFREAL
                CALL RFLU_ComputeExactFlowGaussianPulse(global,x,ro,uo,po, &
                                                        c,g,L,A,d,u,v,w,p)
                x = x - L/2.0_RFREAL
              END IF ! x

              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p

            ELSE
              IF ( x > -L/2.0_RFREAL .AND. x < L/2.0_RFREAL ) THEN
                x = x + L/2.0_RFREAL
                CALL RFLU_ComputeExactFlowGaussianPulse(global,x,ro,uo,po, &
                                                        c,g,L,A,d,u,v,w,p)
                x = x - L/2.0_RFREAL
              END IF ! x

              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
            END IF ! solverType
          END DO ! icg

! ------------------------------------------------------------------------------
!       Generic multiphase wall jet
! ------------------------------------------------------------------------------

        CASE ( "gmpjet" ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)

            IF ( (    x  < 0.00_RFREAL) .AND. & 
                 (ABS(y) < 0.55_RFREAL) .AND. & 
                 (ABS(z) < 0.55_RFREAL) ) THEN 
              d = 34.7418238572_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 5.0E+6_RFREAL
            ELSE 
              d = 1.20903522664_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            END IF ! x
            
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)            
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg     	  

! ------------------------------------------------------------------------------
!       Gradient test 
! ------------------------------------------------------------------------------

! ----- Linear function --------------------------------------------------------

        CASE ( "gtlin" ) 
          xMin = MINVAL(pGrid%cofg(XCOORD,1:pGrid%nCellsTot))        
          yMin = MINVAL(pGrid%cofg(YCOORD,1:pGrid%nCellsTot))
          zMin = MINVAL(pGrid%cofg(ZCOORD,1:pGrid%nCellsTot))
        
          CALL RFLU_SetExactFlowLinear(xMin,yMin,zMin,1,dMin,gx,gy,gz)
          CALL RFLU_SetExactFlowLinear(xMin,yMin,zMin,5,pMin,gx,gy,gz)
                        
          IF ( dMin < 0.0_RFREAL ) THEN 
            dOffs = -2.0_RFREAL*dMin
          ELSE 
            dOffs = 0.0_RFREAL
          END IF ! dMin  
          
          IF ( pMin < 0.0_RFREAL ) THEN 
            pOffs = -2.0_RFREAL*pMin
          ELSE 
            pOffs = 0.0_RFREAL
          END IF ! pMin                        
                                                                      
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)                        
     
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            CALL RFLU_SetExactFlowLinear(x,y,z,1,d,gx,gy,gz)
            CALL RFLU_SetExactFlowLinear(x,y,z,2,u,gx,gy,gz) 
            CALL RFLU_SetExactFlowLinear(x,y,z,3,v,gx,gy,gz)
            CALL RFLU_SetExactFlowLinear(x,y,z,4,w,gx,gy,gz)
            CALL RFLU_SetExactFlowLinear(x,y,z,5,p,gx,gy,gz) 
                                                                           
            d = d + dOffs 
            p = p + pOffs                  
                               
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ----- Trigonometric function -------------------------------------------------

        CASE ( "gttri" ) 
          nx = pRegion%mixtInput%prepRealVal1
          ny = pRegion%mixtInput%prepRealVal2
          nz = pRegion%mixtInput%prepRealVal3            
        
          dOffs = 2.0_RFREAL ! Assume amplitude is unity
          pOffs = 2.0_RFREAL ! Assume amplitude is unity
                                                                      
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)                        
     
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,1,d,gx,gy,gz)
            CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,2,u,gx,gy,gz) 
            CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,3,v,gx,gy,gz)
            CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,4,w,gx,gy,gz)
            CALL RFLU_SetExactFlowTrig(global,nx,ny,nz,x,y,z,5,p,gx,gy,gz)
                                                                            
            d = d + dOffs 
            p = p + pOffs                  
                               
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Jet Flow
! ------------------------------------------------------------------------------

        CASE ( "jet" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,2730, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,2737, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            yp = 8.9E-05_RFREAL + 0.0175_RFREAL*(x - 2.5E-04_RFREAL)

            IF ( (y <= 8.9E-05_RFREAL) .AND. (y >= -8.9E-05_RFREAL) .AND. &
                 (x >= 0.0E+00_RFREAL) .AND. (x <=  2.5E-04_RFREAL) ) THEN

              rl = 880.0_RFREAL
              pl = 1.0_RFREAL
              rg = 11.69_RFREAL
              pg = 0.0_RFREAL
              rv = 0.722_RFREAL
              pv = 0.0_RFREAL
  
              d  = rl*pl + rg*pg + rv*pv
              H  = 8.9E-05_RFREAL
              u  = 300.0*(1.0_RFREAL - (y*y)/(H*H))

              IF ( ABS(y) < 5.5E-05_RFREAL ) THEN
                u = 300.0_RFREAL
              END IF ! ABS(y)

              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 288.0_RFREAL
            ELSE
              rl = 880.0_RFREAL
              pl = 0.0_RFREAL
              rg = 11.69_RFREAL
              pg = 1.0_RFREAL
              rv = 0.722_RFREAL
              pv = 0.0_RFREAL
          
              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 288.0_RFREAL
            END IF ! y
         
            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w
 
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Kieffer jet. NOTE configuration 1 has jet axis along z, configuration 2
!       has jet axis along x.
! ------------------------------------------------------------------------------

        CASE ( "kjet" ) 
          DO icg = 1,pGrid%nCellsTot
            z = pGrid%cofg(ZCOORD,icg)
            
            IF ( z < -0.0254_RFREAL ) THEN 
              d = 7.96_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 7.25E+5_RFREAL
            ELSE 
              d = 1.16871466381_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            END IF ! z        

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)             
          END DO ! icg

        CASE ( "kjet2","kjet2v3","kjet2v4","kjet2v5" ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            
            IF ( x < -0.0127_RFREAL ) THEN 
              d = 7.96_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 7.25E+5_RFREAL
            ELSE 
              d = 1.1220002356_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            END IF ! x        

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)             
          END DO ! icg

        CASE ( "kjet2v3mp","kjet2v4mp","kjet2v5mp" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE
              d = pMixtInput%prepRealVal5
              u = pMixtInput%prepRealVal6
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal7
            END IF ! x

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)

            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Multiphase Shocktube: Water-Air
! ------------------------------------------------------------------------------

        CASE ( "MShock_H2O_Air001" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,2905, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,2912, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq
          
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.7_RFREAL ) THEN
              rl = 1500.0_RFREAL
              pl = 1.0_RFREAL
              rg = 0.0_RFREAL
              pg = 0.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 515.0_RFREAL
            ELSE
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 150.0_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 2.25_RFREAL 
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d 
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

        CASE ( "MShock_H2O_Air002" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,2967, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,2974, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq        
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.7_RFREAL ) THEN
              rl = 1500.0_RFREAL
              pl = 1.0_RFREAL
              rg = 0.0_RFREAL
              pg = 0.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 515.0_RFREAL
            ELSE
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 50.0_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 6.73_RFREAL
            END IF ! x 

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Multiphase Shocktube: Air-Air-He
! ------------------------------------------------------------------------------

        CASE ( "MShock_Air_Air_He001" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,3033, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,3040, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.25_RFREAL ) THEN
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 1.7017_RFREAL 
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 98.956_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 307.13_RFREAL
            ELSE IF ( (x > 0.25_RFREAL) .AND. (x < 0.5_RFREAL) ) THEN
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 1.2763       
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 273.0_RFREAL 
            ELSE
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 0.0_RFREAL
              pg = 0.0_RFREAL
              rv = 0.1760_RFREAL
              pv = 1.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 273.5588_RFREAL
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Moving shock
! ------------------------------------------------------------------------------

        CASE ( "mvsh" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE 
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal6
            END IF ! z
            
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)            
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg     

! ------------------------------------------------------------------------------
!       NSCBC farfield test
! ------------------------------------------------------------------------------

        CASE ( "farf" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.0001_RFREAL
          A_c  = 0.01_RFREAL
          uo_c = 0.1_RFREAL
          uo_c = 0.0_RFREAL

          r1 = 0.006125_RFREAL
          r2 = 0.091875_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            radius = SQRT(x*x+y*y)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            um= uo + A*(SIN(global%pi*(radius-r1)/(r2-r1)))**9.0_RFREAL
            u = um*COS(ATAN2(y,x))
            p = po + (um - uo)*(ro*c)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = um*SIN(ATAN2(y,x))
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       NSCBC tests
! ------------------------------------------------------------------------------

        CASE ( "nscbc1" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.0001_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            u = uo + A*(SIN(global%pi*x/L))**9.0_RFREAL
            p = po + (u - uo)*(ro*c)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = 0.0_RFREAL
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc2" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.01_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            u = uo 
            p = po + ro*c*A*((SIN(global%pi*x/L))**9.0_RFREAL)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = 0.0_RFREAL
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc3" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.0001_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            y = pGrid%cofg(YCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            v = uo + A*(SIN(global%pi*y/L))**9.0_RFREAL
            p = po + (v - uo)*(ro*c)
            d = ro*(p/po)**(1.0_RFREAL/g)
            u = 0.0_RFREAL
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc4" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.01_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            y = pGrid%cofg(YCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            u = 0.0_RFREAL 
            p = po + ro*c*A*((SIN(global%pi*y/L))**9.0_RFREAL)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = uo
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc5" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.0001_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            xx = ((x*x+y*y)**0.5)*COS(ATAN2(y,x)-global%pi/4.0_RFREAL)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
            
            um = uo + A*(SIN(global%pi*xx/L))**9.0_RFREAL
                   
            u = um*COS(global%pi/4.0_RFREAL)
            p = po + (um - uo)*(ro*c)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = um*SIN(global%pi/4.0_RFREAL)
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc6" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.01_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            xx = ((x*x+y*y)**0.5)*COS(ATAN2(y,x)-global%pi/4.0_RFREAL)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
            
            um = uo + A*(SIN(global%pi*xx/L))**9.0_RFREAL
                   
            u = uo*COS(global%pi/4.0_RFREAL)
            p = po + ro*c*A*((SIN(global%pi*xx/L))**9.0_RFREAL)
            d = ro*(p/po)**(1.0_RFREAL/g)
            v = uo*SIN(global%pi/4.0_RFREAL)
            w = 0.0_RFREAL

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

        CASE ( "nscbc7" ) 
          L_l  = 2.0_RFREAL
          A_c  = 0.01_RFREAL
          uo_c = 0.1_RFREAL

          L  = L_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4
          xo = pMixtInput%prepRealVal5

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
      
            x = x + xo
                     
            IF ( x < L ) THEN          
              u = uo 
              p = po + ro*c*A*((SIN(global%pi*x/L))**9.0_RFREAL)
              d = ro*(p/po)**(1.0_RFREAL/g)
              v = 0.0_RFREAL
              w = 0.0_RFREAL
            ELSE
              u = uo 
              p = po
              d = ro
              v = 0.0_RFREAL
              w = 0.0_RFREAL
            END IF ! x

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Shock-tube testing on nscbc
! ------------------------------------------------------------------------------

        CASE ( "nscbc8" ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = 0.0_RFREAL 
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal3
            ELSE 
              d = pMixtInput%prepRealVal4
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal5
            END IF ! x         

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)
                               
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------          
!       Nozzle cavity
! ------------------------------------------------------------------------------

        CASE ( "ncavity" )   
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,3500, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,3507, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
             
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            
            IF ( x <= 3.0E-03_RFREAL ) THEN
              rl = 882.655_RFREAL 
              pl = 1.0_RFREAL
              rg = 24.14836_RFREAL
              pg = 0.0_RFREAL
              rv = 15.547_RFREAL 
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 293.0_RFREAL
            ELSE
              rl = 880.0_RFREAL  
              pl = 0.0_RFREAL
              rg = 24.14836_RFREAL 
              pg = 1.0_RFREAL
              rv = 15.547_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 293_RFREAL
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Proudman-Culick flow. NOTE this problem is two-dimensional and assumed 
!       to lie in the x-y plane, and that the injection boundary is located at 
!       y = -height.
! ------------------------------------------------------------------------------

        CASE ( "onera_c0", "onera_c0_2d_100x50" )
          CALL RFLU_GetParamsHardCodeProudman(dInc,mInj,vInj,pTot)

          height = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            CALL RFLU_ComputeExactFlowProudman(global,x,y,height,dInc,vInj, &
                                               pTot,d,u,v,w,p) 

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w
              pCv(CV_MIXT_PRES,icg) = p

              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p
            ELSE 
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
            END IF ! global%solverType
          END DO ! icg 
          
! ------------------------------------------------------------------------------
!       Culick flow 
! ------------------------------------------------------------------------------

        CASE ( "onera_c0_3d" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)
            z = pGrid%cofg(ZCOORD,icg)

            CALL RFLU_ComputeExactFlowCulick(global,x,y,z,d,u,v,w,p) 

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w
              pCv(CV_MIXT_PRES,icg) = p

              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p
            ELSE 
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
            END IF ! global%solverType
          END DO ! icg            

! ------------------------------------------------------------------------------
!       Pipe acoustics. NOTE the pipe is assumed to have the x-coordinate 
!       running down the axis. 
! ------------------------------------------------------------------------------

        CASE ( "pipeacoust" )
          CALL RFLU_GetParamsHardCodePAcoust(pTot,aTot)
        
          dTot = MixtPerf_D_CGP(aTot,gRef,pTot)        
                                          
          L  = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert)) 
          ro = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))        
                                                        
          im  = MAX(pMixtInput%prepIntVal1,1)
          in  = MAX(pMixtInput%prepIntVal2,1)
          iq  = MAX(pMixtInput%prepIntVal3,1)
          iBc = MAX(MIN(pMixtInput%prepIntVal4,1),0)
                  
          const = MAX(pMixtInput%prepRealVal1,0.0_RFREAL)        
                                                                 
          CALL RFLU_JYZOM(im,iq,dummyReal,etaqm,dummyReal,dummyReal)       
                                                                
          omega = aTot*SQRT((in*global%pi/L)**2 + (etaqm/ro)**2)         
        
          IF ( global%verbLevel > VERBOSE_LOW ) THEN           
            WRITE(STDOUT,'(A,5X,A,1X,I2)'   ) SOLVER_NAME, &
                  'Boundary condition:',iBc          
            WRITE(STDOUT,'(A,5X,A,3(1X,I2))') SOLVER_NAME, &
                  'Mode:',im,in,iq                   
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Total density (kg/m^3):   ',dTot
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Total pressure (N/m^2):   ',pTot            
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Angular frequency (rad/s):',omega
            WRITE(STDOUT,'(A,5X,A,1X,E13.6)') SOLVER_NAME, &
                  'Constant (-):             ',const                  
          END IF ! global%verbLevel                         
        
          IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
            DO icg = 1,pGrid%nCellsTot
              x = pGrid%cofg(XCOORD,icg)
              y = pGrid%cofg(YCOORD,icg)
              z = pGrid%cofg(ZCOORD,icg)

              t = global%currentTime

              CALL RFLU_ComputeExactFlowPAcoust(global,z,y,x,t,L,ro,iBc,im, &
                                                in,iq,etaqm,omega,dTot,pTot, &
                                                aTot,const,d,u,v,w,p) 
                                    
              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w
            
              t = global%currentTime + global%dtImposed/2.0_RFREAL

              CALL RFLU_ComputeExactFlowPAcoust(global,z,y,x,t,L,ro,iBc,im, &
                                                in,iq,etaqm,omega,dTot,pTot, &
                                                aTot,const,d,u,v,w,p) 
                                    
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_PRES,icg) = p

              t = global%currentTime - global%dtImposed/2.0_RFREAL

              CALL RFLU_ComputeExactFlowPAcoust(global,z,y,x,t,L,ro,iBc,im, &
                                                in,iq,etaqm,omega,dTot,pTot, &
                                                aTot,const,d,u,v,w,p) 

              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p
            END DO ! icg    
          ELSE
            DO icg = 1,pGrid%nCellsTot
              x = pGrid%cofg(XCOORD,icg)
              y = pGrid%cofg(YCOORD,icg)
              z = pGrid%cofg(ZCOORD,icg)

              CALL RFLU_ComputeExactFlowPAcoust(global,z,y,x,global%currentTime, &
                                                L,ro,iBc,im,in,iq,etaqm,omega, &
                                                dTot,pTot,aTot,const,d,u,v,w,p) 
                                    
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,gRef,p,u,v,w)               
            END DO ! icg    
          END IF ! solverType

! ------------------------------------------------------------------------------
!       Radial pulse
! ------------------------------------------------------------------------------

        CASE ( "radialpulse" ) 
          A_c  = 0.00001_RFREAL
! TEMPORARY : Have to decide some standard means to initialize
!          uo_c = 0.1_RFREAL
          uo_c = 0.0_RFREAL
  
          L  = pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_c*c
          uo = uo_c*c
  
          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4
  
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            radius = SQRT(x*x+y*y)
            theta  = ATAN2(y,x) 

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE           
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

            u = uo
            p = po
            d = ro
            v = 0.0_RFREAL
            w = 0.0_RFREAL

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              IF ( radius < L ) THEN
                CALL RFLU_ComputeExactFlowRadialPulse(global,radius,theta,ro, &
                                                      uo,po,c,g,L,A,d,u,v,w,p)
              END IF ! radius

              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_PRES,icg) = p

              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p

            ELSE
              IF ( radius < L ) THEN
                CALL RFLU_ComputeExactFlowRadialPulse(global,radius,theta,ro, &
                                                      uo,po,c,g,L,A,d,u,v,w,p)
              END IF ! radius

              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
            END IF ! solverType
          END DO ! icg

! ------------------------------------------------------------------------------
!       Rayleigh problem. Impulsively moving flat plate.
! ------------------------------------------------------------------------------

        CASE ( "rayleighproblem" ) 
          ro  = pMixtInput%prepRealVal1
          uo  = pMixtInput%prepRealVal2
          po  = pMixtInput%prepRealVal3
          muo = global%refVisc
          t   = global%currentTime

          DO icg = 1,pGrid%nCellsTot
            y = pGrid%cofg(YCOORD,icg)

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE           
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

            CALL RFLU_ComputeExactFlowRayleighProblem(global,t,y,ro,uo, &
                                                      po,muo,d,u,v,w,p)

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_PRES,icg) = p

              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p
            ELSE
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
            END IF ! solverType
          END DO ! icg

! ------------------------------------------------------------------------------
!       Ringleb flow. NOTE this problem is two-dimensional and assumed to lie in 
!       the x-y plane and that the exact solution is restricted to gamma = 1.4.
! ------------------------------------------------------------------------------

        CASE ( "ringleb" )
          CALL RFLU_GetParamsHardCodeRingleb(pTot,tTot)

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            CALL RFLU_ComputeExactFlowRingleb(x,y,gcRef,pTot,tTot,d,u,v,w,p)  

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,1.4_RFREAL,p,u,v,w)               
          END DO ! icg    

! ------------------------------------------------------------------------------
!       Shock-bubble interaction (Quirk and Karni 1996)
! ------------------------------------------------------------------------------

        CASE ( "ShockBubble" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,3879, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,3886, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            IF ( ((x-0.085_RFREAL)**2 + y**2) <= 6.25E-04_RFREAL )THEN
              rl = 1000.0_RFREAL
              pl = 0.0_RFREAL
              rg = 1.0_RFREAL
              pg = 0.0_RFREAL
              rv = 0.232127_RFREAL
              pv = 1.0_RFREAL

              d = rl*pl + rg*pg + rv*pv
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              t = 273.0_RFREAL
            ELSE IF ( x < 0.050_RFREAL ) THEN
              rl = 1000.0_RFREAL
              pl = 0.0_RFREAL
              rg = 1.75666_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d = rl*pl + rg*pg + rv*pv
              u = 110.49_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              t = 311.366_RFREAL
            ELSE
              rl = 1000_RFREAL
              pl = 0.0_RFREAL
              rg = 1.2763_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d = rl*pl + rg*pg + rv*pv
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              t = 273.0_RFREAL
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u + v*v + w*w
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCV(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg


! ------------------------------------------------------------------------------
!       Skews diffracting shock
! ------------------------------------------------------------------------------

        CASE ( "skews_ms2p0","skews_ms3p0","skews_ms4p0" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE 
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal6
            END IF ! z
            
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)            
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg           
  
! ------------------------------------------------------------------------------
!       Sommerfeld shock-particle-interaction 
! ------------------------------------------------------------------------------

        CASE ( "somm_spi" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE 
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal6
            END IF ! z
            
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)            
            
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
          END DO ! icg           
  
! ------------------------------------------------------------------------------
!       Sphere shock diffraction
! ------------------------------------------------------------------------------

        CASE ( "sphds" ) 
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE 
              d = pMixtInput%prepRealVal5
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal6
            END IF ! x         
                              
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)                              
                               
            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Potential flow around sphere.
! ------------------------------------------------------------------------------

        CASE ( "sphpotential" ) 
          ro  = pMixtInput%prepRealVal1
          uo  = pMixtInput%prepRealVal2
          po  = pMixtInput%prepRealVal3
          R   = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE           
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

            CALL RFLU_ComputeExactFlowSphPotential(global,x,y,ro,uo,po,R, &
                                                   d,u,v,w,p)

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_PRES,icg) = p

              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p
            ELSE
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
            END IF ! solverType
          END DO ! icg

! ------------------------------------------------------------------------------
!       Shocktubes
! ------------------------------------------------------------------------------

! ----- Generic 1d/2d ----------------------------------------------------------  
  
        CASE ( "stg1d","stg2d" )
          DO icg = 1,pGrid%nCellsTot           
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal4
            ELSE IF ( (x >= pMixtInput%prepRealVal1 ) .AND. &
                      (x  < pMixtInput%prepRealVal11) ) THEN  
              d = pMixtInput%prepRealVal5
              u = pMixtInput%prepRealVal6
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal7
            ELSE
              d = pMixtInput%prepRealVal12
              u = pMixtInput%prepRealVal13
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = pMixtInput%prepRealVal14
            END IF ! x         

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)

            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg 

! ----- Sod case 1 -------------------------------------------------------------  
  
        CASE ( "st_sod1","st_sod1_mp2" )
          DO icg = 1,pGrid%nCellsTot           
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.5_RFREAL ) THEN 
              d = 1.0_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            ELSE 
              d = 0.125_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+4_RFREAL
            END IF ! x        

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg                       
  
! ----- Sod case 2 -------------------------------------------------------------  
  
        CASE ( "st_sod2","st_sod2_mp2" )
          DO icg = 1,pGrid%nCellsTot           
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.0_RFREAL ) THEN 
              d = 1.0_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            ELSE 
              d = 0.01_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+3_RFREAL
            END IF ! x        

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg     
  
! ------------------------------------------------------------------------------
!       Supersonic vortex flow. NOTE this problem is two-dimensional and assumed 
!       to lie in the x-y plane. 
! ------------------------------------------------------------------------------
  
        CASE ( "ssvorth20x5l1"   ,"ssvortp20x5l1"   , &
               "ssvorth20x5l3"   ,"ssvortp20x5l3"   , &
               "ssvorth40x10l1"  ,"ssvortp40x10l1"  , & 
               "ssvorth40x10l3"  ,"ssvortp40x10l3"  , & 
               "ssvorth80x20l1"  ,"ssvortp80x20l1"  , & 
               "ssvorth80x20l3"  ,"ssvortp80x20l3"  , & 
               "ssvorth160x40l1" ,"ssvortp160x40l1" , &
               "ssvorth160x40l3" ,"ssvortp160x40l3" , &
               "ssvorth320x80l1" ,"ssvortp320x80l1" , &
               "ssvorth320x80l3" ,"ssvortp320x80l3" , &
               "ssvorth640x160l1","ssvortp640x160l1", &
               "ssvorth640x160l3","ssvortp640x160l3" )
          CALL RFLU_GetParamsHardCodeSsVortex(ri,Mi,pTot,tTot)     

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            CALL RFLU_ComputeExactFlowSsVortex(x,y,gRef,gcRef,ri,Mi,pTot, & 
                                               tTot,d,u,v,w,p)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,gRef,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Taylor vortex problem 
! ------------------------------------------------------------------------------

        CASE ( "taylorvortex" ) 
          refL  = global%refLength
          refNu = global%refVisc/global%refDensity
          refU  = global%refVelocity   
          refD  = global%refDensity
          refP  = global%refPressure

          pi = global%pi
          d  = refD
          L  = pMixtInput%prepRealVal1 

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              g  = global%refGamma
            ELSE           
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)

              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)
            END IF ! solverType

            IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
              t = global%currentTime
              CALL RFLU_ComputeExactFlowTaylorVortex(t,pi,x,y,L,refL,refNu, &
                                                     refU,refD,refP,u,v,w,p)
              pCv(CV_MIXT_XVEL,icg) = u
              pCv(CV_MIXT_YVEL,icg) = v
              pCv(CV_MIXT_ZVEL,icg) = w

              t = global%currentTime + global%dtImposed/2.0_RFREAL
              CALL RFLU_ComputeExactFlowTaylorVortex(t,pi,x,y,L,refL,refNu, &
                                                     refU,refD,refP,u,v,w,p)
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_PRES,icg) = p

              t = global%currentTime - global%dtImposed
              CALL RFLU_ComputeExactFlowTaylorVortex(t,pi,x,y,L,refL,refNu, &
                                                     refU,refD,refP,u,v,w,p)
              pCvOld(CV_MIXT_DENS,icg) = d
              pCvOld(CV_MIXT_PRES,icg) = p
            ELSE
              t = global%currentTime
              CALL RFLU_ComputeExactFlowTaylorVortex(t,pi,x,y,L,refL,refNu, &
                                                     refU,refD,refP,u,v,w,p)
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
            END IF ! solverType
          END DO ! icg

! ------------------------------------------------------------------------------
!       Multiphase Riemann problem: Two rarefactions (Toro case 1) 
! ------------------------------------------------------------------------------

        CASE ( "Two_Rarefaction" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,4313, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,4320, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.5_RFREAL ) THEN
              rl = 1000_RFREAL
              pl = 0.1_RFREAL
              rg = 1.2342_RFREAL
              pg = 0.9_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = -200.0_RFREAL             
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 273.0_RFREAL
            ELSE
              rl = 1000_RFREAL
              pl = 0.1_RFREAL
              rg = 1.2342_RFREAL
              pg = 0.9_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 200.0_RFREAL             
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 273.0_RFREAL
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCV(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Simple volcano model
! ------------------------------------------------------------------------------

        CASE ( "volcmod2dv3" )
          DO icg = 1,pGrid%nCellsTot
            y = pGrid%cofg(YCOORD,icg)

            IF ( y < -500.0_RFREAL ) THEN
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 12.5E+6_RFREAL
              t = 600.0_RFREAL
            ELSE
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 0.87E+5_RFREAL
              t = 288.15_RFREAL
            END IF ! y

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            d = MixtPerf_D_PRT(p,gc,t)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Vortex test
! ------------------------------------------------------------------------------

        CASE ( "vort" ) 

          Rc_l = 0.15_RFREAL
          A_cl = -0.0005_RFREAL
          uo_c = 1.1_RFREAL
!          uo_c = 0.3_RFREAL

          Rc = Rc_l*pMixtInput%prepRealVal1
          c  = pMixtInput%prepRealVal2
          A  = A_cl*c*pMixtInput%prepRealVal1
          uo = uo_c*c

          ro = pMixtInput%prepRealVal3
          po = pMixtInput%prepRealVal4

          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            radius = SQRT(x*x+y*y)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            IF ( radius <= Rc ) THEN
              psi = A*EXP(-(x*x+y*y)/(2.0_RFREAL*Rc*Rc))                           
            ELSE
              psi = 0.0_RFREAL
            END IF ! radius   

            d = ro 
            u = uo - y*psi/(d*Rc*Rc)
            v = x*psi/(d*Rc*Rc)
            w = 0.0_RFREAL
            p = po + (d*A*psi)/(Rc*Rc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Solid-body vortex  
! ------------------------------------------------------------------------------

        CASE ( "vortex_part" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)
            y = pGrid%cofg(YCOORD,icg)

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            d = 1.0_RFREAL
            u = y
            v = -x
            w = 0.0_RFREAL
            p = 1.0E+5_RFREAL  

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Woodward-Colella ramp. NOTE the initial condition is one-dimensional 
!       and assumed to lie along the x-axis.
! ------------------------------------------------------------------------------

        CASE ( "wcramp","wcrampsc","wcrampsm" )        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < -0.1_RFREAL ) THEN 
              d = 8.0_RFREAL
              u = 2608.87906904_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 116.5E+5_RFREAL
            ELSE 
              d = 1.4_RFREAL
              u = 0.0_RFREAL
              v = 0.0_RFREAL
              w = 0.0_RFREAL
              p = 1.0E+5_RFREAL
            END IF ! x        

            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)
        
            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)                   
          END DO ! icg
          
! ------------------------------------------------------------------------------
!       Gas-liquid mixture: Sod shocktube
! ------------------------------------------------------------------------------

        CASE ( "2DShock001" )
          IF ( pRegion%mixtInput%gasModel /= GAS_MODEL_MIXT_GASLIQ ) THEN 
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,4530, & 
                           'Case initialization only valid with gas-liq model.')
          END IF ! pRegion%mixtInput%gasModel  
          
          IF ( pRegion%specInput%nSpecies /= 2 ) THEN 
            WRITE(errorString,'(A,1X,I2)') 'Should be:', &
                                           pRegion%specInput%nSpecies
            CALL ErrorStop(global,ERR_SPEC_NSPEC_INVALID,4537, &
                           TRIM(errorString))
          END IF ! pRegion%specInput%nSpecies        
        
          rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht, &
                                 rGas)
          rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht, &
                                 rVap)
          cvl  = global%refCvLiq           
        
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < 0.5_RFREAL ) THEN
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 1.0_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 348.43_RFREAL          
            ELSE
              rl = 0.0_RFREAL
              pl = 0.0_RFREAL
              rg = 0.125_RFREAL
              pg = 1.0_RFREAL
              rv = 0.0_RFREAL
              pv = 0.0_RFREAL

              d  = rl*pl + rg*pg + rv*pv
              u  = 0.0_RFREAL
              v  = 0.0_RFREAL
              w  = 0.0_RFREAL
              t  = 278.7_RFREAL 
            END IF ! x

            Cvm = (rl*pl*cvl + rg*pg*cvg + rv*pv*cvv)/d
            Vm2 = u*u+v*v+w*w

            pCv(CV_MIXT_DENS,icg) = d
            pCv(CV_MIXT_XMOM,icg) = d*u
            pCv(CV_MIXT_YMOM,icg) = d*v
            pCv(CV_MIXT_ZMOM,icg) = d*w
            pCv(CV_MIXT_ENER,icg) = d*MixtGasLiq_Eo_CvmTVm2(Cvm,t,Vm2)
          END DO ! icg

! ------------------------------------------------------------------------------
!       Default - must be due to input error
! ------------------------------------------------------------------------------

        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,4595)  
      END SELECT ! global%casename

! ==============================================================================
!   Default
! ==============================================================================  
    
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,4603) 
  END SELECT ! pMixtInput%fluidModel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Initializing flow field from hard code done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitFlowHardCode

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitFlowHardCode.F90,v $
! Revision 1.7  2016/02/05 20:01:55  fred
! Fixing Incorrectly indexed array
!
! Revision 1.6  2016/02/04 22:13:12  fred
! Correcting input value variable name
!
! Revision 1.5  2016/02/04 22:07:31  fred
! Adding JWL EOS iterative capabilities for cylindrical detonation problem
!
! Revision 1.4  2015/07/27 04:45:42  brollin
! 1) Corrected bug in RFLUCONV where global%gridFormat was used instead of global%gridSrcFormat
! 2) Implemented new subroutine for shock tube problems (Shktb)
!
! Revision 1.3  2015/03/16 12:23:55  neal
! Fixed CTINT initialization that was wrong
!
! Revision 1.2  2015/03/11 02:42:50  neal
! Added contact interface select case option (ctint)
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.7  2009/09/28 14:22:03  mparmar
! Added cylpotential, rayleighproblem, and sphpotential
!
! Revision 1.6  2009/07/08 19:12:24  mparmar
! Added radialpulse case
!
! Revision 1.5  2008/12/06 08:43:55  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:07  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/12/03 18:32:55  mparmar
! Changed initialization of gaussian pulse to get sinusoidal pulse
!
! Revision 1.2  2007/11/28 23:05:39  mparmar
! Added acoustic,channelacoust,gaussianpulse & taylorvortex
!
! Revision 1.1  2007/04/09 18:55:41  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.36  2007/04/05 01:13:16  haselbac
! Added stg1d, modified code to allow 2nd interface
!
! Revision 1.35  2007/02/27 13:15:00  haselbac
! Added init for stg1d
!
! Revision 1.34  2007/02/16 20:01:05  haselbac
! Added init for somm_spi case
!
! Revision 1.33  2006/08/19 19:44:09  haselbac
! Significant clean-up and cosmetic changes
!
! Revision 1.32  2006/08/19 15:41:11  mparmar
! Added initializations for vortexNSCBC, farf, nscbc[1-8], vort
!
! Revision 1.31  2006/05/06 16:50:56  haselbac
! Added MP versions of kjet2 cases
!
! Revision 1.30  2006/04/20 20:57:35  haselbac
! Added kjet2v5 case
!
! Revision 1.29  2006/04/19 19:26:23  haselbac
! Fixed order of CASE statements
!
! Revision 1.28  2006/04/13 18:10:12  haselbac
! Added treatment of Culick flow
!
! Revision 1.27  2006/03/30 20:53:04  haselbac
! Changed ShockBubble hard-code, cosmetics
!
! Revision 1.26  2006/03/28 00:36:36  haselbac
! Added kjet2v4 case
!
! Revision 1.25  2006/03/26 20:22:21  haselbac
! Added cases for GL model
!
! Revision 1.24  2006/03/08 23:39:19  haselbac
! Added kjet2v3 case
!
! Revision 1.23  2006/01/06 22:16:49  haselbac
! Added routines to init linear and trig fn for grad testing
!
! Revision 1.22  2005/12/07 17:58:19  fnajjar
! Bug fix for gas constant defs in vortex_part
!
! Revision 1.21  2005/12/07 16:58:37  haselbac
! Added init for vortex_part
!
! Revision 1.20  2005/11/11 16:58:08  haselbac
! Added init for generic 2d shocktube
!
! Revision 1.19  2005/11/10 16:58:33  haselbac
! Bug fix for gmpjet
!
! Revision 1.18  2005/11/10 02:42:20  haselbac
! Added support for variable props, new cases
!
! Revision 1.17  2005/10/09 15:12:17  haselbac
! Added 2d C0 case
!
! Revision 1.16  2005/09/28 22:54:26  mparmar
! Changed cylds init to use prepRealValx
!
! Revision 1.15  2005/09/13 21:38:15  haselbac
! Changed sphds init to use prepRealValx
!
! Revision 1.14  2005/08/25 16:22:36  haselbac
! Bug fix in ds v6 init
!
! Revision 1.13  2005/08/25 03:42:48  haselbac
! Added v6 to ds cases
!
! Revision 1.12  2005/08/24 01:37:04  haselbac
! Added v5 for ds cases
!
! Revision 1.11  2005/08/21 16:03:06  haselbac
! Added v4 to ds cases
!
! Revision 1.10  2005/08/17 23:00:00  haselbac
! Added v3 versions of diffracting shock cases
!
! Revision 1.9  2005/08/08 00:12:46  haselbac
! Modified ds case initialization, now use prepRealValx
!
! Revision 1.8  2005/08/06 00:56:42  haselbac
! Added v2 to ds case names
!
! Revision 1.7  2005/07/05 19:47:58  mparmar
! Added init for diffracting shock over sphere
!
! Revision 1.6  2005/07/04 14:58:37  haselbac
! Added init for mvsh case
!
! Revision 1.5  2005/07/01 16:21:36  haselbac
! Added init for kjet2
!
! Revision 1.4  2005/06/14 00:57:45  haselbac
! Changed init for Skews case
!
! Revision 1.3  2005/04/29 12:51:24  haselbac
! Adapted to changes in interface of RFLU_ComputeExactFlowPAcoust
!
! Revision 1.2  2005/04/20 14:45:37  haselbac
! Extended init of pipeacoust case
!
! Revision 1.1  2005/04/15 15:08:16  haselbac
! Initial revision
!
! ******************************************************************************

