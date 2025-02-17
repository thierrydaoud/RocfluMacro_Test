










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
! Purpose: Check for posivity of variables
!
! Description: None.
!
! Input: 
!   pRegion     Region data
!
! Output: None.
!
! Notes: 
!   1. Compute and check pressure here because, it being computed later 
!      through a call to MixtureProperties, a negative value would not be
!      detected and could cause negative pressure and temperature.
!
! ******************************************************************************
!
! $Id: RFLU_CheckPositivity.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CheckPositivity(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModParameters
  USE ModMPI
  
  USE ModInterfaces, ONLY: MixtPerf_G_CpR, &
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_R_M, &
                           MixtPerf_T_DPR, &
                           RFLU_PrintLocInfo
  USE RFLU_ModJWL ! Subbu 
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

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER, PARAMETER :: MAX_NEGATIVE_LOCS = 10
  INTEGER :: icg,indCp,indMol,indVFracE,nLocs
  INTEGER :: loc(MAX_NEGATIVE_LOCS,MIN_VAL:MAX_VAL)
  REAL(RFREAL) :: Eo,gamma,p,rgas,rho,rrho,t,u,v,Vm2,w,VFracE
  REAL(RFREAL) :: x,y,z
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CheckPositivity.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CheckPositivity',"../libflu/RFLU_CheckPositivity.F90")

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

  nLocs = 0

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol 

  indVFracE = 0


! ******************************************************************************
! Loop over cells and check for positivity
! ******************************************************************************

  cv => pRegion%mixt%cv
  dv => pRegion%mixt%dv
  gv => pRegion%mixt%gv

  DO icg = 1,pGrid%nCells
    x = pGrid%cofg(XCOORD,icg)
    y = pGrid%cofg(YCOORD,icg)
    z = pGrid%cofg(ZCOORD,icg)
    rho  = cv(CV_MIXT_DENS,icg)
    rrho = 1.0_RFREAL/rho
    u    = rrho*cv(CV_MIXT_XMOM,icg)
    v    = rrho*cv(CV_MIXT_YMOM,icg)
    w    = rrho*cv(CV_MIXT_ZMOM,icg)        
    Eo   = rrho*cv(CV_MIXT_ENER,icg)

    ! Subbu - Correction    
    rgas  = MixtPerf_R_M(gv(GV_MIXT_MOL,icg*indMol))
    gamma = MixtPerf_G_CpR(gv(GV_MIXT_CP,icg*indCp),rgas)
    Vm2   = u*u + v*v + w*w !Fred - for Channeling fix try...
   
    ! rho = rho/(1.0_RFREAL - pRegion%plag%vFracE(1,indVFracE*icg))

    ! p = MixtPerf_P_DEoGVm2(rho,Eo,gamma,Vm2)
    ! t = MixtPerf_T_DPR(rho,p,rgas)
      p  = dv(DV_MIXT_PRES,icg)
      t  = dv(DV_MIXT_TEMP,icg)
    ! Subbu - End correction

    ! Subbu - Check for negative pressure
    ! TLJ - added some extra print statements 01/17/2025
    IF ( p <= 0.0_RFREAL ) THEN
      !p  = dv(DV_MIXT_PRES,icg)
      !Eo = p/((gamma-1.0_RFREAL)*rho) + 0.5_RFREAL*Vm2
      !t  = MixtPerf_T_DPR(rho,p,rgas)
      !cv(CV_MIXT_ENER,icg) = rho*Eo

      WRITE(*,'(A,I8,A)') "  ========  Fixed -ve pressure in cell",icg,"  ========  "
      
      IF (pRegion%mixtInput%gasModel == GAS_MODEL_TCPERF) THEN
       WRITE(*,'(A,I8,F10.4,F14.4,F14.4,F10.4,A)') &
               "  ========  Ideal Gas - Before fix values:{cell,rho,e,P,T}" &
               ,icg,rho,Eo,p,t,"  ========  "
       Eo = 1.0e4_RFREAL + 0.5_RFREAL*Vm2 !Energy floor
       p = MixtPerf_P_DEoGVm2(rho,Eo,gamma,Vm2)
       t = MixtPerf_T_DPR(rho,p,rgas) !Recompute P,T
       dv(DV_MIXT_PRES,icg) = p
       dv(DV_MIXT_TEMP,icg) = t  !Reset Dep Var Array
       cv(CV_MIXT_ENER,icg) = rho*Eo !Reset Total energy array
       WRITE(*,'(A,I8,F10.4,F14.4,F14.4,F10.4,A)') &
               "  ========  Ideal Gas - Fix attempted:{cell,rho,e,P,T}" &
               ,icg,rho,Eo,p,t,"  ========  "
       WRITE(*,'(A,2x,3F10.4)') " Location {x,y,z}",x,y,z
      ELSE
       STOP
      END IF 
        !Fred - testing placing an energy floor here for Channeling
        !study IG strong shock cases...2/13/20
    END IF ! p <= 0.0_RFREAL
    ! Subbu - End check for negative pressure

    IF ( (rho <= 0.0_RFREAL) .OR. (p <= 0.0_RFREAL) ) THEN
      nLocs = nLocs + 1   

      IF ( nLocs == 1 ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
              'Negative positive-definite variables detected!'
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Module: Mixture.'        
              
        IF ( global%flowType == FLOW_UNSTEADY ) THEN 
          WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME,'Current time:', &
                                              global%currentTime              
        ELSE 
          WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME, &
                                         'Current iteration number:', &
                                         global%currentIter           
        END IF ! global%flowType                 
                                            
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal 
        WRITE(STDOUT,'(A,6X,A,6(1X,A))') SOLVER_NAME,'#', &
                                         '   Density   ', &
                                         '  x-velocity ', &
                                         '  y-velocity ', &
                                         '  z-velocity ', &
                                         '   Pressure  ', &
                                         ' Temperature '       
      END IF ! nLocs

      IF ( nLocs <= MAX_NEGATIVE_LOCS ) THEN 
        WRITE(STDOUT,'(A,4X,I3,7(1X,E13.6))') SOLVER_NAME,nLocs, & 
                                              rho,u,v,w,p,t,VFracE
        loc(nLocs,MIN_VAL:MAX_VAL) = icg                                   
      END IF ! nLocs
    END IF ! cv    
  END DO ! icg

! ******************************************************************************
! Write out message and call error handling routine
! ******************************************************************************

  IF ( nLocs > 0 ) THEN 
    IF ( nLocs > MAX_NEGATIVE_LOCS ) THEN 
       WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, & 
             'Only wrote the first',MAX_NEGATIVE_LOCS,'of',nLocs, & 
             'cells with negative positive-definite variables.'    
      CALL RFLU_PrintLocInfo(pRegion,loc,MAX_NEGATIVE_LOCS, & 
                             LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
    ELSE 
      CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, & 
                             LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
    END IF ! nLocs
    
    CALL ErrorStop(global,ERR_NEGATIVE_POSDEF,268)   
  END IF ! nLocs

! ******************************************************************************
! End
! ******************************************************************************


  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CheckPositivity

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckPositivity.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:48  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.12  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.11  2006/03/26 20:21:25  haselbac
! Removed dv declaration and definition
!
! Revision 1.10  2005/07/07 22:44:00  haselbac
! Added profiling calls, cosmetics
!
! Revision 1.9  2004/01/29 22:56:27  haselbac
! Some clean-up and cosmetic changes
!
! Revision 1.8  2004/01/11 02:06:38  jiao
! Eliminated some redundant trailing spaces that made some lines too long.
! This changed was needed to compile with NAG F90 compiler.
!
! Revision 1.7  2003/12/04 03:33:53  haselbac
! Clean-up and minor changes
!
! Revision 1.6  2003/07/22 02:12:58  haselbac
! Cosmetics only
!
! Revision 1.5  2003/06/04 22:00:37  haselbac
! Rewrite to be more informative
!
! Revision 1.4  2003/05/13 23:47:21  haselbac
! Added check for pressure, cosmetics
!
! Revision 1.3  2003/04/10 14:38:24  haselbac
! Changed interface statement
!
! Revision 1.2  2003/03/15 16:48:09  haselbac
! Cosmetics, adapted function calls
!
! Revision 1.1  2002/10/27 18:55:24  haselbac
! Initial revision
!
! ******************************************************************************

