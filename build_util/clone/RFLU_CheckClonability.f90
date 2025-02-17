










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
! Purpose: Check clonability.
!
! Description: None.
!
! Input:
!   pRegion     	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_CheckClonability.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CheckClonability(pRegion)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI
  
  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================
   
  TYPE(t_region), POINTER :: pRegion
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: cntr,errorFlag,iPatch,iPatchRelated
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch,pPatchRelated

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CheckClonability.F90,v $ $Revision: 1.1.1.1 $'

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_CheckClonability',"../../utilities/clone/RFLU_CheckClonability.F90")

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking clonability...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal       
  END IF ! global%verbLevel

  pGrid => pRegion%grid  

! ******************************************************************************
! Check that have two linear periodic boundaries
! ******************************************************************************

  cntr = 0

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    IF ( pPatch%bcType == BC_PERIODIC ) THEN 
      cntr = cntr + 1
      
      iPatchRelated =  pPatch%iPatchRelated      
      pPatchRelated => pRegion%patches(iPatchRelated)
      
      IF ( pPatchRelated%iPatchRelated /= iPatch ) THEN 
        CALL ErrorStop(global,ERR_PERIODIC_INCONSISTENT,142)
      END IF ! pPatchRelated%iPatchRelated
      
      IF ( pPatch%angleRelated /= 0.0_RFREAL ) THEN 
        CALL ErrorStop(global,ERR_CLONABILITY,146)
      END IF ! pPatch%angleRelated
    END IF ! pPatch%bcType
  END DO ! iPatch

  IF ( cntr /= 2 ) THEN 
    CALL ErrorStop(global,ERR_CLONABILITY,152)
  END IF ! cntr

! ******************************************************************************
! Check that have two borders
! ******************************************************************************

  IF ( pGrid%nBorders /= 2 ) THEN 
    CALL ErrorStop(global,ERR_NBORDERS_INVALID,160)  
  END IF ! pGrid%nBorders

! ******************************************************************************
! End
! ******************************************************************************
 
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking clonability done.'
  END IF ! global%verbLevel 
 
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CheckClonability


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckClonability.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:54  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/08/07 17:13:59  haselbac
! Initial revision
!
! ******************************************************************************

