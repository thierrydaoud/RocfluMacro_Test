










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
!******************************************************************************
!
! Purpose: Scale grid speeds.
!
! Description: None. 
!
! Input: 
!   region     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Grid speeds need to be scaled because they are derived from an 
!      expression equivalent to a forward Euler discretization, but the actual
!      solution is evolved using the Runge-Kutta scheme. As a result, the 
!      right-hand side differs by a factor involving the Runge-Kutta 
!      coefficients. Its the task of this routine to scale the grid speeds to
!      make up for this discrepancy. 
!   2. Note that at the last Runge-Kutta stage, the coefficient is negative!  
!
!******************************************************************************
!
! $Id: ScaleGridSpeeds.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ScaleGridSpeeds( region )

  USE ModDataTypes
  USE ModGrid,       ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY  : t_patch
  USE ModGlobal, ONLY    : t_global
  USE ModError
  USE ModParameters

  USE RFLU_ModGrid

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  INTEGER      :: ifc, iPatch, irk, irkStep
  REAL(RFREAL) :: scaleFactor, term
  REAL(RFREAL) :: ark(5), grk(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: gs
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  global => region%global

  CALL RegisterFunction(global,'ScaleGridSpeeds',"../libfloflu/ScaleGridSpeeds.F90")
  
! *****************************************************************************
! Set pointers and variables
! *****************************************************************************

  irkStep = region%irkStep

  ark(:) = region%mixtInput%ark(:)
  grk(:) = region%mixtInput%grk(:)

! *****************************************************************************
! Scale grid speeds
! *****************************************************************************

! =============================================================================
! Determine scaling factor
! ============================================================================= 

  IF ( irkStep > 1 .AND. irkStep < global%nrkSteps ) THEN 
    scaleFactor = ark(irkStep-1)/ark(irkStep) 
  ELSE IF ( irkStep == 1 ) THEN
    scaleFactor = 1.0_RFREAL/ark(irkStep)
  ELSE IF ( irkStep == global%nrkSteps ) THEN
    term = 0.0_RFREAL

    DO irk = 1,global%nrkSteps-1
      term = term + grk(irk)/ark(irk)
    END DO ! irk

    scaleFactor = (1.0_RFREAL/ark(global%nrkSteps) - term)/ark(irkStep-1)    
  ELSE ! Defensive programming
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,154)
  END IF ! irkStep

! =============================================================================
! Interior faces
! ============================================================================= 

  gs => region%grid%gs
 
  DO ifc = 1,region%grid%nFaces 
    gs(ifc) = scaleFactor*gs(ifc)
  END DO ! ifc

! =============================================================================
! Patch faces
! ============================================================================= 

  DO iPatch = 1,region%grid%nPatches
    pPatch => region%patches(iPatch)
    gs     => pPatch%gs

    DO ifc = 1,pPatch%nBFaces
      gs(ifc) = scaleFactor*gs(ifc)
    END DO ! ifc
  END DO ! iPatch

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ScaleGridSpeeds

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ScaleGridSpeeds.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:51:17  haselbac
! Initial revision after changing case
!
! Revision 1.8  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.4  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.3  2003/03/15 16:35:29  haselbac
! Changed loop limit
!
! Revision 1.2  2002/12/02 20:10:41  jblazek
! Moved RFLU_ModGrid inside ifdef 1 in ScaleGridSpeeds.
!
! Revision 1.1  2002/11/08 21:55:48  haselbac
! Initial revision
!
!******************************************************************************

