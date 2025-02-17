










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
! Purpose: Determine whether to write data to files.
!
! Description: None.
!
! Input:
!   global                      Pointer to global data
!
! Output: 
!   RFLU_DecideSmallWrite = .TRUE.   If should write to files
!   RFLU_DecideSmallWrite = .FALSE.  If should not write to files
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_DecideSmallWrite.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!******************************************************************************

!begin BBR

LOGICAL FUNCTION RFLU_DecideSmallWrite(global)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
      
  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_global), POINTER :: global

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  LOGICAL :: logical1,logical2

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_DecideSmallWrite.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction(global,'RFLU_DecideSmallWrite',"../libflu/RFLU_DecideSmallWrite.F90")

! *****************************************************************************
! Initialize
! *****************************************************************************

  RFLU_DecideSmallWrite = .FALSE.
  
! *****************************************************************************
! Determine whether should print to screen
! *****************************************************************************

! =============================================================================
! Unsteady flow
! =============================================================================

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    
    logical1 = ABS(global%halftimeSinceWrite-0.5_RFREAL*global%writeTime) & 
             < 0.1_RFREAL*global%dtMin
    logical2 = (global%halftimeSinceWrite > 0.5_RFREAL*global%writeTime)
 
    IF ( logical1 .OR. logical2 ) THEN
      RFLU_DecideSmallWrite = .TRUE.
    END IF ! logical1

! =============================================================================
! Steady flow
! =============================================================================

  ELSE    
    RFLU_DecideSmallWrite = (MOD(global%currentIter,global%writeIter/2) == 0)
  END IF ! global%flowType

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_DecideSmallWrite

! end BBR

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_DecideSmallWrite.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:35  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:52  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.3  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.2  2004/01/29 22:56:29  haselbac
! Removed setting of timeSinceWrite to make routine more usable
!
! Revision 1.1  2003/10/29 21:36:58  haselbac
! Initial revision
!
!******************************************************************************

