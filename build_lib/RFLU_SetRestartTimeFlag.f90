










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
! Purpose: Set restart time and flag.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!
! Output: None.
!
! Notes: 
!   1. Cannot be set in routine which reads restart information file because 
!      in coupled runs that file is not called.
!
! ******************************************************************************
!
! $Id: RFLU_SetRestartTimeFlag.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetRestartTimeFlag(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModMPI
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetRestartTimeFlag.F90,v $ $Revision: 1.1.1.1 $'

  CALL RegisterFunction(global,'RFLU_SetRestartTimeFlag',"../libflu/RFLU_SetRestartTimeFlag.F90")

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Setting restart time and flag...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set restart time and flag
! ******************************************************************************

  IF ( global%flowType == FLOW_STEADY ) THEN ! steady flow
    global%restartIter = global%currentIter

    IF ( global%restartIter == 0 ) THEN 
      global%restartFromScratch = .TRUE.
    ELSE 
      global%restartFromScratch = .FALSE.
    END IF ! global%restartIter
  ELSE ! unsteady flow
    global%restartTime = global%currentTime

    IF ( global%restartTime == 0.0_RFREAL ) THEN 
      global%restartFromScratch = .TRUE.
    ELSE 
      global%restartFromScratch = .FALSE.
    END IF ! global%restartIter     
  END IF ! global%flowType

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Setting restart time and flag done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetRestartTimeFlag

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetRestartTimeFlag.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:36  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.1  2004/06/25 20:04:01  haselbac
! Initial revision
!
! ******************************************************************************

