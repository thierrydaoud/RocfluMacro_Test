










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
! Purpose: Defines the interaction Drag.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: modifies region%inrtInput%inrts
!
! Notes:
!
!   Whereas INRT_Initialize sets up everything generic to all interactions,
!   this routine sets up things specific to this interaction.
!
!   In particular, this is where the designer
!
!   (a) gives values for all the interactions parameters (or leaves them
!       as default values given in INRT_Initialize).  These parameters are
!       described in ModInteract.
!
!   (b) defines the Edges of the interaction.  Here the type of an Edge is
!       given (i.e., whether it transports mass, momentum, or energy), as well
!       as the indices of the Nodes at either end.  For each momentum Edge
!       he must call the routine INRT_FinishMomentumEdge.
!
!******************************************************************************
!
! $Id: INRT_DefineDrag.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_DefineDrag( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract
  USE ModError
  USE INRT_ModParameters

  USE INRT_ModInterfaces, ONLY : INRT_AllocateAuxillary,INRT_FinishMomentumEdge
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: nEdges,nPlag01

  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_inrt_edge),     POINTER :: edge
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_DefineDrag.F90,v $ $Revision: 1.1.1.1 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_DefineDrag',"../rocinteract/INRT_DefineDrag.F90" )

! begin -----------------------------------------------------------------------

! allocate memory for edges, switches and data

  input => region%inrtInput
  inrt  => input%inrts(INRT_TYPE_DRAG)

  inrt%name = "Drag"

  nPlag01 = 1
  IF (input%nPlag == 0) nPlag01 = 0

  nEdges = nPlag01 * INRT_DRAG_NEDGES

  CALL INRT_AllocateAuxillary(global,inrt,nEdges, &
    INRT_SWI_DRAG_TOTAL,INRT_DAT_DRAG_TOTAL)

! set parameters for this interaction:  see comment in ModInteract

! ** for this interaction, all parameters keep their default values **

! define Edges

  IF (nEdges > 0) THEN

    edge => inrt%edges(INRT_DRAG_L_XMOM_G)

    edge%tEdge    = INRT_EDGE_MOME
    edge%iNode(1) = input%indPlagJoint
    edge%iNode(2) = input%indMixt

! - the Insulate Token for this Edge is placed on the Lagrangian particle
! - end of the Edge, i.e., end number 1.

    CALL INRT_FinishMomentumEdge(global,inrt,INRT_DRAG_L_XMOM_G,1)

  END IF ! nEdges

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_DefineDrag

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_DefineDrag.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:49  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:01  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:50:11  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:14  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 21:56:22  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.3  2003/03/24 23:30:52  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************

