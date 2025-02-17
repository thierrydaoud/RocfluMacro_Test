










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
! Purpose: Compute integral 4 of optimal LES approach.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to current region
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_ComputeIntegral4OLES.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ComputeIntegral4OLES(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_ComputeDCUHREInfo

  USE RFLU_ModOLES

  IMPLICIT NONE

! --- parameters      

  TYPE(t_region), POINTER :: pRegion

! --- local variables

  INTEGER :: c1g,c2g,errorFlag,i,iInt,ic1l,ic2l,ifc,ifcp,iFun,iFunNZ,key,l, & 
             loopCounter,m,n,nCells,nFun,nFunNZ,restartFlag,vLoc
  REAL(RFREAL) :: normFact,normFactTerm
  
  TYPE(t_grid), POINTER :: pGrid    
  TYPE(t_global), POINTER :: global           

! ==============================================================================
!     Start
! ==============================================================================

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeIntegral4OLES',"../rocflu/RFLU_ComputeIntegral4OLES.F90")

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Computing integral 4...'        
  END IF ! global%verbLevel      

! ==============================================================================
! Set grid pointer
! ==============================================================================

  pGrid => pRegion%grid

! ==============================================================================
! Set various quantities
! ==============================================================================

! ------------------------------------------------------------------------------
! Set DCUHRE parameters
! ------------------------------------------------------------------------------

  nDim   =  8
  nFun   = 27
  nFunNZ =  7
  key    =  0 

  errAbsReq = 10.0_RFREAL*EPSILON(1.0_RFREAL)
  errRelReq = 5.0E-3_RFREAL 

! ------------------------------------------------------------------------------
! Compute DCUHRE information and allocate work array
! ------------------------------------------------------------------------------

  maxCalls = MAX_CALLS_LIMIT

  CALL RFLU_ComputeDCUHREInfo(global,nDim,nFunNZ,key,maxCalls,workArraySize)
  CALL RFLU_AllocateDCUHREArrays(global,nDim,nFunNZ,nFun)

  ALLOCATE(integral(nFun),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,150,'integral')
  END IF ! global%error

  ALLOCATE(workArray(workArraySize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,156,'workArray')
  END IF ! global%error

! ------------------------------------------------------------------------------
! Set normalization factor
! ------------------------------------------------------------------------------

  normFactTerm = 1.0_RFREAL/(pGrid%deltaOLES**8)

! ==============================================================================
! Loop over protoype faces
! ==============================================================================

  DO ifcp = 1,3    
    ifc = pGrid%fp2fOLES(ifcp)

! ------------------------------------------------------------------------------
!   Get face limits from normal face list
! ------------------------------------------------------------------------------        

    c1g = pGrid%f2c(1,ifc)

    dummy = MAXLOC(ABS(pGrid%fn(1:3,ifc)))
    nzLoc = dummy(1) ! Used in RFLU_DefineCorrelation43x (next two also) 
    nzVal = pGrid%fc(nzLoc,ifc) 
    nzSgn = NINT(pGrid%fn(nzLoc,ifc))

    n = nDim-2

    DO m = 1,3
      IF ( m /= nzLoc ) THEN
        n = n + 1

        lowLim(n) = pGrid%intLimOLES(INT_LIM_LOW,m,c1g) 
        uppLim(n) = pGrid%intLimOLES(INT_LIM_UPP,m,c1g)
      END IF ! m 
    END DO ! m

! ------------------------------------------------------------------------------
!   Set non-zero integral components (NOTE depends on nzLoc)
! ------------------------------------------------------------------------------        

    CALL RFLU_SetMapFunNZ2FunCorr43(nFunNZ)

! ------------------------------------------------------------------------------
!   Loop over cells
! ------------------------------------------------------------------------------        

    nCells = SIZE(pGrid%fsOLES,1)

    DO ic1l = 1,nCells
      c1g = pGrid%fsOLES(ic1l,ifc)

      lowLim(1) = pGrid%intLimOLES(INT_LIM_LOW,XCOORD,c1g) 
      lowLim(2) = pGrid%intLimOLES(INT_LIM_LOW,YCOORD,c1g)
      lowLim(3) = pGrid%intLimOLES(INT_LIM_LOW,ZCOORD,c1g)                

      uppLim(1) = pGrid%intLimOLES(INT_LIM_UPP,XCOORD,c1g)
      uppLim(2) = pGrid%intLimOLES(INT_LIM_UPP,YCOORD,c1g)
      uppLim(3) = pGrid%intLimOLES(INT_LIM_UPP,ZCOORD,c1g)              

      DO ic2l = 1,nCells
        c2g = pGrid%fsOLES(ic2l,ifc)

        lowLim(4) = pGrid%intLimOLES(INT_LIM_LOW,XCOORD,c2g)
        lowLim(5) = pGrid%intLimOLES(INT_LIM_LOW,YCOORD,c2g)
        lowLim(6) = pGrid%intLimOLES(INT_LIM_LOW,ZCOORD,c2g)              

        uppLim(4) = pGrid%intLimOLES(INT_LIM_UPP,XCOORD,c2g)
        uppLim(5) = pGrid%intLimOLES(INT_LIM_UPP,YCOORD,c2g)
        uppLim(6) = pGrid%intLimOLES(INT_LIM_UPP,ZCOORD,c2g)

! ----- Loop over integrals ----------------------------------------------------          

        DO iInt = 1,3    
          errorFlag   = 0
          maxCalls    = MAX_CALLS_START
          restartFlag = 0 

          integral(:)   = 0.0_RFREAL
          integralNZ(:) = 0.0_RFREAL

          DO loopCounter = 1,DCUHRE_LOOP_LIMIT 
            IF ( iInt == 1 ) THEN 
              CALL DCUHRE(nDim,nFunNZ,lowLim,uppLim,MIN_CALLS,maxCalls, & 
                          RFLU_DefineCorrelation430,errAbsReq,errRelReq, & 
                          key,workArraySize,restartFlag,integralNZ, & 
                          errAbsEst,nEval,errorFlag,workArray)        
            ELSE IF ( iInt == 2 ) THEN 
              CALL DCUHRE(nDim,nFunNZ,lowLim,uppLim,MIN_CALLS,maxCalls, & 
                          RFLU_DefineCorrelation431,errAbsReq,errRelReq, & 
                          key,workArraySize,restartFlag,integralNZ, & 
                          errAbsEst,nEval,errorFlag,workArray)
            ELSE IF ( iInt == 3 ) THEN 
              CALL DCUHRE(nDim,nFunNZ,lowLim,uppLim,MIN_CALLS,maxCalls, & 
                          RFLU_DefineCorrelation432,errAbsReq,errRelReq, & 
                          key,workArraySize,restartFlag,integralNZ, & 
                          errAbsEst,nEval,errorFlag,workArray)        
            END IF ! iInt

            IF ( errorFlag == 0 .OR. loopCounter == DCUHRE_LOOP_LIMIT ) THEN
              EXIT
            ELSE IF ( errorFlag == 1 ) THEN 
              restartFlag = 1
              maxCalls    = MAX_CALLS_FACTOR*maxCalls

              CALL RFLU_ComputeDCUHREInfo(global,nDim,nFunNZ,key,maxCalls, & 
                                          workArraySizeNew)

              IF ( workArraySizeNew > workArraySize ) THEN
                EXIT       
              END IF ! workArraySizeNew
            ELSE 
              CALL ErrorStop(global,ERR_DCUHRE_OUTPUT,269)                            
            END IF ! errorFlag
          END DO ! loopCounter

! ------- Normalize integral ---------------------------------------------------

          IF ( iInt == 1 ) THEN 
            normFact = normFactTerm
          ELSE IF ( iInt == 2 ) THEN 
            normFact = normFactTerm*CONST_KOLMOGOROV/ & 
                       (6.0_RFREAL*pGrid%rhoOLES(ifc)**(2.0_RFREAL/3.0_RFREAL))
          ELSE IF ( iInt == 3 ) THEN 
            normFact = normFactTerm*CONST_KOLMOGOROV*CONST_KOLMOGOROV/ & 
                       (36.0_RFREAL*pGrid%rhoOLES(ifc)**(4.0_RFREAL/3.0_RFREAL))
          END IF ! iInt

          DO iFunNZ = 1,nFunNZ 
            iFun = mapFunNZ2FunCorr43(iFunNZ)
            integral(iFun) = normFact*integralNZ(iFunNZ)
          END DO ! iFunNZ

! ------- Store integral in array ----------------------------------------------

          DO iFun = 1,nFun              
            CALL RFLU_MapL2IJK(iFun,l,m,i)            

            vLoc = RFLU_GetI4PosOLES(l,m,ic1l,ic2l,nCells)
 
            IF ( iInt == 1 ) THEN
              pGrid%int40OLES(i,ifcp,vLoc) = integral(iFun) 
            ELSE IF ( iInt == 2 ) THEN 
              pGrid%int41OLES(i,ifcp,vLoc) = integral(iFun)             
            ELSE 
              pGrid%int42OLES(i,ifcp,vLoc) = integral(iFun)             
            END IF ! iInt
          END DO ! iFun     
        END DO ! iInt     

      END DO ! ic2l
    END DO ! ic1l       
  END DO ! ifcp




! ==============================================================================
! Deallocate arrays
! ==============================================================================

  CALL RFLU_DeallocateDCUHREArrays(global)

  DEALLOCATE(integral,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,363,'integral')
  END IF ! global%error

  DEALLOCATE(workArray,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,369,'workArray')
  END IF ! global%error

! ==============================================================================
! End
! ============================================================================== 

  CALL DeregisterFunction(global)     

END SUBROUTINE RFLU_ComputeIntegral4OLES

!*******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeIntegral4OLES.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.3  2008/12/06 08:43:47  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:49:57  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.6  2003/05/16 02:27:45  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.5  2003/03/15 18:32:52  haselbac
! Added KIND qualifyer to NINT, added footer
!
!*******************************************************************************

