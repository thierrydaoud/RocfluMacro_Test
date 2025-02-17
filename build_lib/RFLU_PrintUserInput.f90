










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
! Purpose: Write out user input for checking purposes.
!
! Description: None.
!
! Input:
!   regions        Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_PrintUserInput.F90,v 1.2 2015/12/18 22:55:07 rahul Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PrintUserInput(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMixture, ONLY: t_mixt_input



  USE ModInterfacesSpecies, ONLY: SPEC_PrintUserInput


  USE ModInterfacesInteract, ONLY: INRT_PrintUserInput, &
                                   INRT_PrintMaterialInput

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iReg
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PrintUserInput.F90,v $ $Revision: 1.2 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_PrintUserInput',"../libflu/RFLU_PrintUserInput.F90")

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pMixtInput => regions(1)%mixtInput

  iReg = 1

! ******************************************************************************
! Echo user input
! ******************************************************************************

  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Echoing user input...'

! ==============================================================================
! Region-independent variables
! ==============================================================================

! ------------------------------------------------------------------------------
! Formats
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Formats:'
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Grid:',global%gridFormat
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Solution:',global%solutFormat
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Grid source:',global%gridSource

! ------------------------------------------------------------------------------
! Program burn
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Program burn:'
  WRITE(STDOUT,'(A,5X,A,1X,L1)') SOLVER_NAME,'PBA Flag:',global%pbaFlag

! ------------------------------------------------------------------------------
! Reference values
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')          SOLVER_NAME,'Reference values:'
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Velocity:',global%refVelocity
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Pressure:',global%refPressure
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Density:',global%refDensity
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME, &
                                    'Specific heat at constant pressure:', &
                                    global%refCp
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Ratio of specific heats:', &
                                    global%refGamma
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Length:',global%refLength
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Reynolds number:', &
                                    global%refReNum
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Viscosity:',global%refVisc
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Laminar Prandtl number:', &
                                    global%prLam
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Turbulent Prandtl number:', &
                                    global%prTurb
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Laminar Schmidt number:', &
                                    global%scnLam
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Turbulent Schmidt number:', &
                                    global%scnTurb
                                    
  IF ( pMixtInput%gasModel == GAS_MODEL_MIXT_GASLIQ ) THEN
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'BetaPLiq:',global%refBetaPLiq
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'BetaTLiq:',global%refBetaTLiq
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'CvLiq:',global%refCvLiq
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'DensityLiq:', & 
                                      global%refDensityLiq
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'PressLiq:',global%refPressLiq
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'TempLiq:',global%refTempLiq    
  END IF ! pMixtInput%gasModel                                    

! ------------------------------------------------------------------------------
! Forces
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Forces:'
  WRITE(STDOUT,'(A,5X,A,1X,L1)') SOLVER_NAME,'Flag:',global%forceFlag
  WRITE(STDOUT,'(A,5X,A,1X,L1)') SOLVER_NAME,'Patch coefficient flag:', & 
                                 global%patchCoeffFlag 
                                    
  IF ( global%forceFlag .EQV. .TRUE. ) THEN                                     
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Reference area:', &
                                      global%forceRefArea
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Reference length:', &
                                      global%forceRefLength
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Reference x-coordinate:', &
                                      global%forceRefXCoord
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Reference y-coordinate:', &
                                      global%forceRefYCoord
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Reference z-coordinate:', &
                                      global%forceRefZCoord                                 
  END IF ! global%forceFlag

! ------------------------------------------------------------------------------
! Acceleration
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Acceleration:'
  WRITE(STDOUT,'(A,5X,A,1X,L1)') SOLVER_NAME,'Flag:',global%accelOn 

  IF ( global%accelOn .EQV. .TRUE. ) THEN 
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'x-component:',global%accelX
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'y-component:',global%accelY
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'z-component:',global%accelZ          
  END IF ! global%accelOn
    
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'gravity:',global%gravity

! ------------------------------------------------------------------------------
! Absorbing boundary Condition ABC
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'ABC:'  
  WRITE(STDOUT,'(A,5X,A,1X,L1)') SOLVER_NAME,'ABC Flag:',global%abcFlag

  IF ( global%abcFlag .EQV. .TRUE. ) THEN
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Ref Density:',   &
                                      global%abcDens
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Ref VelocityX:', &
                                      global%abcVelX
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Ref VelocityY:', &
                                      global%abcVelY
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Ref VelocityZ:', &
                                      global%abcVelZ
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Ref Pressure:',  &
                                      global%abcPress

    WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'ABC Kind:',global%abcKind

    IF (global%abcKind == 0 ) THEN
      WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Sponge:',global%abcSponge
      WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'ABC Distrib:', &
                                     global%abcDistrib
      WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'ABC Order:',global%abcOrder

      SELECT CASE( global%abcOrder )
        CASE (0)
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Coeff0:', &
                                            global%abcCoeff0
        CASE (1)
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Coeff0:', &
                                            global%abcCoeff0
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Coeff1:', &
                                            global%abcCoeff1
        CASE (2)
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Coeff0:', &
                                            global%abcCoeff0
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Coeff1:', &
                                            global%abcCoeff1
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Coeff2:', &
                                            global%abcCoeff2
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,286)
      END SELECT ! pPatch%abcOrder

      WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'ABC Domain:',global%abcDomain

      SELECT CASE( global%abcDomain )
        CASE (0)
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Domain Left:', &
                                            global%abcLeft
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Domain Right:', &
                                            global%abcRight
        CASE (1)
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Domain Left:', &
                                            global%abcLeft
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Domain Right:', &
                                            global%abcRight
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Domain Bottom:', &
                                            global%abcBottom
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Domain Top:', &
                                            global%abcTop
        CASE (2)
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Domain Radius:', &
                                            global%abcRadius
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Domain CenterX:', &
                                            global%abcCenterX
          WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Domain CenterY:', &
                                            global%abcCenterY
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,314)
      END SELECT ! pPatch%abcDomain

    END IF ! global%abcKind
  END IF ! global%abcFlag

! ------------------------------------------------------------------------------
! Ghost Fluid Method GFM
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'GFM:'  
  WRITE(STDOUT,'(A,5X,A,1X,L1)') SOLVER_NAME,'GFM Flag:',global%gfmFlag

  IF ( global%gfmFlag .EQV. .TRUE. ) THEN
    WRITE(STDOUT,'(A,5X,A,1X,I5)') SOLVER_NAME,'No. of Rigids:',   &
                                      global%gfmNRigids
    WRITE(STDOUT,'(A,5X,A,1X,I5)') SOLVER_NAME,'No. of Fluids:',   &
                                      global%gfmNFluids
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Left x-coordinates:  ',   &
                                      global%gfmGridX1
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Right x-coordinates: ',   &
                                      global%gfmGridX2
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Bottom y-coordinates:',   &
                                      global%gfmGridY1
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Top y-coordinates:   ',   &
                                      global%gfmGridY2
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Back z-coordinates:  ',   &
                                      global%gfmGridZ1
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Front z-coordinates: ',   &
                                      global%gfmGridZ2
    WRITE(STDOUT,'(A,5X,A,1X,I10)') SOLVER_NAME,'No. of gridpoints in x-direction:',   &
                                      global%gfmGridNx
    WRITE(STDOUT,'(A,5X,A,1X,I10)') SOLVER_NAME,'No. of gridpoints in y-direction:',   &
                                      global%gfmGridNy
    WRITE(STDOUT,'(A,5X,A,1X,I10)') SOLVER_NAME,'No. of gridpoints in z-direction:',   &
                                      global%gfmGridNz
  END IF ! global%gfmFlag

! -----------------------------------------------------------------------------
! Moving Frame
! -----------------------------------------------------------------------------
  
  WRITE(STDOUT,'(A,3X,A)')          SOLVER_NAME,'Moving Frame:'  
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Flag:', &
                                    global%mvFrameFlag
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Acc flag:', &
                                    global%mvfAccFlag
  WRITE(STDOUT,'(A,5X,A,1X,I3)')    SOLVER_NAME,'Acc flag type:', &
                                    global%mvfAccType
  WRITE(STDOUT,'(A,5X,A,1X,I3)')    SOLVER_NAME,'Global patch index:', &
                                    global%iPatchGlobalMvFrame
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Time to start applying acc:', &
                                    global%mvfAccTs
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Time to end applying acc:', &
                                    global%mvfAccTe
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Forced particle acc x:', &
                                    global%mvfAccX
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Forced particle acc y:', &
                                    global%mvfAccY
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Forced particle acc z:', &
                                    global%mvfAccZ
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Forced particle acc, frequency:', &
                                    global%mvfOmega
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Forced particle acc, phase:', &
                                    global%mvfPhase
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Particle mass:', &
                                    global%mvfMass
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Initial location x:', &
                                    global%mvfLocInitX
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Initial location y:', &
                                    global%mvfLocInitY
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Initial location z:', &
                                    global%mvfLocInitZ
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Initial velocity x:', &
                                    global%mvfVelInitX
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Initial velocity y:', &
                                    global%mvfVelInitY
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Initial velocity z:', &
                                    global%mvfVelInitZ

! ------------------------------------------------------------------------------
! Time stepping
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Time stepping:'
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Solver type:',global%solverType
  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Hypre Solver:', &
                                               global%hypreSolver
  END IF ! global%solverType

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    WRITE(STDOUT,'(A,5X,A)')          SOLVER_NAME,'Flow type: Unsteady'
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Timestep:',global%dtImposed
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Current time:', &
                                      global%timeStamp
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Maximum time:', &
                                      global%maxTime
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Run time:', &
                                      global%runTime
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Writing interval:', &
                                      global%writeTime
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Printing interval:', &
                                      global%printTime
    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
      WRITE(STDOUT,'(A,5X,A,1X,I8)')    SOLVER_NAME,'Max PC iteration:', &
                                        global%maxIterPC
      WRITE(STDOUT,'(A,5X,A,1X,I8)')    SOLVER_NAME,'Max Hypre iteration:', &
                                        global%maxIterHypre
      WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Tolerance PC:', &
                                        global%tolPC
      WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Tolerance Hypre:', &
                                        global%tolHypre
    ELSE
      WRITE(STDOUT,'(A,5X,A,1X,I2)')    SOLVER_NAME,'Runge-Kutta scheme:', &
                                        global%rkScheme           
    END IF ! global%solverType
  ELSE
    WRITE(STDOUT,'(A,5X,A)')          SOLVER_NAME,'Flow type: Steady'
    WRITE(STDOUT,'(A,5X,A,1X,I8)')    SOLVER_NAME,'Current iteration:', &
                                      global%currentIter
    WRITE(STDOUT,'(A,5X,A,1X,I8)')    SOLVER_NAME,'Maximum iteration:', &
                                      global%maxIter
    WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Residual tolerance:', &
                                      global%resTol
    WRITE(STDOUT,'(A,5X,A,1X,I8)')    SOLVER_NAME,'Writing interval:', &
                                      global%writeIter
    WRITE(STDOUT,'(A,5X,A,1X,I8)')    SOLVER_NAME,'Printing interval:', &
                                      global%printIter
  END IF ! global%flowType

! ------------------------------------------------------------------------------
! Time zooming
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Time Zooming:'
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Min Plane:',global%tzMinPlane
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Max Plane:',global%tzMaxPlane
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Propellant Density:',global%tzRhos
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Burn A:',global%tzA
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Burn N:',global%tzN
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Throat Radius:',global%tzThroatRad
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Chamber Length:',global%tzLenChamb
  IF(global%tzSubNoz) THEN
     WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, 'Submerged Nozzle: Yes'
     WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Nozzle Inlet:',global%tzNozInlet
     WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Nozzle Radius:',global%tzNozRad
  ELSE
     WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, 'Submerged Nozzle: No'
  ENDIF
  IF(global%tzCoordLong == XCOORD) THEN
     WRITE(STDOUT,'(A,5X,A)')          SOLVER_NAME,'Rocket Axis: X'
  ELSE IF (global%tzCoordLong == YCOORD) THEN
     WRITE(STDOUT,'(A,5X,A)')          SOLVER_NAME,'Rocket Axis: Y'
  ELSE
     WRITE(STDOUT,'(A,5X,A)')          SOLVER_NAME,'Rocket Axis: Z'
  ENDIF

! ------------------------------------------------------------------------------
! Case Constraints
! ------------------------------------------------------------------------------

  IF(global%cnstrCaseRad > 0.0) THEN
     WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Case Constraints:'
     WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Case Radius:',global%cnstrCaseRad
     WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Head End:',global%cnstrHeadEnd
     WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Aft End:',global%cnstrAftEnd
     WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Ellipse Major:',global%cnstrEllipsL
     WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Ellipse Minor:',global%cnstrEllipsT
     WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Sub Nozzle Rad:',global%cnstrNozY
     WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Tolerance 1:',global%cnstrTol1
     WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Tolerance 2:',global%cnstrTol2
     IF(global%cnstrCoordL == XCOORD) THEN
        WRITE(STDOUT,'(A,5X,A)')          SOLVER_NAME,'Rocket Axis: X'
     ELSE IF (global%cnstrCoordL == YCOORD) THEN
        WRITE(STDOUT,'(A,5X,A)')          SOLVER_NAME,'Rocket Axis: Y'
     ELSE
        WRITE(STDOUT,'(A,5X,A)')          SOLVER_NAME,'Rocket Axis: Z'
     ENDIF
  ELSE
     WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'No case constraints:'     
  ENDIF

! ------------------------------------------------------------------------------
! Random number generator values
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Random number generator values:'
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'SeedOffset:',&
                                 global%randSeedOffset
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'SeedType:',&
                                 global%randSeedType


! ------------------------------------------------------------------------------
! Flow specific info
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
! Preprocessor
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Preprocessor:'
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Partition mode:', &
                                 global%prepPartMode

! ------------------------------------------------------------------------------
! Postprocessor
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)')          SOLVER_NAME,'Postprocessor:'
  WRITE(STDOUT,'(A,5X,A,1X,I2)')    SOLVER_NAME,'Output format:', &
                                    global%postOutputFormat
  WRITE(STDOUT,'(A,5X,A,1X,I2)')    SOLVER_NAME,'Plot type:', &
                                    global%postPlotType
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Merge flag:', &
                                    global%postMergeFlag
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Write merged files flag:', &
                                    global%postWriteMergeFlag                               
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Plot volume flag:', &
                                    global%postPlotVolFlag
  WRITE(STDOUT,'(A,5X,A,1X,I2)')    SOLVER_NAME,'Interpolation type:', &
                                    global%postInterpType
  WRITE(STDOUT,'(A,5X,A,1X,I2)')    SOLVER_NAME,'Interpolation order:', &
                                    global%postInterpOrder
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Special cell flag:', &
                                    global%postSpecFlag
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Extraction flag:', &
                                    global%postExtractFlag
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Mixture cv flag:', &
                                    global%postPlotMixtCvFlag
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Mixture dv flag:', &
                                    global%postPlotMixtDvFlag
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Mixture gv flag:', &
                                    global%postPlotMixtGvFlag
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Discontinuity flag:', &
                                    global%postDiscFlag
  WRITE(STDOUT,'(A,5X,A,1X,I3)')    SOLVER_NAME,'Schlieren type:', &
                                    global%postSchType
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Schlieren exponent:', &
                                    global%postSchExp
  WRITE(STDOUT,'(A,5X,A,1X,I3)')    SOLVER_NAME,'Number of fringes:', &
                                    global%postNFringes
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Vorticity flag:', &
                                    global%postVortFlag                                    
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Vortex core flag:', &
                                    global%postVortCoreFlag                                    
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Error flag:', &
                                    global%postCompErrFlag 
  WRITE(STDOUT,'(A,5X,A,1X,I3)')    SOLVER_NAME,'Number of servers:', &
                                    global%postNServers
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Patch plotting flag:', &
                                    global%postPlotPatchFlag
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Lagrangian to Eulerian '// &
                                    'flag:',global%postLag2EulFlag   
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Gradient flag:', &
                                    global%postGradFlag
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Particle fraction:', &
                                    global%postPlagFrac
  WRITE(STDOUT,'(A,5X,A,1X,L1)')    SOLVER_NAME,'Virtual zones flag:', &
                                    global%postZoneVirtFlag

! ------------------------------------------------------------------------------
! Misc global variables
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Misc section:'
  WRITE(STDOUT,'(A,5X,A,3X,E12.5)') SOLVER_NAME,'nBVertEstCorrFactor:', &
                                    global%nBVertEstCorrFactor

! ------------------------------------------------------------------------------
! Material input
! ------------------------------------------------------------------------------

  CALL INRT_PrintMaterialInput(global)

! ==============================================================================
! Region-dependent variables. NOTE these are not really region-dependent,
! just stored that way for commonality with Rocflo, so only write out for
! first region.
! ==============================================================================

! ------------------------------------------------------------------------------
! Flow model and grid motion
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Flow type:'

  IF ( pMixtInput%flowModel == FLOW_EULER ) THEN
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Flow model: Euler'
  ELSE
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Flow model: Navier-Stokes'
  END IF ! pMixtInput%flowModel

  IF ( pMixtInput%moveGrid .EQV. .TRUE. ) THEN
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Grid motion: Active'
  ELSE
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Grid motion: Not active'
  END IF ! pMixtInput%moveGrid

! ------------------------------------------------------------------------------
! Mixture
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Mixture:'
  WRITE(STDOUT,'(A,5X,A,1X,L1)') SOLVER_NAME,'Frozen flag:', &
                                 pMixtInput%frozenFlag
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Gas model:',pMixtInput%gasModel
  
! ------------------------------------------------------------------------------
! Viscosity model
! ------------------------------------------------------------------------------

  IF ( pMixtInput%computeTv ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Viscosity model:'

    IF ( pMixtInput%viscModel == VISC_SUTHR ) THEN
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Viscosity model: Sutherland'
      WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Reference viscosity:', &
                                        pMixtInput%refVisc
      WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Reference visc temperature:', &
                                        pMixtInput%refViscTemp
      WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Sutherland coefficient:', &
                                        pMixtInput%suthCoef
    ELSE IF ( pMixtInput%viscModel == VISC_FIXED ) THEN
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Viscosity model: Fixed'
      WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Viscosity value:', &
                                        pMixtInput%refVisc
    ELSE IF ( pMixtInput%viscModel == VISC_ANTIB ) THEN
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Viscosity model: Antibes'
      WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Reference viscosity:', &
                                        pMixtInput%refVisc
    ELSE IF ( pMixtInput%viscModel == VISC_WILKE_SUTH ) THEN
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Viscosity model: Wilke'
    END IF ! pMixtInput%viscModel
  END IF ! pMixtInput%computeTv

! ------------------------------------------------------------------------------
! Numerics
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Numerics:'

  IF ( pMixtInput%timeScheme == TST_HYB5RK ) THEN
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
          'Time discretization: Explicit multistage'
  ELSE IF ( pMixtInput%timeScheme == TST_STD4RK ) THEN
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
          'Time discretization: Explicit classical Runge-Kutta'
  END IF ! pMixtInput%timeScheme

  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'CFL number:',pMixtInput%cfl

  IF ( pMixtInput%spaceDiscr == DISCR_UPW_ROE ) THEN
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
          'Inviscid flux function: Roe'
  ELSE IF ( pMixtInput%spaceDiscr == DISCR_UPW_HLLC ) THEN
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
          'Inviscid flux function: HLLC'
  ELSE IF ( pMixtInput%spaceDiscr == DISCR_UPW_AUSMPLUS ) THEN
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
          'Inviscid flux function: AUSM+'          
  ELSE IF ( pMixtInput%spaceDiscr == DISCR_UPW_AUSMPLUSUP ) THEN
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME, &
          'Inviscid flux function: AUSM+up'          
  END IF ! pMixtInput%spaceDiscr

  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Spatial order of accuracy:', &
                                 pMixtInput%spaceOrder
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Reconstruction:', &
                                 pMixtInput%reconst
  WRITE(STDOUT,'(A,5X,A)')       SOLVER_NAME,'Constrained reconstruction:'
  WRITE(STDOUT,'(A,7X,A,1X,I2,1X,E12.5)') SOLVER_NAME,'Cells:', & 
                                          pMixtInput%cReconstCells, &
                                          pMixtInput%cReconstCellsWeight
  WRITE(STDOUT,'(A,7X,A,1X,I2,1X,E12.5)') SOLVER_NAME,'Faces:', & 
                                          pMixtInput%cReconstFaces, &
                                          pMixtInput%cReconstFacesWeight
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Entropy fix coefficient:', &
                                    pMixtInput%epsentr
  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'Dissipation factor:', &
                                    pMixtInput%dissFact                                    
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Dimensionality:', &
                                 pMixtInput%dimens 
  WRITE(STDOUT,'(A,5X,A,1X,L1)') SOLVER_NAME,'Axisymmetry flag:', &
                                 pMixtInput%axiFlag
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Cell-stencil dimensionality:', &
                                 pMixtInput%stencilDimensCells  
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Face-stencil dimensionality:', &
                                 pMixtInput%stencilDimensFaces
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME, &
                                 'Boundary-face-stencil dimensionality:', &
                                 pMixtInput%stencilDimensBFaces                 
  WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME, &
                                 'Boundary-face-stencil space Order', &
                                 pMixtInput%spaceOrderBFaces 

  WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'In-cell test tolerance:', &
                                    pMixtInput%tolerICT 

! ------------------------------------------------------------------------------
! Grid motion
! ------------------------------------------------------------------------------

  IF ( pMixtInput%moveGrid .EQV. .TRUE. ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Grid motion:'
    WRITE(STDOUT,'(A,5X,A,2X,I2)')    SOLVER_NAME,'Type:', &
                                      pMixtInput%moveGridType
                                      
    IF ( pMixtInput%moveGridType /= MOVEGRID_TYPE_GENX ) THEN                             
      WRITE(STDOUT,'(A,5X,A,1X,I2)')    SOLVER_NAME,'NIter:', &
                                        pMixtInput%moveGridNIter
      WRITE(STDOUT,'(A,5X,A,1X,E12.5)') SOLVER_NAME,'SFact:', &
                                        pMixtInput%moveGridSFact
    END IF ! pMixtInput%moveGridType                                    
  END IF ! pMixtInput%moveGrid

! ------------------------------------------------------------------------------
! Initialization
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Initialization:'
  WRITE(STDOUT,'(A,5X,A,2X,I2)')    SOLVER_NAME,'Flag:     ', &
                                    global%initFlowFlag
  WRITE(STDOUT,'(A,5X,A,2X,I12)')   SOLVER_NAME,'Integer 1:', &
                                    pMixtInput%prepIntVal1 
  WRITE(STDOUT,'(A,5X,A,2X,I12)')   SOLVER_NAME,'Integer 2:', &
                                    pMixtInput%prepIntVal2 
  WRITE(STDOUT,'(A,5X,A,2X,I12)')   SOLVER_NAME,'Integer 3:', &
                                    pMixtInput%prepIntVal3   
  WRITE(STDOUT,'(A,5X,A,2X,I12)')   SOLVER_NAME,'Integer 4:', &
                                    pMixtInput%prepIntVal4 
  WRITE(STDOUT,'(A,5X,A,2X,I12)')   SOLVER_NAME,'Integer 5:', &
                                    pMixtInput%prepIntVal5 
  WRITE(STDOUT,'(A,5X,A,2X,I12)')   SOLVER_NAME,'Integer 6:', &
                                    pMixtInput%prepIntVal6                                                                                                        
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  1:  ', &
                                    pMixtInput%prepRealVal1 
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  2:  ', &
                                    pMixtInput%prepRealVal2 
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  3:  ', &
                                    pMixtInput%prepRealVal3  
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  4:  ', &
                                    pMixtInput%prepRealVal4
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  5:  ', &
                                    pMixtInput%prepRealVal5 
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  6:  ', &
                                    pMixtInput%prepRealVal6 
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  7:  ', &
                                    pMixtInput%prepRealVal7 
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  8:  ', &
                                    pMixtInput%prepRealVal8
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  9:  ', &
                                    pMixtInput%prepRealVal9
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real 10:  ', &
                                    pMixtInput%prepRealVal10
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  11:  ', &
                                    pMixtInput%prepRealVal11
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  12:  ', &
                                    pMixtInput%prepRealVal12
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  13:  ', &
                                    pMixtInput%prepRealVal13
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  14:  ', &
                                    pMixtInput%prepRealVal14
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  15:  ', &
                                    pMixtInput%prepRealVal15
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  16:  ', &
                                    pMixtInput%prepRealVal16
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  17:  ', &
                                    pMixtInput%prepRealVal17
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  18:  ', &
                                    pMixtInput%prepRealVal18
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  19:  ', &
                                    pMixtInput%prepRealVal19
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real 20:  ', &
                                    pMixtInput%prepRealVal20
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  21:  ', &
                                    pMixtInput%prepRealVal21
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  22:  ', &
                                    pMixtInput%prepRealVal22
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  23:  ', &
                                    pMixtInput%prepRealVal23
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  24:  ', &
                                    pMixtInput%prepRealVal24
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  25:  ', &
                                    pMixtInput%prepRealVal25
  WRITE(STDOUT,'(A,5X,A,2X,E12.5)') SOLVER_NAME,'Real  26:  ', &
                                    pMixtInput%prepRealVal26


! ------------------------------------------------------------------------------
! Multi-physics modules
! ------------------------------------------------------------------------------

  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Multi-physics modules:'

! Turbulence -------------------------------------------------------------------

  IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
    IF ( pMixtInput%turbModel /= TURB_MODEL_NONE ) THEN
    ELSE
      WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Turbulence module: Not active'
    END IF ! MixtInput%turbModel
  END IF ! pMixtInput%flowModel

! Species ----------------------------------------------------------------------

  IF ( global%specUsed .EQV. .TRUE. ) THEN
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Species module: Active'
    CALL SPEC_PrintUserInput(regions(iReg))
  ELSE
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Species module: Not active'
  END IF ! global%specUsed

! Particles --------------------------------------------------------------------

  IF ( global%plagUsed .EQV. .TRUE. ) THEN
  ELSE
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Particle module: Not active'
  END IF ! global%plagUsed

! Radiation --------------------------------------------------------------------

  IF ( pMixtInput%radiUsed .EQV. .TRUE. ) THEN
  ELSE
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Radiation module: Not active'
  END IF ! pMixtInput%radiUsed

! Interactions -----------------------------------------------------------------

  IF ( global%inrtUsed .EQV. .TRUE. ) THEN
    CALL INRT_PrintUserInput( regions(iReg) )
  ENDIF


! ******************************************************************************
! End
! ******************************************************************************

  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Echoing user input done.'

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PrintUserInput

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PrintUserInput.F90,v $
! Revision 1.2  2015/12/18 22:55:07  rahul
! Added write statement for AUSM+up scheme.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.16  2010/03/14 23:26:04  mparmar
! Added hypreSolver and conditional printing for Hypre related variables
!
! Revision 1.15  2009/08/28 18:29:45  mtcampbe
! RocfluMP integration with Rocstar and some makefile tweaks.  To build
! Rocstar with new Rocflu:
! make ROCFLU=RocfluMP
! To build Rocstar with the new RocfluND:
! make ROCFLU=RocfluMP HYPRE=/the/hypre/install/path
!
! Revision 1.14  2009/07/09 20:44:00  mparmar
! Added printing of viscModel for VISC_WILKE_SUTH
!
! Revision 1.13  2009/07/08 19:11:31  mparmar
! Added printing of absorbing layer data
!
! Revision 1.12  2008/12/06 08:43:35  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:16:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2008/05/29 01:35:18  mparmar
! Added printing of new variables for moving reference frame
!
! Revision 1.9  2008/03/27 12:10:00  haselbac
! Added init for axiFlag
!
! Revision 1.8  2008/01/07 14:04:13  haselbac
! Added printing of postZoneVirtFlag
!
! Revision 1.7  2007/11/28 23:05:02  mparmar
! Added SOLV_IMPLICIT_HM parameters, nBVertEstCorrFactor and refViscTemp
!
! Revision 1.6  2007/09/04 13:16:08  haselbac
! Added printing of new plotting variables
!
! Revision 1.5  2007/08/07 17:17:05  haselbac
! Cosmetics only
!
! Revision 1.4  2007/06/18 17:45:59  mparmar
! Added printing of moving reference frame data
!
! Revision 1.3  2007/06/14 01:48:38  haselbac
! Added printing of global%runTime
!
! Revision 1.2  2007/04/12 17:54:23  haselbac
! Added printing of postPlagFrac
!
! Revision 1.1  2007/04/09 18:48:53  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:49  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.55  2006/10/20 21:26:52  mparmar
! Added printing of gravity
!
! Revision 1.54  2006/08/19 15:38:46  mparmar
! Added printing of mixtInput%spaceOrderBFaces
!
! Revision 1.53  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.52  2006/04/07 14:42:45  haselbac
! Added printing of stencilDimens params
!
! Revision 1.51  2006/03/26 20:21:39  haselbac
! Added printing of GL ref params, cosmetics
!
! Revision 1.50  2006/01/06 22:07:13  haselbac
! Added printing of postGradFlag and stencilDimens
!
! Revision 1.49  2005/12/25 15:27:56  haselbac
! Added printing for constrained reconstruction
!
! Revision 1.48  2005/12/24 21:26:26  haselbac
! Added printing of ICT tolerance
!
! Revision 1.47  2005/12/22 19:49:12  gzheng
! fixed a wrong format using L to print an integer, changed to 'I3'
!
! Revision 1.46  2005/12/10 23:28:26  haselbac
! Removed geom post variables
!
! Revision 1.45  2005/12/10 16:58:24  haselbac
! Added printing of postLag2EulFlag
!
! Revision 1.44  2005/12/01 17:09:53  fnajjar
! Added appropriate initialization, checking and printing of random seed type
!
! Revision 1.43  2005/11/17 14:37:53  haselbac
! Added printing of new prepRealVal vars
!
! Revision 1.42  2005/11/10 22:22:05  fnajjar
! ACH: Added frozenFlag
!
! Revision 1.41  2005/11/10 02:02:29  haselbac
! Added printing of acceleration components
!
! Revision 1.40  2005/10/31 19:25:48  haselbac
! Added printing of gasModel
!
! Revision 1.39  2005/10/28 19:17:29  haselbac
! Added printing of postPlotPatchFlag
!
! Revision 1.38  2005/10/27 18:56:51  haselbac
! Added printing of constr variable
!
! Revision 1.37  2005/10/05 20:03:09  haselbac
! Added printing of post output format and nservers
!
! Revision 1.36  2005/08/24 01:35:58  haselbac
! Fixed bug in writing solverType
!
! Revision 1.35  2005/08/10 00:32:08  haselbac
! Added printing of postVortCoreFlag
!
! Revision 1.34  2005/08/09 00:55:01  haselbac
! Added printing of patchCoeffFlag, postWriteMergeFlag, postCompErrFlag
!
! Revision 1.33  2005/08/03 18:19:24  hdewey2
! Added printing solverType
!
! Revision 1.32  2005/07/25 12:21:00  haselbac
! Added postVortFlag, changed format of postSchExp
!
! Revision 1.31  2005/07/14 21:58:45  haselbac
! Added output for AUSM flux function
!
! Revision 1.30  2005/07/11 19:24:13  mparmar
! Added printing of reconst option
!
! Revision 1.29  2005/07/05 19:26:58  haselbac
! Added printing of postSchType and postSchExp
!
! Revision 1.28  2005/05/01 14:19:04  haselbac
! Added printing of postDiscFlag and postNFringes
!
! Revision 1.27  2005/04/20 14:39:39  haselbac
! Added printing of additional int and read vals
!
! Revision 1.26  2005/03/22 03:33:14  haselbac
! Added printing of prep init helper variables
!
! Revision 1.25  2005/03/09 14:53:46  haselbac
! Added printing of dimensionality
!
! Revision 1.24  2004/11/17 16:28:11  haselbac
! Added printing of rkScheme
!
! Revision 1.23  2004/10/26 15:16:35  haselbac
! Added printing of postExtractFlag
!
! Revision 1.22  2004/10/19 19:25:14  haselbac
! Removed printing of surfDiverFlag
!
! Revision 1.21  2004/07/28 15:29:19  jferry
! created global variable for spec use
!
! Revision 1.20  2004/07/23 22:43:15  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.19  2004/07/21 14:54:49  haselbac
! Added printing of postInterpType
!
! Revision 1.18  2004/07/08 02:17:43  haselbac
! Added printing of dissFact
!
! Revision 1.17  2004/06/17 23:04:09  wasistho
! added PERI_PrintUserInput
!
! Revision 1.16  2004/06/16 20:00:22  haselbac
! Added force variables, clean-up, cosmetics
!
! Revision 1.15  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.14  2004/03/02 21:49:21  jferry
! Added inrtUsed flag to mixture data structure
!
! Revision 1.13  2004/02/02 22:49:40  haselbac
! Added printing of input for materials and particles
!
! Revision 1.12  2004/01/29 22:56:36  haselbac
! Cosmetic changes only
!
! Revision 1.11  2003/12/04 03:23:59  haselbac
! Cosmetic changes only
!
! Revision 1.10  2003/11/25 21:02:50  haselbac
! Added support for rocspecies, some clean-up
!
! Revision 1.9  2003/08/07 15:31:20  haselbac
! Added vars and changed some var names
!
! Revision 1.8  2003/07/22 01:57:51  haselbac
! Added writing of global%postInterpOrder
!
! Revision 1.7  2003/05/16 22:06:15  haselbac
! Fixed bug: LOGICALs should be written out using L format
!
! Revision 1.6  2003/05/05 18:40:15  haselbac
! Added printing of plotType, changed pltVolFlag to plotVolFlag
!
! Revision 1.5  2003/04/29 21:48:32  haselbac
! Added printing of surfDiverFlag
!
! Revision 1.4  2003/04/28 22:40:46  haselbac
! Added output of post and prep flags
!
! Revision 1.3  2003/04/10 23:29:59  fnajjar
! Added printouts for viscosity models
!
! Revision 1.2  2003/03/31 16:12:39  haselbac
! Cosmetics, added printing of grid-motion type
!
! Revision 1.1  2003/01/28 15:53:32  haselbac
! Initial revision, moved from rocflu
!
! Revision 1.7  2002/11/04 22:12:16  haselbac
! Bug fix: Added ifdef STATS
!
! Revision 1.6  2002/11/02 02:04:34  wasistho
! Added TURB statistics
!
! Revision 1.5  2002/09/17 22:51:23  jferry
! Removed Fast Eulerian particle type
!
! Revision 1.4  2002/09/09 15:51:56  haselbac
! global now under region and deleted debug output
!
! Revision 1.3  2002/06/17 13:34:12  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.2  2002/05/04 17:11:25  haselbac
! Cosmetic changes
!
! Revision 1.1  2002/03/26 19:25:09  haselbac
! Initial revision
!
! ******************************************************************************

