!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for viscous unsteady force with history kernel
!
! Mei-Adrian history kernel
!
! Copied from rocintereact/
!   INRT_CalcDragUnsteady_AMImplicit.F90
!   INRT_CalcDragUnsteady_AMExplicit.F90
!
! The number of time steps kept for the history
!   kernel is set in libpicl/ppiclF/source/PPICLF_STD.h
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_VU_Rocflu(i,iStage,fvux,fvuy,fvuz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag
      integer*4 i, iStage, iT
      real*8 fvux,fvuy,fvuz
      real*8 time,fH,factor,A,B,kernelVU

!
! Code:
!
      fvux = 0.0d0
      fvuy = 0.0d0
      fvuz = 0.0d0
      iT   = 1
      time = 0.0d0

      fH     = 0.75d0 + .105d0*reyL
      factor = 3.0d0*rpi*rnu*dp*fac

      if (ppiclf_nTimeBH > 1) then
         do iT = 2,ppiclf_nTimeBH-1
            time = ppiclf_timeBH(iT)

            A  = (4.0d0*rpi*time*rnu/dp**2)**(.25d0)
            B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >                 (0.5d0*dp*rnu*fH**3))**(.5d0)

            kernelVU = factor*(A+B)**(-2)

            fvux = fvux + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JX,iT,i)
     >                      -ppiclf_drudtPlag(PPICLF_JX,iT,i) )
            fvuy = fvuy + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JY,iT,i)
     >                      -ppiclf_drudtPlag(PPICLF_JY,iT,i) )
            fvuz = fvuz + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JZ,iT,i)
     >                      -ppiclf_drudtPlag(PPICLF_JZ,iT,i) )
         enddo

         iT = ppiclf_nTimeBH
         time = ppiclf_timeBH(iT)

         A  = (4.0d0*rpi*time*rnu/dp**2)**(.25d0)
         B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >                 (0.5d0*dp*rnu*fH**3))**(.5d0)

         kernelVU = 0.5d0*factor*(A+B)**(-2)

         fvux = fvux + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JX,iT,i)
     >                   -ppiclf_drudtPlag(PPICLF_JX,iT,i) )
         fvuy = fvuy + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JY,iT,i)
     >                   -ppiclf_drudtPlag(PPICLF_JY,iT,i) )
         fvuz = fvuz + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JZ,iT,i)
     >                   -ppiclf_drudtPlag(PPICLF_JZ,iT,i) )
      endif


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for viscous unsteady force with history kernel
!
! Mei-Adrian history kernel
!
! Using Hinsberg second-order method for integrating
!   the history integral
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_VU_Hinsberg(i,iStage,fvux,fvuy,fvuz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag
      integer*4 i, iStage, iT
      real*8 fvux,fvuy,fvuz
      real*8 time,fH,factor,A,B,kernelVU

!
! Code:
!
      fvux = 0.0d0
      fvuy = 0.0d0
      fvuz = 0.0d0
      iT   = 1
      time = 0.0d0

      fH     = 0.75d0 + .105d0*reyL
      factor = 3.0d0*rpi*rnu*dp*fac

      if (ppiclf_nTimeBH > 1) then
         do iT = 2,ppiclf_nTimeBH-1
            time = ppiclf_timeBH(iT)

            A  = (4.0d0*rpi*time*rnu/dp**2)**(.25d0)
            B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >                 (0.5d0*dp*rnu*fH**3))**(.5d0)

            kernelVU = factor*(A+B)**(-2)

            fvux = fvux + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JX,iT,i)
     >                      -ppiclf_drudtPlag(PPICLF_JX,iT,i) )
            fvuy = fvuy + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JY,iT,i)
     >                      -ppiclf_drudtPlag(PPICLF_JY,iT,i) )
            fvuz = fvuz + kernelVU*
     >                   ( ppiclf_drudtMixt(PPICLF_JZ,iT,i)
     >                      -ppiclf_drudtPlag(PPICLF_JZ,iT,i) )
         enddo

         iT = ppiclf_nTimeBH
         time = ppiclf_timeBH(iT)

         A  = (4.0d0*rpi*time*rnu/dp**2)**(.25d0)
         B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >                 (0.5d0*dp*rnu*fH**3))**(.5d0)

         kernelVU = 0.5d0*factor*(A+B)**(-2)

         fvux = fvux + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JX,iT,i)
     >                   -ppiclf_drudtPlag(PPICLF_JX,iT,i) )
         fvuy = fvuy + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JY,iT,i)
     >                   -ppiclf_drudtPlag(PPICLF_JY,iT,i) )
         fvuz = fvuz + kernelVU*
     >                ( ppiclf_drudtMixt(PPICLF_JZ,iT,i)
     >                   -ppiclf_drudtPlag(PPICLF_JZ,iT,i) )
      endif


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for viscous unsteady force with history kernel
!
! Mei-Adrian history kernel
!
! Using Hinsberg second-order method for integrating
!   the history integral
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_ShiftUnsteadyData
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag
      integer*4 i, iT

!
! Code:
!
      do i=1,ppiclf_npart
         do iT = ppiclf_nUnsteadyData,2,-1
            ppiclf_drudtMixt(PPICLF_JX,iT,i) = 
     >                      ppiclf_drudtMixt(PPICLF_JX,iT-1,i)
            ppiclf_drudtMixt(PPICLF_JY,iT,i) = 
     >                      ppiclf_drudtMixt(PPICLF_JY,iT-1,i)
            ppiclf_drudtMixt(PPICLF_JZ,iT,i) = 
     >                      ppiclf_drudtMixt(PPICLF_JZ,iT-1,i)

            ppiclf_drudtPlag(PPICLF_JX,iT,i) = 
     >                      ppiclf_drudtPlag(PPICLF_JX,iT-1,i)
            ppiclf_drudtPlag(PPICLF_JY,iT,i) = 
     >                      ppiclf_drudtPlag(PPICLF_JY,iT-1,i)
            ppiclf_drudtPlag(PPICLF_JZ,iT,i) = 
     >                      ppiclf_drudtPlag(PPICLF_JZ,iT-1,i)
         enddo
      enddo


      if (ppiclf_nTimeBH < ppiclf_nUnsteadyData) then
            ppiclf_nTimeBH = ppiclf_nTimeBH + 1
      endif

      do iT = ppiclf_nTimeBH,2,-1
            ppiclf_timeBH(it) = ppiclf_timeBH(iT-1) + ppiclf_dt
      enddo


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Update arrays for Viscous Unsteady Force
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_UpdatePlag(i)
!
      implicit none
!
      include "PPICLF"
!
      integer*4 i
      real*8 SDrho

!
! Code:
!
      SDrho = ppiclf_rprop(PPICLF_R_JRHSR,i)
     >         + ppiclf_y(PPICLF_JVX,i) * ppiclf_rprop(PPICLF_R_JPGCX,i)
     >         + ppiclf_y(PPICLF_JVY,i) * ppiclf_rprop(PPICLF_R_JPGCY,i)
     >         + ppiclf_y(PPICLF_JVZ,i) * ppiclf_rprop(PPICLF_R_JPGCZ,i)

      ppiclf_drudtMixt(PPICLF_JX,1,i) =
     >            ppiclf_rprop(PPICLF_R_JSDRX,i)
      ppiclf_drudtMixt(PPICLF_JY,1,i) =
     >            ppiclf_rprop(PPICLF_R_JSDRY,i)
      ppiclf_drudtMixt(PPICLF_JZ,1,i) =
     >            ppiclf_rprop(PPICLF_R_JSDRZ,i)

      ppiclf_drudtPlag(PPICLF_JX,1,i) =
     >            ppiclf_y(PPICLF_JVX,i)*SDrho
     >            + rhof*ppiclf_ydot(PPICLF_JVX,i)
      ppiclf_drudtPlag(PPICLF_JY,1,i) =
     >            ppiclf_y(PPICLF_JVY,i)*SDrho
     >            + rhof*ppiclf_ydot(PPICLF_JVY,i)
      ppiclf_drudtPlag(PPICLF_JZ,1,i) =
     >            ppiclf_y(PPICLF_JVZ,i)*SDrho
     >            + rhof*ppiclf_ydot(PPICLF_JVZ,i)


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Reset arrays for Viscous Unsteady Force
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_prop2plag(pp)
!
      implicit none
!
      include "PPICLF"
!
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag
      integer*4 i,k,ic,iT
      integer*4 pp
!
! Code:
!
      do i=1,ppiclf_npart
         k = 0
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_drudtMixt(ic,iT,i) = ppiclf_rprop3(k,i)
         enddo
         enddo
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_drudtPlag(ic,iT,i) = ppiclf_rprop3(k,i)
         enddo
         enddo
      enddo


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Reset arrays for Viscous Unsteady Force
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_plag2prop(pp)
!
      implicit none
!
      include "PPICLF"
!
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag
      integer*4 i,k,ic,iT
      integer*4 pp
!
! Code:
!
      do i=1,ppiclf_npart
         k = 0
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_drudtMixt(ic,iT,i)
         enddo
         enddo
         do ic = 1,3
         do iT = 1, ppiclf_nUnsteadyData
            k = k+1
            ppiclf_rprop3(k,i) = ppiclf_drudtPlag(ic,iT,i)
         enddo
         enddo
      enddo


      return
      end
