#include "PPICLF_USER.h"
#include "PPICLF_STD.h"
!-----------------------------------------------------------------------
!
! Sam - Code for solving ydot = F(y,t)
!
! Set internal flags
!   These flags allows the user to turn on or
!      off various force terms, or to select
!      between different versions of each force model
!   To turn off all particle motion, set stationary = 0
!   To turn off a force set the flag to zero
!   For two-way coupling set collisional_flag = 0
!   For four-way coupling set collisional_flag = 1
!   To turn off feedback force set feedback_flag = 0
!   To use fluctuations turn corresponding flag = 1
!
!
!   stationary = 0 if 1, do not move particles but do
!         calculate drag forces; feedback force can also
!         be turned on
!
!   qs_flag = 2  ! none = 0; Parmar = 1; Osnes = 2
!   am_flag = 1  ! Parmar = 1; Briney = 2
!   pg_flag = 1
!   collisional_flag = 1 ! two way coupled = 0 ; four-way = 1
!
!   heattransfer_flag = 1
!
!   ViscousUnsteady_flag = 0 no viscous unsteady drag
!                        = 1 history kernal for visc. unsteady drag
!
!   feedback_flag = 1
!
!   qs_fluct_flag = 1 ! None = 0 ; Lattanzi = 1 ; Osnes = 2 
!
!
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_SetYdot
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
      real*8 :: ppiclf_rcp_part
      common /RFLU_ppiclf_misc/ ppiclf_rcp_part

      real*8 fqsx, fqsy, fqsz
      real*8 fqsforce
      real*8 fqs_fluct(3)
      real*8 famx, famy, famz 
      real*8 fdpdx, fdpdy, fdpdz
      real*8 fcx, fcy, fcz
      real*8 fbx, fby, fbz 
      real*8 fvux, fvuy, fvuz

      real*8 beta,cd

      real*8 factor, One, OneThird, rcp_fluid,
     >   rmass_add

      real*8 gkern
  
!-----------------------------------------------------------------------
      ! Thierry - 06/27/204 - added mass variables declaration
      integer*4 j, l
      real*8 SDrho
!-----------------------------------------------------------------------

      integer*4 i, n, ic, k

! Needed for heat transfer
      real*8 qq, rmass_therm, temp

! Needed for angular velocity
      real*8 taux, tauy, tauz, rmass_omega
      real*8 tau

! Finite Diff Material derivative Variables
      integer*4 nstage, istage
      integer*4 icallb
      save      icallb
      data      icallb /0/
      integer*4 idebug
      save      idebug
      data      idebug /0/

! Print Data to file
      LOGICAL I_EXIST 
      Character(LEN=25) str 
      integer*4 f_dump
      save      f_dump  
      data      f_dump /1/

      logical exist_file
!
!-----------------------------------------------------------------------
!    
      ! Avery added 10/10/2024 for subbin nearest neighbor search
      
      INTEGER*4 SBin_map( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >       + 1) - 1), (ppiclf_npart+ppiclf_npart_gp))
      INTEGER*4  SBin_counter( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >       + 1) - 1))
      INTEGER*4 i_Bin(3), n_SBin(3), tot_SBin

      INTEGER*4 nsubbin_size
      INTEGER*4 nbin_total
! 
!-----------------------------------------------------------------------
!
!
! Code:
!
      icallb = icallb + 1
      nstage = 3
      istage = mod(icallb,nstage)
      if (istage .eq. 0) istage = 3  

      ! Count every iStage=1 for debug output
      if (iStage .eq. 1) idebug = idebug + 1

      ! Print dt and time every time step
      if (ppiclf_nid==0) then
      if (istage .eq. 1) then
        write(6,'(a,2x,2(1pe14.6),2x,i3)') '*** PPICLF dt, time = ',
     >      ppiclf_dt,ppiclf_time
      endif
      endif

      rpi        = 4.0*atan(1.0)
      rcp_part   = ppiclf_rcp_part
      rpr        = 0.72
      rcp_fluid  = 1006.0

      fac = ppiclf_rk3ark(iStage)*PPICLF_DT
      if (1==2) then
         if (ppiclf_nid==0) print*,'dt,fac=',
     >      istage,ppiclf_dt,fac,
     >      stationary, qs_flag, am_flag, pg_flag,
     >      collisional_flag, heattransfer_flag, feedback_flag,
     >      qs_fluct_flag, ppiclf_debug, ppiclf_nTimeBH,
     >      ppiclf_nUnsteadyData
      endif

      One = 1.d0
      OneThird = 1.0d0/3.0d0

      ! Set initial max values
      phimax    = 0.d0

      fqsx_max  = 0.d0
      fqsy_max  = 0.d0
      fqsz_max  = 0.d0
      famx_max  = 0.d0
      famy_max  = 0.d0
      famz_max  = 0.d0
      fdpdx_max = 0.d0
      fdpdy_max = 0.d0
      fdpdz_max = 0.d0
      fcx_max   = 0.d0
      fcy_max   = 0.d0
      fcz_max   = 0.d0
      fvux_max  = 0.d0
      fvuy_max  = 0.d0
      fvuz_max  = 0.d0
      qq_max    = 0.d0

      fqsx_fluct_max = 0.d0
      fqsy_fluct_max = 0.d0
      fqsz_fluct_max = 0.d0
      fqsx_total_max = 0.d0
      fqsy_total_max = 0.d0
      fqsz_total_max = 0.d0

      fqs_mag = 0.0
      fam_mag = 0.0
      fdp_mag = 0.0
      fc_mag  = 0.0

      umean_max = 0.d0
      vmean_max = 0.d0
      wmean_max = 0.d0

!
!-----------------------------------------------------------------------
!
! Reapply axi-sym collision correction
! Right now hard coding smallest radius  
!     do i=1,ppiclf_npart
!        ppiclf_rprop(PPICLF_R_JDPe,i) = 
!     > (0.00005/ppiclf_rprop(PPICLF_R_JSPT,i))
!     > * ppiclf_rprop(PPICLF_R_JDP,i)  
        !if (ppiclf_npart .gt. 0) then
        !if ((i .eq. 1) .or. (i .eq. ppiclf_npart)) then
        !  write(*,*) "i,JSPT",i,ppiclf_rprop(PPICLF_R_JSPT,i)       
        !endif
        !endif
!      end do 
!
!-----------------------------------------------------------------------
!
! Reset arrays for Viscous Unsteady Force
!
      if (ViscousUnsteady_flag>=1) then
         call ppiclf_user_prop2plag(ppiclf_nUnsteadyData)
      endif
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
! Avery added 10/10/2024 - Map particles to subbins if collisional force, 
! Briney Added Mass force, or QS fluctation force is flagged
!
      !nearest neighbor search is used for col, am_flag 2, qs_fluct
      if (sbNearest_flag == 1) then

         if ((am_flag==2).or.(collisional_flag>=1)
     >          .or.(qs_fluct_flag>=1)) then

            call ppiclf_user_subbinMap(i_Bin, n_SBin, tot_SBin 
     >                               ,SBin_counter ,SBin_map)

         endif ! Collisions, QS Fluct, or Brinery AM flags on
      
         ! Print out relevant information about subbin
         if (ppiclf_nid==0) then
         if (iStage==1) then

         nbin_total = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)
         nsubbin_size =
     >     (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     >     (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     >     (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >       + 1) 

         if (ppiclf_time .EQ. 0.0) then
         write(*,*) 'Subbin Method used!'
         write(6,*) 'SUBBIN ', 
     >     ppiclf_time,
     >     ppiclf_bins_dx(1:3),
     >     nsubbin_size,
     >     tot_SBin,n_SBin(1:3),
     >     ppiclf_npart,ppiclf_npart_gp,
     >     nsubbin_size*(ppiclf_npart+ppiclf_npart_gp),
     >     ' GB: ',nsubbin_size*
     >             (ppiclf_npart+ppiclf_npart_gp)*4/1e9

         endif ! end ppiclf_time = 0

         if (ppiclf_debug==2) write(7001,*)
     >     ppiclf_time,
     >     ppiclf_bins_dx(1:3),
     >     nsubbin_size,
     >     tot_SBin,n_SBin(1:3),
     >     ppiclf_npart,ppiclf_npart_gp,
     >     nsubbin_size*(ppiclf_npart+ppiclf_npart_gp),
     >     nsubbin_size*(ppiclf_npart+ppiclf_npart_gp)*4/1e9 
         ! last entry in GB; assuming 4 bytes for integer*4

         endif ! end iStage = 1
         endif ! end ppiclf_nid = 0

      endif ! end sbNearest_flag = 1

!
!-----------------------------------------------------------------------
!
!
! Evaluate ydot, the rhs of the equations of motion
! for the particles
!

      do i=1,ppiclf_npart

         ! Choose viscosity law
         if (rmu_flag==rmu_fixed_param) then
            ! Constant viscosity law
            rmu = rmu_ref
         elseif (rmu_flag==rmu_suth_param) then
            ! Sutherland law
            temp    = ppiclf_rprop(PPICLF_R_JT,i)
            rmu     = rmu_ref*sqrt(temp/tref)
     >                   *(1.0+suth/tref)/(1.0+suth/temp)
         else
             call ppiclf_exittr('Unknown viscosity law$', 0.0d0, 0)
         endif
         rkappa = rcp_fluid*rmu/rpr


         ! Useful values
         rmass  = ppiclf_rprop(PPICLF_R_JVOLP,i)
     >              *ppiclf_rprop(PPICLF_R_JRHOP,i)
         vx     = ppiclf_rprop(PPICLF_R_JUX,i) - ppiclf_y(PPICLF_JVX,i)
         vy     = ppiclf_rprop(PPICLF_R_JUY,i) - ppiclf_y(PPICLF_JVY,i)
         vz     = ppiclf_rprop(PPICLF_R_JUZ,i) - ppiclf_y(PPICLF_JVZ,i)
         vmag   = sqrt(vx*vx + vy*vy + vz*vz)
         rhof   = ppiclf_rprop(PPICLF_R_JRHOF,i)
         dp     = ppiclf_rprop(PPICLF_R_JDP,i)
         rep    = vmag*dp*rhof/rmu
         rphip  = ppiclf_rprop(PPICLF_R_JPHIP,i)
         rphif  = 1.0-ppiclf_rprop(PPICLF_R_JPHIP,i)
         asndf  = ppiclf_rprop(PPICLF_R_JCS,i)
         rmachp = vmag/asndf
         rhop   = ppiclf_rprop(PPICLF_R_JRHOP,i)

         ! TLJ - redefined rprop(PPICLF_R_JSPT,i) to be the particle
         !   velocity magnitude for plotting purposes - 01/03/2025
         ppiclf_rprop(PPICLF_R_JSPT,i) = sqrt(
     >       ppiclf_y(PPICLF_JVX,i)**2 +
     >       ppiclf_y(PPICLF_JVY,i)**2 +
     >       ppiclf_y(PPICLF_JVZ,i)**2)

         rep = max(0.1,rep)

         ! Redefine volume fractions
         ! Need to make sure phi_p + phi_f = 1
         rphip = ppiclf_rprop(PPICLF_R_JPHIP,i)
         rphip = min(rphip,0.62)
         rphif = 1.0-rphip

         ! TLJ: Needed for viscous unsteady force
         rhoMixt = rhof/(1.0d0-rphip)
         reyL = dp*vmag*rhoMixt/rmu
         rnu = rmu/rhoMixt

         phimax = max(phimax,abs(rphip))

         if (ppiclf_debug==2 .and. ppiclf_nid==0) then
            if (iStage==3) then
               if (i==1) then
                  write(7010,*) i,ppiclf_time,rmass,vmag,rhof,dp,
     >             rep,rphip,rphif,rmachp,rhop,rhoMixt,reyL,
     >             rmu,rnu,rkappa
               endif
               if (i==ppiclf_npart) then
                  write(7011,*) i,ppiclf_time,rmass,vmag,rhof,dp,
     >             rep,rphip,rphif,rmachp,rhop,rhoMixt,reyL,
     >             rmu,rnu,rkappa
               endif
            endif
         endif

         ! Zero out for each particle i
         famx = 0.0d0; famy = 0.0d0; famz = 0.0d0; rmass_add = 0.0d0;
         fqsx = 0.0d0; fqsy = 0.0d0; fqsz = 0.0d0; beta = 0.0d0;
         fqs_fluct = 0.0d0
         fdpdx = 0.0d0; fdpdy = 0.0d0; fdpdz = 0.0d0;
         fcx = 0.0d0; fcy = 0.0d0; fcz = 0.0d0;
         taux = 0.0d0; tauy = 0.0d0; tauz = 0.0d0;
         fvux = 0.0d0; fvuy = 0.0d0; fvuz = 0.0d0;
         qq=0.0d0
         upmean = 0.0; vpmean = 0.0; wpmean = 0.0;
         u2pmean = 0.0; v2pmean = 0.0; w2pmean = 0.0;

!
! Step 1a: New Added-Mass model of Briney
!
         ! 07/15/2024 - If am_flag = 2, then we need to call
         !   the Unary term before we make any calls to nearest
         !   neighbor
         ! 06/05/2024 - Thierry - For each particle i, initialize
         ! variables to be used in nearest neighbors to zero
         ! before looping over particle j (j neq i)
         ! Briney Added Mass flag
         if (am_flag == 2) then 
            ! zero variables initially
            nneighbors = 0
            do j=1,3
               Fam(j) = 0.0
               Wdot_neighbor_mean(j) = 0.0
            enddo

            ! 07/14/24 - Thierry - If Briney Algorithm flag and fluct_flag
            !   are ON -> evaluate added-mass unary term before evaluating
            !   neighbor-induced acceleration in EvalNearestNeighbor

            call ppiclf_user_AM_Briney_Unary(i,iStage,
     >           famx,famy,famz,rmass_add)

         endif ! end am_flag = 2

!
! Step 1b: Call NearestNeighbor if particles i and j interact
!
         if ((am_flag==2).or.(collisional_flag>=1)
     >          .or.(qs_fluct_flag>=1)) then

         if ((qs_fluct_flag>=1) .and. (vmag .gt. 1.d-8)) then
            ! Compute mean for particle i
            !    add neighbor particle j afterward
            ! Box filter is used if qs_fluct_filter_flag=0
            !   The box filter used here is a simple cube centered
            !     at particle i with half-width dist2 (see
            !     ppiclf_user_EvalNearestNeighbor.f for definition)
            !   We use a simple arithmetic mean
            !   phipmean is not used
            ! Gaussian filter is used if qs_fluct_filter_flag=1
            !   We use the value of the Gaussian times the volume
            !     of the particle to get the filtered particle volume

            if (qs_fluct_filter_flag==0) then
               ! box filter
               phipmean = ppiclf_rprop(PPICLF_R_JVOLP,i)
               upmean   = ppiclf_y(PPICLF_JVX,i)
               vpmean   = ppiclf_y(PPICLF_JVY,i)
               wpmean   = ppiclf_y(PPICLF_JVZ,i)
               u2pmean  = upmean**2
               icpmean  = 1
            else if (qs_fluct_filter_flag==1) then
               ! gaussian kernel
               ! r = 0
               gkern = sqrt(rpi*ppiclf_filter**2/
     >                (4.0d0*log(2.0d0)))**(-ppiclf_ndim)
               phipmean = gkern*ppiclf_rprop(PPICLF_R_JVOLP,i)
               upmean   = gkern*ppiclf_y(PPICLF_JVX,i)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               vpmean   = gkern*ppiclf_y(PPICLF_JVY,i)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               wpmean   = gkern*ppiclf_y(PPICLF_JVZ,i)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               u2pmean  = gkern*(ppiclf_y(PPICLF_JVX,i)**2)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               v2pmean  = gkern*(ppiclf_y(PPICLF_JVY,i)**2)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               w2pmean  = gkern*(ppiclf_y(PPICLF_JVZ,i)**2)*
     >                    ppiclf_rprop(PPICLF_R_JVOLP,i)
               icpmean = 1
            end if
         end if

         ! add neighbors
         IF ( sbNearest_flag .EQ. 1) THEN
            CALL ppiclf_solve_NearestNeighborSB(
     >           i,tot_SBin,SBin_counter,SBin_map,n_SBin,i_Bin)
         ELSE
             CALL ppiclf_solve_NearestNeighbor(i)
         END IF

         end if ! end Step 1b; nearestneighbor

!
! Step 2: Force component quasi-steady
!
         if (qs_flag==1) call ppiclf_user_QS_Parmar(i,beta,cd)
         if (qs_flag==2) call ppiclf_user_QS_Osnes (i,beta,cd)
         fqsx = beta*vx
         fqsy = beta*vy
         fqsz = beta*vz

         fqsx_max = max(fqsx_max,abs(fqsx))
         fqsy_max = max(fqsy_max,abs(fqsy))
         fqsz_max = max(fqsz_max,abs(fqsz))

         fqs_mag  = max(fqs_mag,sqrt(fqsx*fqsx+fqsy*fqsy+fqsz*fqsz))

!
! Step 3: Force fluctuation for quasi-steady force
!
         ! Note: QS fluctuations needs nearest neighbors,
         !   and is called above in Step 1b
         if (qs_fluct_flag==1) then
            call ppiclf_user_QS_fluct_Lattanzi(i,iStage,fqs_fluct)
         elseif (qs_fluct_flag==2) then
            call ppiclf_user_QS_fluct_Osnes(i,iStage,fqs_fluct)
         endif

         ! Add fluctuation part to quasi-steady force
         fqsx = fqsx + fqs_fluct(1)
         fqsy = fqsy + fqs_fluct(2)
         fqsz = fqsz + fqs_fluct(3)

         ! Store quasi-steady fluctuating force
         ppiclf_rprop(PPICLF_R_FLUCTFX,i) = fqs_fluct(1)
         ppiclf_rprop(PPICLF_R_FLUCTFY,i) = fqs_fluct(2)
         ppiclf_rprop(PPICLF_R_FLUCTFZ,i) = fqs_fluct(3)


         fqsx_fluct_max = max(fqsx_fluct_max, abs(fqs_fluct(1)))
         fqsy_fluct_max = max(fqsy_fluct_max, abs(fqs_fluct(2)))
         fqsz_fluct_max = max(fqsz_fluct_max, abs(fqs_fluct(3)))

         fqsx_total_max = max(fqsx_total_max, abs(fqsx))
         fqsy_total_max = max(fqsy_total_max, abs(fqsy))
         fqsz_total_max = max(fqsz_total_max, abs(fqsz))

         umean_max = max(umean_max, abs(upmean))
         vmean_max = max(vmean_max, abs(vpmean))
         wmean_max = max(wmean_max, abs(wpmean))

!
! Step 4: Force component added mass
!
         if (am_flag == 1) then 
            call ppiclf_user_AM_Parmar(i,iStage,
     >           famx,famy,famz,rmass_add)

!-----------------------------------------------------------------------
!Thierry - Added Mass code continues here
         
         elseif (am_flag == 2) then 

            ! Thierry - binary_model.f90 evaluates the terms
            ! in the folllowing order:
            !   (1) Unary Term
            !   (2) Evaluates Neighbor Acceleration
            !   (3) Binary Term
            ! We replicate that here by calling them in the same order
            ! Unary and Binary calculations are now under 
            !   two separate subroutines
            ! Thierry - need to make sure NearestNeighbor is called
            !    if fluct_flag = 0 (ie, no QS fluctuations)

            ! Binary subroutine only valid when number of neighbors .gt. 0
            if (nneighbors .gt. 0) then
               call ppiclf_user_AM_Briney_Binary(i,iStage,
     >              famx,famy,famz,rmass_add)
            else
            ! if particle has no neighbors, need to multiply added mass forces
            ! by volume, as this is taken care of in Binary subroutine
               famx = famx*ppiclf_rprop(PPICLF_R_JVOLP,i)
               famy = famy*ppiclf_rprop(PPICLF_R_JVOLP,i)
               famz = famz*ppiclf_rprop(PPICLF_R_JVOLP,i)

            endif

         else
           famx = 0.0
           famy = 0.0
           famz = 0.0 
           !call ppiclf_exittr('Unknown Added-Mass Law$', 0.0d0, 0)
         endif

!-----------------------------------------------------------------------

         famx_max = max(famx_max,abs(famx))
         famy_max = max(famy_max,abs(famy))
         famz_max = max(famz_max,abs(famz))
         fam_mag =  max(fam_mag,sqrt(famx*famx+famy*famy+famz*famz))

!
! Step 5: Force component pressure gradient
!
         if (pg_flag == 1) then
            fdpdx = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >               ppiclf_rprop(PPICLF_R_JDPDX,i)
            fdpdy = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >               ppiclf_rprop(PPICLF_R_JDPDY,i)
            fdpdz = -ppiclf_rprop(PPICLF_R_JVOLP,i)*
     >               ppiclf_rprop(PPICLF_R_JDPDZ,i)
         endif ! end pg_flag = 1

         fdpdx_max = max(fdpdx_max,abs(fdpdx))
         fdpdy_max = max(fdpdy_max,abs(fdpdy))
         fdpdz_max = max(fdpdz_max,abs(fdpdz))
         fdp_mag =  max(fdp_mag,sqrt(fdpdx*fdpdx+fdpdy*fdpdy
     >                  +fdpdz*fdpdz))

!
! Step 6: Force component collisional force, ie, particle-particle
!
         if (collisional_flag >= 1) then
            ! Collision force:
            !  A discrete numerical model for granular assemblies
            !  - Cundall and Strack (1979)
            !  - Geotechnique

            ! Sam - STILL NEED TO VALIDATE COLLISION FORCE
            ! Sam - Step 1b already calls nearest neighbor
            
            fcx  = ppiclf_ydotc(PPICLF_JVX,i)
            fcy  = ppiclf_ydotc(PPICLF_JVY,i)
            fcz  = ppiclf_ydotc(PPICLF_JVZ,i) 

         endif ! collisional_flag >= 1

         fcx_max = max(fcx_max, abs(fcx))
         fcy_max = max(fcy_max, abs(fcy))
         fcz_max = max(fcz_max, abs(fcz))
         fc_mag =  max(fc_mag,sqrt(fcx*fcx+fcy*fcy+fcz*fcz))

!
! Step 7: Viscous unsteady force with history kernel
!
         if (ViscousUnsteady_flag==1) then
            call ppiclf_user_VU_Rocflu(i,iStage,fvux,fvuy,fvuz)
         elseif (ViscousUnsteady_flag==2) then
            call ppiclf_user_VU_Hinsberg(i,iStage,fvux,fvuy,fvuz)
         endif

         fvux_max = max(fvux_max, abs(fvux))
         fvuy_max = max(fvuy_max, abs(fvuy))
         fvuz_max = max(fvuz_max, abs(fvuz))

!
! Step 8: Heat transfer model
!
         rmass_therm = rmass*rcp_part

         if (heattransfer_flag >= 1) then
            call ppiclf_user_HT_driver(i,qq,rmass_therm)
         endif ! heattransfer_flag >= 1

         qq_max = max(qq_max, abs(qq))

!
! Step 9: Angular velocity model
!
         rmass_omega = rmass*dp*dp/10.0d0

         if (collisional_flag >= 2) then
            taux  = ppiclf_ydotc(PPICLF_JOX,i)
            tauy  = ppiclf_ydotc(PPICLF_JOY,i)
            tauz  = ppiclf_ydotc(PPICLF_JOZ,i) 
            call ppiclf_user_Torque_driver(i,iStage,taux,tauy,tauz)
         endif ! collisional_flag >= 2

         tau = sqrt(taux*taux + tauy*tauy + tauz*tauz)
         tau_max = max(tau_max, abs(tau))

!
! Step 10: Set ydot for all PPICLF_SLN number of equations
!
         ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
         ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
         ppiclf_ydot(PPICLF_JZ, i) = ppiclf_y(PPICLF_JVZ,i)
         ppiclf_ydot(PPICLF_JVX,i) = (fqsx+famx+fdpdx+fcx+fvux)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JVY,i) = (fqsy+famy+fdpdy+fcy+fvuy)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JVZ,i) = (fqsz+famz+fdpdz+fcz+fvuz)/
     >                               (rmass+rmass_add)
         ppiclf_ydot(PPICLF_JT,i)  = qq/rmass_therm
         ppiclf_ydot(PPICLF_JOX,i) = taux/rmass_omega
         ppiclf_ydot(PPICLF_JOY,i) = tauy/rmass_omega
         ppiclf_ydot(PPICLF_JOZ,i) = tauz/rmass_omega

!
! Update and Shift data for viscous unsteady case
!
         if (ViscousUnsteady_flag>=1) then
            call ppiclf_user_UpdatePlag(i)
         endif

!
! Step 11: Feed Back force to the gas phase
!
         ! Comment: ydotc represented the collisional force in the
         !    particle eqautions above. Here, we over-write the
         !    ydotc vectors for the feedback force used in Rocflu.
         !    Note that Rocflu uses a negative of the RHS, and
         !    so ppiclf must respect this odd convention.
         !
         ! Project work done by hydrodynamic forces:
         !   Inter-phase heat transfer and energy coupling in turbulent 
         !   dispersed multiphase flows
         !   - Ling et al. (2016)
         !   - Physics of Fluids
         ! See also for more details
         !   Explosive dispersal of particles in high speed environments
         !   - Durant et al. (2022)
         !   - Journal of Applied Physics

         if (feedback_flag==0) then
            ppiclf_ydotc(PPICLF_JVX,i) = 0.0d0 
            ppiclf_ydotc(PPICLF_JVY,i) = 0.0d0 
            ppiclf_ydotc(PPICLF_JVZ,i) = 0.0d0 
            ppiclf_ydotc(PPICLF_JT,i)  = 0.0d0
         endif

         if (feedback_flag==1) then
            ! Momentum equations feedback terms
            ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_rprop(PPICLF_R_JSPL,i) *
     >         (ppiclf_ydot(PPICLF_JVX,i)*rmass - fcx)
            ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_rprop(PPICLF_R_JSPL,i) *
     >         (ppiclf_ydot(PPICLF_JVY,i)*rmass - fcy)
            ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_rprop(PPICLF_R_JSPL,i) *
     >         (ppiclf_ydot(PPICLF_JVZ,i)*rmass - fcz)

            ! Energy equation feedback term
            ppiclf_ydotc(PPICLF_JT,i) = ppiclf_rprop(PPICLF_R_JSPL,i) *
     >         ( ppiclf_ydotc(PPICLF_JVX,i)*ppiclf_y(PPICLF_JVX,i) + 
     >           ppiclf_ydotc(PPICLF_JVY,i)*ppiclf_y(PPICLF_JVY,i) + 
     >           ppiclf_ydotc(PPICLF_JVZ,i)*ppiclf_y(PPICLF_JVZ,i) +
     >           qq )
         endif 

!
! Step 12: If stationary, don't move particles. Feedback can still be on
! though.
!
         if (stationary .gt. 0) then
            if (stationary==1) then
               ppiclf_ydot(PPICLF_JX ,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JY ,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JZ, i)  = 0.0d0
               ppiclf_ydot(PPICLF_JVX,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JVY,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JVZ,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JT,i)   = 0.0d0
               ppiclf_ydot(PPICLF_JOX,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JOY,i)  = 0.0d0
               ppiclf_ydot(PPICLF_JOZ,i)  = 0.0d0
            else
               call ppiclf_exittr('Unknown stationary flag$', 0.0d0, 0)
            endif
         elseif(stationary .lt. 0) then
            call ppiclf_user_unit_tests(i,iStage,famx,famy,famz)
            ppiclf_ydot(PPICLF_JX ,i) = ppiclf_y(PPICLF_JVX,i)
            ppiclf_ydot(PPICLF_JY ,i) = ppiclf_y(PPICLF_JVY,i)
            ppiclf_ydot(PPICLF_JZ, i) = ppiclf_y(PPICLF_JVZ,i)
            ppiclf_ydot(PPICLF_JVX,i) = ppiclf_ydot(PPICLF_JVX,i)+famx
            ppiclf_ydot(PPICLF_JVY,i) = ppiclf_ydot(PPICLF_JVY,i)+famy
            ppiclf_ydot(PPICLF_JVZ,i) = ppiclf_ydot(PPICLF_JVZ,i)+famz
            ppiclf_ydot(PPICLF_JT,i)  = 0.0d0
            ppiclf_ydot(PPICLF_JOX,i) = 0.0d0
            ppiclf_ydot(PPICLF_JOY,i) = 0.0d0
            ppiclf_ydot(PPICLF_JOZ,i) = 0.0d0
         endif


      enddo ! do i=1,ppiclf_npart


!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!06/05/2024 - Thierry - Store density-weighted acceleration
!
      ! Briney Added Mass flag
      if (am_flag==2) then 
         do i=1,ppiclf_npart
     
            ! Substantial derivative of density - how rocflu does it 
            SDrho = ppiclf_rprop(PPICLF_R_JRHSR,i)
     >         + ppiclf_y(PPICLF_JVX,i) * ppiclf_rprop(PPICLF_R_JPGCX,i)
     >         + ppiclf_y(PPICLF_JVY,i) * ppiclf_rprop(PPICLF_R_JPGCY,i)
     >         + ppiclf_y(PPICLF_JVZ,i) * ppiclf_rprop(PPICLF_R_JPGCZ,i)
      
            ! Fluid density
            rhof   = ppiclf_rprop(PPICLF_R_JRHOF,i)

            ! X-acceleration
            ppiclf_rprop(PPICLF_R_WDOTX,i) = 
     >                  ppiclf_rprop(PPICLF_R_JSDRX,i) 
     >                 -(rhof*ppiclf_ydot(PPICLF_JVX,i)) 
     >                 -(ppiclf_y(PPICLF_JVX,i)*SDrho)
          
            ! Y-acceleration
            ppiclf_rprop(PPICLF_R_WDOTY,i) = 
     >                  ppiclf_rprop(PPICLF_R_JSDRY,i) 
     >                 -(rhof*ppiclf_ydot(PPICLF_JVY,i))
     >                 -(ppiclf_y(PPICLF_JVY,i)*SDrho)
          
            ! Z-acceleration
            ppiclf_rprop(PPICLF_R_WDOTZ,i) = 
     >                  ppiclf_rprop(PPICLF_R_JSDRZ,i) 
     >                 -(rhof*ppiclf_ydot(PPICLF_JVZ,i))
     >                 -(ppiclf_y(PPICLF_JVZ,i)*SDrho)

            ! write out for debug
            if (ppiclf_debug==2) then
            if (ppiclf_nid==0 .and. iStage==1) then
            if (mod(idebug,10)==0) then
               if (i<=3) then
                  write(7020+i,*) i, ppiclf_time, rhof,
     >                ppiclf_rprop(PPICLF_R_JSDRX,i),                   
     >                ppiclf_rprop(PPICLF_R_JSDRY,i), 
     >                ppiclf_rprop(PPICLF_R_JSDRZ,i),
     >                ppiclf_ydot(PPICLF_JVX,i),
     >                ppiclf_ydot(PPICLF_JVY,i),
     >                ppiclf_ydot(PPICLF_JVZ,i),
     >                ppiclf_y(PPICLF_JVX,i),
     >                ppiclf_y(PPICLF_JVY,i),
     >                ppiclf_y(PPICLF_JVZ,i)

                  write(7030+i,*) i, ppiclf_time, 
     >              ppiclf_rprop(PPICLF_R_WDOTX:PPICLF_R_WDOTZ,i)
               endif
            endif
            endif
            endif

         enddo
      endif

!
! ----------------------------------------------------------------------
!

      ! Use ppiclf ALLREDUCE to compute values across processors
      ! Note that ALLREDUCE uses MPI_BARRIER, which is cpu expensive
      ! Print out every 10th iStage=1 counts
      if (ppiclf_debug   .ge. 1) then
      if (iStage         .eq. 1) then
      if (mod(idebug,10) .eq. 0) then

         call ppiclf_user_debug

      endif
      endif
      endif

!
! ----------------------------------------------------------------------
!
!

      !
      ! Reset arrays for Viscous Unsteady Force
      !
      if (ViscousUnsteady_flag>=1) then
         if (iStage==3) call ppiclf_user_ShiftUnsteadyData
         call ppiclf_user_plag2prop(ppiclf_nUnsteadyData)
      endif


! ----------------------------------------------------------------------

      return
      end
!-----------------------------------------------------------------------
!
! Created Oct. 18, 2024
!
! Subroutine to map both real and ghost particles to subbins
! for nearest neighbor search
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ppiclf_user_subbinMap(i_Bin, n_SBin, tot_SBin,
     >                                  SBin_counter, SBin_map)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4  SBin_map( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >       + 1) - 1), (ppiclf_npart+ppiclf_npart_gp))
      INTEGER*4  SBin_counter( 0 : (
     > (FLOOR((ppiclf_bins_dx(1)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(2)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3))
     >        + 1) *
     > (FLOOR((ppiclf_bins_dx(3)+2*ppiclf_d2chk(3))/ppiclf_d2chk(3)) 
     >       + 1) - 1))
      INTEGER*4  i_Bin(3), n_SBin(3), tot_SBin
     >          
!
! Internal:
!
      REAL*8    xp(3), bin_xMin(3)
      INTEGER*4 temp_SBin, i_SBin(3),ibinTemp, i, j, k, l 

!
! Code:
!
        ! Determine ppiclf bin in each dimension for this processor
        ! All real particles are in the same bin.  Look at 1st r particle
        DO l = 1,3
            i_Bin(l) = FLOOR((ppiclf_y(l,1) - ppiclf_binb(2*l-1))
     >                 /ppiclf_bins_dx(l))
            bin_xMin(l) = ppiclf_binb(2*l-1)+i_Bin(l)*ppiclf_bins_dx(l) 
        END DO 
        ! Determine the number of subbins in each dimension
        DO l = 1,3
          IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
            n_SBin(l) = FLOOR((ppiclf_bins_dx(l)+2*ppiclf_d2chk(3))
     >                       /ppiclf_d2chk(3))
          ELSE
            n_SBin(l) = 0
          END IF
        END DO

        ! Determine total number of subbins
        tot_SBin = (n_SBin(1)+1)*(n_SBin(2)+1)
        IF (ppiclf_ndim .EQ. 3) THEN
          tot_SBin = tot_SBin*(n_SBin(3)+1)
        END IF
        ! Assign Subbin counters to 0
        SBin_counter = 0

        ! Map each real particle to a subbin
        DO i = 1 , ppiclf_npart
           DO l = 1,3
              IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
                 xp(l) = ppiclf_y(l,i)
              ELSE
                 xp(l) = 0.0
              END IF
           END DO 
           ! Determine subbin
           DO l = 1,3
              IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
                  i_SBin(l) = FLOOR((xp(l) - (bin_xMin(l) 
     >            - ppiclf_d2chk(3)))/ppiclf_d2chk(3)) 
              ELSE
                 i_SBin(l) = 0
              END IF
           END DO
           temp_SBin = i_SBin(1) + n_SBin(1)*i_SBin(2) +
     >                n_SBin(1)*n_SBin(2)*i_SBin(3)
           SBin_counter(temp_SBin) = SBin_counter(temp_SBin) + 1
           SBin_map(temp_SBin,SBin_counter(temp_SBin)) = i
        END DO ! real particle loop


        ! Map each ghost particle to a subbin
        DO i = 1 , ppiclf_npart_gp
          DO l = 1,3
            IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
              xp(l) = ppiclf_rprop_gp(l,i)
            ELSE
              xp(l) = 0.0
            END IF
          END DO
          ! Only map ghost particles within one neighborwidth
          ! from bin edge to subbins. All others are outside
          ! of collision search distance.
          IF (xp(1) .GT. (bin_xMin(1)-ppiclf_d2chk(3))
     >  .AND. xp(2) .GT. (bin_xMin(2)-ppiclf_d2chk(3))
     >  .AND. xp(3) .GT. (bin_xMin(3)-ppiclf_d2chk(3))
     >  .AND. xp(1) .LT. (bin_xMin(1)+ppiclf_bins_dx(1)+ppiclf_d2chk(3))
     >  .AND. xp(2) .LT. (bin_xMin(2)+ppiclf_bins_dx(2)+ppiclf_d2chk(3))
     >  .AND. xp(3) .LT. (bin_xMin(3)+ppiclf_bins_dx(3)+ppiclf_d2chk(3))
     >        ) THEN
            ! Determine subbin
            DO l = 1,3
              IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
                i_SBin(l) = FLOOR((xp(l) - (bin_xMin(l) 
     >          - ppiclf_d2chk(3)))/ppiclf_d2chk(3)) 
              ELSE
                i_SBin(l) = 0
              END IF
            END DO
            temp_SBin = i_SBin(1) + n_SBin(1)*i_SBin(2) +
     >                 n_SBin(1)*n_SBin(2)*i_SBin(3)
            SBin_counter(temp_SBin) = SBin_counter(temp_SBin) + 1
            ! negative in subbin map means it is ghost particle
            SBin_map(temp_SBin,SBin_counter(temp_SBin)) = -i
          END IF 
        END DO ! gp loop

      RETURN
      END 
!-----------------------------------------------------------------------
!
! Created Spet. 14, 2024
!
! Subroutine for unit test problems
!      
! stationary = -1: azimuthal velocity only
!            = -2: radial velocity only
!            = -3: radial + azimuthal velocity
!
! NOTE: Time step ppiclf_dt is hardcoded in ppiclf_solve_IntegrateRK3s_Rocflu
!       subroutine to be ppiclf_dt = 5.0000000000000004E-008      
!-----------------------------------------------------------------------
      subroutine ppiclf_user_unit_tests(i,iStage,famx,famy,famz)
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

      integer*4 i, iStage
      real*8 famx, famy, famz
!
! Code:
!
      if (stationary==-1) then
         call ppiclf_user_unit01(i,iStage,famx,famy,famz)
      elseif (stationary==-2) then
         call ppiclf_user_unit02(i,iStage,famx,famy,famz)
      elseif (stationary==-3) then
         call ppiclf_user_unit03(i,iStage,famx,famy,famz)
      endif

      return
      end
!      
!-----------------------------------------------------------------------
!
! Azimuthal(theta) velocity only
!
      subroutine ppiclf_user_unit01(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
      real*8 pi, omg00, omg02, tau

!
! Code:
!
      pi = acos(-1.0d0)
      tau = 0.1d-3 ! time for 2pi revolution. Can be max simulation time
      omg00 = (2.0d0*pi)/tau
      omg02 = omg00*omg00

      if(ppiclf_time .eq. 0.0) then
        ppiclf_y(PPICLF_JVX,i) = -omg00*ppiclf_y(PPICLF_JY,i)
        ppiclf_y(PPICLF_JVY,i) =  omg00*ppiclf_y(PPICLF_JX,i)
        ppiclf_y(PPICLF_JVZ,i) =  0.0d0
      endif
      
      famx = -omg02*ppiclf_y(PPICLF_JX,i)
      famy = -omg02*ppiclf_y(PPICLF_JY,i)
      famz = 0.0d0

      return
      end
!
!-----------------------------------------------------------------------
!
! Radial velocity only
!
      subroutine ppiclf_user_unit02(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
      real*8 drdt00, drdt02, theta
!
! Code:
!
      ! dx/dt = dr/dt * cos theta
      ! dy/dt = dr/dt * sin theta
      drdt00 = 79.928163327808903
      drdt02 = drdt00*drdt00
      
      ! angle between particle and +ve x-axis
      theta = atan2(ppiclf_y(PPICLF_JY,i), ppiclf_y(PPICLF_JX,i))

      if(ppiclf_time .eq. 0.0) then
        ppiclf_y(PPICLF_JVX,i) =  -drdt00*cos(theta)
        ppiclf_y(PPICLF_JVY,i) =  -drdt00*sin(theta)
        ppiclf_y(PPICLF_JVZ,i) =  0.0d0
      endif
      
      famx = -drdt02*cos(theta)
      famy = -drdt02*sin(theta)
      famz = 0.0d0


      return
      end
!
!-----------------------------------------------------------------------
!
! Radial + Azimuthal velocity
!
      subroutine ppiclf_user_unit03(i,iStage,famx,famy,famz)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i, iStage
      real*8 famx, famy, famz
      real*8 pi, tau, omg00, omg02, drdt00, theta
!
! Code:
!
      ! dx/dt = dr/dt * cos theta
      ! dy/dt = dr/dt * sin theta
      
      pi = acos(-1.0d0)

      tau = 0.1d-3
      omg00 = (2.0d0*pi)/tau
      omg02 = omg00*omg00
      drdt00 = 79.928163327808903
      
      ! angle between particle and +ve x-axis
      theta = atan2(ppiclf_y(PPICLF_JY,i), ppiclf_y(PPICLF_JX,i))

      if(ppiclf_time .eq. 0.0) then
        ppiclf_y(PPICLF_JVX,i) = -omg00*ppiclf_y(PPICLF_JY,i)
     >                           -drdt00*cos(theta)
        ppiclf_y(PPICLF_JVY,i) =  omg00*ppiclf_y(PPICLF_JX,i) 
     >                           -drdt00*sin(theta)
        ppiclf_y(PPICLF_JVZ,i) =  0.0d0
      endif
      
      famx = 2.0d0*drdt00*omg00*sin(theta) 
     >       -omg02*ppiclf_y(PPICLF_JX,i)
      famy = -2.0d0*drdt00*omg00*cos(theta) 
     >       -omg02*ppiclf_y(PPICLF_JY,i)
      famz = 0.0d0

      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for quasi-steady force
!
! Quasi-steady force (Re_p and Ma_p corrections):
!   Improved Drag Correlation for Spheres and Application 
!   to Shock-Tube Experiments 
!   - Parmar et al. (2010)
!   - AIAA Journal
!
! Quasi-steady force (phi corrections):
!   The Added Mass, Basset, and Viscous Drag Coefficients 
!   in Nondilute Bubbly Liquids Undergoing Small-Amplitude 
!   Oscillatory Motion
!   - Sangani et al. (1991)
!   - Phys. Fluids A
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_Parmar(i,beta,cd)
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

      integer*4 i
      real*8 gamma,mp,phi,re
      real*8 rcd1,rmacr,rcd_mcr,rcd_std,rmach_rat,rcd_M1,
     >   rcd_M2,C1,C2,C3,f1M,f2M,f3M,lrep,factor,cd,beta

!
! Code:
!
      gamma = 1.4d0
      mp  = dmax1(rmachp,0.01d0)
      phi = dmax1(rphip,0.0001d0)
      re  = dmax1(rep,0.1d0)

      if(re .lt. 1E-14) then
         rcd1 = 1.0
      else 
         rmacr= 0.6 ! Critical rmachp no.
         rcd_mcr = (1.+0.15*re**(0.684)) + 
     >                (re/24.0)*(0.513/(1.+483./re**(0.669)))
       if (mp .le. rmacr) then
          rcd_std = (1.+0.15*re**(0.687)) + 
     >                (re/24.0)*(0.42/(1.+42500./re**(1.16)))
          rmach_rat = mp/rmacr
          rcd1 = rcd_std + (rcd_mcr - rcd_std)*rmach_rat
       else if (mp .le. 1.0) then
         rcd_M1 = (1.0+0.118*re**0.813) +
     >                (re/24.0)*0.69/(1.0+3550.0/re**.793)
         C1 =  6.48
         C2 =  9.28
         C3 = 12.21
         f1M = -1.884 +8.422*mp -13.70*mp**2 +8.162*mp**3
         f2M = -2.228 +10.35*mp -16.96*mp**2 +9.840*mp**3
         f3M =  4.362 -16.91*mp +19.84*mp**2 -6.296*mp**3
         lrep = log(re)
         factor = f1M*(lrep-C2)*(lrep-C3)/((C1-C2)*(C1-C3))
     >              +f2M*(lrep-C1)*(lrep-C3)/((C2-C1)*(C2-C3))
     >              +f3M*(lrep-C1)*(lrep-C2)/((C3-C1)*(C3-C2)) 
         rcd1 = rcd_mcr + (rcd_M1-rcd_mcr)*factor
       else if (mp .lt. 1.75) then
         rcd_M1 = (1.0+0.118*re**0.813) +
     >              (re/24.0)*0.69/(1.0+3550.0/re**.793)
         rcd_M2 = (1.0+0.107*re**0.867) +
     >              (re/24.0)*0.646/(1.0+861.0/re**.634)
         C1 =  6.48
         C2 =  8.93
         C3 = 12.21
         f1M = -2.963 +4.392*mp -1.169*mp**2 -0.027*mp**3
     >             -0.233*exp((1.0-mp)/0.011)
         f2M = -6.617 +12.11*mp -6.501*mp**2 +1.182*mp**3
     >             -0.174*exp((1.0-mp)/0.010)
         f3M = -5.866 +11.57*mp -6.665*mp**2 +1.312*mp**3
     >             -0.350*exp((1.0-mp)/0.012)
         lrep = log(re)
         factor = f1M*(lrep-C2)*(lrep-C3)/((C1-C2)*(C1-C3))
     >              +f2M*(lrep-C1)*(lrep-C3)/((C2-C1)*(C2-C3))
     >              +f3M*(lrep-C1)*(lrep-C2)/((C3-C1)*(C3-C2)) 
         rcd1 = rcd_M1 + (rcd_M2-rcd_M1)*factor
       else
         rcd1 = (1.0+0.107*re**0.867) +
     >                (re/24.0)*0.646/(1.0+861.0/re**.634)
       end if ! mp
      endif    ! re

      cd = (24.0/re)*rcd1*(1.+2.*phi)/((1.-phi)**3)

      beta = rcd1*3.0*rpi*rmu*dp

      beta = beta*(1.+2.*phi)/((1.-phi)**3)


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
! Modified March 5, 2024
!
! Subroutine for quasi-steady force
!
! QS Force calculated as a function of Re, Ma and phi
!
! Use Osnes etal (2023) correlations
! A.N. Osnes, M. Vartdal, M. Khalloufi, 
!    J. Capecelatro, and S. Balachandar.
! Comprehensive quasi-steady force correlations for compressible flow
!    through random particle suspensions.
! International Journal of Multiphase Flow, Vol. 165, 104485, (2023).
! doi: https://doi.org/10.1016/j.imultiphaseflow.2023.104485.
!
! E. Loth, J.T. Daspit, M. Jeong, T. Nagata, and T. Nonomura.
! Supersonic and hypersonic drag coefficients for a sphere.
! AIAA Journal, Vol. 59(8), pp. 3261-3274, (2021).
! doi: https://doi.org/10.2514/1.J060153.
!
! NOTE: Re<45 Rarified formula of Loth et al has been redefined by Balachandar
! to avoid singularity as Ma -> 0.
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_Osnes(i,beta,cd)
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

      integer*4 i
      real*8 gamma,mp,phi,re,Knp,fKn,CD1,s,JM,CD2,
     >   cd_loth,CM,GM,HM,b1,b2,b3,cd,beta
      real*8 sgby2, JMt

!
! Code:
!
      gamma = 1.4d0
      mp  = dmax1(rmachp,0.01d0)
      phi = dmax1(rphip,0.0001d0)
      re  = dmax1(rep,0.1d0)

      ! Loth's correlation
      if (re .le. 45.0) then
         ! Rarefied-dominated regime
         Knp = sqrt(0.5d0 * rpi * gamma) * mp / re
         if (Knp > 0.01) then
            fKn = 1.0d0 / (1.0d0 
     >            + Knp*(2.514d0 + 0.8d0*exp(-0.55d0/Knp)))
         else
            fKn = 1.0d0 / (1.0d0 
     >            + Knp*(2.514d0 + 0.8d0*exp(-0.55d0/0.01)))

         end if
         CD1 = (24.0/re)*(1.0d0 + 0.15d0* re**(0.687d0)) * fKn
         s = mp * sqrt(0.5d0 * gamma)
         sgby2 = sqrt(0.5d0 * gamma)
         if (mp <= 1) then
            !JMt = 2.26d0*(mp**4) - 0.1d0*(mp**3) + 0.14d0*mp
            JMt = 2.26d0*(mp**4) + 0.14d0*mp
         else
            JMt = 1.6d0*(mp**4) + 0.25d0*(mp**3) 
     >            + 0.11d0*(mp**2) + 0.44d0*mp
         end if
!
! Reformulated version of Loth et al. to avoid singularity at mp = 0
!
         CD2 = (1.0d0 + 2.0d0*(s**2)) * exp(-s**2) * mp
     >          / ((sgby2**3)*sqrt(rpi)) 
     >          + (4.0d0*(s**4) + 4.0d0*(s**2) - 1.0d0) 
     >          * erf(s) / (2.0d0*(sgby2**4)) 
     >          + (2.0d0*(mp**3) / (3.0d0 * sgby2)) * sqrt(rpi)

         CD2 = CD2 / (1.0d0 + (((CD2/JMt) - 1.0d0) * sqrt(re/45.0d0)))
         cd_loth = CD1 / (1.0d0 + (mp**4)) 
     >          +  CD2 / (1.0d0 + (mp**4))
      else
         ! Compression-dominated regime
         ! TLJ: coefficients tweaked to get continuous values
         !      on the two branches at the critical points
         if (mp < 1.5d0) then
            CM = 1.65d0 + 0.65d0*tanh(4d0*mp - 3.4d0)
         else
            !CM = 2.18d0 - 0.13d0*tanh(0.9d0*mp - 2.7d0)
            CM = 2.18d0 - 0.12913149918318745d0*tanh(0.9d0*mp - 2.7d0)
         end if
         if (mp < 0.8) then
            GM = 166.0d0*(mp**3) + 3.29d0*(mp**2) - 10.9d0*mp + 20.d0
         else
            !GM = 5.0d0 + 40.d0*(mp**(-3))
            GM = 5.0d0 + 47.809331200000017d0*(mp**(-3))
         end if
         if (mp < 1) then
            HM = 0.0239d0*(mp**3) + 0.212d0*(mp**2) 
     >           - 0.074d0*mp + 1.d0
         else
            !HM =   0.93d0 + 1.0d0 / (3.5d0 + (mp**5))
            HM =   0.93967777777777772d0 + 1.0d0 / (3.5d0 + (mp**5))
         end if

         cd_loth = (24.0/re)*(1 + 0.15 * (re**(0.687)))*HM + 
     >      0.42*CM/(1+42500/re**(1.16*CM) + GM/sqrt(re))

      end if

      b1 = 5.81*phi/((1.0-phi)**2) + 
     >     0.48*(phi**(1.d0/3.d0))/((1.0-phi)**3)

      b2 = ((1.0-phi)**2)*(phi**3)*
     >     re*(0.95+0.61*(phi**3)/((1.0-phi)*2))

      b3 = dmin1(sqrt(20.0d0*mp),1.0d0)*
     >     (5.65*phi-22.0*(phi**2)+23.4*(phi**3))*
     >     (1+tanh((mp-(0.65-0.24*phi))/0.35))

      cd = cd_loth/(1.0-phi) + b3 + (24.0/re)*(1.0-phi)*(b1+b2)

      beta = 3.0*rpi*rmu*dp*(re/24.0)*cd


      return
      end
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Quasi-steady force fluctuations
! Osnes, Vartdal, Khalloufi, Capecelatro (2023)
!   Comprehensive quasi-steady force correlations
!   for compressible flow through random particle
!   suspensions.
!   International Journal of Multiphase Flows,
!   Vo. 165, 104485.
! Lattanzi, Tavanashad, Subramaniam, Capecelatro (2022)
!   Stochastic model for the hydrodynamic force in
!   Euler-Lagrange silumations of particle-laden flows.
!   Physical Review Fluids, Vol. 7, 014301.
! Note: To compute the granular temperature, we assume
!   the velocity fluctuations are uncorrelated.
! Note: The means are computed using a box filter with an
!   adaptive filter width.
! Compute mean using box filter for langevin model - not for fedback
!
! The mean is calcuated according to Lattanzi etal,
!   Physical Review Fluids, 2022.
!
! Sam - TODO: either couple the projection filter to the
! fluctuation filter or completely decouple them.
! Right now they use the same filter width and are assumed to
! be the same. 
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_fluct_Lattanzi(i,iStage,fqs_fluct)
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
      integer*4 i, iStage
      real*8 fqs_fluct(3)

      real*8 aSDE,bq,bSDE,chi,denum,dW1,dW2,dW3,fq,Fs,gkern,
     >   sigF,tF_inv,theta,upflct,vpflct,wpflct,Z1,Z2,Z3
      real*8 TwoPi

!
! Code:
!
      TwoPi = 2.0d0*acos(-1.0d0)

      if (qs_fluct_filter_flag==0) then
         denum = max(dfloat(icpmean),1.d0)  ! for arithmetic mean
      else if (qs_fluct_filter_flag==1) then
         denum = max(phipmean,1.0e-6)  ! for volume mean
      endif
      upmean = upmean / denum
      vpmean = vpmean / denum
      wpmean = wpmean / denum
      u2pmean = u2pmean / denum
      v2pmean = v2pmean / denum
      w2pmean = w2pmean / denum

      if (ppiclf_debug==2) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,"(2x,E16.8,i5,16(1x,F13.8))") ppiclf_time,
     >                denum,upmean,vpmean,wpmean
         endif
      endif

      ! Lattenzi is valid only for incompressible flows,
      !   so here we use the compressible correction of Osnes
      !   Equations (9-12) of Osnes paper, with eqn (9) corrected
      !   Note that sigD has units of force; N = Pa-m^2
      fq = 6.52*rphip - 22.56*(rphip**2) + 49.90*(rphip**3)
      Fs = 3.0*rpi*rmu*dp*(1.0+0.15*((rep*(1.0-rphip))**0.687))
     >          *(1.0-rphip)*vmag
      bq = min(sqrt(20.0*rmachp),1.0)*0.55*(rphip**0.7) 
     >          *(1.0+tanh((rmachp-0.5)/0.2))
      sigF = (fq + bq)*Fs

      chi = (1.0+2.50*rphip+4.51*(rphip**2)+4.52*(rphip**3))
     >         /((1.0-(rphip/0.64)**3)**0.68)

      fqs_fluct = 0.0d0

! Particle velocity fluctuation and Granular Temperature
! Need particle velocity mean
! Though the theory assumes granular temperature to be an 
!    average over neighboring particles, here it is approximated 
!    as that of the chosen particle - Comment 3/6/24
!    This is now fixed - Comment 4/12/24
!
      ! Particle velocity fluctuation
      upflct = ppiclf_y(PPICLF_JVX,i) - upmean
      vpflct = ppiclf_y(PPICLF_JVY,i) - vpmean
      wpflct = ppiclf_y(PPICLF_JVZ,i) - wpmean

      ! Granular temperature
      ! theta = (upflct*upflct + vpflct*vpflct + wpflct*wpflct)/3.0
      ! This is averaged over neighboring particles
      theta  = ((u2pmean + v2pmean + w2pmean) - 
     >          (upmean**2 + vpmean**2 + wpmean**2))/3.0d0

      ! 11/21/24 - Thierry - prevent NaN variables
      if(theta.le.1.d-12) then
        theta = 0.0d0
      endif

      tF_inv = (24.0*rphip*chi)/dp * sqrt(theta/rpi)

      aSDE = tF_inv
      bSDE = sigF*sqrt(2.0*tF_inv)

      call RANDOM_NUMBER(UnifRnd)

      Z1 = sqrt(-2.0d0*log(UnifRnd(1)))*cos(TwoPi*UnifRnd(2))
      Z2 = sqrt(-2.0d0*log(UnifRnd(3)))*cos(TwoPi*UnifRnd(4))
      Z3 = sqrt(-2.0d0*log(UnifRnd(5)))*cos(TwoPi*UnifRnd(6))

      dW1 = sqrt(fac)*Z1
      dW2 = sqrt(fac)*Z2
      dW3 = sqrt(fac)*Z3

      fqs_fluct(1) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFX,i)
     >           + bSDE*dW1
      fqs_fluct(2) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFY,i)
     >           + bSDE*dW2
      fqs_fluct(3) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFZ,i)
     >           + bSDE*dW3


      if (ppiclf_debug==2 .and. (iStage==1 .and. ppiclf_nid==0)) then
         if (ppiclf_time.gt.2.d-8) then
         if (i<=10) then
            write(7350+(i-1)*1,*) i,ppiclf_time,             ! 0-1
     >         rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,  ! 2-10
     >         rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac, ! 11-18
     >         vx,vy,vz,ppiclf_dt,                           ! 19-22
     >         ppiclf_npart,ppiclf_n_bins(1:3),              ! 23-26
     >         ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3), ! 27
     >         ppiclf_binb(1:6),                             ! 28-33
     >         upmean,vpmean,wpmean,phipmean,                ! 34-37
     >         ppiclf_y(PPICLF_JVX:PPICLF_JVZ,i),            ! 38-40
     >         upflct,vpflct,wpflct,icpmean,                 ! 41-44
     >         fq,Fs,bq,theta,chi,tF_inv,                    ! 45-50
     >         aSDE,bSDE,sigF,                               ! 51-53
     >         fqs_fluct(1:3),                               ! 54-56
     >         Z1,Z2,Z3,dW1,dW2,dW3,                         ! 57-62
     >         ppiclf_np
         endif
         endif
      endif


      return
      end
!
!
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024 - T.L. Jackson
! Modified 3/6/24 - Balachandar
!
! Quasi-steady force fluctuations
! Osnes, Vartdal, Khalloufi, Capecelatro (2023)
!   Comprehensive quasi-steady force correlations
!   for compressible flow through random particle
!   suspensions.
!   International Journal of Multiphase Flows,
!   Vo. 165, 104485.
! Lattanzi, Tavanashad, Subramaniam, Capecelatro (2022)
!   Stochastic model for the hydrodynamic force in
!   Euler-Lagrange silumations of particle-laden flows.
!   Physical Review Fluids, Vol. 7, 014301.
! Note: To compute the granular temperature, we assume
!   the velocity fluctuations are uncorrelated.
! Note: The means are computed using a box filter with an
!   adaptive filter width.
! Compute mean using box filter for langevin model - not for fedback
!
! The mean is calcuated according to Lattanzi etal,
!   Physical Review Fluids, 2022.
!
! Sam - TODO: either couple the projection filter to the
! fluctuation filter or completely decouple them.
! Right now they use the same filter width and are assumed to
! be the same. 
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_QS_fluct_Osnes(i,iStage,fqs_fluct)
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

      integer*4 i, iStage
      real*8 fqs_fluct(3)

      real*8 aSDE,bq,chi,denum,dW1,dW2,fq,Fs,gkern,
     >   sigD,tF_inv,theta,upflct,vpflct,wpflct,Z1,Z2
      real*8 TwoPi
      real*8 bSDE_CD, bSDE_CL,CD_frac,CD_prime
      real*8 sigT,sigCT
      real*8 sigmoid_cf, f_CF
      real*8 avec(3)
      real*8 bvec(3)
      real*8 cvec(3)
      real*8 dvec(3)
      real*8 cosrand,sinrand
      real*8 eunit(3)

!
! Code:
!
      TwoPi = 2.0d0*acos(-1.0d0)

      if (qs_fluct_filter_flag==0) then
         denum = max(dfloat(icpmean),1.d0)  ! for arithmetic mean
      else if (qs_fluct_filter_flag==1) then
         denum = max(phipmean,1.0e-6)  ! for volume mean
      endif
      upmean = upmean / denum
      vpmean = vpmean / denum
      wpmean = wpmean / denum
      u2pmean = u2pmean / denum
      v2pmean = v2pmean / denum
      w2pmean = w2pmean / denum

      if (ppiclf_debug==2) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,*) 'FLUC1 ',
     >          ppiclf_time,denum,upmean,vpmean,wpmean,
     >          abs(upmean),abs(vpmean),abs(wpmean),
     >          u2pmean,v2pmean,w2pmean
         endif
      endif

      ! Equations (9-12) of Osnes paper, with eqn (9) corrected
      ! Note that sigD has units of force; N = Pa-m^2
      fq = 6.52*rphip - 22.56*(rphip**2) + 49.90*(rphip**3)
      Fs = 3.0*rpi*rmu*dp*(1.0+0.15*((rep*(1.0-rphip))**0.687))
     >          *(1.0-rphip)*vmag
      bq = min(sqrt(20.0*rmachp),1.0)*0.55*(rphip**0.7) 
     >          *(1.0+tanh((rmachp-0.5)/0.2))
      sigD = (fq + bq)*Fs

! Particle velocity fluctuation and Granular Temperature
! Need particle velocity mean
! Though the theory assumes granular temperature to be an 
!    average over neighboring particles, here it is approximated 
!    as that of the chosen particle - Comment 3/6/24
!    This is now fixed - Comment 4/12/24
!
      upflct = ppiclf_y(PPICLF_JVX,i) - upmean
      vpflct = ppiclf_y(PPICLF_JVY,i) - vpmean
      wpflct = ppiclf_y(PPICLF_JVZ,i) - wpmean

      ! Granular temperature
      ! theta = (upflct*upflct + vpflct*vpflct + wpflct*wpflct)/3.0
      ! This is averaged over neighboring particles
      theta  = ((u2pmean + v2pmean + w2pmean) - 
     >          (upmean**2 + vpmean**2 + wpmean**2))/3.0d0

      ! 11/21/24 - Thierry - prevent NaN variables
      if(theta.le.1.d-12) then
        theta = 0.0d0
      endif

      chi = (1.0 + 2.50*rphip + 4.51*(rphip**2) + 4.52*(rphip**3))
     >         /((1.0-(rphip/0.64)**3)**0.68)

      tF_inv = (24.0*rphip*chi/dp) * sqrt(theta/rpi)

      aSDE = tF_inv
      bSDE_CD = sigD*sqrt(2.0*tF_inv)  ! Modified 3/6/24

! Fluctuating perpendicular force
! Compare CD' (units of force) against sigma_CD (units of force) 
!    to determine which of 5 bins to use for sigCT
!
! Added 3/6/24 
! Modified 3/14/24 
!
      avec = [vx,vy,vz]/max(1.d-8,vmag)

      CD_prime = ppiclf_rprop(PPICLF_R_FLUCTFX,i)*avec(1) +
     >           ppiclf_rprop(PPICLF_R_FLUCTFY,i)*avec(2) +
     >           ppiclf_rprop(PPICLF_R_FLUCTFZ,i)*avec(3)
      CD_frac  = CD_prime/sigD

      ! 11/21/24 - Thierry - prevent NaN variables
      if(CD_prime.eq.0.0 .and. sigD.eq.0.0) then
        CD_frac = 0.0d0
      endif

      ! Thierry Daoud - Updated June 2, 2024
      sigmoid_cf = 1.0 / (1.0 + exp(-CD_frac))
      f_CF = 0.39356905*sigmoid_cf + 0.43758848
      sigT  = f_CF*sigD
      bSDE_CL = sigT*sqrt(2.0*tF_inv)

      if (ppiclf_debug==2) then
      if (i<=4) then
         if (ppiclf_nid==0 .and. iStage==1) then
            write(6,*) 'FLUC1 ',i,
     >        CD_prime,CD_frac,sigD,theta,bSDE_CD,bSDE_CL
         endif
      endif
      endif


! Calculate the three orthogonal unit vectors
! The first vector (avec) is vx/vmag, vy/vmag, and vz/vmag
! The second (bvec) is constructued by taking cross-product with eunit
!   Note: if avec is in dir. of e_x=(1,0,0), use e_y=(0,1,0) to get e_z
!       : if avec is in dir. of e_y=(0,1,0), use e_z=(0,0,1) to get e_x
!       : if avec is in dir. of e_z=(0,0,1), use e_x=(1,0,0) to get e_y
! The third (cvec) is cross product of the first two
! written 3/6/24
!
      eunit = [1,0,0]
      if (abs(avec(2))+abs(avec(3)) <= 1.d-8) then
         eunit = [0,1,0]
      elseif (abs(avec(1))+abs(avec(3)) <= 1.d-8) then
         eunit = [0,0,1]
      endif

      bvec(1) = avec(2)*eunit(3) - avec(3)*eunit(2)
      bvec(2) = avec(3)*eunit(1) - avec(1)*eunit(3)
      bvec(3) = avec(1)*eunit(2) - avec(2)*eunit(1)
      denum   = max(1.d-8,sqrt(bvec(1)**2 + bvec(2)**2 + bvec(3)**2))
      bvec    = bvec / denum

      cvec(1) = avec(2)*bvec(3) - avec(3)*bvec(2)
      cvec(2) = avec(3)*bvec(1) - avec(1)*bvec(3)
      cvec(3) = avec(1)*bvec(2) - avec(2)*bvec(1)
      denum   = max(1.d-8,sqrt(cvec(1)**2 + cvec(2)**2 + cvec(3)**2))
      cvec    = cvec / denum

      call RANDOM_NUMBER(UnifRnd)

      Z1 = sqrt(-2.0d0*log(UnifRnd(1)))*cos(TwoPi*UnifRnd(2))
      Z2 = sqrt(-2.0d0*log(UnifRnd(3)))*cos(TwoPi*UnifRnd(4))

      dW1 = sqrt(fac)*Z1
      dW2 = sqrt(fac)*Z2

! Random mixture of bvec and cvec - make sure the new one is a unit vector
! Added 3/6/24
!
      cosrand = cos(TwoPi*UnifRnd(5))
      sinrand = sin(TwoPi*UnifRnd(5)) 
      dvec(1) = bvec(1)*cosrand + cvec(1)*sinrand
      dvec(2) = bvec(2)*cosrand + cvec(2)*sinrand
      dvec(3) = bvec(3)*cosrand + cvec(3)*sinrand
      denum   = max(1.d-8,sqrt(dvec(1)**2 + dvec(2)**2 + dvec(3)**2))
      dvec    = dvec/denum

      fqs_fluct(1) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFX,i)
     >           + bSDE_CD*dW1*avec(1) + bSDE_CL*dW2*dvec(1)
      fqs_fluct(2) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFY,i)
     >           + bSDE_CD*dW1*avec(2) + bSDE_CL*dW2*dvec(2)
      fqs_fluct(3) = 
     >           (1.0-aSDE*fac)*ppiclf_rprop(PPICLF_R_FLUCTFZ,i)
     >           + bSDE_CD*dW1*avec(3) + bSDE_CL*dW2*dvec(3)


      if (ppiclf_debug==2 .and. (iStage==1 .and. ppiclf_nid==0)) then
         if (ppiclf_time.gt.2.d-8) then
         if (i<=10) then
            write(7350+(i-1)*1,*) i,ppiclf_time,             ! 0-1
     >         rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,  ! 2-10
     >         rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac, ! 11-18
     >         vx,vy,vz,ppiclf_dt,                           ! 19-22
     >         ppiclf_npart,ppiclf_n_bins(1:3),              ! 23-26
     >         ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3), ! 27
     >         ppiclf_binb(1:6),                             ! 28-33
     >         upmean,vpmean,wpmean,phipmean,                ! 34-37
     >         ppiclf_y(PPICLF_JVX:PPICLF_JVZ,i),            ! 38-40
     >         upflct,vpflct,wpflct,icpmean,                 ! 41-44
     >         fq,Fs,bq,theta,chi,tF_inv,                    ! 45-50
     >         aSDE,bSDE_CD,bSDE_CL,                         ! 51-53
     >         sigD,sigT,                                    ! 54-55
     >         CD_prime,CD_frac,sigmoid_cf,f_CF,             ! 56-59
     >         eunit,avec,bvec,                              ! 60-68
     >         cvec,dvec,rpi,                                ! 69-75
     >         fqs_fluct(1:3),                               ! 76-78
     >         Z1,Z2,dW1,dW2,                                ! 79-82
     >         ppiclf_np,                                    ! 83
     >         ppiclf_y(PPICLF_JX:PPICLF_JZ,i)               ! 84-86
         endif
         endif
      endif


      return
      end
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for added mass
!   also called the quasi-unsteady force,
!   or the inviscid unsteady force in case of the Euler equations
!      
! ADDED from rocflu
! Added mass force (phi corrections):
!   On the dispersed two-phase flow in the laminar flow regime
!   - Zuber (1964), Chem. Engng Sci.
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_AM_Parmar(i,iStage,
     >                 famx,famy,famz,rmass_add)
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
      integer*4 i, iStage
      real*8 famx, famy, famz, rmass_add
      real*8 rcd_am
      real*8 SDrho

!
! Code:
!
      ! Mach number correction
      if (rmachp .gt. 0.6) then        
         rcd_am = 1.0 + 1.8*((0.6)**2)  + 7.6*((0.6)**4)
      else
         rcd_am = 1.0 + 1.8*(rmachp**2) + 7.6*(rmachp**4)
      endif
      rcd_am = rcd_am * 0.5

      ! Volume fraction correction
      rcd_am = rcd_am*(1.0+2.0*rphip)

      rmass_add = rhof*ppiclf_rprop(PPICLF_R_JVOLP,i)*rcd_am

      !NEW Added mass, using how rocflu does it
      !1st Derivative, substantial how rocflu does it
      SDrho = ppiclf_rprop(PPICLF_R_JRHSR,i) 
     >      + ppiclf_y(PPICLF_JVX,i) * ppiclf_rprop(PPICLF_R_JPGCX,i)
     >      + ppiclf_y(PPICLF_JVY,i) * ppiclf_rprop(PPICLF_R_JPGCY,i)
     >      + ppiclf_y(PPICLF_JVZ,i) * ppiclf_rprop(PPICLF_R_JPGCZ,i)

      famx = rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i) *
     >   (ppiclf_rprop(PPICLF_R_JSDRX,i)-(ppiclf_y(PPICLF_JVX,i)*SDrho))

      famy = rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i) *
     >   (ppiclf_rprop(PPICLF_R_JSDRY,i)-(ppiclf_y(PPICLF_JVY,i)*SDrho))

      famz = rcd_am*ppiclf_rprop(PPICLF_R_JVOLP,i) *
     >   (ppiclf_rprop(PPICLF_R_JSDRZ,i)-(ppiclf_y(PPICLF_JVZ,i)*SDrho))


      return
      end
!
!
!-----------------------------------------------------------------------
!
! Created May 20, 2024
!
! Subroutine for added mass
!   also called the quasi-unsteady force,
!   or the inviscid unsteady force in case of the Euler equations
!      
! Implementing Added Mass Algorithm from S.Briney (2024)
!  
! n       = number of points
! alpha   = volume fraction
! rad     = particle radius
! d       = center-to-center distance
! rmax    = center-to-center max neighbor distance
! R       = resistance matrix (output)
! dr_max  = max interaction distance between particles considered 
! poins   = 3xn array of points x, y, z. Initialized as points(3,n)
!
! correction_analytical_always = 
!    if true, always use the analytical distant neighbor correction
!    if false, use numerical when dr_max/rad < 3.49
! 
! This subroutine corresponds to the flag:  am_flag = 2 
!
!
! IMPORTANT NOTE: THIS SUBROUTINE IS CURRENTLY ONLY 
!    VALID FOR MONODISPERSE PARTICLE BEDS 
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!                          Unary part 
!-----------------------------------------------------------------------
!
!
!
      subroutine ppiclf_user_AM_Briney_Unary(i,iStage,
     >                 famx,famy,famz,rmass_add)
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
      integer i, j, k, l, n, jj
      integer*4 iStage
      real*8 rad
      real*8 famx, famy, famz, rmass_add
      real*8 rcd_am
      real*8 SDrho

!
! Code:
!
      ! Thierry - need to decide if we want to keep 
      !    1) Mach no. correction  - yes, according to Bala
      !    2) Vol frac. correction - no, this is in binary term
      !

      ! Mach number correction
      if (rmachp .gt. 0.6) then
         rcd_am = 1.0 + 1.8*((0.6)**2)  + 7.6*((0.6)**4)
      else
         rcd_am = 1.0 + 1.8*(rmachp**2) + 7.6*(rmachp**4)
      endif
      rcd_am = rcd_am * 0.5

      ! Volume fraction correction - In the Briney model the
      !    volume fraction correction is in the binary term,
      !    not the unary term
      !rcd_am = rcd_am*(1.0+2.0*rphip)

      rmass_add = rhof*ppiclf_rprop(PPICLF_R_JVOLP,i)*rcd_am

      !NEW Added mass, using how rocflu does it
      !1st Derivative, substantial how rocflu does it
      SDrho = ppiclf_rprop(PPICLF_R_JRHSR,i)
     >      + ppiclf_y(PPICLF_JVX,i) * ppiclf_rprop(PPICLF_R_JPGCX,i)
     >      + ppiclf_y(PPICLF_JVY,i) * ppiclf_rprop(PPICLF_R_JPGCY,i)
     >      + ppiclf_y(PPICLF_JVZ,i) * ppiclf_rprop(PPICLF_R_JPGCZ,i)

      ! Take care of volume in Binary subroutine
      famx = rcd_am*
     >   (ppiclf_rprop(PPICLF_R_JSDRX,i)-(ppiclf_y(PPICLF_JVX,i)*SDrho))

      famy = rcd_am*
     >   (ppiclf_rprop(PPICLF_R_JSDRY,i)-(ppiclf_y(PPICLF_JVY,i)*SDrho))

      famz = rcd_am*
     >   (ppiclf_rprop(PPICLF_R_JSDRZ,i)-(ppiclf_y(PPICLF_JVZ,i)*SDrho))

      Fam(1) = famx
      Fam(2) = famy
      Fam(3) = famz

      if (ppiclf_debug==2) then
      if (ppiclf_nid .eq. 0 .or. ppiclf_np == 1) then
      if (i<=5 .and. iStage==1) then  
         open(unit=7051,file='fort.7051',position='append')
         open(unit=7052,file='fort.7052',position='append')
         open(unit=7053,file='fort.7053',position='append')
         open(unit=7054,file='fort.7054',position='append')
         open(unit=7055,file='fort.7055',position='append')
         write(7050+i,*) i, ppiclf_nid, ppiclf_np, ppiclf_time, ! 0-3
     >      ppiclf_rprop(PPICLF_R_JRHSR,i),        ! 4
     >      ppiclf_rprop(PPICLF_R_JPGCX,i),        ! 5
     >      ppiclf_rprop(PPICLF_R_JPGCY,i),        ! 6
     >      ppiclf_rprop(PPICLF_R_JPGCZ,i),        ! 7
     >      SDrho,                                 ! 8
     >      rhof, rphip, rmachp, SDrho,            ! 9-12
     >      ppiclf_rprop(PPICLF_R_WDOTX,i),        ! 13
     >      ppiclf_rprop(PPICLF_R_WDOTY,i),        ! 14
     >      ppiclf_rprop(PPICLF_R_WDOTZ,i),        ! 15
     >      Wdot_neighbor_mean(1:3), nneighbors,   ! 16-19
     >      famx, famy, famz,                      ! 20-22
     >      ppiclf_rprop(PPICLF_R_JVOLP,i)*Fam(1), ! 23
     >      ppiclf_rprop(PPICLF_R_JVOLP,i)*Fam(2), ! 24
     >      ppiclf_rprop(PPICLF_R_JVOLP,i)*Fam(3)  ! 25
         flush(7051)
         flush(7052)
         flush(7053)
         flush(7054)
         flush(7055)
      end if  
      end if  
      end if  

      return
      end
!
!-----------------------------------------------------------------------
!
!
! IMPORTANT NOTE: THIS SUBROUTINE IS CURRENTLY ONLY 
!    VALID FOR MONODISPERSE PARTICLE BEDS 
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!                          Binary part 
!-----------------------------------------------------------------------
!
!
      subroutine ppiclf_user_AM_Briney_Binary(i,iStage,
     >                 famx,famy,famz,rmass_add)
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
      integer i, j, k, l, n, jj
      integer*4 iStage
      real*8 rad
      real*8 dr_max
      real*8 IA, II
      real*8 famx, famy, famz, rmass_add
      real*8 alpha

! Declare functions
      real*8 IA_analytical, IA_numerical, II_analytical, II_numerical
!
! Code:
!
      ! particle radius
      rad = ppiclf_rprop(PPICLF_R_JDP,i) * 0.5d0
      
      ! ppiclf_d2chk(3) is neighbor width - user defined
      dr_max = ppiclf_d2chk(3)

      ! In the example program binary_model.f90, the nearest
      ! neighbor calculations are done here at this point to
      ! calculate wdot. However, these have been previously 
      ! calculated in ppiclf_user.f in a separate do loop
      ! over 1 <= i <= ppiclf_npart

      ! Model only valid for local volume fraction
      ! less than 0.4, so we limit it here without
      ! over riding rphip
      alpha = min(0.4, rphip) ! limit alpha to mitigate misuse

      ! we need to make the numerical integration more efficient
      ! do it in a table and save it 
      if (dr_max / rad >= 3.49) then
         IA = IA_analytical(dr_max, rad, alpha) ! self acceleration
         II = II_analytical(dr_max, rad, alpha) ! neighbor acceleration (induced)
      else
         if (ppiclf_nid==0 .and. iStage==1) then
             print*, "***WARNING*** - NUMERICAL 
     >               FUNCTIONS USED IN ADDED MASS"
         endif
         IA = IA_numerical(dr_max, rad, alpha)
         II = II_numerical(dr_max, rad, alpha)
      end if

      do j=1,3
         jj = PPICLF_R_WDOTX + (j-1)
         Fam(j) = Fam(j) + IA*ppiclf_rprop(jj, i) ! added mass
         Fam(j) = Fam(j) + II*Wdot_neighbor_mean(j) 
     >                               / nneighbors ! induced added mass
      end do

      ! multiply by volume before adding unary term
      ! doing so here implies that the particle volume is
      ! the same for all particles; i.e., monodisperse packs
      do j=1,3
         Fam(j) = ppiclf_rprop(PPICLF_R_JVOLP,i)*Fam(j) 
      end do

      famx = Fam(1)
      famy = Fam(2)
      famz = Fam(3)
        
      if (ppiclf_debug==2) then
      if (ppiclf_nid .eq. 0 .or. ppiclf_np == 1) then
      if (i<=5 .and. iStage==1) then
         open(unit=7061,file='fort.7061',position='append')
         open(unit=7062,file='fort.7062',position='append')
         open(unit=7063,file='fort.7063',position='append')
         open(unit=7064,file='fort.7064',position='append')
         open(unit=7065,file='fort.7065',position='append')
         write(7060+i,*) i,iStage, ppiclf_time,     ! 0-2
     >      rphip, rmachp,                          ! 3-4
     >      IA, II,                                 ! 5-6
     >      ppiclf_rprop(PPICLF_R_WDOTX,i),         ! 7
     >      ppiclf_rprop(PPICLF_R_WDOTY,i),         ! 8
     >      ppiclf_rprop(PPICLF_R_WDOTZ,i),         ! 9
     >      Wdot_neighbor_mean(1:3), nneighbors,    ! 10-13
     >      famx, famy, famz,                       ! 14-16
     >      Fam(1), Fam(2), Fam(3)                  ! 17-19
         flush(7061)
         flush(7062)
         flush(7063)
         flush(7064)
         flush(7065)
      end if
      end if
      end if
        

      return
      end
!------------------------------------------------------------------------
!
! Created May 20, 2024
!
! Subroutine for added mass
!   also called the quasi-unsteady force,
!   or the inviscid unsteady force in case of the Euler equations
!
! Contains functions:
!    real*8 function B11_11(d, alpha)
!    real*8 function B11_22(d, alpha)
!    real*8 function B12_11(d, alpha)
!    real*8 function B12_22(d, alpha)
!    real*8 function IA_analytical(rmax, rad, alpha)
!    real*8 function II_analytical(rmax, rad, alpha)
!    real*8 function rdf_analytical(r, alpha)
!    real*8 function IA_numerical_integrand(d, rad, alpha)
!    real*8 function II_numerical_integrand(d, rad, alpha)
!    real*8 function IA_numerical(rmax, rad, alpha)
!    real*8 function II_numerical(rmax, rad, alpha)
!    
! Implementing Added Mass Algorithm from S.Briney (2024)
!  
! n       = number of points
! alpha   = volume fraction
! rad     = particle radius
! d       = center-to-center distance
! rmax    = center-to-center max neighbor distance
! R       = resistance matrix (output)
! x       = x_2 - x_1
! y       = y_2 - y_1
! z       = z_2 - z_1
! dr_max  = max interaction distance between particles considered 
! poins   = 3xn array of points x, y, z. Initialized as points(3,n)
!
! correction_analytical_always 
!    = if true, always use the analytical distant neighbor correction
!    > if false, use numerical when dr_max/rad < 3.49
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Functions for calculating the four curves that define the binary model
! These the the vol fraction-corrected binary added mass terms
!   taken from Appendix C of Briney et at (2024).
!
! parallel added mass
      real*8 function B11_11(d, alpha)
      implicit none
       
      ! input
      real*8 d, alpha
      
      ! local vars
      ! asymptotic solution (Beguin et al. 2016)
      ! empirical higher order terms
      ! empirical volume fraction correction
      real*8 asym, hot, alpha_corr
      
      ! Beguin et al. 2016
      asym = 0.5*(3.0/64.0 * (2.0/d)**6 + 9.0/256.0*(2.0/d)**8) 

      hot = 0.059/((d-1.9098)*(d-0.4782)**2 * (d**3))
      alpha_corr = (alpha*alpha - 0.0902*alpha) 
     >                 * 70.7731/((d+2.0936)**6)
      
      B11_11 = asym + hot + alpha_corr
      
      return
      end

! perpendicular added mass
      real*8 function B11_22(d, alpha)
      implicit none
      
      ! input
      real*8 d, alpha
      
      ! local vars
      ! hot combined with alpha correction in this case
      ! See comment in calc_B11_11
      real*8 asym, hot
      real*8 A, B
          
      ! Beguin et al. 2016
      asym = 0.5*(3.0/256.0 * (2.0/d)**6 + 3.0/256.0 
     >               *(2.0/d)**8) 
      
      A = 0.0003 + 0.0262*alpha*alpha - 0.0012*alpha
      B = 1.3127 + 1.0401*alpha*alpha - 1.2519*alpha
      hot = A / ((d-B)**6)
      
      B11_22 = asym + hot
      
      return
      end
      
! parallel induced added mass
      real*8 function B12_11(d, alpha)
      implicit none
         
      ! input
      real*8 d, alpha
      
      ! local vars
      real*8 asym, hot, alpha_corr
      
      ! Beguin et al. 2016
      asym = 0.5*(-3.0/8.0*(2.0/d)**3 - 3.0/512.0
     >               *(2.0/d)**9) 

      hot = -0.0006/((d-1.5428)**5)
      alpha_corr = -0.7913*alpha
     >               *exp(-(0.9801 - 0.1075*alpha)*d)
      
      B12_11 = asym + hot + alpha_corr

      return
      end
      
! perpendicular induced added mass
      real*8 function B12_22(d, alpha)
      implicit none
          
      ! input
      real*8 d, alpha
      
      ! local vars
      real*8 asym, hot
          
      ! Beguin et al. 2016
      asym = 0.5*(3.0/16.0 * (2.0/d)**3 + 3.0/4096.0 * (2.0/d)**9) 
          
      ! includes alpha correction
      hot = (-0.1985 + 16.7372*alpha)/(d**6)
     >           + (1.3907 - 48.2604*alpha)/(d**8) 
      
      B12_22 = asym + hot
      
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Analytical correction functions for distant neighbors. 
!     Assumes g(r) = 1. Good for maxr >= 3.5.
!      
!     These can be found in Appendix D of Briney etal (2024).
!
!     added mass
!     rmax = center-to-center maximum neighbor distance
!     rad = particle radius
!     alpha = volume fraction
      real*8 function IA_analytical(rmax, rad, alpha)
      implicit none
      
      ! input
      real*8 rmax, rad, alpha
      
      ! local vars
      real*8 term1, term2, term3, c, numerator, B11, B22, maxr
      
      maxr = rmax / rad
      
      ! B11
      term1 = (5.0*maxr**2 + 9.0)/(10.0*maxr**5)
          
      term2 = ((0.00720826 - 0.0150737 * maxr) * log(maxr - 1.9098)
     >         - 0.120023 * maxr * log(maxr - 0.4782) 
     >         + 0.057395 * log(maxr - 0.4782) 
     >         + (0.135097 * maxr - 0.0646033) * log(maxr) - 0.0861828)
     >                                           /(maxr - 0.4782)
      
      term3 = (alpha**2 - 0.0902*alpha) 
     >     * (0.333333 * maxr**2 + 0.348933 * maxr + 0.146105)
     >     /(maxr + 2.0936)**5
      
      B11 = term1 + term2 + term3
      
      ! B22
      term1 = (5.0*maxr**2 + 12.0)/(40.0*maxr**5)
      
      numerator = 0.0262*alpha**2 - 0.0012*alpha + 0.0003
      c = -1.0401*alpha**2 + 1.2519*alpha - 1.3127
      term2 = numerator * (10.0*maxr**2 + 5.0*maxr*c + c**2) 
     >             / (30.0*(maxr+c)**5)
      
      B22 = term1 + term2
      
      IA_analytical = alpha*(B11 + 2.0*B22)
          
!     write(1,*) IA_analytical, rmax, rad, alpha,
!    >                maxr, term1, term2, term3, B11, B22

      return
      end

      ! induced added mass
      real*8 function II_analytical(rmax, rad, alpha)
      implicit none
      
      ! input
      real*8 rmax, rad, alpha
      
      ! local vars
      real*8 term1, term2, term3, B11, B22, c, maxr, prefactor
      
      maxr = rmax / rad ! scale the filter width
      
      ! B11
      term1 = -1.0/(4.0*maxr**6) 
      
      term2 = -6.0e-4 * (0.5*maxr**2 - 0.516067*maxr + 0.199744)
     >             /(1.5482 - maxr)**4
      
      prefactor = -0.7913*alpha
      c = 0.9801 - 0.1075*alpha
      term3 = prefactor * (maxr *c *(maxr *c + 2.0) + 2.0) 
     >                       * exp((-maxr *c))/c**3
      
      B11 = term1 + term2 + term3
      
      ! B22
      term1 = 1.0 / (32.0 * maxr**6)
      
      term2 = (-0.1985 + 16.7372*alpha) / (3.0*maxr**3)
      
      term3 = (1.3907 - 48.2604*alpha) / (5.0*maxr**5)
      
      B22 = term1 + term2 + term3
      
      II_analytical = alpha*(B11 + 2.0*B22)
          
!          write(2,*) II_analytical, rmax, rad, alpha,
!     >                maxr, term1, term2, term3, B11, B22
      
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     Numerical distant neighbor corrections
!      
!     Uses the rdf function
!      
!     r = dist between particles / diameter 
!          - normalization different than my convention!
!     phi = particle volume fraction
!     Trokhymchuk et al. 2005. The Journal of Chemical Physics.
!     https://doi.org/10.1063/1.1979488
!     erratum: https://doi.org/10.1063/1.2188941

      real*8 function rdf_analytical(r, alpha)
      implicit none
      
      ! inputs
      real*8 r, alpha
      
      ! local vars
      real*8  eta, d, mu, alpha0, beta0, denom_gamma, 
     >              gamma, omega, kappa
      real*8  alpha1, beta, rstar, gm, gsig, Brad, Arad, delta, 
     >               Crad, g2
      
      eta = alpha
      
      d   = (2.0*eta*(eta*eta - 3.0*eta - 3.0
     >     + sqrt(3.0*(eta**4-2.0*(eta**3)
     >     +eta**2+6.0*eta+3.0))))**(1.0/3.0)
      
      !mu    = (2.0*eta/(1.0-eta))*(-1.0-(d/(2.0*eta))-(eta/d))
      mu    = (2.0*eta/(1.0-eta))*(-1.0-(d/(2.0*eta))+(eta/d)) ! see erratum
      
      alpha0 = (2.0*eta/(1.0-eta))
     >         *(-1.0+(d/(4.0*eta))-(eta/(2.0*d)))
      
      beta0 = (2.0*eta/(1.0-eta))*sqrt(3.0)
     >         *(-(d/(4.0*eta))-(eta/(2.0*d)))
      
      denom_gamma = (alpha0**2 + beta0**2 - 2.0*mu*alpha0)
     >         *(1.0 + 0.5*eta) - mu*(1.0 + 2.0*eta) ! erratum
      
      gamma = atan(-(1.0/beta0)*((alpha0*(alpha0**2+beta0**2)
     >         -mu*(alpha0**2-beta0**2))*(1.0+0.5*eta)
     >         +(alpha0**2+beta0**2-mu*alpha0)*(1.0+2.0*eta))
     >         /(denom_gamma)) ! erratum
      
      !gamma = atan(-(1.0/beta0)*((alpha0*(alpha0**2+beta0**2)
      ! > -mu*(alpha0**2+beta0**2))*(1.0+0.5*eta)+
      ! > (alpha0**2+beta0**2-mu*alpha0)*(1.0+2.0*eta)));
      
      omega = -0.682*exp(-24.697*eta)+4.72+4.45*eta
      
      kappa = 4.674*exp(-3.935*eta)+3.536*exp(-56.27*eta)
      
      alpha1 = 44.554+79.868*eta+116.432*eta*eta-44.652*exp(2.0*eta)
      
      beta  = -5.022+5.857*eta+5.089*exp(-4.0*eta)
      
      rstar = 2.0116-1.0647*eta+0.0538*eta*eta
      
      gm    = 1.0286-0.6095*eta+3.5781*(eta**2)-21.3651*(eta**3)
     >             +42.6344*(eta**4)-33.8485*(eta**5)
      
      gsig  = (1.0/(4.0*eta))*(((1.0+eta+(eta**2)-(2.0/3.0)*(eta**3)
     >             -(2.0/3.0)*(eta**4))/((1.0-eta)**3))-1.0)
      
      Brad = (gm-(gsig/rstar)*exp(mu*(rstar - 1.0)))
     >           *rstar/(cos(beta*(rstar-1.0)+gamma)
     >           *exp(alpha1*(rstar-1.0))
     >           -cos(gamma)*exp(mu*(rstar-1.0)))
      
      Arad = gsig - Brad*cos(gamma)
      
      delta = -omega*rstar - atan((kappa*rstar+1.0)/(omega*rstar))
      
      Crad = rstar*(gm-1.0)*exp(kappa*rstar)/cos(omega*rstar+delta)
      
      g2 = (Arad/r)*exp(mu*(r-1.0)) + (Brad/r)
     >         *cos(beta*(r-1.0)+gamma)*exp(alpha1*(r-1.0))
      
      if (r > rstar) then
         g2 = 1.0 +(Crad/r)*cos(omega*r+delta)*exp(-kappa*r)
      else if (r < 1) then
         g2 = 0.0
      end if
      
      rdf_analytical = g2
      
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     d = center to center distance
!     rad = radius
!     alpha = volume fraction

      real*8 function IA_numerical_integrand(d, rad, alpha)
      implicit none
      
      ! input
      real*8 d, rad, alpha
      
      ! local vars
      real*8 g, f11, f22, rho, pi, w

      ! declare funcitons used
      real*8 rdf_analytical, B11_11, B11_22
      
      pi = 4.0*ATAN(1.0)
      
      g = rdf_analytical((d/(2.0*rad)), alpha)
      
      f11 = B11_11(d/rad, alpha)
      f22 = B11_22(d/rad, alpha)
      
      rho = alpha / (4.0/3.0*pi) ! number density
      w = 1.0/3.0 * rho * 4.0*pi*(d/rad)**2
      
      IA_numerical_integrand = w*g*(f11 + 2.0*f22)

      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     d = center to center distance
!     rad = radius
!     alpha = volume fraction

      real*8 function II_numerical_integrand(d, rad, alpha)
      implicit none
      
      ! input
      real*8 d, rad, alpha
      
      ! local vars
      real*8 g, f11, f22, rho, pi, w
          
      ! declare funcitons used
      real*8 rdf_analytical, B12_11, B12_22
      
      pi = 4.0*ATAN(1.0)
      
      g = rdf_analytical((d/(2.0*rad)), alpha)
      
      f11 = B12_11(d/rad, alpha)
      f22 = B12_22(d/rad, alpha)
      
      rho = alpha / (4.0/3.0*pi) ! number density
      w = 1.0/3.0 * rho * 4.0*pi*(d/rad)**2
      
      II_numerical_integrand = w*g*(f11 + 2.0*f22)

      return
      end
!
!-----------------------------------------------------------------------
!
      
      real*8 function IA_numerical(rmax, rad, alpha)
      implicit none
      
      ! input
      real*8 rmax, rad, alpha
      
      ! local vars
      real*8 maxr, dr, r, integral, f, coef(2)
      integer npts, i, j
          
      ! declare funcitons used
      real*8 IA_numerical_integrand
      
      npts = 1001 ! must be odd
      
      maxr = rad*20.0 ! essentially infinity
      
      coef(1) = 4.0
      coef(2) = 2.0
      
      if (maxr > rmax) then
         dr = (maxr - rmax) / (npts + 1)
     
         integral = 0.0
         r = rmax
      
         ! Simpson's rule
         ! endpoints first
         f = IA_numerical_integrand(r, rad, alpha)
         integral = integral + f
      
         f = IA_numerical_integrand(maxr, rad, alpha)
         integral = integral + f
      
         ! middle points
         do i=2,npts-1
            r = r + dr
            f = IA_numerical_integrand(r, rad, alpha)
      
            j = mod(i, 2) + 1
            integral = integral + coef(j)*f ! 4, 2, 4, 2, ...
         end do
      
         IA_numerical = dr * integral / 3.0
      
      else
         IA_numerical = 0.0
      end if
      
      return
      end
      
!
!-----------------------------------------------------------------------
!
      
      real*8 function II_numerical(rmax, rad, alpha)
      implicit none
      
      ! input
      real*8 rmax, rad, alpha
      
      ! local vars
      real*8 maxr, dr, r, integral, f, coef(2)
      integer npts, i, j
      
      ! declare funcitons used
      real*8 II_numerical_integrand
          
      npts = 1001 ! must be odd
      
      maxr = rad*20.0 ! essentially infinity
      
      coef(1) = 4.0
      coef(2) = 2.0
      
      if (maxr > rmax) then
         dr = (maxr - rmax) / (npts + 1)
      
         integral = 0.0
         r = rmax
      
         ! Simpson's rule
         ! endpoints first
         f = II_numerical_integrand(r, rad, alpha)
         integral = integral + f
      
         f = II_numerical_integrand(maxr, rad, alpha)
         integral = integral + f
      
         ! middle points
         do i=2,npts-1
            r = r + dr
            f = II_numerical_integrand(r, rad, alpha)
      
            j = mod(i, 2) + 1
            integral = integral + coef(j)*f ! 4, 2, 4, 2, ...
         end do
      
         II_numerical = dr * integral / 3.0
      
      else
         II_numerical = 0.0
      end if
      
      return 
      end
!------------------------------------------------------------------------
!
! Created May 20, 2024
!
! Subroutine for added mass
!   also called the quasi-unsteady force,
!   or the inviscid unsteady force in case of the Euler equations
!
! Contains subroutines:
!     rotation_matrix(x, y, z, Q)
!     resistance_pair(x, y, z, alpha, rad, R)
!
! Implementing Added Mass Algorithm from S.Briney (2024)
!  
! n       = number of points
! alpha   = volume fraction
! rad     = particle radius
! d       = center-to-center distance
! rmax    = center-to-center max neighbor distance
! R       = resistance matrix (output)
! x       = x_2 - x_1
! y       = y_2 - y_1
! z       = z_2 - z_1
! dr_max  = max interaction distance between particles considered 
! poins   = 3xn array of points x, y, z. Initialized as points(3,n)
!
! correction_analytical_always 
!    = if true, always use the analytical distant neighbor correction
!    > if false, use numerical when dr_max/rad < 3.49
!
!
!-----------------------------------------------------------------------
!
! Calculates the rotation matrix Q to align the coordinate 
!    system such that the point (x, y, z) is on the x-axis

      subroutine rotation_matrix(x, y, z, Q)
       
      ! inputs
      real*8 x, y, z
       
      ! output rotation matrix
      real*8 Q(3, 3)
       
      ! local vars
      real*8 gamma, beta
           
      gamma =  atan2(y, x)
      beta  = -atan2(z, sqrt(x*x + y*y))
       
      Q(1, 1) = cos(beta)*cos(gamma)
      Q(1, 2) = cos(beta)*sin(gamma)
      Q(1, 3) = -sin(beta)
       
      Q(2, 1) = -sin(gamma)
      Q(2, 2) = cos(gamma)
      Q(2, 3) = 0.0
       
      Q(3, 1) = sin(beta)*cos(gamma)
      Q(3, 2) = sin(beta)*sin(gamma)
      Q(3, 3) = cos(beta)
       
      return
      end
!
!-----------------------------------------------------------------------
!
      ! This is the one of the functions the user should typically call.
      ! Returns the resistance matrix for a 2 particle system 
      !    of arbitrary orientation
      ! x = x2 - x1, y = y2 - y1, z = z2 - z1 (relative position)
      ! rad = particle radius
      ! R = resistance matrix (output)

      subroutine resistance_pair(x, y, z, alpha, rad, R)
      implicit none
       
      ! input: x, y, z of second particle relative to the second particle
      ! x = x2 - x1, etc.
      ! rad = particle radius (monodisperse)
      real*8 x, y, z, alpha, rad
       
      ! output: resistance matrix
      real*8 R(6, 6)
       
      ! local vars
      real*8, dimension(3, 3) :: Q, Qt, B11, B12
       
      real*8 dist
       
      ! declare functions
      real*8 B11_11, B11_22, B12_11, B12_22
       
      ! get rotation matrix Q
      call rotation_matrix(x, y, z, Q)
       
      ! normalize the distance by the particle radius
      dist = sqrt(x*x + y*y + z*z) / rad 
       
      ! initially set coefficients to 0
      B11(1:3, 1:3) = 0.0
      B12(1:3, 1:3) = 0.0
      
      ! self acceleration
      B11(1, 1) = B11_11(dist, alpha)
      B11(2, 2) = B11_22(dist, alpha)
      B11(3, 3) = B11(2, 2)
       
      ! neighbor acceleration
      B12(1, 1) = B12_11(dist, alpha)
      B12(2, 2) = B12_22(dist, alpha)
      B12(3, 3) = B12(2, 2)
       
      ! rotate
      Qt  = transpose(Q)
      B11 = matmul(Qt, matmul(B11, Q))
      B12 = matmul(Qt, matmul(B12, Q))
            
      ! store resitance matrix
      R(1:3, 1:3) = B11
      R(1:3, 4:6) = B12
       
      R(4:6, 1:3) = B12
      R(4:6, 4:6) = B11
          
      ! write(7059,*) x, y, z, alpha, rad, dist 
      !
      ! write(7060,*) B11(1,1:3), B11(2,1:3), B11(3,1:3)
      ! write(7061,*) B12(1,1:3), B12(2,1:3), B12(3,1:3)
      ! write(7062,*) C11(1,1:3), C11(2,1:3), C11(3,1:3)
      ! write(7063,*) C12(1,1:3), C12(2,1:3), C12(3,1:3)
      ! write(7064,*) Qt(1,1:3), Qt(2,1:3), Qt(3,1:3)
      !
      ! write(7069,*) R(1,1:6)
      ! write(7069,*) R(2,1:6)
      ! write(7069,*) R(3,1:6)
      ! write(7069,*) R(4,1:6)
      ! write(7069,*) R(5,1:6)
      ! write(7069,*) R(6,1:6)
      ! write(7069,*) "------------------------------------"

      return
      end
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer models for non-burning
!    particles
! May update if we include particle combustion
!
! if heattransfer_flag = 0  ignore heat transfer
!                      = 1  Stokes
!                      = 2  Ranz-Marshall (1952)
!                      = 3  Gunn (1977)
!                      = 4  Fox (1978)
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_HT_driver(i,qq,rmass_therm)
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

      integer*4 i
      real*8 qq, rmass_therm

!
! Code:
!
      if (heattransfer_flag == 1) then
         call HTModel_Stokes(i,qq,rmass_therm)
      elseif (heattransfer_flag == 2) then
         call HTModel_RM(i,qq,rmass_therm)
      elseif (heattransfer_flag == 3) then
         call HTModel_Gunn(i,qq,rmass_therm)
      elseif (heattransfer_flag == 4) then
         call HTModel_Fox(i,qq,rmass_therm)
      else
         call ppiclf_exittr('Unknown heat transfer model$', 0.0d0, 0)
      endif


      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Stokes limit
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Stokes(i,qq,rmass_therm)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!

      integer*4 i
      real*8 qq, rmass_therm
      real*8 OneThird
      real*8 Nuss

!
! Code:
!
      OneThird = 1.0d0/3.0d0

      rmass_therm = rmass*rcp_part

      qq = rpi*rkappa*dp*(ppiclf_rprop(PPICLF_R_JT,i) -
     >                          ppiclf_y(PPICLF_JT,i) )

      ! define Nusselt number
      Nuss = 2.0d0

      qq = qq*Nuss

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Ranz-Marshall
!    define Nusselt number Nu = Nu(Pr,Re)
!
! (1) Ling et al (2012) 
!     "Interaction of a planar shock wave with a dense particle
!     curtain: Modeling and experiments"
!     Physics of Fluids, Vol. 24, 113301
! (2) Durant et al (2022) 
!     "Explosive dispersal of particles in high speed environments"
!     Journal of Applied Physics, Vol. 132, 184902
!
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_RM(i,qq,rmass_therm)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!

      integer*4 i
      real*8 qq, rmass_therm
      real*8 OneThird
      real*8 Nuss

!
! Code:
!
      OneThird = 1.0d0/3.0d0

      rmass_therm = rmass*rcp_part

      qq = rpi*rkappa*dp*(ppiclf_rprop(PPICLF_R_JT,i) -
     >                          ppiclf_y(PPICLF_JT,i) )

      ! define Nusselt number Nu = Nu(Pr,Re)
      Nuss = 2.0d0+0.6d0*(rep**0.5d0)*(rpr**OneThird)

      qq = qq*Nuss

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Gunn
!    define Nusselt number Nu = Nu(Pr,Re,phi)
!
! (1)  Gunn (1978)
!        "Transfer of heat or mass to particles in
!        fixed and fluidised beds"
!        Int. J. Heat Mass Transfer, Vol. 21, pp. 467-476
! (2)  Houim and Oran (2016)
!        "A multiphase model for compressible granular-gaseous
!        flows: formulation and initial tests"
!        J. Fluid Mech., Vol. 789, pp. 166-220
! (3) Boniou and Fox (2023)
!        "Shock-particle-curtain-interaction study with
!        a hyperbolic two-fluid model: Effect of particle
!        force models"
!        Int. Journal of Mutiphase Flow, Vol. 169, 104591
!
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Gunn(i,qq,rmass_therm)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!

      integer*4 i
      real*8 qq, rmass_therm
      real*8 OneThird
      real*8 vg
      real*8 Nuss

!
! Code:
!
      OneThird = 1.0d0/3.0d0

      rmass_therm = rmass*rcp_part

      qq = rpi*rkappa*dp*(ppiclf_rprop(PPICLF_R_JT,i) -
     >                          ppiclf_y(PPICLF_JT,i) )

      ! define bed voidage = ratio of free volume avaliable
      ! for flow to the total volume of bed; aka, volume
      ! fraction of the gas phase
      ! vg = 1.0 - rphip == rphif
      vg = rphif

      ! define Nusselt number Nu = Nu(Pr,Re,phi)
      Nuss = (7.0d0-10.0*vg+5.0*vg*vg)
     >           *(1.0+0.7d0*(rep**0.2d0)*(rpr**OneThird))
     >     + (1.33d0-2.4d0*vg+1.2d0*vg*vg)*(rep**0.7d0)*(rpr**OneThird)

      qq = qq*Nuss

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for quasi-steady heat transfer
! Quasi-steady heat transfer: Fox
!    define Nusselt number Nu = Nu(Pr,Re,M)
!
! (1)  Fox, Rackett, and Nicholls (1978)
!        "Shock wave ignition of magnesium powders"
!        in Proc. 11th Int. Symp. Shock Tubes and Waves.
!        University of Washington Press, Seattle, WA, pp. 262-268
! (2)  Ling, Haselbacher, and Balachandar (2011)
!        "Importance of unsteady contributions to force and heating
!        for particles in compressible flows. Part I: Modeling and
!        anaylsis for shock-particle interaction"
!        Int. Journal of Multiphase Flow, Vol. 37, pp. 1026-1044
!
!
!-----------------------------------------------------------------------
!
      subroutine HTModel_Fox(i,qq,rmass_therm)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!

      integer*4 i
      real*8 qq, rmass_therm
      real*8 OneThird
      real*8 Nuss

!
! Code:
!
      OneThird = 1.0d0/3.0d0

      rmass_therm = rmass*rcp_part

      qq = rpi*rkappa*dp*(ppiclf_rprop(PPICLF_R_JT,i) -
     >                          ppiclf_y(PPICLF_JT,i) )

      ! define Nusselt number Nu = Nu(Pr,Re,M)
      Nuss = 2.0d0*exp(-rmachp)/(1.0d0+17.0d0*rmachp/rep)
     >     + 0.495d0*(rpr**OneThird)*(rep**0.55d0)*
     >       ((1.0d0+0.5d0*exp(-17.0d0*rmachp/rep))/1.5d0)

      qq = qq*Nuss

      return
      end
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for computing the torque terms on the
!    RHS of the angular velocity equations
!
!
! if collisional_flag = 1  F = Fn
!                     = 2  F = Fn + Ft + Tt
!                     = 3  F = Fn + Ft + Tt + Th + Tr
! where Tt = collisional torque
!       Th = hydrodynamic torque
!       Tr = rolling torque
!
! Note that the collisional and rolling torques are due to
!     particle-particle interactions and thus are evaulated
!     in ppiclf_user_EvalNearestNeighbor. Therefore, only the
!     hydrodynamic torque is left to be calculated.
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_Torque_driver(i,iStage,taux,tauy,tauz)
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

      integer*4 i,iStage
      real*8 taux, tauy, tauz
      real*8 taux_hydro, tauy_hydro, tauz_hydro
      real*8 rmass_local

!
! Code:
!
      taux_hydro = 0.0d0
      tauy_hydro = 0.0d0
      tauz_hydro = 0.0d0

      if (collisional_flag == 3) then
         call Torque_Hydro(i,taux_hydro,tauy_hydro,tauz_hydro)
      endif

      taux = taux + taux_hydro
      tauy = tauy + tauy_hydro
      tauz = tauz + tauz_hydro

      return
      end
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Created June 17, 2024
!
! Subroutine for hydrodynamic torque
!
!
!-----------------------------------------------------------------------
!
      subroutine Torque_Hydro(i,taux_hydro,tauy_hydro,tauz_hydro)
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
      real*8 taux_hydro, tauy_hydro, tauz_hydro
      real*8 omgrx, omgry, omgrz, omgr_mag
      real*8 Ct1, Ct2, Ct3, Ct
      real*8 reyr, beta, rIp, factor

!
! Code:
!
      ! Compute relative angular velocity components
      !    and magnitude
      omgrx = 0.5d0*ppiclf_rprop(PPICLF_R_JXVOR,i)
     >          - ppiclf_y(PPICLF_JOX,i)
      omgry = 0.5d0*ppiclf_rprop(PPICLF_R_JYVOR,i)
     >          - ppiclf_y(PPICLF_JOY,i)
      omgrz = 0.5d0*ppiclf_rprop(PPICLF_R_JZVOR,i)
     >          - ppiclf_y(PPICLF_JOZ,i)
      omgr_mag = sqrt(omgrx*omgrx + omgry*omgry + omgrz*omgrz)

      ! Particle rotational Reynolds number
      reyr = rhof*dp*dp*omgr_mag/(4.0d0*rmu)
      reyr = max(reyr,0.001d0)

      ! Compute the hydrodynamic torque parameter Ct=Ct(Re_r)
      if (reyr < 1) then
         Ct1 = 0.0d0
         Ct2 = 16.0d0*rpi
         Ct3 = 0.0d0
      elseif (reyr < 10) then
         Ct1 = 0.0d0
         Ct2 = 16.0d0*rpi
         Ct3 = 0.0418d0
      elseif (reyr < 20) then
         Ct1 = 5.32d0
         Ct2 = 37.2d0
         Ct3 = 0.0d0
      elseif (reyr < 50) then
         Ct1 = 6.44d0
         Ct2 = 32.2d0
         Ct3 = 0.0d0
      elseif (reyr < 100) then
         Ct1 = 6.45d0
         Ct2 = 32.1d0
         Ct3 = 0.0d0
      else
         call ppiclf_exittr('Re rotational too large$', reyr, 0)
      endif

      Ct = Ct1/sqrt(reyr) + Ct2/Reyr + Ct3*reyr

      ! Now compute hydrodynamic torque components
      beta = rhop/rhof
      rIp  = rmass*dp*dp/10.0d0
      factor = rIp*60.0d0*Ct*omgr_mag/(64.0d0*rpi*beta)

      taux_hydro = factor*omgrx
      tauy_hydro = factor*omgry
      tauz_hydro = factor*omgrz


      return
      end
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
!-----------------------------------------------------------------------
      subroutine ppiclf_user_EvalNearestNeighbor
     >                                        (i,j,yi,rpropi,yj,rpropj)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 i
      integer*4 j
      real*8 yi    (PPICLF_LRS)    
      real*8 rpropi(PPICLF_LRP)
      real*8 yj    (PPICLF_LRS)    
      real*8 rpropj(PPICLF_LRP)
!
! Internal:
!
      real*8 ksp,erest
      !common /ucollision/ ksp,erest

      real*8 rpi2, rthresh, rxdiff, rydiff, rzdiff, rdiff, rm1, rm2,
     >       rmult, eta, rbot, rn_12x, rn_12y, rn_12z, rdelta12,
     >       rv12_mag, rv12_mage, rksp_max, rnmag, rksp_wall, rextra,
     >       JDP_i,JDP_j   
!
      rpi2  =  9.869604401089358d0
      ksp = 100.0 !7.09E+3 !3.05E+7 !500000 !40000 !21000
      erest = 0.7
      ! other particles
      if (j .ne. 0) then
         !Added spload and radius factor
!         rthresh  = 0.5d0*(rpropi(PPICLF_R_JDP) + rpropj(PPICLF_R_JDP))

!        JDP_i =(rpropi(PPICLF_R_JDPi)/rpropi(PPICLF_R_JSPT))**(1.0/3.0)
!        JDP_i = JDP_i * rpropi(PPICLF_R_JDP)  
!        JDP_j =(rpropj(PPICLF_R_JDPi)/rpropj(PPICLF_R_JSPT))**(1.0/3.0)
!        JDP_j = JDP_j * rpropj(PPICLF_R_JDP)
!        rthresh  = 0.5d0*(JDP_i + JDP_j)
         rthresh  = 0.5d0*(rpropi(PPICLF_R_JDP) + rpropj(PPICLF_R_JDP))

         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = rxdiff**2 + rydiff**2
     >          +rzdiff**2
         rdiff = sqrt(rdiff)
         
         if (rdiff .gt. rthresh) return
         
         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         rm2 = rpropj(PPICLF_R_JRHOP)*rpropj(PPICLF_R_JVOLP)
         
         rmult = 1.0d0/sqrt(1.0d0/rm1+1.0d0/rm2)
         eta   = 2.0d0*sqrt(ksp)*log(erest)/sqrt(log(erest)**2+rpi2)
     >           *rmult
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = (yj(PPICLF_JVX)-yi(PPICLF_JVX))*rn_12x
     >            + (yj(PPICLF_JVY)-yi(PPICLF_JVY))*rn_12y
     >            + (yj(PPICLF_JVZ)-yi(PPICLF_JVZ))*rn_12z

         rv12_mage = rv12_mag*eta
         rksp_max  = ksp*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
     >                              + rnmag*rn_12y
         ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
     >                              + rnmag*rn_12z

      ! boundaries
      elseif (j .eq. 0) then

         rksp_wall = ksp!10000

         ! give a bit larger collision threshold for walls
         rextra   = 0.0d0
         ! add sploading and radius factor 
         rthresh  = (0.5d0+rextra)*rpropi(PPICLF_R_JDP)
         
         rxdiff = yj(PPICLF_JX) - yi(PPICLF_JX)
         rydiff = yj(PPICLF_JY) - yi(PPICLF_JY)
         rzdiff = yj(PPICLF_JZ) - yi(PPICLF_JZ)
         
         rdiff = rxdiff**2 + rydiff**2
     >          +rzdiff**2
         rdiff = sqrt(rdiff)
         
         if (rdiff .gt. rthresh) return

         rm1 = rpropi(PPICLF_R_JRHOP)*rpropi(PPICLF_R_JVOLP)
         
         rmult = sqrt(rm1)
         eta   = 2.0d0*sqrt(rksp_wall)*log(erest)
     >           /sqrt(log(erest)**2+rpi2)*rmult
         
         rbot = 1.0d0/rdiff
         rn_12x = rxdiff*rbot
         rn_12y = rydiff*rbot
         rn_12z = rzdiff*rbot
         
         rdelta12 = rthresh - rdiff
         
         rv12_mag = -yi(PPICLF_JVX)*rn_12x
     >              -yi(PPICLF_JVY)*rn_12y
     >              -yi(PPICLF_JVZ)*rn_12z

         rv12_mage = rv12_mag*eta
         rksp_max  = rksp_wall*rdelta12
         rnmag     = -rksp_max - rv12_mage
         
         ppiclf_ydotc(PPICLF_JVX,i) = ppiclf_ydotc(PPICLF_JVX,i)
     >                              + rnmag*rn_12x
         ppiclf_ydotc(PPICLF_JVY,i) = ppiclf_ydotc(PPICLF_JVY,i)
     >                              + rnmag*rn_12y
         ppiclf_ydotc(PPICLF_JVZ,i) = ppiclf_ydotc(PPICLF_JVZ,i)
     >                              + rnmag*rn_12z
        
       !write(*,*) "Wall NEAR",i,ppiclf_ydotc(PPICLF_JVY,i)  
      endif

      return
      end
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for projection map
!
!
! This subroutine maps certain ppiclf terms onto the fluid mesh
!
!    This routine is needed for
!       (i)  feedback_flag = 1 (rocpicl/PICL_TEMP_Runge.F90)
!       (ii) # PROBE is in the *.inp file (libfloflu/WriteProbe.F90)

!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_MapProjPart(map,y,ydot,ydotc,rprop)
!
      implicit none
      include "PPICLF"
!
! Input:
!
      real*8 y    (PPICLF_LRS)
      real*8 ydot (PPICLF_LRS)
      real*8 ydotc(PPICLF_LRS)
      real*8 rprop(PPICLF_LRP)
!
! Output:
!
      real*8 map  (PPICLF_LRP_PRO)
!
! Code:
!
!
      !Add here sploading factor         
      map(PPICLF_P_JPHIP) = rprop(PPICLF_R_JVOLP)*rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JFX)   = ydotc(PPICLF_JVX) * rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JFY)   = ydotc(PPICLF_JVY) * rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JFZ)   = ydotc(PPICLF_JVZ) * rprop(PPICLF_R_JSPL)
      map(PPICLF_P_JE)    = ydotc(PPICLF_JT)  * rprop(PPICLF_R_JSPL)

      ! TLJ - modified 12/21/2024
      map(PPICLF_P_JPHIPD) = rprop(PPICLF_R_JVOLP)*rprop(PPICLF_R_JRHOP)
      map(PPICLF_P_JPHIPU) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JVX)
      map(PPICLF_P_JPHIPV) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JVY)
      map(PPICLF_P_JPHIPW) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JVZ)
      map(PPICLF_P_JPHIPT) = rprop(PPICLF_R_JVOLP)*y(PPICLF_JT)


      return
      end
!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine to set user-defined values at time t=0
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_InitZero
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
      integer*4 i,j,k

!
! Code:
!
      ppiclf_TimeBH = 0.0d0



      return
      end

!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine for output if ppiclf_debug=1
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_debug
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
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
      integer*4 i, n, ic, k, iStage

! Needed for allreduce
      integer*4 ngop
      parameter(ngop = 31)
      real*8 xin(ngop),wout(ngop)

! Needed for viscous unsteady
      integer*4 iT,ii
      real*8 time, factor, A, B, fH, kernelVU
      real*8 FVUoutput
      real*8 ppiclf_npart_sum
      integer*4 npart_tot, npart_max, npart_min
      integer*4 ppiclf_iglsum,ppiclf_iglmax,ppiclf_iglmin
      external  ppiclf_iglsum,ppiclf_iglmax,ppiclf_iglmin

      iStage = 1 ! for internal use only

!
! Code:
!
      ! Use ppiclf ALLREDUCE to compute values across processors
      ! Note that ALLREDUCE uses MPI_BARRIER, which is cpu expensive
      ! Print out every 10th iStage=1 counts

         !xin(1) = dfloat(ppiclf_npart)
         !call ppiclf_gop(xin, wout, '+  ', 1)
         !ppiclf_npart_sum = wout(1)
         npart_tot = ppiclf_iglsum(PPICLF_NPART,1)
         npart_max = ppiclf_iglmax(PPICLF_NPART,1)
         npart_min = ppiclf_iglmin(PPICLF_NPART,1)


         xin=(/phimax,
     >         fqsx_max,fqsy_max,fqsz_max,
     >         famx_max,famy_max,famz_max, 
     >         fdpdx_max,fdpdy_max,fdpdz_max, 
     >         fcx_max,fcy_max,fcz_max,
     >         umean_max,vmean_max,wmean_max,
     >         fqs_mag,fam_mag,fdp_mag,fc_mag,
     >         fqsx_fluct_max,fqsy_fluct_max,fqsz_fluct_max,
     >         fqsx_total_max,fqsy_total_max,fqsz_total_max,
     >         fvux_max,fvuy_max,fvuz_max,
     >         qq_max,tau_max/)
         call ppiclf_gop(xin, wout, 'M  ', ngop)
         phimax     = wout(1)
         fqsx_max   = wout(2)
         fqsy_max   = wout(3)
         fqsz_max   = wout(4)
         famx_max   = wout(5)
         famy_max   = wout(6)
         famz_max   = wout(7)
         fdpdx_max  = wout(8)
         fdpdy_max  = wout(9)
         fdpdz_max  = wout(10)
         fcx_max    = wout(11)
         fcy_max    = wout(12)
         fcz_max    = wout(13)
         umean_max  = wout(14)
         vmean_max  = wout(15)
         wmean_max  = wout(16)
         fqs_mag    = wout(17)
         fam_mag    = wout(18)
         fdp_mag    = wout(19)
         fc_mag     = wout(20)
         fqsx_fluct_max = wout(21)
         fqsy_fluct_max = wout(22)
         fqsz_fluct_max = wout(23)
         fqsx_total_max = wout(24)
         fqsy_total_max = wout(25)
         fqsz_total_max = wout(26)
         fvux_max   = wout(27)
         fvuy_max   = wout(28)
         fvuz_max   = wout(29)
         qq_max     = wout(30)
         tau_max    = wout(31)

         ! Sam - logging for debugging purposes
         ! TLJ - below is a mess I created, need to clean up
         if (ppiclf_nid.eq.0) then

         goto 500
            fH     = 0.75d0 + .105d0*reyL
            factor = 3.0d0*rpi*rnu*dp*fac
            FVUoutput = 0.0
            if (ppiclf_nTimeBH>1) then
               do iT = 2,ppiclf_nTimeBH-1
                  time = ppiclf_timeBH(iT)
                  A  = (4.0d0*rpi*time*rnu/dp**2)**(.25d0)
                  B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >              (0.5d0*dp*rnu*fH**3))**(.5d0)
                  kernelVU = factor*(A+B)**(-2)
                  FVUoutput = FVUoutput + kernelVU*
     >               (ppiclf_drudtMixt(1,iT,1)-ppiclf_drudtPlag(1,iT,1))
                  if (abs(FVUoutput) < 1.d-20) FVUoutput = 0.0d0
                  !if (abs(FVUoutput(iT)) < 1.d-20) FVUoutput(iT) = 0.0d0
               enddo
               iT = ppiclf_nTimeBH
               time = ppiclf_timeBH(iT)
               A  = (4.0d0*rpi*time*rnu/dp**2)**(.25d0)
               B  = (0.5d0*rpi*(vmag**3)*(time**2)/ 
     >              (0.5d0*dp*rnu*fH**3))**(.5d0)
               kernelVU = 0.5*factor*(A+B)**(-2)
               FVUoutput = FVUoutput + kernelVU*
     >             (ppiclf_drudtMixt(1,iT,1)-ppiclf_drudtPlag(1,iT,1))
            endif

            WRITE(7225,"(700(1x,E14.6))") ppiclf_time, 
     >        ((ppiclf_drudtMixt(1,i,1)-ppiclf_drudtPlag(1,i,1))
     >        ,i=1,6)
            WRITE(7227,"(2i5,2x,700(2x,E14.6))") iT,iStage,ppiclf_time,
     >        time,FVUoutput,A,B,kernelVU
            WRITE(7229,"(i5,2x,i5,2x,800(1x,E14.6))") ppiclf_nTimeBH,
     >          ppiclf_nUnsteadyData,ppiclf_dt,
     >          ppiclf_time,ppiclf_TimeBH(1:6)

 500        continue

            WRITE(7226,"(700(1x,E14.6))") ppiclf_time,
     >        ((ppiclf_drudtMixt(1,i,1)-ppiclf_drudtPlag(1,i,1))
     >        ,i=1,ppiclf_nUnsteadyData)
            WRITE(7228,"(70(1x,E14.6))") ppiclf_time,
     >        ((ppiclf_drudtMixt(3,i,1)-ppiclf_drudtPlag(3,i,1))
     >        ,i=1,ppiclf_nUnsteadyData)
            WRITE(7230,"(7(1x,E23.16))") ppiclf_time, ppiclf_y(1:6, 1)
            WRITE(7231,"(28(1x,E23.16))") ppiclf_time,phimax,
     >             fqsx_max, fqsy_max, fqsz_max,
     >             famx_max, famy_max, famz_max,
     >             fdpdx_max, fdpdy_max, fdpdz_max,
     >             fcx_max, fcy_max, fcz_max,
     >             qq_max,tau_max,
     >             fqsx_total_max,fqsy_total_max,fqsz_total_max,
     >             fvux_max, fvuy_max, fvuz_max
            WRITE(7232,"(26(1x,F13.8))") ppiclf_time,
     >             umean_max,vmean_max,wmean_max
            WRITE(7233,"(i5,2x,28(1x,E23.16))")
     >             ppiclf_nid,ppiclf_dt,ppiclf_time,
     >             fac, phimax,
     >             fqsx_fluct_max, fqsy_fluct_max, fqsz_fluct_max
            WRITE(7234,*) ppiclf_nid,istage,PPICLF_LRS ,PPICLF_LPART,
     >             PPICLF_NPART,ppiclf_time,
     >             ppiclf_rprop(PPICLF_R_FLUCTFX:PPICLF_R_FLUCTFZ,1),
     >             ppiclf_ydotc(PPICLF_JVX:PPICLF_JT,1)
            WRITE(7235,"(26(1x,F13.8))") ppiclf_time,
     >             fqs_mag,fam_mag,fdp_mag,fc_mag
            WRITE(7236,"(26(1x,F13.8))") ppiclf_time,
     >             fcx_max, fcy_max, fcz_max
            WRITE(7237,"(26(1x,F13.8))") ppiclf_time,UnifRnd
            WRITE(7240,"(26(1x,F13.8))") ppiclf_time,
     >             fqsx_max, fqsy_max, fqsz_max,
     >             fqsx_fluct_max, fqsy_fluct_max, fqsz_fluct_max,
     >             fqsx_total_max,fqsy_total_max,fqsz_total_max
            WRITE(7241,"(5(1x,E23.16))") ppiclf_time,
     >             phipmean, upmean, vpmean, wpmean
            WRITE(7243,"(1x,E23.16,5(2x,I8))") ppiclf_time,
     >             npart_tot,npart_max,npart_min
         endif



      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_comm_InitMPI(comm,id,np)
     > bind(C, name="ppiclc_comm_InitMPI")
#else
      subroutine ppiclf_comm_InitMPI(comm,id,np)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input: 
!
      integer*4 comm
      integer*4 id
      integer*4 np
!
      if (PPICLF_LINIT .or. PPICLF_LFILT .or. PPICLF_OVERLAP)
     >   call ppiclf_exittr('InitMPI must be called first$',0.0d0,0)

      ppiclf_comm = comm
      ppiclf_nid  = id
      ppiclf_np   = np

      call ppiclf_prints('   *Begin InitCrystal$')
         call ppiclf_comm_InitCrystal
      call ppiclf_prints('    End InitCrystal$')

      PPICLF_LCOMM = .true.

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_InitCrystal
!
      implicit none
!
      include "PPICLF"
!
      call pfgslib_crystal_setup(ppiclf_cr_hndl,ppiclf_comm,ppiclf_np)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >                            nfacegp, nedgegp, ncornergp
      integer*4 exit_1_array(3), exit_2_array(3), finished(3)
      integer*4 ix, iy, iz, iperiodicx, iperiodicy, iperiodicz, 
     >          npt_total, j, i, idum, jdum, kdum, total_bin, 
     >          sum_value, count
      real*8 xmin, ymin, zmin, xmax, ymax, zmax, rduml, rdumr, rthresh,
     >       rmiddle, rdiff
      logical exit_1, exit_2
      integer*4 ppiclf_iglsum
      external ppiclf_iglsum
      real*8 ppiclf_glmin,ppiclf_glmax,ppiclf_glsum
      external ppiclf_glmin,ppiclf_glmax,ppiclf_glsum
!

! face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/
     >                 -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                 -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim.eq. 3)
     >iz = 3

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)
         
      ! TLJ this line is not necessary 12/21/2024
      ppiclf_d2chk(1) = max(ppiclf_d2chk(2),ppiclf_d2chk(3))


      ! binning requires > 1 global particle. This takes care of 
      ! single particle case
      npt_total = ppiclf_iglsum(ppiclf_npart,1)
c     if (npt_total .eq. 1) then
      if (.not. ppiclf_lproj .and. .not. ppiclf_lsubsubbin) 
     >   ppiclf_d2chk(1) = 1E-16

      !if (ppiclf_nid==0) print*,'Bins: ', ppiclf_time, ppiclf_d2chk

      ! compute binb
      xmin = 1E10
      ymin = 1E10
      zmin = 1E10
      xmax = -1E10
      ymax = -1E10
      zmax = -1E10
      do i=1,ppiclf_npart
         rduml = ppiclf_y(ix,i) - ppiclf_d2chk(1)
         rdumr = ppiclf_y(ix,i) + ppiclf_d2chk(1)
         if (rduml .lt. xmin) xmin = rduml
         if (rdumr .gt. xmax) xmax = rdumr

         rduml = ppiclf_y(iy,i) - ppiclf_d2chk(1)
         rdumr = ppiclf_y(iy,i) + ppiclf_d2chk(1)
         if (rduml .lt. ymin) ymin = rduml
         if (rdumr .gt. ymax) ymax = rdumr

         if (ppiclf_ndim .eq. 3) then
            rduml = ppiclf_y(iz,i) - ppiclf_d2chk(1)
            rdumr = ppiclf_y(iz,i) + ppiclf_d2chk(1)
            if (rduml .lt. zmin) zmin = rduml
            if (rdumr .gt. zmax) zmax = rdumr
         endif
      enddo

      ppiclf_binb(1) = ppiclf_glmin(xmin,1)
      ppiclf_binb(2) = ppiclf_glmax(xmax,1)
      ppiclf_binb(3) = ppiclf_glmin(ymin,1)
      ppiclf_binb(4) = ppiclf_glmax(ymax,1)
      ppiclf_binb(5) = 0.0d0
      ppiclf_binb(6) = 0.0d0
      if(ppiclf_ndim .gt. 2) ppiclf_binb(5) = ppiclf_glmin(zmin,1)
      if(ppiclf_ndim .gt. 2) ppiclf_binb(6) = ppiclf_glmax(zmax,1)

      if (npt_total .gt. 0) then
      do i=1,ppiclf_ndim
         if (ppiclf_bins_balance(i) .eq. 1) then
            rmiddle = 0.0
            do j=1,ppiclf_npart
               rmiddle = rmiddle + ppiclf_y(i,j)
            enddo
            rmiddle = ppiclf_glsum(rmiddle,1)
            rmiddle = rmiddle/npt_total

            rdiff =  max(abs(rmiddle-ppiclf_binb(2*(i-1)+1)),
     >                   abs(ppiclf_binb(2*(i-1)+2)-rmiddle))
            ppiclf_binb(2*(i-1)+1) = rmiddle - rdiff
            ppiclf_binb(2*(i-1)+2) = rmiddle + rdiff
         endif
      enddo
      endif

      ! Thierry - we comment this out to prevent periodic
      !           algorithm to overwrite bin boundaries

!      if (ang_case==111) then
!      if (ppiclf_xdrange(2,1) .lt. ppiclf_binb(2) .or.
!     >    ppiclf_xdrange(1,1) .gt. ppiclf_binb(1) .or. 
!     >    iperiodicx .eq. 0) then
!         ppiclf_binb(1) = ppiclf_xdrange(1,1)
!         ppiclf_binb(2) = ppiclf_xdrange(2,1)
!      endif
!
!      if (ppiclf_xdrange(2,2) .lt. ppiclf_binb(4) .or.
!     >    ppiclf_xdrange(1,2) .gt. ppiclf_binb(3) .or.
!     >    iperiodicy .eq. 0) then
!         ppiclf_binb(3) = ppiclf_xdrange(1,2)
!         ppiclf_binb(4) = ppiclf_xdrange(2,2)
!      endif
!      
!      endif ! ang_case

      ! Thierry - we make the bins in z-direction as big as the fluid mesh
      !           this is also needed for the bin calculation
      if (ppiclf_ndim .gt. 2) then
      if (ppiclf_xdrange(2,3) .lt. ppiclf_binb(6) .or.
     >    ppiclf_xdrange(1,3) .gt. ppiclf_binb(5) .or. 
     >    iperiodicz .eq. 0) then
         ppiclf_binb(5) = ppiclf_xdrange(1,3)
         ppiclf_binb(6) = ppiclf_xdrange(2,3)
      endif ! ndim
      endif ! xdrange

      if (npt_total .lt. 1) return

      finished(1) = 0
      finished(2) = 0
      finished(3) = 0
      total_bin = 1 

      do i=1,ppiclf_ndim
         finished(i) = 0
         exit_1_array(i) = ppiclf_bins_set(i)
         exit_2_array(i) = 0
         if (ppiclf_bins_set(i) .ne. 1) ppiclf_n_bins(i) = 1
         ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                        ppiclf_binb(2*(i-1)+1)  ) / 
     >                       ppiclf_n_bins(i)
         ! Make sure exit_2 is not violated by user input
         if (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1)) then
            do while (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1))
               ppiclf_n_bins(i) = max(1, ppiclf_n_bins(i)-1)
               ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                              ppiclf_binb(2*(i-1)+1)  ) / 
     >                             ppiclf_n_bins(i)
         WRITE(*,*) "Inf. loop in CreateBin", i, 
     >              ppiclf_bins_dx(i), ppiclf_d2chk(1)
         call ppiclf_exittr('Inf. loop in CreateBin$',0.0,0)
            enddo
         endif
         total_bin = total_bin*ppiclf_n_bins(i)
      enddo

      ! Make sure exit_1 is not violated by user input
      count = 0
      do while (total_bin > ppiclf_np)
          count = count + 1;
          i = modulo((ppiclf_ndim-1)+count,ppiclf_ndim)+1
          ppiclf_n_bins(i) = max(ppiclf_n_bins(i)-1,1)
          ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                         ppiclf_binb(2*(i-1)+1)  ) / 
     >                        ppiclf_n_bins(i)
          total_bin = 1
          do j=1,ppiclf_ndim
             total_bin = total_bin*ppiclf_n_bins(j)
          enddo
          if (total_bin .le. ppiclf_np) exit
       enddo

       exit_1 = .false.
       exit_2 = .false.

       do while (.not. exit_1 .and. .not. exit_2)
          do i=1,ppiclf_ndim
             if (exit_1_array(i) .eq. 0) then
                ppiclf_n_bins(i) = ppiclf_n_bins(i) + 1
                ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                               ppiclf_binb(2*(i-1)+1)  ) / 
     >                              ppiclf_n_bins(i)

                ! Check conditions
                ! exit_1
                total_bin = 1
                do j=1,ppiclf_ndim
                   total_bin = total_bin*ppiclf_n_bins(j)
                enddo
                if (total_bin .gt. ppiclf_np) then
                   ! two exit arrays aren't necessary for now, but
                   ! to make sure exit_2 doesn't slip through, we
                   ! set both for now
                   exit_1_array(i) = 1
                   exit_2_array(i) = 1
                   ppiclf_n_bins(i) = ppiclf_n_bins(i) - 1
                   ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                                  ppiclf_binb(2*(i-1)+1)  ) / 
     >                                  ppiclf_n_bins(i)
                   exit
                endif
                
                ! exit_2
                if (ppiclf_bins_dx(i) .lt. ppiclf_d2chk(1)) then
                   ! two exit arrays aren't necessary for now, but
                   ! to make sure exit_2 doesn't slip through, we
                   ! set both for now
                   exit_1_array(i) = 1
                   exit_2_array(i) = 1
                   ppiclf_n_bins(i) = ppiclf_n_bins(i) - 1
                   ppiclf_bins_dx(i) = (ppiclf_binb(2*(i-1)+2) -
     >                                  ppiclf_binb(2*(i-1)+1)  ) / 
     >                                  ppiclf_n_bins(i)
                   exit
                endif
             endif
          enddo

          ! full exit_1
          sum_value = 0
          do i=1,ppiclf_ndim
             sum_value = sum_value + exit_1_array(i)
          enddo
          if (sum_value .eq. ppiclf_ndim) then
             exit_1 = .true.
          endif

          ! full exit_2
          sum_value = 0
          do i=1,ppiclf_ndim
             sum_value = sum_value + exit_2_array(i)
          enddo
          if (sum_value .eq. ppiclf_ndim) then
             exit_2 = .true.
          endif
       enddo

! -------------------------------------------------------
! SETUP 3D BACKGROUND GRID PARAMETERS FOR GHOST PARTICLES
! -------------------------------------------------------
      ! Check for too small bins 
      rthresh = 1E-12
      total_bin = 1
      do i=1,ppiclf_ndim
         total_bin = total_bin*ppiclf_n_bins(i)
         if (ppiclf_bins_dx(i) .lt. rthresh) ppiclf_bins_dx(i) = 1.0
      enddo

!     current box coordinates
      if (ppiclf_nid .le. total_bin-1) then
         idum = modulo(ppiclf_nid,ppiclf_n_bins(1))
         jdum = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
         kdum = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
         if (ppiclf_ndim .lt. 3) kdum = 0
         ppiclf_binx(1,1) = ppiclf_binb(1) + idum    *ppiclf_bins_dx(1)
         ppiclf_binx(2,1) = ppiclf_binb(1) + (idum+1)*ppiclf_bins_dx(1)
         ppiclf_biny(1,1) = ppiclf_binb(3) + jdum    *ppiclf_bins_dx(2)
         ppiclf_biny(2,1) = ppiclf_binb(3) + (jdum+1)*ppiclf_bins_dx(2)
         ppiclf_binz(1,1) = 0.0d0
         ppiclf_binz(2,1) = 0.0d0
         if (ppiclf_ndim .gt. 2) then
            ppiclf_binz(1,1) = ppiclf_binb(5)+kdum    *ppiclf_bins_dx(3)
            ppiclf_binz(2,1) = ppiclf_binb(5)+(kdum+1)*ppiclf_bins_dx(3)
         endif
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateSubBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 nbin, idum, jdum, kdum, ndumx, ndumy, itmp, jtmp, ktmp,
     >          i, j, k
!

      nbin = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)

c     current box coordinates
      if (ppiclf_nid .le. nbin-1) then
         idum = modulo(ppiclf_nid,ppiclf_n_bins(1))
         jdum = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
         kdum = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
         if (ppiclf_ndim .lt. 3) kdum = 0
         ! interior grid of each bin
         ! +1 for making mesh smaller and +1 since these are vertice counts
         ppiclf_bx = floor(ppiclf_bins_dx(1)/ppiclf_filter) + 1 + 1
         ppiclf_by = floor(ppiclf_bins_dx(2)/ppiclf_filter) + 1 + 1
         ppiclf_bz = 1
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_bz = floor(ppiclf_bins_dx(3)/ppiclf_filter) + 1 + 1

         ppiclf_bx = ppiclf_bx*ppiclf_ngrids
         ppiclf_by = ppiclf_by*ppiclf_ngrids
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_bz = ppiclf_bz*ppiclf_ngrids

         if (ppiclf_bx .gt. PPICLF_BX1)
     >      call ppiclf_exittr('Increase PPICLF_BX1$',0.,ppiclf_bx)
         if (ppiclf_by .gt. PPICLF_BY1)
     >      call ppiclf_exittr('Increase PPICLF_BY1$',0.,ppiclf_by)
         if (ppiclf_bz .gt. PPICLF_BZ1)
     >      call ppiclf_exittr('Increase PPICLF_BZ1$',0.,ppiclf_bz)

         ppiclf_rdx = ppiclf_bins_dx(1)/(ppiclf_bx-1)
         ppiclf_rdy = ppiclf_bins_dx(2)/(ppiclf_by-1)
         ppiclf_rdz = 0
         if (ppiclf_ndim .gt. 2) 
     >      ppiclf_rdz = ppiclf_bins_dx(3)/(ppiclf_bz-1)

         ndumx = ppiclf_n_bins(1)*(ppiclf_bx-1) + 1
         ndumy = ppiclf_n_bins(2)*(ppiclf_by-1) + 1
    
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            ppiclf_grid_x(i,j,k) = sngl(ppiclf_binx(1,1) +
     >                                  (i-1)*ppiclf_rdx)
            ppiclf_grid_y(i,j,k) = sngl(ppiclf_biny(1,1) +
     >                                  (j-1)*ppiclf_rdy)
            ppiclf_grid_z(i,j,k) = sngl(ppiclf_binz(1,1) +
     >                                  (k-1)*ppiclf_rdz)

            itmp = idum*(ppiclf_bx-1) + (i-1)
            jtmp = jdum*(ppiclf_by-1) + (j-1)
            ktmp = kdum*(ppiclf_bz-1) + (k-1)
    
            ppiclf_grid_i(i,j,k)  = itmp + ndumx*jtmp + ndumx*ndumy*ktmp

         enddo
         enddo
         enddo

      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_MapOverlapMesh
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Internal:
!
      integer*4 icalld
      save      icalld
      data      icalld /0/
      integer*4 nkey(2), i, j, k,l, ie, iee, ii, jj, kk, ndum, nrank,
     >          nl, nii, njj, nrr, ilow, jlow, klow, nxyz, il,
     >          ihigh, jhigh, khigh, ierr
      real*8 rxval, ryval, rzval
      logical partl
      real*8 ppiclf_vlmin, ppiclf_vlmax
      external ppiclf_vlmin, ppiclf_vlmax

      ! Sam - for ghost cells
      real*8 xmin(3), xmax(3), xminb(3), xmaxb(3)
      integer*4 nsendg, iig(3), iin(3), iing(3)

      integer*4 ix, iy, iz, ixLow, ixHigh, iyLow,
     >          iyHigh, izLow, izHigh
      ! Avery - for largest cell size
      real*8 EleSizei(3), MaxPoint(3), MinPoint(3) 
     
      ppiclf_neltb = 0 !counts number of Rocflu elements on this processor
                       !that are within one of the ppiclf bins
      DO ie=1,ppiclf_nee
      ! Avery added - find cell max x,y,z lengths
        DO l=1,3
          MaxPoint(l) = -1000000.0d0
          MinPoint(l) =  1000000.0d0 
          EleSizei(l) =  0.0d0
          DO k=1,PPICLF_LEZ
            DO j=1,PPICLF_LEY
              DO i=1,PPICLF_LEX
                IF (ppiclf_xm1bs(i,j,k,l,ie) .GT. MaxPoint(l)) 
     >              MaxPoint(l) = ppiclf_xm1bs(i,j,k,l,ie)
                IF (ppiclf_xm1bs(i,j,k,l,ie) .LT. MinPoint(l)) 
     >              MinPoint(l) = ppiclf_xm1bs(i,j,k,l,ie)
              END DO !i
            END DO !j
          END DO !k
          EleSizei(l) = 1.1*(MaxPoint(l) - MinPoint(l))
        END DO !l 
      ! Avery - end
      do k=1,PPICLF_LEZ
      do j=1,PPICLF_LEY
      do i=1,PPICLF_LEX
         ! Grid positions without additional length
         rxval = ppiclf_xm1bs(i,j,k,1,ie)
         ryval = ppiclf_xm1bs(i,j,k,2,ie)
         rzval = 0.0d0
         if(ppiclf_ndim.gt.2) rzval = ppiclf_xm1bs(i,j,k,3,ie)
       
         ! Exits if cell is outside of all bin boundaries
         if (rxval .gt. ppiclf_binb(2)) goto 1255
         if (rxval .lt. ppiclf_binb(1)) goto 1255
         if (ryval .gt. ppiclf_binb(4)) goto 1255
         if (ryval .lt. ppiclf_binb(3)) goto 1255
         if (ppiclf_ndim.gt.2 .and. rzval .gt. ppiclf_binb(6)) 
     >      goto 1255
         if (ppiclf_ndim.gt.2 .and. rzval .lt. ppiclf_binb(5))
     >      goto 1255
 
         ! Determining what bin the cell is in
         ii    = floor((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
         jj    = floor((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
         kk    = floor((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3))

         ! Default is Do loop with ix=iy=iz=2 for no additional length
         ixLow =2
         ixHigh=2
         iyLow =2
         iyHigh=2
         izLow =2
         izHigh=2

         ! These series of if statements check if near bin boundary 
         ! Default is for no additional bin checks to be applied (do loop set to 2,2)
         ! Add or subtract cell distance with Do loop if cell is near bin boundary
         
         if (floor((rxval + EleSizei(1) -ppiclf_binb(1))
     >   /ppiclf_bins_dx(1)) .NE. ii) then
         ixHigh = 3
         endif

         if (floor((rxval - EleSizei(1) -ppiclf_binb(1))
     >   /ppiclf_bins_dx(1)) .NE. ii) then
         ixLow = 1
         endif

         if (floor((ryval + EleSizei(2) -ppiclf_binb(3))
     >   /ppiclf_bins_dx(2)) .NE. jj) then
         iyHigh = 3
         endif

         if (floor((ryval - EleSizei(2) -ppiclf_binb(3))
     >   /ppiclf_bins_dx(2)) .NE. jj) then
         iyLow = 1
         endif

         if (ppiclf_ndim .gt. 2 .and. floor((rzval + EleSizei(3)
     >   -ppiclf_binb(5))/ppiclf_bins_dx(3)) .NE. kk) then
         izHigh = 3
         endif

         if (ppiclf_ndim .gt. 2 .and. floor((rzval - EleSizei(3)
     >   -ppiclf_binb(5))/ppiclf_bins_dx(3)) .NE. kk) then
         izLow = 1
         endif

      do ix=ixLow,ixHigh
      do iy=iyLow,iyHigh
      do iz=izLow,izHigh
         
         ! Changes r value by element size if near bin
         rxval = ppiclf_xm1bs(i,j,k,1,ie) + (ix-2)*EleSizei(1)
         ryval = ppiclf_xm1bs(i,j,k,2,ie) + (iy-2)*EleSizei(2)
         rzval = 0.0d0
         if(ppiclf_ndim.gt.2) rzval = ppiclf_xm1bs(i,j,k,3,ie)
     >           + (iz-2)*EleSizei(3)

         ! Finds correct bin indicies for cell
         ii    = floor((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
         jj    = floor((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
         kk    = floor((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
         if (ppiclf_ndim.lt.3) kk = 0
         if (ii .eq. ppiclf_n_bins(1)) ii = ppiclf_n_bins(1) - 1
         if (jj .eq. ppiclf_n_bins(2)) jj = ppiclf_n_bins(2) - 1
         if (kk .eq. ppiclf_n_bins(3)) kk = ppiclf_n_bins(3) - 1
         if (ii .eq. -1) ii = 0
         if (jj .eq. -1) jj = 0
         if (kk .eq. -1) kk = 0

         ! Calculates processor rank
         ndum  = ii + ppiclf_n_bins(1)*jj + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
         nrank = ndum

         if (ii .lt. 0 .or. ii .gt. ppiclf_n_bins(1)-1) goto 1233
         if (jj .lt. 0 .or. jj .gt. ppiclf_n_bins(2)-1) goto 1233
         if (kk .lt. 0 .or. kk .gt. ppiclf_n_bins(3)-1) goto 1233

         ppiclf_neltb = ppiclf_neltb + 1
         if(ppiclf_neltb .gt. PPICLF_LEE) then
           call ppiclf_exittr('Increase PPICLF_LEE$',0.0d0,ppiclf_neltb)
         endif

         ppiclf_er_map(1,ppiclf_neltb) = ie
         ppiclf_er_map(2,ppiclf_neltb) = ppiclf_nid
         ppiclf_er_map(3,ppiclf_neltb) = ndum
         ppiclf_er_map(4,ppiclf_neltb) = nrank
         ppiclf_er_map(5,ppiclf_neltb) = nrank
         ppiclf_er_map(6,ppiclf_neltb) = nrank

         if (ppiclf_neltb .gt. 1) then
         do il=1,ppiclf_neltb-1
            if (ppiclf_er_map(1,il) .eq. ie) then
            if (ppiclf_er_map(4,il) .eq. nrank) then
               ppiclf_neltb = ppiclf_neltb - 1
               goto 1233
            endif
            endif
         enddo
         endif
 1233 continue
      enddo !iz
      enddo !iy
      enddo !ix
 1255 continue ! When a cell is outside the bin boundary
      enddo !k
      enddo !i
      enddo !j
      enddo !ie

      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_neltb
       iee = ppiclf_er_map(1,ie)
       call ppiclf_copy(ppiclf_xm1b(1,1,1,1,ie)
     >                 ,ppiclf_xm1bs(1,1,1,1,iee),nxyz)
       call ppiclf_copy(ppiclf_xm1b(1,1,1,2,ie)
     >                 ,ppiclf_xm1bs(1,1,1,2,iee),nxyz)
       call ppiclf_copy(ppiclf_xm1b(1,1,1,3,ie)
     >                 ,ppiclf_xm1bs(1,1,1,3,iee),nxyz)
      enddo

      ppiclf_neltbb = ppiclf_neltb
      do ie=1,ppiclf_neltbb
         call ppiclf_icopy(ppiclf_er_maps(1,ie),ppiclf_er_map(1,ie)
     >             ,PPICLF_LRMAX)
      enddo


      nl   = 0
      nii  = PPICLF_LRMAX
      njj  = 6
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      nrr  = nxyz*3
      nkey(1) = 2
      nkey(2) = 1
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl,ppiclf_neltb
     >       ,PPICLF_LEE,ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,njj)
      call pfgslib_crystal_tuple_sort    (ppiclf_cr_hndl,ppiclf_neltb
     >       ,ppiclf_er_map,nii,partl,nl,ppiclf_xm1b,nrr,nkey,2)


      do ie=1,ppiclf_neltb
      do k=1,PPICLF_LEZ
      do j=1,PPICLF_LEY
      do i=1,PPICLF_LEX
         rxval = ppiclf_xm1b(i,j,k,1,ie)
         ryval = ppiclf_xm1b(i,j,k,2,ie)
         rzval = 0.0d0
         if(ppiclf_ndim.gt.2) rzval = ppiclf_xm1b(i,j,k,3,ie)
         
         ii    = floor((rxval-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
         jj    = floor((ryval-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
         kk    = floor((rzval-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
         if (ppiclf_ndim.eq.2) kk = 0
          if (ii .eq. ppiclf_n_bins(1)) ii = ppiclf_n_bins(1) - 1
          if (jj .eq. ppiclf_n_bins(2)) jj = ppiclf_n_bins(2) - 1
          if (kk .eq. ppiclf_n_bins(3)) kk = ppiclf_n_bins(3) - 1
          if (ii .eq. -1) ii = 0
          if (jj .eq. -1) jj = 0
          if (kk .eq. -1) kk = 0
          ndum  = ii + ppiclf_n_bins(1)*jj + 
     >                 ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk

         ppiclf_modgp(i,j,k,ie,1) = ii
         ppiclf_modgp(i,j,k,ie,2) = jj
         ppiclf_modgp(i,j,k,ie,3) = kk
         ppiclf_modgp(i,j,k,ie,4) = ndum
   
      enddo
      enddo
      enddo
      enddo

      do ie=1,ppiclf_neltb
         ppiclf_xerange(1,1,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,1,ie),nxyz)
         ppiclf_xerange(2,1,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,1,ie),nxyz)
         ppiclf_xerange(1,2,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,2,ie),nxyz)
         ppiclf_xerange(2,2,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,2,ie),nxyz)
         ppiclf_xerange(1,3,ie) = 
     >      ppiclf_vlmin(ppiclf_xm1b(1,1,1,3,ie),nxyz)
         ppiclf_xerange(2,3,ie) = 
     >      ppiclf_vlmax(ppiclf_xm1b(1,1,1,3,ie),nxyz)

         ilow  = 
     >     floor((ppiclf_xerange(1,1,ie) - ppiclf_binb(1))/
     >                                             ppiclf_bins_dx(1))
         ihigh = 
     >     floor((ppiclf_xerange(2,1,ie) - ppiclf_binb(1))/
     >                                             ppiclf_bins_dx(1))
         jlow  = 
     >     floor((ppiclf_xerange(1,2,ie) - ppiclf_binb(3))/
     >                                             ppiclf_bins_dx(2))
         jhigh = 
     >     floor((ppiclf_xerange(2,2,ie) - ppiclf_binb(3))/
     >                                             ppiclf_bins_dx(2))
         klow  = 
     >     floor((ppiclf_xerange(1,3,ie) - ppiclf_binb(5))/
     >                                             ppiclf_bins_dx(3))
         khigh = 
     >     floor((ppiclf_xerange(2,3,ie) - ppiclf_binb(5))/
     >                                             ppiclf_bins_dx(3))
         if (ppiclf_ndim.lt.3) then
            klow = 0
            khigh = 0
         endif

         ppiclf_el_map(1,ie) = ilow  + ppiclf_n_bins(1)*jlow  
     >                         + ppiclf_n_bins(1)*ppiclf_n_bins(2)*klow
         ppiclf_el_map(2,ie) = ihigh + ppiclf_n_bins(1)*jhigh 
     >                         + ppiclf_n_bins(1)*ppiclf_n_bins(2)*khigh
         ppiclf_el_map(3,ie) = ilow
         ppiclf_el_map(4,ie) = ihigh
         ppiclf_el_map(5,ie) = jlow
         ppiclf_el_map(6,ie) = jhigh
         ppiclf_el_map(7,ie) = klow
         ppiclf_el_map(8,ie) = khigh
      enddo

      if (icalld .eq. 0) then 

         icalld = icalld + 1

         call ppiclf_prints('   *Begin mpi_comm_split$')
            call mpi_comm_split(ppiclf_comm
     >                         ,ppiclf_nid
     >                         ,0
     >                         ,ppiclf_comm_nid
     >                         ,ierr)
         call ppiclf_prints('    End mpi_comm_split$')

         ! TLJ commented out recursive loop
         !call ppiclf_prints('   *Begin InitSolve$')
         !   call ppiclf_solve_InitSolve
         !call ppiclf_prints('    End InitSolve$')

         call ppiclf_io_OutputDiagGrid
      endif

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_comm_InitOverlapMesh(ncell,lx1,ly1,lz1,
     >                                       xgrid,ygrid,zgrid)
     > bind(C, name="ppiclc_comm_InitOverlapMesh")
#else
      subroutine ppiclf_comm_InitOverlapMesh(ncell,lx1,ly1,lz1,
     >                                       xgrid,ygrid,zgrid)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 ncell
      integer*4 lx1
      integer*4 ly1
      integer*4 lz1
      real*8    xgrid(*)
      real*8    ygrid(*)
      real*8    zgrid(*)
!
! External:
!
      integer*4 nxyz, i, j, ie
      integer*4 k, jj, icont
!
      ppiclf_overlap = .true.

      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitOverlap$',0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitOverlap$'
     >                  ,0.0d0,0)

      if (ncell .gt. PPICLF_LEE .or. ncell .lt. 0) 
     >   call ppiclf_exittr('Increase LEE in InitOverlap$',0.0d0,ncell)
      if (lx1 .ne. PPICLF_LEX) 
     >   call ppiclf_exittr('LX1 != LEX in InitOverlap$',0.0d0,ncell)
      if (ly1 .ne. PPICLF_LEY)
     >   call ppiclf_exittr('LY1 != LEY in InitOverlap$',0.0d0,ncell)
      if (lz1 .ne. PPICLF_LEZ)
     >   call ppiclf_exittr('LZ1 != LEZ in InitOverlap$',0.0d0,ncell)

      ppiclf_nee = ncell
      nxyz       = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ

      do ie=1,ppiclf_nee
         ! TLJ changing loop structure
         !do i=1,nxyz
         !   j = (ie-1)*nxyz + i
         !   ppiclf_xm1bs(i,1,1,1,ie) = xgrid(j)
         !   ppiclf_xm1bs(i,1,1,2,ie) = ygrid(j)
         !   ppiclf_xm1bs(i,1,1,3,ie) = zgrid(j)
         !enddo
         icont = 0
         do k=1,PPICLF_LEZ
         do j=1,PPICLF_LEY
         do i=1,PPICLF_LEX
            icont = icont + 1
            jj = (ie-1)*nxyz + icont
            ppiclf_xm1bs(i,j,k,1,ie) = xgrid(jj)
            ppiclf_xm1bs(i,j,k,2,ie) = ygrid(jj)
            ppiclf_xm1bs(i,j,k,3,ie) = zgrid(jj)
         enddo
         enddo
         enddo
      enddo
      
      call ppiclf_solve_InitSolve

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_FindParticle
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 ix, iy, iz, i, ii, jj, kk, ndum, nrank
!
      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim.eq.3)
     >iz = 3

      do i=1,ppiclf_npart
         ! check if particles are greater or less than binb bounds....
         ii  = floor((ppiclf_y(ix,i)-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
         jj  = floor((ppiclf_y(iy,i)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
         kk  = floor((ppiclf_y(iz,i)-ppiclf_binb(5))/ppiclf_bins_dx(3)) 
         if (ppiclf_ndim .lt. 3) kk = 0
         ndum  = ii + ppiclf_n_bins(1)*jj + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kk
         nrank = ndum

         ppiclf_iprop(8,i)  = ii
         ppiclf_iprop(9,i)  = jj
         ppiclf_iprop(10,i) = kk
         ppiclf_iprop(11,i) = ndum

         ppiclf_iprop(3,i)  = nrank ! where particle is actually moved
         ppiclf_iprop(4,i)  = nrank ! where particle is actually moved
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_comm_MoveParticle
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      logical partl    
      integer*4 lrf
      parameter(lrf = PPICLF_LRS*4 + PPICLF_LRP + PPICLF_LRP2
     >       + PPICLF_LRP3)
      real*8 rwork(lrf,PPICLF_LPART)
      integer*4 i, ic, j0
!

      do i=1,ppiclf_npart
         ic = 1
         call ppiclf_copy(rwork(ic,i),ppiclf_y(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_y1((i-1)*PPICLF_LRS+1)
     >                   ,PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_ydot(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_ydotc(1,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop(1,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop2(1,i),PPICLF_LRP2)
         ic = ic + PPICLF_LRP2
         call ppiclf_copy(rwork(ic,i),ppiclf_rprop3(1,i),PPICLF_LRP3)
      enddo

      j0 = 4
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >                                  ,ppiclf_npart,PPICLF_LPART
     >                                  ,ppiclf_iprop,PPICLF_LIP
     >                                  ,partl,0
     >                                  ,rwork,lrf
     >                                  ,j0)

      if (ppiclf_npart .gt. PPICLF_LPART .or. ppiclf_npart .lt. 0)
     >   call ppiclf_exittr('Increase LPART$',0.0d0,ppiclf_npart)

      do i=1,ppiclf_npart
         ic = 1
         call ppiclf_copy(ppiclf_y(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_y1((i-1)*PPICLF_LRS+1),rwork(ic,i)
     >                   ,PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_ydot(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_ydotc(1,i),rwork(ic,i),PPICLF_LRS)
         ic = ic + PPICLF_LRS
         call ppiclf_copy(ppiclf_rprop(1,i),rwork(ic,i),PPICLF_LRP)
         ic = ic + PPICLF_LRP
         call ppiclf_copy(ppiclf_rprop2(1,i),rwork(ic,i),PPICLF_LRP2)
         ic = ic + PPICLF_LRP2
         call ppiclf_copy(ppiclf_rprop3(1,i),rwork(ic,i),PPICLF_LRP3)
      enddo
        
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_comm_CreateGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 xdlen,ydlen,zdlen,rxdrng(3),rxnew(3), rfac, rxval, ryval,
     >       rzval, rxl, ryl, rzl, rxr, ryr, rzr, distchk, dist
      integer*4 iadd(3),gpsave(27)
      real*8 map(PPICLF_LRP_PRO)
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >           nfacegp, nedgegp, ncornergp, iperiodicx, iperiodicy,
     >           iperiodicz, jx, jy, jz, ip, idum, iip, jjp, kkp, ii1,
     >           jj1, kk1, iig, jjg, kkg, iflgx, iflgy, iflgz,
     >           isave, iflgsum, ndumn, nrank, ibctype, i, ifc, ist, j,
     >           k
!

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      jx    = 1
      jy    = 2
      jz    = 3

      ! Thierry - we do not assign the bins to be as big as
      !           the periodic domain in x/y directions anymore. only in z. 
      
      !xdlen = ppiclf_binb(2) - ppiclf_binb(1) ! when bins = periodic domain
      !ydlen = ppiclf_binb(4) - ppiclf_binb(3) ! when bins = periodic domain
      
      ! Thierry - this works whether the bins are as big as periodic domain, or not.
      xdlen = ppiclf_xdrange(2,1) - ppiclf_xdrange(1,1)
      ydlen = ppiclf_xdrange(2,2) - ppiclf_xdrange(1,2)
      
      zdlen = -1.
      if (ppiclf_ndim .gt. 2) 
!     >   zdlen = ppiclf_binb(6) - ppiclf_binb(5)
      ! Thierry - this works whether the bins are as big as periodic domain, or not.
     >   zdlen = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
      if (iperiodicx .ne. 0) xdlen = -1
      if (iperiodicy .ne. 0) ydlen = -1
      if (iperiodicz .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

      ppiclf_npart_gp = 0

      rfac = 1.0d0

      do ip=1,ppiclf_npart

         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))

c        idum = 1
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 2
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 3
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)

         idum = 0
         do j=1,PPICLF_LRS
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip)
         enddo
         idum = PPICLF_LRS
         do j=1,PPICLF_LRP
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip)
         enddo
         idum = PPICLF_LRS+PPICLF_LRP
         do j=1,PPICLF_LRP_PRO
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = map(j)
         enddo

         rxval = ppiclf_cp_map(1,ip)
         ryval = ppiclf_cp_map(2,ip)
         rzval = 0.0d0
         if (ppiclf_ndim .gt. 2) rzval = ppiclf_cp_map(3,ip)

         iip    = ppiclf_iprop(8,ip)
         jjp    = ppiclf_iprop(9,ip)
         kkp    = ppiclf_iprop(10,ip)

         rxl = ppiclf_binb(1) + ppiclf_bins_dx(1)*iip
         rxr = rxl + ppiclf_bins_dx(1)
         ryl = ppiclf_binb(3) + ppiclf_bins_dx(2)*jjp
         ryr = ryl + ppiclf_bins_dx(2)
         rzl = 0.0d0
         rzr = 0.0d0
         if (ppiclf_ndim .gt. 2) then
            rzl = ppiclf_binb(5) + ppiclf_bins_dx(3)*kkp
            rzr = rzl + ppiclf_bins_dx(3)
         endif

         isave = 0

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) 
            jj1 = jjp + el_face_num(ist+2)
            kk1 = kkp + el_face_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_n_bins(1))
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_n_bins(2))
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_n_bins(3))
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)

            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_n_bins(1))
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_n_bins(2))
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_n_bins(3))
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)

            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)
            if (dist .gt. distchk) cycle

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
               iflgx = 1
               iig =modulo(iig,ppiclf_n_bins(1))
               if (iperiodicx .ne. 0) cycle
            endif
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
               iflgy = 1
               jjg =modulo(jjg,ppiclf_n_bins(2))
               if (iperiodicy .ne. 0) cycle
            endif
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
               iflgz = 1  
               kkg =modulo(kkg,ppiclf_n_bins(3))
               if (iperiodicz .ne. 0) cycle
            endif

            iflgsum = iflgx + iflgy + iflgz
            ndumn = iig + ppiclf_n_bins(1)*jjg 
     >                  + ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
            nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
                 
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval
       
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1

            call ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1)
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2)
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3)

            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  333 continue
         enddo

      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_AngularCreateGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8 xdlen,ydlen,zdlen,rxdrng(3),rxnew(3), rfac, rxval, ryval,
     >       rzval, rxl, ryl, rzl, rxr, ryr, rzr, distchk, dist
      integer*4 iadd(3),gpsave(27)
      real*8 map(PPICLF_LRP_PRO)
      integer*4  el_face_num(18),el_edge_num(36),el_corner_num(24),
     >           nfacegp, nedgegp, ncornergp, iperiodicx, iperiodicy,
     >           iperiodicz, jx, jy, jz, ip, idum, iip, jjp, kkp, ii1,
     >           jj1, kk1, iig, jjg, kkg, iflgx, iflgy, iflgz,
     >           isave, iflgsum, ndumn, nrank, ibctype, i, ifc, ist, j,
     >           k
      ! 08/27/24 - Thierry - added for angular periodicty starts here
      real*8 alpha
      integer*4 xrank, yrank, zrank
      ! 08/27/24 - Thierry - added for angular periodicty ends here
      ! 09/26/24 - Thierry - added for angular periodicty starts here
      real*8 dist1, dist2
      ! 09/26/24 - Thierry - added for angular periodicty ends here
!

c     face, edge, and corner number, x,y,z are all inline, so stride=3
      el_face_num = (/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /)
      el_edge_num = (/ -1,-1,0 , 1,-1,0, 1,1,0 , -1,1,0 ,
     >                  0,-1,-1, 1,0,-1, 0,1,-1, -1,0,-1,
     >                  0,-1,1 , 1,0,1 , 0,1,1 , -1,0,1  /)
      el_corner_num = (/ -1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
     >                   -1,-1,1,  1,-1,1,  1,1,1,  -1,1,1 /)

      nfacegp   = 4  ! number of faces
      nedgegp   = 4  ! number of edges
      ncornergp = 0  ! number of corners

      if (ppiclf_ndim .gt. 2) then
         nfacegp   = 6  ! number of faces
         nedgegp   = 12 ! number of edges
         ncornergp = 8  ! number of corners
      endif

      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)

! ------------------------
c CREATING GHOST PARTICLES
! ------------------------
      jx    = 1
      jy    = 2
      jz    = 3

      ! Thierry - we dont use xdlen and ydlen in this algorithm. no need to modify them.
      xdlen = ppiclf_binb(2) - ppiclf_binb(1)
      ydlen = ppiclf_binb(4) - ppiclf_binb(3)
      zdlen = -1.
      if (ppiclf_ndim .gt. 2) 
!     >   zdlen = ppiclf_binb(6) - ppiclf_binb(5)
     >   zdlen = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
      if (iperiodicx .ne. 0) xdlen = -1
      if (iperiodicy .ne. 0) ydlen = -1
      if (iperiodicz .ne. 0) zdlen = -1

      rxdrng(1) = xdlen
      rxdrng(2) = ydlen
      rxdrng(3) = zdlen

      ppiclf_npart_gp = 0

      rfac = 1.0d0

      do ip=1,ppiclf_npart

         call ppiclf_user_MapProjPart(map,ppiclf_y(1,ip)
     >         ,ppiclf_ydot(1,ip),ppiclf_ydotc(1,ip),ppiclf_rprop(1,ip))

c        idum = 1
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 2
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)
c        idum = 3
c        ppiclf_cp_map(idum,ip) = ppiclf_y(idum,ip)

         idum = 0
         do j=1,PPICLF_LRS
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_y(j,ip) ! ppiclf_y(PPICLF_JX/ JY/ JZ/ JVX/ JVY/ JVZ/ JT, ip)
         enddo
         idum = PPICLF_LRS
         do j=1,PPICLF_LRP
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = ppiclf_rprop(j,ip) ! ppiclf_rprop(PPICLF_R_JRHOP/ R_JRHOF/ .../ R_WDOTZ, ip)
         enddo
         idum = PPICLF_LRS+PPICLF_LRP
         do j=1,PPICLF_LRP_PRO
            idum = idum + 1
            ppiclf_cp_map(idum,ip) = map(j) ! map(PPICLF_P_JPHIP/ JFX/ .../ JPHIPW) - these are found in ppiclf_user_MapProjPart
         enddo

         rxval = ppiclf_cp_map(1,ip) ! ppiclf_y(PPICLF_JX,ip)
         ryval = ppiclf_cp_map(2,ip) ! ppiclf_y(PPICLF_JY,ip)
         rzval = 0.0d0
         if (ppiclf_ndim .gt. 2) rzval = ppiclf_cp_map(3,ip) ! ppiclf_y(PPICLF_JZ,ip)

         iip    = ppiclf_iprop(8,ip) ! ith coordinate of bin
         jjp    = ppiclf_iprop(9,ip) ! jth coordinate of bin
         kkp    = ppiclf_iprop(10,ip) ! kth coordinate of bin

         rxl = ppiclf_binb(1) + ppiclf_bins_dx(1)*iip ! min x of bin
         rxr = rxl + ppiclf_bins_dx(1)                ! max x of bin
         ryl = ppiclf_binb(3) + ppiclf_bins_dx(2)*jjp
         ryr = ryl + ppiclf_bins_dx(2)
         rzl = 0.0d0
         rzr = 0.0d0
         if (ppiclf_ndim .gt. 2) then
            rzl = ppiclf_binb(5) + ppiclf_bins_dx(3)*kkp
            rzr = rzl + ppiclf_bins_dx(3)
         endif

         isave = 0

         ! faces
         do ifc=1,nfacegp
            ist = (ifc-1)*3
            ii1 = iip + el_face_num(ist+1) 
            jj1 = jjp + el_face_num(ist+2)
            kk1 = kkp + el_face_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)

            if (ang_case==1) then  ! for wedge geometry

               ! Thierry - I dont think it's code efficient to call this subroutine
               !           for every particle, every ghost face, at every time step
               !           I'm wondering if it's better if we make the plane values 
               !           as global values that are initialized in the beginning 
            
               call ppiclf_solve_InitAngularPlane(ip,
     >                                 ang_per_rin  , ang_per_rout  ,
     >                                 ang_per_angle, ang_per_xangle,
     >                                 dist1, dist2)
               if ((dist .gt. distchk).and.(dist1.gt.distchk)
     >           .and.(dist2.gt.distchk)) cycle
            else
               if (dist .gt. distchk) cycle
            endif

            iflgx = 0
            iflgy = 0
            iflgz = 0
!-----------------------------------------------------------------------
            ! 08/27/24 - Thierry - modification for angular periodicty starts here

               ! angle between particle and x-axis
                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
     >                        ppiclf_y(PPICLF_JX,ip))
                

                call ppiclf_solve_InvokeAngularPeriodic(ip, 
     >                                                  ang_per_flag,
     >                                                  alpha,         
     >                                                  ang_per_angle,  
     >                                                  ang_per_xangle, 
     >                                                  0)

              ! Thierry - this is how FindParticle implements it
              ! need to find a way to make the code deal with negative xrot values

            xrank = iig ; yrank=jjg; zrank = kkg
            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
            
            ! Angular periodicity check in x- and y-directions
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
              iflgx = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicx .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
              iflgy = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicy .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            ! Linear periodicity check in z-direction
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
              iflgz = 1
              kkg =modulo(kkg,ppiclf_n_bins(3))
              if (iperiodicz .ne. 0) cycle
              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
              ! rxdrng(3) = -1.0  if not periodic in Z
              if (rxdrng(3) .gt. 0) then 
                if (iadd(3) .ge. ppiclf_n_bins(3)) then ! particle leaving from max z-face
                  rxnew(3) = rxnew(3) - rxdrng(3)
                elseif (iadd(3) .lt. 0) then ! particle leaving from min z-face
                  rxnew(3) = rxnew(3) + rxdrng(3)
                end if ! iadd
              end if ! rxrdrng
            else ! z-linear periodicity not applicable
              kkg = zrank
            end if ! kkg
            
            iflgsum = iflgx + iflgy + iflgz
            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
             nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!-----------------------------------------------------------------------

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 111
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz
            
            rxnew(1) = xrot(1)
            rxnew(2) = xrot(2)
            ppiclf_cp_map(4,ip) = vrot(1)
            ppiclf_cp_map(5,ip) = vrot(2)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
            !           as this is now taken care of when anticipating where the particle might be when calling
            !           ppiclf_comm_InvokeAngularPeriodic
            !           we only need to assign xr and vr to ppiclf_rprop_gp

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ
            
            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  111 continue
         enddo

         ! edges
         do ifc=1,nedgegp
            ist = (ifc-1)*3
            ii1 = iip + el_edge_num(ist+1) 
            jj1 = jjp + el_edge_num(ist+2)
            kk1 = kkp + el_edge_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)

            if (ang_case==1) then  ! for wedge geometry

               call ppiclf_solve_InitAngularPlane(ip,
     >                                 ang_per_rin  , ang_per_rout  ,
     >                                 ang_per_angle, ang_per_xangle,
     >                                 dist1, dist2)
               if ((dist .gt. distchk).and.(dist1.gt.distchk)
     >           .and.(dist2.gt.distchk)) cycle
            else
               if (dist .gt. distchk) cycle
            endif

            iflgx = 0
            iflgy = 0
            iflgz = 0
            ! periodic if out of domain - add some ifsss
!-----------------------------------------------------------------------
            ! 08/27/24 - Thierry - modification for angular periodicty starts here

               ! angle between particle and x-axis
                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
     >                        ppiclf_y(PPICLF_JX,ip))
                

                call ppiclf_solve_InvokeAngularPeriodic(ip, 
     >                                                  ang_per_flag,
     >                                                  alpha,         
     >                                                  ang_per_angle,  
     >                                                  ang_per_xangle, 
     >                                                  0)

              ! Thierry - this is how FindParticle implements it
              ! need to find a way to make the code deal with negative xrot values

            xrank = iig ; yrank=jjg; zrank = kkg
            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
            
            ! Angular periodicity check in x- and y-directions
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
              iflgx = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicx .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
              iflgy = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicy .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            ! Linear periodicity check in z-direction
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
              iflgz = 1
              kkg =modulo(kkg,ppiclf_n_bins(3))
              if (iperiodicz .ne. 0) cycle
              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
              ! rxdrng(3) = -1.0  if not periodic in Z
              if (rxdrng(3) .gt. 0) then ! particle leaving from max z-face
                if (iadd(3) .ge. ppiclf_n_bins(3)) then
                  rxnew(3) = rxnew(3) - rxdrng(3)
                elseif (iadd(3) .lt. 0) then
                  rxnew(3) = rxnew(3) + rxdrng(3)
                end if ! iadd
              end if ! rxrdrng
            else ! z-linear periodicity not applicable
              kkg = zrank
            end if ! kkg
            
            iflgsum = iflgx + iflgy + iflgz
            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
             nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!-----------------------------------------------------------------------

            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 222
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz

            rxnew(1) = xrot(1)
            rxnew(2) = xrot(2)
            ppiclf_cp_map(4,ip) = vrot(1)
            ppiclf_cp_map(5,ip) = vrot(2)
                 
            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
            !           as this is now taken care of when anticipating where the particle might be when calling
            !           ppiclf_comm_InvokeAngularPeriodic
            !           we only need to assign xr and vr to ppiclf_rprop_gp

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ
            
            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  222 continue
         enddo

         ! corners
         do ifc=1,ncornergp
            ist = (ifc-1)*3
            ii1 = iip + el_corner_num(ist+1) 
            jj1 = jjp + el_corner_num(ist+2)
            kk1 = kkp + el_corner_num(ist+3)

            iig = ii1
            jjg = jj1
            kkg = kk1

            distchk = 0.0d0
            dist = 0.0d0
            if (ii1-iip .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (ii1-iip .lt. 0) dist = dist +(rxval - rxl)**2
               if (ii1-iip .gt. 0) dist = dist +(rxval - rxr)**2
            endif
            if (jj1-jjp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (jj1-jjp .lt. 0) dist = dist +(ryval - ryl)**2
               if (jj1-jjp .gt. 0) dist = dist +(ryval - ryr)**2
            endif
            if (ppiclf_ndim .gt. 2) then
            if (kk1-kkp .ne. 0) then
               distchk = distchk + (rfac*ppiclf_d2chk(1))**2
               if (kk1-kkp .lt. 0) dist = dist +(rzval - rzl)**2
               if (kk1-kkp .gt. 0) dist = dist +(rzval - rzr)**2
            endif
            endif
            distchk = sqrt(distchk)
            dist = sqrt(dist)

            if (ang_case==1) then  ! for wedge geometry
            
               call ppiclf_solve_InitAngularPlane(ip,
     >                                 ang_per_rin  , ang_per_rout  ,
     >                                 ang_per_angle, ang_per_xangle,
     >                                 dist1, dist2)
               if ((dist .gt. distchk).and.(dist1.gt.distchk)
     >           .and.(dist2.gt.distchk)) cycle
            else
               if (dist .gt. distchk) cycle
            endif

            iflgx = 0
            iflgy = 0
            iflgz = 0

!-----------------------------------------------------------------------
            ! 08/27/24 - Thierry - modification for angular periodicty starts here

               ! angle between particle and x-axis
                alpha = atan2(ppiclf_y(PPICLF_JY,ip), 
     >                        ppiclf_y(PPICLF_JX,ip))
                

                call ppiclf_solve_InvokeAngularPeriodic(ip, 
     >                                                  ang_per_flag,
     >                                                  alpha,         
     >                                                  ang_per_angle,  
     >                                                  ang_per_xangle, 
     >                                                  0)

              ! Thierry - this is how FindParticle implements it
              ! need to find a way to make the code deal with negative xrot values

            xrank = iig ; yrank=jjg; zrank = kkg
            ! Thierry - previously placed before the CheckPeriodicBC call, had to move them for the periodic check
            iadd(1) = ii1
            iadd(2) = jj1
            iadd(3) = kk1
            rxnew(1) = rxval
            rxnew(2) = ryval
            rxnew(3) = rzval ! z-coordinate does not change when angular periodicity is invoked
            
            ! Angular periodicity check in x- and y-directions
            if (iig .lt. 0 .or. iig .gt. ppiclf_n_bins(1)-1) then
              iflgx = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicx .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            if (jjg .lt. 0 .or. jjg .gt. ppiclf_n_bins(2)-1) then
              iflgy = 1
              call ppiclf_comm_CheckAngularBC(xrank,yrank,zrank)
              if (iperiodicy .ne. 0) cycle
              iig = xrank
              jjg = yrank
            end if
            
            ! Linear periodicity check in z-direction
            if (kkg .lt. 0 .or. kkg .gt. ppiclf_n_bins(3)-1) then
              iflgz = 1
              kkg =modulo(kkg,ppiclf_n_bins(3))
              if (iperiodicz .ne. 0) cycle
              ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
              ! rxdrng(3) = -1.0  if not periodic in Z
              if (rxdrng(3) .gt. 0) then ! particle leaving from max z-face
                if (iadd(3) .ge. ppiclf_n_bins(3)) then
                  rxnew(3) = rxnew(3) - rxdrng(3)
                elseif (iadd(3) .lt. 0) then
                  rxnew(3) = rxnew(3) + rxdrng(3)
                end if ! iadd
              end if ! rxrdrng
            else ! z-linear periodicity not applicable
              kkg = zrank
            end if ! kkg
            
            iflgsum = iflgx + iflgy + iflgz
            ndumn  = iig + ppiclf_n_bins(1)*jjg + 
     >                ppiclf_n_bins(1)*ppiclf_n_bins(2)*kkg
             nrank = ndumn

            if (nrank .eq. ppiclf_nid .and. iflgsum .eq. 0) cycle

            ! 08/27/24 - Thierry - modification for angular periodicty ends here
!-----------------------------------------------------------------------
            do i=1,isave
               if (gpsave(i) .eq. nrank .and. iflgsum .eq.0) goto 333
            enddo
            isave = isave + 1
            gpsave(isave) = nrank

            ibctype = iflgx+iflgy+iflgz

            rxnew(1) = xrot(1)
            rxnew(2) = xrot(2)
            ppiclf_cp_map(4,ip) = vrot(1)
            ppiclf_cp_map(5,ip) = vrot(2)

            ppiclf_npart_gp = ppiclf_npart_gp + 1
            ppiclf_iprop_gp(1,ppiclf_npart_gp) = nrank
            ppiclf_iprop_gp(2,ppiclf_npart_gp) = iig
            ppiclf_iprop_gp(3,ppiclf_npart_gp) = jjg
            ppiclf_iprop_gp(4,ppiclf_npart_gp) = kkg
            ppiclf_iprop_gp(5,ppiclf_npart_gp) = ndumn

            ! Thierry - we don't need ppiclf_comm_CheckPeriodicBC anymore for the angular periodic ghost algorithm
            !           as this is now taken care of when anticipating where the particle might be when calling
            !           ppiclf_comm_InvokeAngularPeriodic
            !           we only need to assign xr and vr to ppiclf_rprop_gp

            ppiclf_rprop_gp(1,ppiclf_npart_gp) = rxnew(1) ! ppiclf_y(PPICLF_JX, ip) for the periodic ghost particle
            ppiclf_rprop_gp(2,ppiclf_npart_gp) = rxnew(2) ! JY
            ppiclf_rprop_gp(3,ppiclf_npart_gp) = rxnew(3) ! JZ

            do k=4,PPICLF_LRP_GP
               ppiclf_rprop_gp(k,ppiclf_npart_gp) = ppiclf_cp_map(k,ip)
            enddo
  333 continue
         enddo

      enddo ! ip 

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_CheckPeriodicBC(rxnew,rxdrng,iadd)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      real*8 rxdrng(3)
      integer*4 iadd(3)
!
! Input/Output:
!
      real*8 rxnew(3)
!
      ! rxdrng(1) = ppiclf_xdrange(2,1) - ppiclf_xdrange(1,1)
      ! rxdrng(1) = -1.0  if not periodic in X
      ! particle leaving from max x periodic face
      if (rxdrng(1) .gt. 0 ) then
      if (iadd(1) .ge. ppiclf_n_bins(1)) then
         rxnew(1) = rxnew(1) - rxdrng(1)
         goto 123
      endif
      endif
      ! particle leaving from min x periodic face
      if (rxdrng(1) .gt. 0 ) then
      if (iadd(1) .lt. 0) then
         rxnew(1) = rxnew(1) + rxdrng(1)
         goto 123
      endif
      endif

  123 continue    
      ! rxdrng(2) = ppiclf_xdrange(2,2) - ppiclf_xdrange(1,2)
      ! rxdrng(2) = -1.0  if not periodic in Y
      ! particle leaving from max y periodic face
      if (rxdrng(2) .gt. 0 ) then
      if (iadd(2) .ge. ppiclf_n_bins(2)) then
         rxnew(2) = rxnew(2) - rxdrng(2)
         goto 124
      endif
      endif
      if (rxdrng(2) .gt. 0 ) then
      ! particle leaving from min y periodic face
      if (iadd(2) .lt. 0) then
         rxnew(2) = rxnew(2) + rxdrng(2)
         goto 124
      endif
      endif
  124 continue

      if (ppiclf_ndim .gt. 2) then
        ! rxdrng(3) = ppiclf_xdrange(2,3) - ppiclf_xdrange(1,3)
        ! rxdrng(3) = -1.0  if not periodic in Z
      ! particle leaving from max z periodic face
         if (rxdrng(3) .gt. 0 ) then
         if (iadd(3) .ge. ppiclf_n_bins(3)) then
            rxnew(3) = rxnew(3) - rxdrng(3)
            goto 125
         endif
         endif
      ! particle leaving from min z periodic face
         if (rxdrng(3) .gt. 0 ) then
         if (iadd(3) .lt. 0) then
            rxnew(3) = rxnew(3) + rxdrng(3)
            goto 125
         endif
         endif
      endif
  125 continue

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_CheckAngularBC(xrank, yrank, zrank)
!
      implicit none
!
      include "PPICLF"
!
! Local:
!
      integer*4 xrank, yrank, zrank
!
! Output:
!

      SELECT CASE (ang_case)
        CASE(1) ! general wedge ; 0 <= angle < 90
!          print*, "Wedge CheckAngularBC"
          xrank  = floor((xrot(1)-ppiclf_binb(1))/ppiclf_bins_dx(1)) 
          yrank  = floor((xrot(2)-ppiclf_binb(3))/ppiclf_bins_dx(2)) 
          zrank  = floor((xrot(3)-ppiclf_binb(5))/ppiclf_bins_dx(3))

        CASE(2) ! quarter cylinder ; angle = 90
!          print*, "Quarter Cylinder CheckAngularBC"
          xrank  = floor((abs(xrot(1))-ppiclf_binb(1))
     >                    /ppiclf_bins_dx(1)) 
          yrank  = floor((abs(xrot(2))-ppiclf_binb(3))
     >                   /ppiclf_bins_dx(2)) 
          zrank  = floor((xrot(3)-ppiclf_binb(5))/ppiclf_bins_dx(3))

        CASE(3) ! half cylinder ; angle = 180
          print*, "Half Cylinder CheckAngularBC"

        CASE DEFAULT
            call ppiclf_exittr('Invalid Ghost Rotational Case!$',
     >       0.0d0 ,ppiclf_nid)
          END SELECT

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_comm_MoveGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      logical partl         
!
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl
     >                                  ,ppiclf_npart_gp,PPICLF_LPART_GP
     >                                  ,ppiclf_iprop_gp,PPICLF_LIP_GP
     >                                  ,partl,0
     >                                  ,ppiclf_rprop_gp,PPICLF_LRP_GP
     >                                  ,1)

      return
      end
c----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine ppiclf_gop( x, w, op, n)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      real*8 x(n), w(n)
      character*3 op
      integer*4 n
!
! Internal:
!
      integer*4 i, ie
!
      if (op.eq.'+  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_sum,ppiclf_comm,ie)
      elseif (op.EQ.'M  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_max,ppiclf_comm,ie)
      elseif (op.EQ.'m  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_min,ppiclf_comm,ie)
      elseif (op.EQ.'*  ') then
      call mpi_allreduce
     >        (x,w,n,MPI_DOUBLE_PRECISION,mpi_prod,ppiclf_comm,ie)
      endif

      do i=1,n
         x(i) = w(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_igop( x, w, op, n)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      integer*4 x(n), w(n)
      character*3 op
      integer*4 n
!
! Internal:
!
      integer*4 i, ierr
!
      if     (op.eq.'+  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_sum ,ppiclf_comm,ierr)
      elseif (op.EQ.'M  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_max ,ppiclf_comm,ierr)
      elseif (op.EQ.'m  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_min ,ppiclf_comm,ierr)
      elseif (op.EQ.'*  ') then
        call MPI_Allreduce (x,w,n,mpi_integer,mpi_prod,ppiclf_comm,ierr)
      endif

      do i=1,n
         x(i) = w(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer*4 function ppiclf_iglsum(a,n)
! 
      implicit none
! 
! Input:
! 
      integer*4 a(1)
      integer*4 n
! 
! Internal:
! 
      integer*4 tsum
      integer*4 tmp(1),work(1)
      integer*4 i
!
      tsum= 0
      do i=1,n
         tsum=tsum+a(i)
      enddo
      tmp(1)=tsum
      call ppiclf_igop(tmp,work,'+  ',1)
      ppiclf_iglsum=tmp(1)
      return
      end
C-----------------------------------------------------------------------
      real*8 function ppiclf_glsum (x,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A
      real*8 x(n),tsum
      real*8 tmp(1),work(1)
      integer*4 i,n
!
      TSUM = 0.0d0
      DO 100 I=1,N
         TSUM = TSUM+X(I)
 100  CONTINUE
      TMP(1)=TSUM
      CALL ppiclf_GOP(TMP,WORK,'+  ',1)
      ppiclf_GLSUM = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      real*8 function ppiclf_glmax(a,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A
      REAL*8 A(n),tmax
      real*8 TMP(1),WORK(1)
      integer*4 i,n
!
      TMAX=-99.0e20
      DO 100 I=1,N
         TMAX=MAX(TMAX,A(I))
  100 CONTINUE
      TMP(1)=TMAX
      CALL ppiclf_GOP(TMP,WORK,'M  ',1)
      ppiclf_GLMAX=TMP(1)
      return
      END
c-----------------------------------------------------------------------
      integer*4 function ppiclf_iglmax(a,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A
      integer*4 a(n),tmax
      integer*4 tmp(1),work(1)
      integer*4 i,n
!
      tmax= -999999999
      do i=1,n
         tmax=max(tmax,a(i))
      enddo
      tmp(1)=tmax
      call ppiclf_igop(tmp,work,'M  ',1)
      ppiclf_iglmax=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      real*8 function ppiclf_glmin(a,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A
      REAL*8 A(n),tmin
      real*8 TMP(1),WORK(1)
      integer*4 i,n
!
      TMIN=99.0e20
      DO 100 I=1,N
         TMIN=MIN(TMIN,A(I))
  100 CONTINUE
      TMP(1)=TMIN
      CALL ppiclf_GOP(TMP,WORK,'m  ',1)
      ppiclf_GLMIN = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      integer*4 function ppiclf_iglmin(a,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of a
      integer*4 a(n),tmin
      integer*4 tmp(1),work(1)
      integer*4 i, n
!
      tmin=  999999999
      do i=1,n
         tmin=min(tmin,a(i))
      enddo
      tmp(1)=tmin
      call ppiclf_igop(tmp,work,'m  ',1)
      ppiclf_iglmin=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      real*8 function ppiclf_vlmin(vec,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of VEC
      REAL*8 VEC(n),tmin
      integer*4 i, n
!
      TMIN = 99.0E20
      DO 100 I=1,N
         TMIN = MIN(TMIN,VEC(I))
 100  CONTINUE
      ppiclf_VLMIN = TMIN
      return
      END
c-----------------------------------------------------------------------
      real*8 function ppiclf_vlmax(vec,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of VEC
      REAL*8 VEC(n),tmax
      integer*4 i, n
!
      TMAX =-99.0E20
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      ppiclf_VLMAX = TMAX
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_copy(a,b,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A and B
      real*8 a(n),b(n)
      integer*4 i,n
!

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_icopy(a,b,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed dimension of A and B
      INTEGER*4 A(n), B(n)
      integer*4 i,n
!
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_chcopy(a,b,n)
! 
      implicit none
! 
! Vars:
! 
      ! TLJ changed A and B dimenions
      CHARACTER*1 A(n), B(n)
      integer*4 i,n
!
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_exittr(stringi,rdata,idata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
! 
! Vars:
! 
      character*1 stringi(132)
      character*1 stringo(132)
      character*25 s25
      integer*4 ilen, ierr, k, idata
      integer*4 ppiclf_indx1
      real*8 rdata
      external ppiclf_indx1
!
      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s25,25) rdata,idata
   25 format(1x,1p1e14.6,i10)
      call ppiclf_chcopy(stringo(ilen),s25,25)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+24)
    1 format('PPICLF: ERROR ',132a1)

c     call mpi_finalize (ierr)
      call mpi_abort(ppiclf_comm, 1, ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_printsri(stringi,rdata,idata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character*1 stringi(132)
      character*1 stringo(132)
      character*25 s25
      integer*4 ilen, idata, k, ierr
      integer*4 ppiclf_indx1
      real*8 rdata
      external ppiclf_indx1
!

      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s25,25) rdata,idata
   25 format(1x,1p1e14.6,i10)
      call ppiclf_chcopy(stringo(ilen),s25,25)

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+24)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_printsi(stringi,idata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character*1 stringi(132)
      character*1 stringo(132)
      character*10 s10
      integer*4 ilen, idata, k, ierr
      integer*4 ppiclf_indx1
      external ppiclf_indx1
!
      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s10,10) idata
   10 format(1x,i9)
      call ppiclf_chcopy(stringo(ilen),s10,10)

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+9)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_printsr(stringi,rdata)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character*1 stringi(132)
      character*1 stringo(132)
      character*15 s15
      integer*4 ilen, k, ierr
      integer*4 ppiclf_indx1
      real*8 rdata
      external ppiclf_indx1
!
      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')
      write(s15,15) rdata
   15 format(1x,1p1e14.6)
      call ppiclf_chcopy(stringo(ilen),s15,15)

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen+14)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_prints(stringi)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character*1 stringi(132)
      character*1 stringo(132)
      integer*4 ilen, k, ierr
      integer*4 ppiclf_indx1
      external ppiclf_indx1
!
      call ppiclf_blank(stringo,132)
      call ppiclf_chcopy(stringo,stringi,132)
      ilen = ppiclf_indx1(stringo,'$')

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid.eq.0) write(6,1) (stringo(k),k=1,ilen-1)
    1 format('PPICLF: ',132a1)

      call mpi_barrier(ppiclf_comm,ierr)

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE PPICLF_BLANK(A,N)
! 
      implicit none
! 
! Vars:
!
      ! TLJ changed dimension of A
      CHARACTER*1 A(N)
      CHARACTER*1 BLNK
      SAVE        BLNK
      DATA        BLNK /' '/
      integer*4 i,n
!
C
      DO 10 I=1,N
         A(I)=BLNK
   10 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      INTEGER*4 FUNCTION PPICLF_INDX1(S1,S2)
! 
      implicit none
! 
! Vars:
!
      CHARACTER*1 S1(132),S2
      integer*4 n1, i
!
      N1=132
      PPICLF_INDX1=0
      IF (N1.LT.1) return
C
      DO 100 I=1,N1
         IF (S1(I).EQ.S2) THEN
            PPICLF_INDX1=I
            return
         ENDIF
  100 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      character*132 FUNCTION PPICLF_CHSTR(S1,indx1)
! 
      implicit none
! 
! Vars:
!
      ! TLJ modified, but not sure why I had to
      CHARACTER S1
      INTEGER indx1
!
      PPICLF_CHSTR = S1(1:indx1)
      return
      END
c-----------------------------------------------------------------------
      subroutine ppiclf_byte_open_mpi(fnamei,mpi_fh,ifro,ierr)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      character fnamei*(*)
      logical ifro
      CHARACTER*1 BLNK
      DATA BLNK/' '/
      character*132 fname
      character*1   fname1(132)
      equivalence  (fname1,fname)
      integer*4 imode, ierr, mpi_fh
!
      imode = MPI_MODE_WRONLY+MPI_MODE_CREATE
      if(ifro) then
        imode = MPI_MODE_RDONLY 
      endif

      call MPI_file_open(ppiclf_comm,fnamei,imode,
     &                   MPI_INFO_NULL,mpi_fh,ierr)

      return
      end
C--------------------------------------------------------------------------
      subroutine ppiclf_byte_read_mpi(buf,icount,mpi_fh,ierr)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      real*4 buf(1)          ! buffer
      integer*4 iout, icount, mpi_fh, ierr
!

      iout = icount ! icount is in 4-byte words
      call MPI_file_read_all(mpi_fh,buf,iout,MPI_REAL,
     &                       MPI_STATUS_IGNORE,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_write_mpi(buf,icount,iorank,mpi_fh,ierr)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      real*4 buf(1)          ! buffer
      integer*4 icount, iorank, mpi_fh, ierr, iout
!

      iout = icount ! icount is in 4-byte words
      if(iorank.ge.0 .and. ppiclf_nid.ne.iorank) iout = 0
      call MPI_file_write_all(mpi_fh,buf,iout,MPI_REAL,
     &                        MPI_STATUS_IGNORE,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_close_mpi(mpi_fh,ierr)
! 
      implicit none
! 
      include 'mpif.h'
!
! Vars:
!
      integer*4 mpi_fh, ierr
!

      call MPI_file_close(mpi_fh,ierr)

      return
      end
c--------------------------------------------------------------------------
      subroutine ppiclf_byte_set_view(ioff_in,mpi_fh)
! 
      implicit none
! 
      include 'mpif.h'
!
! Vars:
!
      integer*8 ioff_in
      integer*4 mpi_fh, ierr
!
      call MPI_file_set_view(mpi_fh,ioff_in,MPI_BYTE,MPI_BYTE,
     &                       'native',MPI_INFO_NULL,ierr)

      return
      end
C--------------------------------------------------------------------------
      subroutine ppiclf_bcast(buf,len)
! 
      implicit none
! 
      include "PPICLF"
      include 'mpif.h'
!
! Vars:
!
      real*4 buf(1)
      integer*4 len, ierr
!

      call mpi_bcast (buf,len,mpi_byte,0,ppiclf_comm,ierr)

      return
      end
C--------------------------------------------------------------------------
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_io_ReadParticleVTU(filein1,istartoutin,
     > npart,dp_max)
     > bind(C, name="ppiclc_io_ReadParticleVTU")
#else
      subroutine ppiclf_io_ReadParticleVTU(filein1,istartoutin,
     > npart,dp_max)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      character*1 filein1(132)
!
! Internal:
!
      real*4  rout_pos(3      *PPICLF_LPART) 
     >       ,rout_sln(PPICLF_LRS*PPICLF_LPART)
     >       ,rout_lrp(PPICLF_LRP*PPICLF_LPART)
     >       ,rout_lip(3      *PPICLF_LPART)
      real*8 dp_max
      character*1 dum_read
      character*132 filein2
      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip,stride_len
      integer*4 vtu, isize, jx, jy, jz, ivtu_size, ifound, i, npt_total,
     >          npart_min, npmax, npart, ndiff, iorank, icount_pos,
     >          icount_sln, icount_lrp, icount_lip, j, ic_pos, ic_sln,
     >          ic_lrp, ic_lip, pth, ierr, indx1
      integer*4 ppiclf_indx1
      external ppiclf_indx1
      character*132 PPICLF_CHSTR
      EXTERNAL PPICLF_CHSTR
      ! Sam - this modifies the interface for ppiclC. Will now need to
      ! include NULL for optional arguments
      integer*4, optional :: istartoutin
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!
      call ppiclf_prints(' *Begin ReadParticleVTU$')
      
      call ppiclf_solve_InitZero
      PPICLF_RESTART = .true.

      if (present(istartoutin)) then
        istartout = istartoutin
      else
        istartout = 0
      end if


      indx1 = ppiclf_indx1(filein1,'.')
      indx1 = indx1 + 3 ! v (1) t (2) u (3)
      filein2 = ppiclf_chstr(filein1(1:indx1),indx1)
      if (ppiclf_nid == 0) then
         print*,'   * ParticleVTU filein2 ',trim(filein2)
      endif

      isize = 4
      jx    = 1
      jy    = 2
      jz    = 1
      if (ppiclf_ndim .eq. 3)
     >jz    = 3

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=trim(filein2)
     >    ,access='stream',form="unformatted")

      ivtu_size = -1
      ifound = 0
      do i=1,1000000
      read(vtu) dum_read
      if (dum_read == '_') ifound = ifound + 1
      if (ifound .eq. 2) then
         ivtu_size = i
         exit
      endif
      enddo
      read(vtu) npt_total
      close(vtu)
      npt_total = npt_total/isize/3
      endif

      call ppiclf_bcast(npt_total,isize)
      call ppiclf_bcast(ivtu_size,isize)


      npart_min = npt_total/ppiclf_np+1
      if (npart_min*ppiclf_np .gt. npt_total) npart_min = npart_min-1

      if (npt_total .gt. PPICLF_LPART*ppiclf_np) 
     >   call ppiclf_exittr('Increase LPART to at least$',0.0d0
     >    ,npart_min)


      npmax = min(npt_total/PPICLF_LPART+1,ppiclf_np)
      stride_len = 0
      if (ppiclf_nid .le. npmax-1 .and. ppiclf_nid. ne. 0) 
     >stride_len = ppiclf_nid*PPICLF_LPART

      npart = PPICLF_LPART
      if (ppiclf_nid .gt. npmax-1) npart = 0

      ndiff = npt_total - (npmax-1)*PPICLF_LPART
      if (ppiclf_nid .eq. npmax-1) npart = ndiff

      iorank = -1

      call ppiclf_byte_open_mpi(trim(filein2),pth,.true.,ierr)

      idisp_pos = ivtu_size + isize*(3*stride_len + 1)
      icount_pos = npart*3   
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_read_mpi(rout_pos,icount_pos,pth,ierr)

      do i=1,PPICLF_LRS
         idisp_sln = ivtu_size + isize*(3*npt_total 
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + i)
         j = 1 + (i-1)*npart
         icount_sln = npart

         call ppiclf_byte_set_view(idisp_sln,pth)
         call ppiclf_byte_read_mpi(rout_sln(j),icount_sln
     >                            ,pth,ierr)
      enddo

      do i=1,PPICLF_LRP
         idisp_lrp = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + i)

         j = 1 + (i-1)*npart
         icount_lrp = npart

         call ppiclf_byte_set_view(idisp_lrp,pth)
         call ppiclf_byte_read_mpi(rout_lrp(j),icount_lrp
     >                            ,pth,ierr)
      enddo

      do i=1,3
         idisp_lip = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + PPICLF_LRP*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + PPICLF_LRP + i)

         j = 1 + (i-1)*npart
         icount_lip = npart

         call ppiclf_byte_set_view(idisp_lip,pth)
         call ppiclf_byte_read_mpi(rout_lip(j),icount_lip
     >                            ,pth,ierr)
      enddo

      call ppiclf_byte_close_mpi(pth,ierr)

      ic_pos = 0
      ic_sln = 0
      ic_lrp = 0
      ic_lip = 0
      do i=1,npart
         ic_pos = ic_pos + 1
         ppiclf_y(jx,i) = rout_pos(ic_pos)
         ic_pos = ic_pos + 1
         ppiclf_y(jy,i) = rout_pos(ic_pos)
         ic_pos = ic_pos + 1
         if (ppiclf_ndim .eq. 3) then
         ppiclf_y(jz,i) = rout_pos(ic_pos)
         endif
      enddo
      do j=1,PPICLF_LRS
      do i=1,npart
         ic_sln = ic_sln + 1
         ppiclf_y(j,i) = rout_sln(ic_sln)
      enddo
      enddo
      do j=1,PPICLF_LRP
      do i=1,npart
         ic_lrp = ic_lrp + 1
         ppiclf_rprop(j,i) = rout_lrp(ic_lrp)
      enddo
      enddo
      do j=5,7
      do i=1,npart
         ic_lip = ic_lip + 1
         ppiclf_iprop(j,i) = int(rout_lip(ic_lip))
      enddo
      enddo

      ppiclf_npart = npart

      ! TLJ - for restarts we need dp_max
      ! For the moment we assume monodispersed packs
      dp_max = ppiclf_rprop(PPICLF_R_JDP,1)
      call ppiclf_printsi('  End ReadParticleVTU$',npt_total)

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_io_ReadWallVTK(filein1)
     > bind(C, name="ppiclc_io_ReadWallVTK")
#else
      subroutine ppiclf_io_ReadWallVTK(filein1)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      character*1 filein1(132)
!
! Internal:
!
      real*8 points(3,4*PPICLF_LWALL)
      integer*4 fid, nmax, i, j, i1, i2, i3, isize, irsize, ierr,
     >          npoints, nwalls, indx1
      character*1000 text
      character*132 filein2
      integer*4 ppiclf_indx1
      external ppiclf_indx1
      character*132 PPICLF_CHSTR
      external PPICLF_CHSTR
!
      !! THROW ERRORS HERE IN FUTURE
      call ppiclf_prints(' *Begin ReadWallVTK$')
      
      indx1 = ppiclf_indx1(filein1,'.')
      indx1 = indx1 + 3 ! v (1) t (2) k (3)
      !filein2 = ppiclf_chstr(filein1(1:indx1))
      filein2 = 'filein.vtk' !ppiclf_chstr(filein1(1:indx1))

      if (ppiclf_nid .eq. 0) then

      fid = 432

      open (unit=fid,file=trim(filein2),action="read")
      
      nmax = 10000
      do i=1,nmax
         read (fid,*,iostat=ierr) text,npoints

         do j=1,npoints
            read(fid,*) points(1,j),points(2,j),points(3,j)
         enddo

         read (fid,*,iostat=ierr) text,nwalls

         do j=1,nwalls
            if (ppiclf_ndim .eq. 2) then
               read(fid,*) i1,i2
    
               i1 = i1 + 1
               i2 = i2 + 1

               call ppiclf_solve_InitWall( 
     >                 (/points(1,i1),points(2,i1)/),
     >                 (/points(1,i2),points(2,i2)/),
     >                 (/points(1,i1),points(2,i1)/))  ! dummy 2d
    
            elseif (ppiclf_ndim .eq. 3) then
               read(fid,*) i1,i2,i3

               i1 = i1 + 1
               i2 = i2 + 1
               i3 = i3 + 1
    
               call ppiclf_solve_InitWall( 
     >                 (/points(1,i1),points(2,i1),points(3,i1)/),
     >                 (/points(1,i2),points(2,i2),points(3,i2)/),
     >                 (/points(1,i3),points(2,i3),points(3,i3)/))

            endif
         enddo
            
         exit
      enddo

      close(fid)

      endif

      isize  = 4
      irsize = 8
      call ppiclf_bcast(ppiclf_nwall, isize)
      call ppiclf_bcast(ppiclf_wall_c,9*PPICLF_LWALL*irsize)
      call ppiclf_bcast(ppiclf_wall_n,4*PPICLF_LWALL*irsize)

      call ppiclf_printsi('  End ReadWallVTK$',nwalls)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteSubBinVTU(filein1)
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      character*3 filein
      character*12 vtufile
      character*6  prostr
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth, nvtx_total, ncll_total
      integer*8 idisp_pos
      integer*8 stride_lenv
      integer*4 icount_pos(PPICLF_BX1, PPICLF_BY1, PPICLF_BZ1)
      real*4 rpoint(3)
      integer*4 iint, nnp, nxx, ndxgpp1, ndygpp1, ndxygpp1, if_sz,
     >          isize, ibin, jbin, kbin, ndumx, ndumy, i, j, k,
     >          ie, ndum, kmax, ii, jj, kk, itmp, jtmp, ktmp, il, ir,
     >          jl, jr, kl, kr, npa, npb, npc, npd, npe, npf, npg, nph,
     >          ioff_dum, itype, ivtu_size, if_pos, ierr, icount_dum,
     >          iorank
!

      call ppiclf_printsi(' *Begin WriteSubBinVTU$',ppiclf_cycle)

      call ppiclf_prints(' *Begin InitSolve$')
         call ppiclf_solve_InitSolve
      call ppiclf_prints('  End InitSolve$')

      call ppiclf_prints(' *Begin CreateSubBin$')
         call ppiclf_comm_CreateSubBin
      call ppiclf_prints('  End CreateSubBin$')

      call ppiclf_prints(' *Begin ProjectParticleSubBin$')
         call ppiclf_solve_ProjectParticleSubBin
      call ppiclf_prints('  End ProjectParticleSubBin$')

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      ncll_total = ppiclf_n_bins(1)*(ppiclf_bx-1)
     >            *ppiclf_n_bins(2)*(ppiclf_by-1)
      if (ppiclf_ndim.eq.3) ncll_total = ncll_total
     >            *ppiclf_n_bins(3)*(ppiclf_bz-1)
      nvtx_total = (ppiclf_n_bins(1)*(ppiclf_bx-1)+1)
     >            *(ppiclf_n_bins(2)*(ppiclf_by-1)+1)
      if (ppiclf_ndim.eq.3) nvtx_total = nvtx_total
     >            *(ppiclf_n_bins(3)*(ppiclf_bz-1)+1)

      ndxgpp1 = ppiclf_n_bins(1) + 1
      ndygpp1 = ppiclf_n_bins(2) + 1
      ndxygpp1 = ndxgpp1*ndygpp1

      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'grd'
      else 
         write(filein,'(A3)') filein1
      endif

      isize  = 4

      ! get which bin this processor holds
      ibin = modulo(ppiclf_nid,ppiclf_n_bins(1))
      jbin = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
      kbin = 0
      if (ppiclf_ndim .eq. 3)
     >kbin = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

! test skip
c     goto 1511

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="ProjectParticleSubBinInt32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') nvtx_total
      write(vtu,'(A)',advance='no') '" NumberOfCells="'
      write(vtu,'(I0)',advance='no') ncll_total
      write(vtu,'(A)',advance='yes') '"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*nvtx_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'
      do ie=1,PPICLF_LRP_PRO
         write(prostr,'(A4,I2.2)') "PRO-",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*nvtx_total + isize
      enddo
c     call ppiclf_io_WriteDataArrayVTU(vtu,"PPR",1,iint)
c     iint = iint + 1*isize*nvtx_total + isize
c     call ppiclf_io_WriteDataArrayVTU(vtu,"NID",1,iint)
c     iint = iint + 1*isize*nvtx_total + isize
      write(vtu,'(A)',advance='yes') '   </PointData> '
      write(vtu,'(A)',advance='yes') '   <CellData>'
      write(vtu,'(A)',advance='yes') '   </CellData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write connectivity here
      ndumx = ppiclf_n_bins(1)*(ppiclf_bx-1) + 1
      ndumy = ppiclf_n_bins(2)*(ppiclf_by-1) + 1
      ndum = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)
      do ie=0,ndum-1
         i = modulo(ie,ppiclf_n_bins(1))
         j = modulo(ie/ppiclf_n_bins(1),ppiclf_n_bins(2))
         k = ie/(ppiclf_n_bins(1)*ppiclf_n_bins(2))

      kmax = 1
      if (ppiclf_ndim .eq. 3) kmax = ppiclf_bz-1

      do kk=1,kmax
      do jj=1,ppiclf_by-1
      do ii=1,ppiclf_bx-1

         itmp = i*(ppiclf_bx-1) + (ii-1)
         jtmp = j*(ppiclf_by-1) + (jj-1)
         ktmp = 0
         if (ppiclf_ndim .eq. 3)
     >   ktmp = k*(ppiclf_bz-1) + (kk-1)

         kl = ktmp
         kr = ktmp+1
         jl = jtmp
         jr = jtmp+1
         il = itmp
         ir = itmp+1

c        ndum = itmp + ndumx*jtmp + ndumx*ndumy*ktmp

         if (ppiclf_ndim .eq. 3) then
            npa = il + ndumx*jl + ndumx*ndumy*kl
            npb = ir + ndumx*jl + ndumx*ndumy*kl
            npc = il + ndumx*jr + ndumx*ndumy*kl
            npd = ir + ndumx*jr + ndumx*ndumy*kl
            npe = il + ndumx*jl + ndumx*ndumy*kr
            npf = ir + ndumx*jl + ndumx*ndumy*kr
            npg = il + ndumx*jr + ndumx*ndumy*kr
            nph = ir + ndumx*jr + ndumx*ndumy*kr

            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0,A)',advance='no')  npd, ' '
            write(vtu,'(I0,A)',advance='no')  npe, ' '
            write(vtu,'(I0,A)',advance='no')  npf, ' '
            write(vtu,'(I0,A)',advance='no')  npg, ' '
            write(vtu,'(I0)'  ,advance='yes') nph
         else
            npa = il + ndumx*jl + ndumx*ndumy*kl
            npb = ir + ndumx*jl + ndumx*ndumy*kl
            npc = il + ndumx*jr + ndumx*ndumy*kl
            npd = ir + ndumx*jr + ndumx*ndumy*kl

            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0)'  ,advance='yes') npd
         endif
      enddo
      enddo
      enddo
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ioff_dum = 4
      if (ppiclf_ndim .eq. 3) ioff_dum = 8
      ! write offsetts here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') ioff_dum*i
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="UInt8" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      itype = 8
      if (ppiclf_ndim .eq. 3) itype = 11
      ! write types here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') itype
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

c1511 continue

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif

      call ppiclf_bcast(ivtu_size, isize)

      iorank = -1

      if_pos = 3*isize*nvtx_total

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write points first
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      do k=1,ppiclf_bz
      do j=1,ppiclf_by
      do i=1,ppiclf_bx
         icount_pos(i,j,k) = 0
      enddo
      enddo
      enddo
      if (ppiclf_nid .le. 
     >      ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ppiclf_ndim .eq. 3) then
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            if (i .ne. ppiclf_bx .and.
     >          j .ne. ppiclf_by .and.
     >          k .ne. ppiclf_bz) then
                  icount_pos(i,j,k) = 3
            endif

            if (i .eq. ppiclf_bx) then
            if (j .ne. ppiclf_by) then
            if (k .ne. ppiclf_bz) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (j .eq. ppiclf_by) then
            if (i .ne. ppiclf_bx) then
            if (k .ne. ppiclf_bz) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (k .eq. ppiclf_bz) then
            if (i .ne. ppiclf_bx) then
            if (j .ne. ppiclf_by) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (j .eq. ppiclf_by) then
            if (k .ne. ppiclf_bz) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (k .eq. ppiclf_bz) then
            if (j .ne. ppiclf_by) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (j .eq. ppiclf_by) then
            if (k .eq. ppiclf_bz) then
            if (i .ne. ppiclf_bx) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (j .eq. ppiclf_by) then
            if (k .eq. ppiclf_bz) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
            if (kbin .eq. ppiclf_n_bins(3)-1) then
               icount_pos(i,j,k) = 3
            endif
            endif
            endif
            endif
            endif
            endif

         enddo
         enddo
         enddo
         ! 2D
         else
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            if (i .ne. ppiclf_bx .and.
     >          j .ne. ppiclf_by) then
                  icount_pos(i,j,1) = 3
            endif

            if (i .eq. ppiclf_bx) then
            if (j .ne. ppiclf_by) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
               icount_pos(i,j,1) = 3
            endif
            endif
            endif

            if (j .eq. ppiclf_by) then
            if (i .ne. ppiclf_bx) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,1) = 3
            endif
            endif
            endif

            if (i .eq. ppiclf_bx) then
            if (j .eq. ppiclf_by) then
            if (ibin .eq. ppiclf_n_bins(1)-1) then
            if (jbin .eq. ppiclf_n_bins(2)-1) then
               icount_pos(i,j,1) = 3
            endif
            endif
            endif
            endif
         enddo
         enddo

         endif

      endif

      do k=1,ppiclf_bz
      do j=1,ppiclf_by
      do i=1,ppiclf_bx
         stride_lenv = ppiclf_grid_i(i,j,k)
         idisp_pos   = ivtu_size + isize*(3*stride_lenv + 1)
         icount_dum  = icount_pos(i,j,k)
         rpoint(1)   = ppiclf_grid_x(i,j,k)
         rpoint(2)   = ppiclf_grid_y(i,j,k)
         rpoint(3)   = 0.0d0
         if (ppiclf_ndim .eq. 3)
     >   rpoint(3)   = ppiclf_grid_z(i,j,k)
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_dum,iorank,pth,ierr)

      enddo
      enddo
      enddo

      call ppiclf_byte_close_mpi(pth,ierr)


      ! projected fields
      do ie=1,PPICLF_LRP_PRO

      if_pos = 1*isize*nvtx_total

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      do k=1,ppiclf_bz
      do j=1,ppiclf_by
      do i=1,ppiclf_bx
           stride_lenv = ppiclf_grid_i(i,j,k)
           idisp_pos   = ivtu_size + isize*(3*nvtx_total ! position fld
     >                   + (ie-1)*nvtx_total ! prev fields
     >                   + 1*stride_lenv    ! this fld
     >                   + 1 + ie)          ! ints
           icount_dum  = icount_pos(i,j,k)/3 ! either zero or 1
           rpoint(1)   = ppiclf_grid_fld(i,j,k,ie)
           call ppiclf_byte_set_view(idisp_pos,pth)
           call ppiclf_byte_write_mpi(rpoint,icount_dum,iorank,pth,ierr)
      enddo
      enddo
      enddo

      call ppiclf_byte_close_mpi(pth,ierr)

      enddo

      call mpi_barrier(ppiclf_comm,ierr)

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi('  End WriteSubBinVTU$',ppiclf_cycle)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteBinVTU(filein1)
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      character*3 filein
      character*12 vtufile
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth, nvtx_total, ncll_total
      integer*8 idisp_pos,idisp_cll,stride_lenv(8),stride_lenc
      integer*4 iint, nnp, nxx, ndxgpp1, ndygpp1, ndxygpp1, if_sz, ibin,
     >          jbin, kbin, il, ir, jl, jr, kl, kr, nbinpa, nbinpb,
     >          nbinpc, nbinpd, nbinpe, nbinpf, nbinpg, nbinph, i, j, k,
     >          ii, npa, npb, npc, npd, npe, npf, npg, nph, 
     >          ioff_dum, itype, iorank, if_cll, if_pos, icount_pos, 
     >          icount_cll, ierr, isize, ivtu_size
      real*4 rpoint(3)
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!

      call ppiclf_printsi(' *Begin WriteBinVTU$',ppiclf_cycle)

      if (icalld1 .eq. 0) icalld1 = istartout

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      nvtx_total = (ppiclf_n_bins(1)+1)*(ppiclf_n_bins(2)+1)
      if (ppiclf_ndim .gt. 2) 
     >    nvtx_total = nvtx_total*(ppiclf_n_bins(3)+1)
      ncll_total = ppiclf_n_bins(1)*ppiclf_n_bins(2)
      if (ppiclf_ndim .gt. 2) ncll_total = ncll_total*ppiclf_n_bins(3)

      ndxgpp1 = ppiclf_n_bins(1) + 1
      ndygpp1 = ppiclf_n_bins(2) + 1
      ndxygpp1 = ndxgpp1*ndygpp1

      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'bin'
      else 
         write(filein,'(A3)') filein1
      endif

      isize  = 4

      ! get which bin this processor holds
      ibin = modulo(ppiclf_nid,ppiclf_n_bins(1))
      jbin = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
      kbin = 0
      if (ppiclf_ndim .eq. 3)
     >kbin = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))

      il = ibin
      ir = ibin+1
      jl = jbin
      jr = jbin+1
      kl = kbin
      kr = kbin
      if (ppiclf_ndim .eq. 3) then
         kl = kbin
         kr = kbin+1
      endif

      nbinpa = il + ndxgpp1*jl + ndxygpp1*kl
      nbinpb = ir + ndxgpp1*jl + ndxygpp1*kl
      nbinpc = il + ndxgpp1*jr + ndxygpp1*kl
      nbinpd = ir + ndxgpp1*jr + ndxygpp1*kl
      if (ppiclf_ndim .eq. 3) then
         nbinpe = il + ndxgpp1*jl + ndxygpp1*kr
         nbinpf = ir + ndxgpp1*jl + ndxygpp1*kr
         nbinpg = il + ndxgpp1*jr + ndxygpp1*kr
         nbinph = ir + ndxgpp1*jr + ndxygpp1*kr
      endif

      stride_lenv(1) = 0
      stride_lenv(2) = 0
      stride_lenv(3) = 0
      stride_lenv(4) = 0
      stride_lenv(5) = 0
      stride_lenv(6) = 0
      stride_lenv(7) = 0
      stride_lenv(8) = 0
 
      stride_lenc = 0
      if (ppiclf_nid .le. 
     >      ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         stride_lenv(1) = nbinpa
         stride_lenv(2) = nbinpb
         stride_lenv(3) = nbinpc
         stride_lenv(4) = nbinpd
         if (ppiclf_ndim .eq. 3) then
            stride_lenv(5) = nbinpe
            stride_lenv(6) = nbinpf
            stride_lenv(7) = nbinpg
            stride_lenv(8) = nbinph
         endif
         
         stride_lenc    = ppiclf_nid
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') nvtx_total
      write(vtu,'(A)',advance='no') '" NumberOfCells="'
      write(vtu,'(I0)',advance='no') ncll_total
      write(vtu,'(A)',advance='yes') '"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*nvtx_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'
      write(vtu,'(A)',advance='yes') '   </PointData> '
      write(vtu,'(A)',advance='yes') '   <CellData>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"PPR",1,iint)
      iint = iint + 1   *isize*ncll_total + isize
      write(vtu,'(A)',advance='yes') '   </CellData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write connectivity here
      do ii=0,ncll_total-1
         i = modulo(ii,ppiclf_n_bins(1))
         j = modulo(ii/ppiclf_n_bins(1),ppiclf_n_bins(2))
         k = 0
         if (ppiclf_ndim .eq. 3)
     >   k = ii/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
          
c     do K=0,ppiclf_n_bins(3)-1
         kl = K
         kr = K
         if (ppiclf_ndim .eq. 3) then
            kl = K
            kr = K+1
         endif
c     do J=0,ppiclf_n_bins(2)-1
         jl = J
         jr = J+1
c     do I=0,ppiclf_n_bins(1)-1
         il = I
         ir = I+1

         if (ppiclf_ndim .eq. 3) then
            npa = il + ndxgpp1*jl + ndxygpp1*kl
            npb = ir + ndxgpp1*jl + ndxygpp1*kl
            npc = il + ndxgpp1*jr + ndxygpp1*kl
            npd = ir + ndxgpp1*jr + ndxygpp1*kl
            npe = il + ndxgpp1*jl + ndxygpp1*kr
            npf = ir + ndxgpp1*jl + ndxygpp1*kr
            npg = il + ndxgpp1*jr + ndxygpp1*kr
            nph = ir + ndxgpp1*jr + ndxygpp1*kr
            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0,A)',advance='no')  npd, ' '
            write(vtu,'(I0,A)',advance='no')  npe, ' '
            write(vtu,'(I0,A)',advance='no')  npf, ' '
            write(vtu,'(I0,A)',advance='no')  npg, ' '
            write(vtu,'(I0)'  ,advance='yes') nph
         else
            npa = il + ndxgpp1*jl + ndxygpp1*kl
            npb = ir + ndxgpp1*jl + ndxygpp1*kl
            npc = il + ndxgpp1*jr + ndxygpp1*kl
            npd = ir + ndxgpp1*jr + ndxygpp1*kl
            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0)'  ,advance='yes') npd
         endif
c     enddo
c     enddo
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write offsetts here
      ioff_dum = 4
      if (ppiclf_ndim .eq. 3) ioff_dum = 8
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') ioff_dum*i
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="UInt8" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      itype = 8
      if (ppiclf_ndim .eq. 3) itype = 11
      ! write types here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') itype
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

c1511 continue

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif

      call ppiclf_bcast(ivtu_size, isize)

      iorank = -1

      if_pos = 3*isize*nvtx_total


      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif


      call mpi_barrier(ppiclf_comm,ierr)

      ! write points first
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      ! point A
      icount_pos = 0
      if (ppiclf_nid .le. 
     >    ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         icount_pos = 3
      endif
      idisp_pos  = ivtu_size + isize*(3*stride_lenv(1) + 1)
      rpoint(1)  = sngl(ppiclf_binx(1,1))
      rpoint(2)  = sngl(ppiclf_biny(1,1))
      rpoint(3)  = 0.0
      if (ppiclf_ndim .eq. 3)
     >rpoint(3)  = sngl(ppiclf_binz(1,1))
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)

      ! 3d
      if (ppiclf_ndim .gt. 2) then

         ! point B
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(2) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
            rpoint(3)  = sngl(ppiclf_binz(1,1))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point C
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(3) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(1,1))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point E
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(5) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point D
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(4) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(1,1))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point F
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(6) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point G
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(7) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point H
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(8) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = sngl(ppiclf_binz(2,1))
         endif
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)


      ! 2d
      else

         ! point B
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(2) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(1,1))
            rpoint(3)  = 0.0
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point C
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(3) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(1,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = 0.0
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point D
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(4) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_binx(2,1))
            rpoint(2)  = sngl(ppiclf_biny(2,1))
            rpoint(3)  = 0.0
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)

      endif

      call ppiclf_byte_close_mpi(pth,ierr)

      call mpi_barrier(ppiclf_comm,ierr)

      if_cll = 1*isize*ncll_total

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_cll
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write points first
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      !
      ! cell values
      idisp_cll = ivtu_size + isize*(3*(nvtx_total) 
     >     + 1*stride_lenc + 2)
      icount_cll = 0
      if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         icount_cll = 1
      endif
      rpoint(1)  = real(nxx)
      call ppiclf_byte_set_view(idisp_cll,pth)
      call ppiclf_byte_write_mpi(rpoint,icount_cll,iorank,pth,ierr)

      call ppiclf_byte_close_mpi(pth,ierr)

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi('  End WriteBinVTU$',ppiclf_cycle)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteParticleVTU(filein1)
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      real*4  rout_pos(3      *PPICLF_LPART) 
     >       ,rout_sln(PPICLF_LRS*PPICLF_LPART)
     >       ,rout_lrp(PPICLF_LRP*PPICLF_LPART)
     >       ,rout_lip(3      *PPICLF_LPART)
      character*3 filein
      character*12 vtufile
      character*6  prostr
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth,prevs(2,ppiclf_np)
      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip,stride_len
      integer*4 iint, nnp, nxx, npt_total, jx, jy, jz, if_sz, isize,
     >          iadd, if_pos, if_sln, if_lrp, if_lip, ic_pos, ic_sln,
     >          ic_lrp, ic_lip, i, j, ie, nps, nglob, nkey, ndum,
     >          icount_pos, icount_sln, icount_lrp, icount_lip, iorank,
     >          ierr, ivtu_size
      integer*4 ppiclf_iglsum
      external ppiclf_iglsum
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!

      call ppiclf_printsi(' *Begin WriteParticleVTU$',ppiclf_cycle)

      if (icalld1 .eq. 0) icalld1 = istartout

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      npt_total = ppiclf_iglsum(nxx,1)

      jx    = 1
      jy    = 2
      jz    = 1
      if (ppiclf_ndim .eq. 3)
     >jz    = 3

      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'par'
      else 
         write(filein,'(A3)') filein1
      endif

! --------------------------------------------------
! COPY PARTICLES TO OUTPUT ARRAY
! --------------------------------------------------

      isize = 4

      iadd = 0
      if_pos = 3*isize*npt_total
      if_sln = 1*isize*npt_total
      if_lrp = 1*isize*npt_total
      if_lip = 1*isize*npt_total

      ic_pos = iadd
      ic_sln = iadd
      ic_lrp = iadd
      ic_lip = iadd
      do i=1,nxx

         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(ppiclf_y(jx,i))
         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(ppiclf_y(jy,i))
         ic_pos = ic_pos + 1
         if (ppiclf_ndim .eq. 3) then
            rout_pos(ic_pos) = sngl(ppiclf_y(jz,i))
         else
            rout_pos(ic_pos) = 0.0
         endif
      enddo
      do j=1,PPICLF_LRS
      do i=1,nxx
         ic_sln = ic_sln + 1
         rout_sln(ic_sln) = sngl(ppiclf_y(j,i))
      enddo
      enddo
      do j=1,PPICLF_LRP
      do i=1,nxx
         ic_lrp = ic_lrp + 1
         rout_lrp(ic_lrp) = sngl(ppiclf_rprop(j,i))
      enddo
      enddo
      do j=5,7
      do i=1,nxx
         ic_lip = ic_lip + 1
         rout_lip(ic_lip) = real(ppiclf_iprop(j,i))
      enddo
      enddo

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      do i=1,nnp
         prevs(1,i) = i-1
         prevs(2,i) = nxx
      enddo

      nps   = 1 ! index of new proc for doing stuff
      nglob = 1 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      ndum = 2
      call pfgslib_crystal_ituple_transfer(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nnp,nps)
      call pfgslib_crystal_ituple_sort(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nglob,nkey)

      stride_len = 0
      if (ppiclf_nid .ne. 0) then
      do i=1,ppiclf_nid
         stride_len = stride_len + prevs(2,i)
      enddo
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') npt_total
      write(vtu,'(A)',advance='yes') '" NumberOfCells="0"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*npt_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'


      do ie=1,PPICLF_LRS
         write(prostr,'(A1,I2.2)') "y",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      do ie=1,PPICLF_LRP
         write(prostr,'(A4,I2.2)') "rprop",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      do ie=1,3
         write(prostr,'(A3,I2.2)') "tag",ie
         call ppiclf_io_WriteDataArrayVTU(vtu,prostr,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      write(vtu,'(A)',advance='yes') '   </PointData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif

      call ppiclf_bcast(ivtu_size, isize)

      ! byte-displacements
      idisp_pos = ivtu_size + isize*(3*stride_len + 1)

      ! how much to write
      icount_pos = 3*nxx
      icount_sln = 1*nxx
      icount_lrp = 1*nxx
      icount_lip = 1*nxx

      iorank = -1

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_write_mpi(rout_pos,icount_pos,iorank,pth,ierr)
      call ppiclf_byte_close_mpi(pth,ierr)

      call mpi_barrier(ppiclf_comm,ierr)

      do i=1,PPICLF_LRS
         idisp_sln = ivtu_size + isize*(3*npt_total 
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + i)

         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_sln
           close(vtu)
         endif
   
         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*ppiclf_npart + 1
   
         ! write
         call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
         call ppiclf_byte_set_view(idisp_sln,pth)
         call ppiclf_byte_write_mpi(rout_sln(j),icount_sln,iorank,pth
     >                             ,ierr)
         call ppiclf_byte_close_mpi(pth,ierr)
      enddo

      do i=1,PPICLF_LRP
         idisp_lrp = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + i)

         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_lrp
           close(vtu)
         endif
   
         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*ppiclf_npart + 1
   
         ! write
         call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
         call ppiclf_byte_set_view(idisp_lrp,pth)
         call ppiclf_byte_write_mpi(rout_lrp(j),icount_lrp,iorank,pth
     >                             ,ierr)
         call ppiclf_byte_close_mpi(pth,ierr)
      enddo

      do i=1,3
         idisp_lip = ivtu_size + isize*(3*npt_total  
     >                         + PPICLF_LRS*npt_total
     >                         + PPICLF_LRP*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + PPICLF_LRS + PPICLF_LRP + i)
         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_lip
           close(vtu)
         endif

         call mpi_barrier(ppiclf_comm,ierr)

         j = (i-1)*ppiclf_npart + 1
   
         ! write
         call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
         call ppiclf_byte_set_view(idisp_lip,pth)
         call ppiclf_byte_write_mpi(rout_lip(j),icount_lip,iorank,pth
     >                             ,ierr)
         call ppiclf_byte_close_mpi(pth,ierr)
      enddo

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi(' *End WriteParticleVTU$',ppiclf_cycle)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteDataArrayVTU(vtu,dataname,ncomp,idist)
!
      implicit none
!
! Input:
!
      integer*4 vtu,ncomp
      integer*4 idist
      character (len = *) dataname
!
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="'
      write(vtu,'(A)',advance='no') dataname
      write(vtu,'(A)',advance='no') '" NumberOfComponents="'
      write(vtu,'(I0)',advance='no') ncomp
      write(vtu,'(A)',advance='no') '" format="append" '
      write(vtu,'(A)',advance='no') 'offset="'
      write(vtu,'(I0)',advance='no') idist
      write(vtu,'(A)',advance='yes') '"/>'

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagAll
!
      implicit none
!
      include "PPICLF"
!
      call ppiclf_prints('*********** PPICLF OUTPUT *****************$')
      call ppiclf_io_OutputDiagGen
      call ppiclf_io_OutputDiagGhost
      if (ppiclf_lsubbin) call ppiclf_io_OutputDiagSubBin
      if (ppiclf_overlap) call ppiclf_io_OutputDiagGrid

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGen
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 npart_max, npart_min, npart_tot, npart_ide, nbin_total
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
      call ppiclf_prints(' *Begin General Info$')
         npart_max = ppiclf_iglmax(ppiclf_npart,1)
         npart_min = ppiclf_iglmin(ppiclf_npart,1)
         npart_tot = ppiclf_iglsum(ppiclf_npart,1)
         npart_ide = npart_tot/ppiclf_np

         nbin_total = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)

      call ppiclf_printsi('  -Cycle                  :$',ppiclf_cycle)
      call ppiclf_printsi('  -Output Freq.           :$',ppiclf_iostep)
      call ppiclf_printsr('  -Time                   :$',ppiclf_time)
      call ppiclf_printsr('  -dt                     :$',ppiclf_dt)
      call ppiclf_printsi('  -Global particles       :$',npart_tot)
      call ppiclf_printsi('  -Local particles (Max)  :$',npart_max)
      call ppiclf_printsi('  -Local particles (Min)  :$',npart_min)
      call ppiclf_printsi('  -Local particles (Ideal):$',npart_ide)
      call ppiclf_printsi('  -Total ranks            :$',ppiclf_np)
      call ppiclf_printsi('  -Problem dimensions     :$',ppiclf_ndim)
      call ppiclf_printsi('  -Integration method     :$',ppiclf_imethod)
      call ppiclf_printsi('  -Number of bins total   :$',nbin_total)
      call ppiclf_printsi('  -Number of bins (x)     :$',
     >                    ppiclf_n_bins(1))
      call ppiclf_printsi('  -Number of bins (y)     :$'
     >                    ,ppiclf_n_bins(2))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsi('  -Number of bins (z)     :$'
     >                    ,ppiclf_n_bins(3))
      call ppiclf_printsr('  -Bin xl coordinate      :$',ppiclf_binb(1))
      call ppiclf_printsr('  -Bin xr coordinate      :$',ppiclf_binb(2))
      call ppiclf_printsr('  -Bin yl coordinate      :$',ppiclf_binb(3))
      call ppiclf_printsr('  -Bin yr coordinate      :$',ppiclf_binb(4))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin zl coordinate      :$',ppiclf_binb(5))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin zr coordinate      :$',ppiclf_binb(6))

      call ppiclf_prints('  End General Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGrid
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 nel_max_orig, nel_min_orig, nel_total_orig, nel_max_map,
     >          nel_min_map, nel_total_map
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
      call ppiclf_prints(' *Begin Grid Info$')

         nel_max_orig   = ppiclf_iglmax(ppiclf_nee,1)
         nel_min_orig   = ppiclf_iglmin(ppiclf_nee,1)
         nel_total_orig = ppiclf_iglsum(ppiclf_nee,1)

         nel_max_map   = ppiclf_iglmax(ppiclf_neltb,1)
         nel_min_map   = ppiclf_iglmin(ppiclf_neltb,1)
         nel_total_map = ppiclf_iglsum(ppiclf_neltb,1)

      call ppiclf_printsi('  -Orig. Global cells     :$',nel_total_orig)
      call ppiclf_printsi('  -Orig. Local cells (Max):$',nel_max_orig)
      call ppiclf_printsi('  -Orig. Local cells (Min):$',nel_min_orig)
      call ppiclf_printsi('  -Map Global cells       :$',nel_total_map)
      call ppiclf_printsi('  -Map Local cells (Max)  :$',nel_max_map)
      call ppiclf_printsi('  -Map Local cells (Min)  :$',nel_min_map)

      call ppiclf_prints('  End Grid Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGhost
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 npart_max, npart_min, npart_tot
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
      call ppiclf_prints(' *Begin Ghost Info$')

         npart_max = ppiclf_iglmax(ppiclf_npart_gp,1)
         npart_min = ppiclf_iglmin(ppiclf_npart_gp,1)
         npart_tot = ppiclf_iglsum(ppiclf_npart_gp,1)

      call ppiclf_printsi('  -Global ghosts          :$',npart_tot)
      call ppiclf_printsi('  -Local ghosts (Max)     :$',npart_max)
      call ppiclf_printsi('  -Local ghosts (Min)     :$',npart_min)

      call ppiclf_prints('  End Ghost Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagSubBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 nbin_total
!

      call ppiclf_prints(' *Begin SubBin Info$')

      nbin_total = ppiclf_bx*ppiclf_by*ppiclf_bz

      call ppiclf_printsi('  -Number of local bin    :$',nbin_total)
      call ppiclf_printsi('  -Number of local bin (x):$',PPICLF_BX)
      call ppiclf_printsi('  -Number of local bin (y):$',PPICLF_BY)
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsi('  -Number of local bin (z):$',PPICLF_BZ)
      call ppiclf_printsr('  -Bin width (x)          :$',ppiclf_rdx)
      call ppiclf_printsr('  -Bin width (y)          :$',ppiclf_rdy)
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin width (z)          :$',ppiclf_rdz)
      call ppiclf_printsr('  -Filter width           :$',ppiclf_filter)
      call ppiclf_printsr('  -Filter cut-off         :$'
     >                                                 ,ppiclf_d2chk(2))
      call ppiclf_printsi('  -SubBins per filter res.:$',ppiclf_ngrids)

      call ppiclf_prints('  End SubBin Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_Initialize(xi1,xpmin,xpmax,
     >        yi1,ypmin,ypmax,zi1,zpmin,zpmax,
     >        ai1,apa,apxa,aprin,aprout)
!
      implicit none
!
      include "PPICLF"
!
      integer*4 xi1, yi1, zi1, ai1
      real*8 xpmin,xpmax,ypmin,ypmax,zpmin,zpmax,
     >       apa,apxa,aprin,aprout
      real*8 pi, angled

      ! Linear X-Periodicity
      x_per_flag = xi1
      if(x_per_flag.eq.1) then
        if (xpmin .ge. xpmax) call ppiclf_exittr('PeriodicX 
     >      must have xmin < xmax$',xpmin,0)
        ppiclf_iperiodic(1) = 0
        x_per_min = xpmin
        x_per_max = xpmax
        ppiclf_xdrange(1,1) = xpmin
        ppiclf_xdrange(2,1) = xpmax
      endif

      ! Linear Y-Periodicity
      y_per_flag = yi1
      if(y_per_flag.eq.1) then
        if (ypmin .ge. ypmax) call ppiclf_exittr('PeriodicY 
     >     must have ymin < ymax$',ypmin,0)
        ppiclf_iperiodic(2) = 0
        y_per_min = ypmin
        y_per_max = ypmax
        ppiclf_xdrange(1,2) = ypmin
        ppiclf_xdrange(2,2) = ypmax
      endif

      ! Linear Z-Periodicity
      z_per_flag = zi1
      if(z_per_flag.eq.1) then
        if (zpmin .ge. zpmax) call ppiclf_exittr('PeriodicZ 
     >     must have zmin < zmax$',zpmin,0)
        ppiclf_iperiodic(3) = 0
        z_per_min = zpmin
        z_per_max = zpmax
        ppiclf_xdrange(1,3) = zpmin
        ppiclf_xdrange(2,3) = zpmax
      endif


      ! Angular Periodicity
      ang_per_flag = ai1
      if(ang_per_flag.eq.1) then
        ppiclf_iperiodic(1) = 0 ! X-Periodicity
        ppiclf_iperiodic(2) = 0 ! Y-Periodicity
        ang_per_angle  = apa
        ang_per_xangle = apxa
        ang_per_rin    = aprin
        ang_per_rout   = aprout
      endif

      ! User cannot initialize X/Y-Periodicity with Angular Periodicity
      if((x_per_flag.eq.1).or.(y_per_flag.eq.1).and.(ang_per_flag.eq.1))
     >   call ppiclf_exittr('PPICLF: Invalid Periodicity choice$',0,0)

      ! Thierry - compute ang_case

      pi = acos(-1.0)
      angled = ang_per_angle * 180.0d0 / pi ! store angle value in degrees

      if (ang_per_flag==0) then
         ang_case = 0 ! standard geometry
      else
         if (angled .lt. 90.0)        ang_case = 1 ! general wedge
         if (nint(angled) .eq. 90.0)  ang_case = 2 ! quarter cylinder
         if (nint(angled) .eq. 180.0) ang_case = 3 ! half cylinder
      endif

      if (ppiclf_nid==0) then
         print*, " "
         print*, " ======================================="
         print*, " "
         print*, "!!! PPICLF Angular Periodicity Initialized !!!!"
         print*, "  Angular periodicity flag =", ang_per_flag
         if (ang_per_flag==0) then
            print*, "  Init Angular- ang_case =", ang_case
         else
            print*, "  Init Angular- angle =", ang_per_angle
            print*, "  Init Angular- angled =", angled
            print*, "  Init Angular- nint(angled) =", nint(angled)
            print*, "  Init Angular- ang_case =", ang_case
         endif
         print*, " "
         print*, " ======================================="
         print*, " "
      endif

      return
      end
!
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_AddParticles(npart,y,rprop)
     > bind(C, name="ppiclc_solve_AddParticles")
#else
      subroutine ppiclf_solve_AddParticles(npart,y,rprop)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input: 
!
      integer*4  npart
      real*8     y(*)
      real*8     rprop(*)
!
! Internal:
!
      integer*4 ppiclf_iglsum,ntotal
      external ppiclf_iglsum
!

      call ppiclf_prints('   *Begin AddParticles$')

      if (ppiclf_npart+npart .gt. PPICLF_LPART .or. npart .lt. 0)
     >   call ppiclf_exittr('Invalid number of particles$',
     >                      0.0,ppiclf_npart+npart)

      call ppiclf_printsi('      -Begin copy particles$',npart)

      ! First, append arrays onto existing arrays
      call ppiclf_copy(ppiclf_y(1,ppiclf_npart+1),
     >                 y,
     >                 npart*PPICLF_LRS)
      call ppiclf_copy(ppiclf_rprop(1,ppiclf_npart+1),
     >                 rprop,
     >                 npart*PPICLF_LRP)
      ppiclf_npart = ppiclf_npart + npart

      call ppiclf_printsi('      -Begin copy particles$',ppiclf_npart)

      if (.not. PPICLF_RESTART) then
         call ppiclf_prints('      -Begin ParticleTag$')
            call ppiclf_solve_SetParticleTag(npart)
         call ppiclf_prints('       End ParticleTag$')
      endif

      if (ppiclf_iglsum(ppiclf_npart,1).gt.0) then
         call ppiclf_prints('      -Begin CreateBin$')
            call ppiclf_comm_CreateBin
         call ppiclf_prints('       End CreateBin$')

         call ppiclf_prints('      -Begin FindParticle$')
            call ppiclf_comm_FindParticle
         call ppiclf_prints('       End FindParticle$')

         call ppiclf_prints('      -Begin MoveParticle$')
            call ppiclf_comm_MoveParticle
         call ppiclf_prints('       End MoveParticle$')

      endif

      call ppiclf_prints('    End AddParticles$')

      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y,
     >                                     rprop,filt2,filt3)
     > bind(C, name="ppiclc_solve_InitParticle")
#else
      subroutine ppiclf_solve_InitParticle(imethod,ndim,iendian,npart,y,
     >                                     rprop,filt2,filt3)
#endif
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input: 
!
      integer*4  imethod
      integer*4  ndim
      integer*4  iendian
      integer*4  npart
      integer*4  ierr
      real*8     y(*)
      real*8     rprop(*)
      real*8     filt2,filt3
!
      call mpi_barrier(ppiclf_comm,ierr)
!
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitParticle$',0.0d0
     >   ,ppiclf_nid)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter must be before InitParticle$',0.0d0
     >                  ,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.0d0
     >                  ,0)

      call ppiclf_prints('*Begin InitParticle$')

         call ppiclf_prints('   *Begin InitParam$')
            call ppiclf_solve_InitParam(imethod,ndim,iendian)
            ppiclf_d2chk(2) = filt2
            ppiclf_d2chk(3) = filt3
            ppiclf_d2chk(1) = max( ppiclf_d2chk(2),ppiclf_d2chk(3) )

            ! TLJ added 12/21/2024
            if (ppiclf_nid==0) then
               print*,'TLJ checking d2chk(2) = ',ppiclf_d2chk(2)
               print*,'TLJ checking d2chk(3) = ',ppiclf_d2chk(3)
               print*,'TLJ checking d2chk(1) = ',ppiclf_d2chk(1)
            endif

         call ppiclf_prints('    End InitParam$')

         if (.not. PPICLF_RESTART) then
            call ppiclf_prints('   *Begin InitZero$')
               call ppiclf_solve_InitZero
            call ppiclf_prints('   *End InitZero$')
            call ppiclf_prints('   *Begin AddParticles$')
               call ppiclf_solve_AddParticles(npart,y,rprop)
            call ppiclf_prints('   *End AddParticles$')

            ! TLJ - 11/23/2024
            ! The write files at t=0 has been moved to a single
            !   location in ppiclf_solve_WriteVTK
            !call ppiclf_prints('   *Begin WriteParticleVTU$')
            !   call ppiclf_io_WriteParticleVTU('')
            !call ppiclf_prints('    End WriteParticleVTU$')
            !call ppiclf_prints('   *Begin WriteBinVTU$')
            !   call ppiclf_io_WriteBinVTU('')
            !call ppiclf_prints('    End WriteBinVTU$')
         endif

      call ppiclf_prints(' End InitParticle$')
!
      call mpi_barrier(ppiclf_comm,ierr)
!

      ! This prints out initial bin information
      call ppiclf_io_OutputDiagGen

      PPICLF_LINIT = .true.

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitParam(imethod,ndim,iendian)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4  imethod
      integer*4  ndim
      integer*4  iendian
!
      if (imethod .eq. 0 .or. imethod .ge. 3 .or. imethod .le. -2)
     >   call ppiclf_exittr('Invalid integration method$',0.0d0,imethod)
      if (ndim .le. 1 .or. ndim .ge. 4)
     >   call ppiclf_exittr('Invalid problem dimension$',0.0d0,ndim)
      if (iendian .lt. 0 .or. iendian .gt. 1)
     >   call ppiclf_exittr('Invalid Endian$',0.0d0,iendian)

      ppiclf_imethod      = imethod
      ppiclf_ndim         = ndim
      ppiclf_iendian      = iendian

      ppiclf_filter       = -1   ! filt default for no projection

      ppiclf_iperiodic(1) = 1    ! periodic in x (== 0) 
      ppiclf_iperiodic(2) = 1    ! periodic in y (==0)
      ppiclf_iperiodic(3) = 1    ! periodic in z (==0)

      ppiclf_cycle  = 0
      ppiclf_iostep = 1
      ppiclf_dt     = 0.0d0
      ppiclf_time   = 0.0d0

      ppiclf_overlap    = .false.
      ppiclf_linit      = .false.
      ppiclf_lfilt      = .false.
      ppiclf_lfiltgauss = .false.
      ppiclf_lfiltbox   = .false.
      ppiclf_lintp      = .false.
      ppiclf_lproj      = .false.
      ppiclf_lsubbin    = .false.
      ppiclf_lsubsubbin = .false.
      if (PPICLF_INTERP .eq. 1)  ppiclf_lintp = .true.
      if (PPICLF_PROJECT .eq. 1) ppiclf_lproj = .true.

      ppiclf_xdrange(1,1) = -1E20
      ppiclf_xdrange(2,1) =  1E20
      ppiclf_xdrange(1,2) = -1E20
      ppiclf_xdrange(2,2) =  1E20
      ppiclf_xdrange(1,3) = -1E20
      ppiclf_xdrange(2,3) =  1E20

      ppiclf_d2chk(1) = 0.0d0
      ppiclf_d2chk(2) = 0.0d0
      ppiclf_d2chk(3) = 0.0d0

      ppiclf_n_bins(1) = 1
      ppiclf_n_bins(2) = 1
      ppiclf_n_bins(3) = 1

      ppiclf_bins_set(1) = 0
      ppiclf_bins_set(2) = 0
      ppiclf_bins_set(3) = 0

      ppiclf_bins_balance(1) = 0
      ppiclf_bins_balance(2) = 0
      ppiclf_bins_balance(3) = 0

      ppiclf_nwall    = 0
      ppiclf_iwallm   = 0

      PPICLF_INT_ICNT = 0

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitNeighborBin(rwidth)
     > bind(C, name="ppiclc_solve_InitNeighborBin")
#else
      subroutine ppiclf_solve_InitNeighborBin(rwidth)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      real*8 rwidth
!
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitNeighborBin$',0.0d0
     >                  ,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitNeighborBin$'
     >                  ,0.0d0,0)

      ppiclf_lsubsubbin = .true.

      ppiclf_d2chk(3) = rwidth

      ! TLJ added 12/21/2024
      if (ppiclf_nid==0) then
         print*,'TLJ checking d2chk(3) = ',ppiclf_d2chk(3)
      endif

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitTargetBins(str,n,balance)
     > bind(C, name="ppiclc_solve_InitTargetBins")
#else
      subroutine ppiclf_solve_InitTargetBins(str,n,balance)
#endif
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      character*1 str
      integer*4 n
      integer*4 balance
!
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitTargetBins$'
     >                   ,0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitTargetBins$'
     >                  ,0.0d0,0)

      if (str == 'x' .or. str == 'X') then 
         ppiclf_n_bins(1) = n
         if (n .gt. 1) ppiclf_bins_set(1) = 1
         ppiclf_bins_balance(1) = balance
      elseif (str == 'y' .or. str == 'Y') then 
         ppiclf_n_bins(2) = n
         if (n .gt. 1) ppiclf_bins_set(2) = 1
         ppiclf_bins_balance(2) = balance
      elseif (str == 'z' .or. str == 'Z') then 
        if (ppiclf_ndim .lt. 3)
     >   call ppiclf_exittr('Dim must be 3 to use InitTargetBins on z$'
     >                   ,0.,ppiclf_ndim)
         ppiclf_n_bins(3) = n
         if (n .gt. 1) ppiclf_bins_set(3) = 1
         ppiclf_bins_balance(3) = balance
      endif

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SetNeighborBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i
!
      do i=1,ppiclf_npart
         ppiclf_nb_r(1,i) = floor((ppiclf_cp_map(1,i)-ppiclf_binb(1))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_r(2,i) = floor((ppiclf_cp_map(2,i)-ppiclf_binb(3))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_r(3,i) = 0
         if (ppiclf_ndim .eq. 3)
     >   ppiclf_nb_r(3,i) = floor((ppiclf_cp_map(3,i)-ppiclf_binb(5))/
     >                             ppiclf_d2chk(3))
      enddo

      do i=1,ppiclf_npart_gp
         ppiclf_nb_g(1,i) = floor((ppiclf_rprop_gp(1,i)-ppiclf_binb(1))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_g(2,i) = floor((ppiclf_rprop_gp(2,i)-ppiclf_binb(3))/
     >                             ppiclf_d2chk(3))
         ppiclf_nb_g(3,i) = 0
         if (ppiclf_ndim .eq. 3)
     >   ppiclf_nb_g(3,i) = floor((ppiclf_rprop_gp(3,i)-ppiclf_binb(5))/
     >                             ppiclf_d2chk(3))
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitZero
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i,j,ic,k,ie
!
      ic = 0
      do i=1,PPICLF_LPART
      do j=1,PPICLF_LRS
         ic = ic + 1
         ppiclf_y    (j,i) = 0.0d0
         ppiclf_ydot (j,i) = 0.0d0
         ppiclf_ydotc(j,i) = 0.0d0
         ppiclf_y1   (ic ) = 0.0d0
      enddo
      do j=1,PPICLF_LRP
         ppiclf_rprop(j,i) = 0.0d0
      enddo
      do j=1,PPICLF_LRP2
         ppiclf_rprop2(j,i) = 0.0d0
      enddo
      do j=1,PPICLF_LRP3
         ppiclf_rprop3(j,i) = 0.0d0
      enddo
      do j=1,PPICLF_LIP
         ppiclf_iprop(j,i) = 0
      enddo
      enddo
      ppiclf_npart = 0

      do ie=1,PPICLF_LEE
      do ic=1,PPICLF_LRP_INT
      do k=1,PPICLF_LEZ
      do j=1,PPICLF_LEY
      do i=1,PPICLF_LEX
        ppiclf_int_fld(i,j,k,ic,ie) = 0.0d0
      enddo
      enddo
      enddo
      enddo
      enddo

      !!call ppiclf_user_InitZero

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_NearestNeighbor(i)
!
      implicit none
!
      include "PPICLF"
! 
! Input:
! 
      integer*4 i
! 
! Internal: 
! 
      real*8 ydum(PPICLF_LRS), rpropdum(PPICLF_LRP)
      real*8 A(3),B(3),C(3),AB(3),AC(3), dist2, xdist2, ydist2,
     >       dist_total
      integer*4 i_iim, i_iip, i_jjm, i_jjp, i_kkm, i_kkp, j, j_ii, j_jj,
     >          j_kk, jp
      real*8 rnx, rny, rnz, area, rpx1, rpy1, rpz1, rpx2, rpy2, rpz2,
     >       rflip, a_sum, rd, rdist, theta, tri_area, rthresh,
     >       ab_dot_ac, ab_mag, ac_mag, zdist2
      integer*4 istride, k, kmax, kp, kkp, kk
! 
      i_iim = ppiclf_nb_r(1,i) - 1
      i_iip = ppiclf_nb_r(1,i) + 1
      i_jjm = ppiclf_nb_r(2,i) - 1
      i_jjp = ppiclf_nb_r(2,i) + 1
      i_kkm = ppiclf_nb_r(3,i) - 1
      i_kkp = ppiclf_nb_r(3,i) + 1

      dist2 = ppiclf_d2chk(3)**2

      do j=1,ppiclf_npart
         if (j .eq. i) cycle

         j_ii = ppiclf_nb_r(1,j)
         j_jj = ppiclf_nb_r(2,j)
         j_kk = ppiclf_nb_r(3,j)

         if (j_ii .gt. i_iip .or. j_ii .lt. i_iim) cycle
         if (j_jj .gt. i_jjp .or. j_jj .lt. i_jjm) cycle
         if (ppiclf_ndim .eq. 3) then
         if (j_kk .gt. i_kkp .or. j_kk .lt. i_kkm) cycle
         endif

         xdist2 = (ppiclf_cp_map(1,i)-ppiclf_cp_map(1,j))**2
         if (xdist2 .gt. dist2) cycle
         ydist2 = (ppiclf_cp_map(2,i)-ppiclf_cp_map(2,j))**2
         if (ydist2 .gt. dist2) cycle
         dist_total = xdist2 + ydist2
         if (ppiclf_ndim .eq. 3) then
         zdist2 = (ppiclf_cp_map(3,i)-ppiclf_cp_map(3,j))**2
         if (zdist2 .gt. dist2) cycle
         dist_total = dist_total+zdist2
         endif
         if (dist_total .gt. dist2) cycle

         call ppiclf_user_EvalNearestNeighbor(i,j,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ppiclf_cp_map(1,j)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,j))

      enddo

      do j=1,ppiclf_npart_gp
         j_ii = ppiclf_nb_g(1,j)
         j_jj = ppiclf_nb_g(2,j)
         j_kk = ppiclf_nb_g(3,j)

         if (j_ii .gt. i_iip .or. j_ii .lt. i_iim) cycle
         if (j_jj .gt. i_jjp .or. j_jj .lt. i_jjm) cycle
         if (ppiclf_ndim .eq. 3) then
         if (j_kk .gt. i_kkp .or. j_kk .lt. i_kkm) cycle
         endif

         xdist2 = (ppiclf_cp_map(1,i)-ppiclf_rprop_gp(1,j))**2
         if (xdist2 .gt. dist2) cycle
         ydist2 = (ppiclf_cp_map(2,i)-ppiclf_rprop_gp(2,j))**2
         if (ydist2 .gt. dist2) cycle
         dist_total = xdist2 + ydist2
         if (ppiclf_ndim .eq. 3) then
         zdist2 = (ppiclf_cp_map(3,i)-ppiclf_rprop_gp(3,j))**2
         if (zdist2 .gt. dist2) cycle
         dist_total = dist_total+zdist2
         endif
         if (dist_total .gt. dist2) cycle

         jp = -1*j
         call ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ppiclf_rprop_gp(1,j)
     >                                 ,ppiclf_rprop_gp(1+PPICLF_LRS,j))

      enddo

      istride = ppiclf_ndim
      do j=1,ppiclf_nwall

         rnx  = ppiclf_wall_n(1,j)
         rny  = ppiclf_wall_n(2,j)
         rnz  = 0.0d0
         area = ppiclf_wall_n(3,j)
         rpx1 = ppiclf_cp_map(1,i)
         rpy1 = ppiclf_cp_map(2,i)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(1,j)
         rpy2 = ppiclf_wall_c(2,j)
         rpz2 = 0.0d0
         rpx2 = rpx2 - rpx1
         rpy2 = rpy2 - rpy1

         if (ppiclf_ndim .eq. 3) then
            rnz  = ppiclf_wall_n(3,j)
            area = ppiclf_wall_n(4,j)
            rpz1 = ppiclf_cp_map(3,i)
            rpz2 = ppiclf_wall_c(3,j)
            rpz2 = rpz2 - rpz1
         endif
    
         rflip = rnx*rpx2 + rny*rpy2 + rnz*rpz2
         if (rflip .gt. 0.0d0) then
            rnx = -1.0d0*rnx
            rny = -1.0d0*rny
            rnz = -1.0d0*rnz
         endif


         a_sum = 0.0d0
         kmax = 2
         if (ppiclf_ndim .eq. 3) kmax = 3
         do k=1,kmax 
            kp = k+1
            if (kp .gt. kmax) kp = kp-kmax ! cycle
            
            kk   = istride*(k-1)
            kkp  = istride*(kp-1)
            rpx1 = ppiclf_wall_c(kk+1,j)
            rpy1 = ppiclf_wall_c(kk+2,j)
            rpz1 = 0.0d0
            rpx2 = ppiclf_wall_c(kkp+1,j)
            rpy2 = ppiclf_wall_c(kkp+2,j)
            rpz2 = 0.0d0

            if (ppiclf_ndim .eq. 3) then
               rpz1 = ppiclf_wall_c(kk+3,j)
               rpz2 = ppiclf_wall_c(kkp+3,j)
            endif

            rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)

            rdist = abs(rnx*ppiclf_cp_map(1,i)+rny*ppiclf_cp_map(2,i)
     >                 +rnz*ppiclf_cp_map(3,i)+rd)
            rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

            ! give a little extra room for walls (2x)
            if (rdist .gt. 2.0d0*ppiclf_d2chk(3)) goto 1511

            ydum(1) = ppiclf_cp_map(1,i) - rdist*rnx
            ydum(2) = ppiclf_cp_map(2,i) - rdist*rny
            ydum(3) = 0.0d0

            A(1) = ydum(1)
            A(2) = ydum(2)
            A(3) = 0.0d0

            B(1) = rpx1
            B(2) = rpy1
            B(3) = 0.0d0

            C(1) = rpx2
            C(2) = rpy2
            C(3) = 0.0d0

            AB(1) = B(1) - A(1)
            AB(2) = B(2) - A(2)
            AB(3) = 0.0d0

            AC(1) = C(1) - A(1)
            AC(2) = C(2) - A(2)
            AC(3) = 0.0d0

            if (ppiclf_ndim .eq. 3) then
               ydum(3) = ppiclf_cp_map(3,i) - rdist*rnz
               A(3) = ydum(3)
               B(3) = rpz1
               C(3) = rpz2
               AB(3) = B(3) - A(3)
               AC(3) = C(3) - A(3)

               AB_DOT_AC = AB(1)*AC(1) + AB(2)*AC(2) + AB(3)*AC(3)
               AB_MAG = sqrt(AB(1)**2 + AB(2)**2 + AB(3)**2)
               AC_MAG = sqrt(AC(1)**2 + AC(2)**2 + AC(3)**2)
               theta  = acos(AB_DOT_AC/(AB_MAG*AC_MAG))
               tri_area = 0.5d0*AB_MAG*AC_MAG*sin(theta)
            elseif (ppiclf_ndim .eq. 2) then
               AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
               tri_area = AB_MAG
            endif
            a_sum = a_sum + tri_area
         enddo

         rthresh = 1.10d0 ! keep it from slipping through crack on edges
         if (a_sum .gt. rthresh*area) cycle

         jp = 0
         call ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ydum
     >                                 ,rpropdum)

 1511 continue
      enddo

      return
      end
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_NearestNeighborSB(i,SBt,SBc,SBm,SBn,iB)
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
! 
! Input:
!
      INTEGER*4 i, SBt, SBn(3), iB(3)  
      INTEGER*4 SBc(0:(SBt-1)),
     >  SBm(0:(SBt-1),(ppiclf_npart+ppiclf_npart_gp))
! 
! Internal: 
! 
      REAL*8 ydum(PPICLF_LRS), rpropdum(PPICLF_LRP), xp(3), bin_xMin(3)
      REAL*8 A(3),B(3),C(3),AB(3),AC(3), dist2, xdist2, ydist2,
     >       dist_total
      INTEGER*4 j, jp, l, iSB, jSB, kSB, loopSB, tempSB, iSBin(3)
      REAL*8 rnx, rny, rnz, area, rpx1, rpy1, rpz1, rpx2, rpy2, rpz2,
     >       rflip, a_sum, rd, rdist, theta, tri_area, rthresh,
     >       ab_dot_ac, ab_mag, ac_mag, zdist2
      INTEGER*4 istride, k, kmax, kp, kkp, kk
! 
      dist2 = ppiclf_d2chk(3)**2
      
      ! find ith particle subbin (tempSB)
      DO l = 1,3
        IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
          xp(l) = ppiclf_y(l,i)
          bin_xMin(l) = ppiclf_binb(2*l-1)+iB(l)*ppiclf_bins_dx(l)
        ELSE
          xp(l) = 0.0
        END IF
      END DO
      DO l = 1,3
        IF (l .LT. 3 .OR. ppiclf_ndim .GT. 2) THEN
          iSBin(l) = FLOOR((xp(l) - (bin_xMin(l)
     >         - ppiclf_d2chk(3)))/ppiclf_d2chk(3))
        ELSE
          iSBin(l) = 0
        END IF
      END DO
      tempSB = iSBin(1) + iSBin(2)*SBn(1) + iSBin(3)*SBn(1)*SBn(2)
      
      ! Loop through real particles
      DO iSB = 1,3     !to look at -1,current,+1 x-dir subbins
        DO jSB = 1,3   !to look at -1,current,+1 x-dir subbins
          DO kSB = 1,3 !to look at -1,current,+1 x-dir subbins
          ! Loops through 27 adjacent subbins
          loopSB = tempSB + (-2+iSB) + (-2+jSB)*SBn(1) 
     >             + (-2+kSB)*SBn(1)*SBn(2)
          IF (loopSB .GT. -1 .AND. loopSB .LT. SBt) THEN
            DO k = 1,SBc(loopSB) 
              j = SBm(loopSB,k)
              IF (j .GT. 0) THEN ! Real particle
                IF (j .EQ. i) CYCLE
                xdist2 = (ppiclf_cp_map(1,i)-ppiclf_cp_map(1,j))**2
                IF (xdist2 .GT. dist2) CYCLE
                ydist2 = (ppiclf_cp_map(2,i)-ppiclf_cp_map(2,j))**2
                IF (ydist2 .GT. dist2) CYCLE
                dist_total = xdist2 + ydist2
                IF (ppiclf_ndim .EQ. 3) THEN
                  zdist2 = (ppiclf_cp_map(3,i)-ppiclf_cp_map(3,j))**2
                  IF (zdist2 .GT. dist2) CYCLE
                  dist_total = dist_total+zdist2
                END IF
                IF (dist_total .GT. dist2) CYCLE
                CALL ppiclf_user_EvalNearestNeighbor(i,j
     >                                   ,ppiclf_cp_map(1,i)
     >                                   ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                   ,ppiclf_cp_map(1,j)
     >                                   ,ppiclf_cp_map(1+PPICLF_LRS,j))
              ELSE IF (j .LT. 0) THEN ! Ghost Particle
                ! Negative was just use for ghost particle indicator
                ! in subbin mapping array. Need to flip sign
                j = - j                 
                xdist2 = (ppiclf_cp_map(1,i)-ppiclf_rprop_gp(1,j))**2
                IF (xdist2 .GT. dist2) CYCLE
                ydist2 = (ppiclf_cp_map(2,i)-ppiclf_rprop_gp(2,j))**2
                IF (ydist2 .GT. dist2) CYCLE
                dist_total = xdist2 + ydist2
                IF (ppiclf_ndim .EQ. 3) THEN
                zdist2 = (ppiclf_cp_map(3,i)-ppiclf_rprop_gp(3,j))**2
                IF (zdist2 .GT. dist2) CYCLE
                dist_total = dist_total+zdist2
                END IF
                IF (dist_total .GT. dist2) CYCLE
                jp = -1*j
                CALL ppiclf_user_EvalNearestNeighbor(i,jp
     >                                 ,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ppiclf_rprop_gp(1,j)
     >                                 ,ppiclf_rprop_gp(1+PPICLF_LRS,j))
              END IF
            END DO !k
          END IF ! if loopSB is valid
        END DO !kSB
      END DO !jSB
      END DO !iSB


      istride = ppiclf_ndim
      do j=1,ppiclf_nwall

         rnx  = ppiclf_wall_n(1,j)
         rny  = ppiclf_wall_n(2,j)
         rnz  = 0.0d0
         area = ppiclf_wall_n(3,j)
         rpx1 = ppiclf_cp_map(1,i)
         rpy1 = ppiclf_cp_map(2,i)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(1,j)
         rpy2 = ppiclf_wall_c(2,j)
         rpz2 = 0.0d0
         rpx2 = rpx2 - rpx1
         rpy2 = rpy2 - rpy1

         if (ppiclf_ndim .eq. 3) then
            rnz  = ppiclf_wall_n(3,j)
            area = ppiclf_wall_n(4,j)
            rpz1 = ppiclf_cp_map(3,i)
            rpz2 = ppiclf_wall_c(3,j)
            rpz2 = rpz2 - rpz1
         endif
    
         rflip = rnx*rpx2 + rny*rpy2 + rnz*rpz2
         if (rflip .gt. 0.0d0) then
            rnx = -1.0d0*rnx
            rny = -1.0d0*rny
            rnz = -1.0d0*rnz
         endif


         a_sum = 0.0d0
         kmax = 2
         if (ppiclf_ndim .eq. 3) kmax = 3
         do k=1,kmax 
            kp = k+1
            if (kp .gt. kmax) kp = kp-kmax ! cycle
            
            kk   = istride*(k-1)
            kkp  = istride*(kp-1)
            rpx1 = ppiclf_wall_c(kk+1,j)
            rpy1 = ppiclf_wall_c(kk+2,j)
            rpz1 = 0.0d0
            rpx2 = ppiclf_wall_c(kkp+1,j)
            rpy2 = ppiclf_wall_c(kkp+2,j)
            rpz2 = 0.0d0

            if (ppiclf_ndim .eq. 3) then
               rpz1 = ppiclf_wall_c(kk+3,j)
               rpz2 = ppiclf_wall_c(kkp+3,j)
            endif

            rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)

            rdist = abs(rnx*ppiclf_cp_map(1,i)+rny*ppiclf_cp_map(2,i)
     >                 +rnz*ppiclf_cp_map(3,i)+rd)
            rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)

            ! give a little extra room for walls (2x)
            if (rdist .gt. 2.0d0*ppiclf_d2chk(3)) goto 1519

            ydum(1) = ppiclf_cp_map(1,i) - rdist*rnx
            ydum(2) = ppiclf_cp_map(2,i) - rdist*rny
            ydum(3) = 0.0d0

            A(1) = ydum(1)
            A(2) = ydum(2)
            A(3) = 0.0d0

            B(1) = rpx1
            B(2) = rpy1
            B(3) = 0.0d0

            C(1) = rpx2
            C(2) = rpy2
            C(3) = 0.0d0

            AB(1) = B(1) - A(1)
            AB(2) = B(2) - A(2)
            AB(3) = 0.0d0

            AC(1) = C(1) - A(1)
            AC(2) = C(2) - A(2)
            AC(3) = 0.0d0

            if (ppiclf_ndim .eq. 3) then
               ydum(3) = ppiclf_cp_map(3,i) - rdist*rnz
               A(3) = ydum(3)
               B(3) = rpz1
               C(3) = rpz2
               AB(3) = B(3) - A(3)
               AC(3) = C(3) - A(3)

               AB_DOT_AC = AB(1)*AC(1) + AB(2)*AC(2) + AB(3)*AC(3)
               AB_MAG = sqrt(AB(1)**2 + AB(2)**2 + AB(3)**2)
               AC_MAG = sqrt(AC(1)**2 + AC(2)**2 + AC(3)**2)
               theta  = acos(AB_DOT_AC/(AB_MAG*AC_MAG))
               tri_area = 0.5d0*AB_MAG*AC_MAG*sin(theta)
            elseif (ppiclf_ndim .eq. 2) then
               AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
               tri_area = AB_MAG
            endif
            a_sum = a_sum + tri_area
         enddo

         rthresh = 1.10d0 ! keep it from slipping through crack on edges
         if (a_sum .gt. rthresh*area) cycle

         jp = 0
         call ppiclf_user_EvalNearestNeighbor(i,jp,ppiclf_cp_map(1,i)
     >                                 ,ppiclf_cp_map(1+PPICLF_LRS,i)
     >                                 ,ydum
     >                                 ,rpropdum)

 1519 continue
      enddo

      return
      end
!-----------------------------------------------------------------------
       subroutine ppiclf_solve_InitWall(xp1,xp2,xp3)
!
      implicit none
!
      include "PPICLF"
! 
! Input:
! 
      real*8 xp1(*)
      real*8 xp2(*)
      real*8 xp3(*)
!
! Internal:
!
      real*8 rpx1, rpy1, rpz1, rpx2, rpy2, rpz2,
     >       a_sum, theta, tri_area, 
     >       ab_dot_ac, ab_mag, ac_mag, rise, run, rmag, 
     >       rpx3, rpy3, rpz3
      integer*4 istride, k, kmax, kp, kkp, kk
      real*8 A(3),B(3),C(3),AB(3),AC(3)
!
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitWall$',0.d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitWall$'
     >                  ,0.d0,0)

      ppiclf_nwall = ppiclf_nwall + 1 

      if (ppiclf_nwall .gt. PPICLF_LWALL)
     >call ppiclf_exittr('Increase LWALL in user file$'
     >                  ,0.d0,ppiclf_nwall)

      istride = ppiclf_ndim
      a_sum = 0.0d0
      kmax = 2
      if (ppiclf_ndim .eq. 3) kmax = 3

      if (ppiclf_ndim .eq. 3) then
         ppiclf_wall_c(1,ppiclf_nwall) = xp1(1)
         ppiclf_wall_c(2,ppiclf_nwall) = xp1(2)
         ppiclf_wall_c(3,ppiclf_nwall) = xp1(3)
         ppiclf_wall_c(4,ppiclf_nwall) = xp2(1)
         ppiclf_wall_c(5,ppiclf_nwall) = xp2(2)
         ppiclf_wall_c(6,ppiclf_nwall) = xp2(3)
         ppiclf_wall_c(7,ppiclf_nwall) = xp3(1)
         ppiclf_wall_c(8,ppiclf_nwall) = xp3(2)
         ppiclf_wall_c(9,ppiclf_nwall) = xp3(3)

         A(1) = (xp1(1) + xp2(1) + xp3(1))/3.0d0
         A(2) = (xp1(2) + xp2(2) + xp3(2))/3.0d0
         A(3) = (xp1(3) + xp2(3) + xp3(3))/3.0d0
      elseif (ppiclf_ndim .eq. 2) then
         ppiclf_wall_c(1,ppiclf_nwall) = xp1(1)
         ppiclf_wall_c(2,ppiclf_nwall) = xp1(2)
         ppiclf_wall_c(3,ppiclf_nwall) = xp2(1)
         ppiclf_wall_c(4,ppiclf_nwall) = xp2(2)

         A(1) = (xp1(1) + xp2(1))/2.0d0
         A(2) = (xp1(2) + xp2(2))/2.0d0
         A(3) = 0.0d0
      endif

      ! compoute area:
      do k=1,kmax 
         kp = k+1
         if (kp .gt. kmax) kp = kp-kmax ! cycle
         
         kk   = istride*(k-1)
         kkp  = istride*(kp-1)
         rpx1 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy1 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(kkp+1,ppiclf_nwall)
         rpy2 = ppiclf_wall_c(kkp+2,ppiclf_nwall)
         rpz2 = 0.0d0

         B(1) = rpx1
         B(2) = rpy1
         B(3) = 0.0d0
        
         C(1) = rpx2
         C(2) = rpy2
         C(3) = 0.0d0
        
         AB(1) = B(1) - A(1)
         AB(2) = B(2) - A(2)
         AB(3) = 0.0d0
        
         AC(1) = C(1) - A(1)
         AC(2) = C(2) - A(2)
         AC(3) = 0.0d0

         if (ppiclf_ndim .eq. 3) then
             rpz1 = ppiclf_wall_c(kk+3,ppiclf_nwall)
             rpz2 = ppiclf_wall_c(kkp+3,ppiclf_nwall)
             B(3) = rpz1
             C(3) = rpz2
             AB(3) = B(3) - A(3)
             AC(3) = C(3) - A(3)
        
             AB_DOT_AC = AB(1)*AC(1) + AB(2)*AC(2) + AB(3)*AC(3)
             AB_MAG = sqrt(AB(1)**2 + AB(2)**2 + AB(3)**2)
             AC_MAG = sqrt(AC(1)**2 + AC(2)**2 + AC(3)**2)
             theta  = acos(AB_DOT_AC/(AB_MAG*AC_MAG))
             tri_area = 0.5d0*AB_MAG*AC_MAG*sin(theta)
         elseif (ppiclf_ndim .eq. 2) then
             AB_MAG = sqrt(AB(1)**2 + AB(2)**2)
             tri_area = AB_MAG
         endif
         a_sum = a_sum + tri_area
      enddo
      
      ppiclf_wall_n(ppiclf_ndim+1,ppiclf_nwall) = a_sum

      ! wall normal:
      if (ppiclf_ndim .eq. 2) then

         rise = xp2(2) - xp1(2)
         run  = xp2(1) - xp1(1)

         rmag = sqrt(rise**2 + run**2)
         rise = rise/rmag
         run  = run/rmag
         
         ppiclf_wall_n(1,ppiclf_nwall) = rise
         ppiclf_wall_n(2,ppiclf_nwall) = -run

      elseif (ppiclf_ndim .eq. 3) then

         k  = 1
         kk = istride*(k-1)
         rpx1 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy1 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz1 = ppiclf_wall_c(kk+3,ppiclf_nwall)
         
         k  = 2
         kk = istride*(k-1)
         rpx2 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy2 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz2 = ppiclf_wall_c(kk+3,ppiclf_nwall)
         
         k  = 3
         kk = istride*(k-1)
         rpx3 = ppiclf_wall_c(kk+1,ppiclf_nwall)
         rpy3 = ppiclf_wall_c(kk+2,ppiclf_nwall)
         rpz3 = ppiclf_wall_c(kk+3,ppiclf_nwall)
    
         A(1) = rpx2 - rpx1
         A(2) = rpy2 - rpy1
         A(3) = rpz2 - rpz1

         B(1) = rpx3 - rpx2
         B(2) = rpy3 - rpy2
         B(3) = rpz3 - rpz2

         ppiclf_wall_n(1,ppiclf_nwall) = A(2)*B(3) - A(3)*B(2)
         ppiclf_wall_n(2,ppiclf_nwall) = A(3)*B(1) - A(1)*B(3)
         ppiclf_wall_n(3,ppiclf_nwall) = A(1)*B(2) - A(2)*B(1)

         rmag = sqrt(ppiclf_wall_n(1,ppiclf_nwall)**2 +
     >               ppiclf_wall_n(2,ppiclf_nwall)**2 +
     >               ppiclf_wall_n(3,ppiclf_nwall)**2)

         ppiclf_wall_n(1,ppiclf_nwall) = ppiclf_wall_n(1,ppiclf_nwall)
     >                                  /rmag
         ppiclf_wall_n(2,ppiclf_nwall) = ppiclf_wall_n(2,ppiclf_nwall)
     >                                  /rmag
         ppiclf_wall_n(3,ppiclf_nwall) = ppiclf_wall_n(3,ppiclf_nwall)
     >                                  /rmag

      endif

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitPeriodicX(xl,xr)
     > bind(C, name="ppiclc_solve_InitPeriodicX")
#else
      subroutine ppiclf_solve_InitPeriodicX(xl,xr)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8 xl
      real*8 xr
! 
      if (xl .ge. xr)
     >call ppiclf_exittr('PeriodicX must have xl < xr$',xl,0)

      ppiclf_iperiodic(1) = 0

      ppiclf_xdrange(1,1) = xl
      ppiclf_xdrange(2,1) = xr

      call ppiclf_solve_InitSolve

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitPeriodicY(yl,yr)
     > bind(C, name="ppiclc_solve_InitPeriodicY")
#else
      subroutine ppiclf_solve_InitPeriodicY(yl,yr)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8 yl
      real*8 yr
! 
      if (yl .ge. yr)
     >call ppiclf_exittr('PeriodicY must have yl < yr$',yl,0)

      ppiclf_iperiodic(2) = 0

      ppiclf_xdrange(1,2) = yl
      ppiclf_xdrange(2,2) = yr

      call ppiclf_solve_InitSolve

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitPeriodicZ(zl,zr)
     > bind(C, name="ppiclc_solve_InitPeriodicZ")
#else
      subroutine ppiclf_solve_InitPeriodicZ(zl,zr)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8 zl
      real*8 zr
! 
      if (zl .ge. zr)
     >call ppiclf_exittr('PeriodicZ must have zl < zr$',zl,0)
      if (ppiclf_ndim .lt. 3)
     >call ppiclf_exittr('Cannot do PeriodicZ if not 3D$',zl,0)

      ppiclf_iperiodic(3) = 0

      ppiclf_xdrange(1,3) = zl
      ppiclf_xdrange(2,3) = zr

      call ppiclf_solve_InitSolve

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitAngularPeriodic(flag,
     >              rin, rout, angle, xangle)
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      ! Thierry -  07/24/24 - modified code begings here
      ! global variables - user input file
      integer*4 flag
      real*8 rin, rout, angle, xangle
      ! local variables
      real*8 pi, angled

        ! Thierry - 07/24/24 - modified code begins here
        ! Implementation of Angular Periodicity
        ! Just like how Rocflu does it in modflu/RFLU_ModRelatedPatches.F90
        ! this is invoked when particle is leaving x-axis or y-axis
       
        ! sign convention for theta is +ve when measured CCW
        ! switch angle sign when particle is leaving from upper face
        if (rin .ge. rout)
     >   call ppiclf_exittr('Angular Per must have rin < rout$',rout,0)

            ppiclf_iperiodic(1) = 0 ! X-periodic
            ppiclf_iperiodic(2) = 0 ! Y-periodic

            SELECT CASE (ang_case)
              CASE (1) ! general wedge ; 0 <= angle < 90
                if (ppiclf_nid.eq.0) print*,"General Wedge Initialized!"
                ppiclf_xdrange(1,1) = rin  ! Min X-periodic face
                ppiclf_xdrange(2,1) = rout ! Max X-periodic face
                ppiclf_xdrange(1,2) = tan(xangle)*rout ! Min Y-periodic face
                ppiclf_xdrange(2,2) = tan(angle - abs(xangle))*rout ! Max Y-periodic face

              CASE (2) ! quarter cylinder ; angle = 90
                if (ppiclf_nid.eq.0)
     >             print*,"Quarter Cylinder Initialized!"
                ppiclf_xdrange(1,1) = rin  
                ppiclf_xdrange(2,1) = rout 
                ppiclf_xdrange(1,2) = tan(xangle)*rout
                ppiclf_xdrange(2,2) = rout 
              
              CASE (3) ! half cylinder ; angle = 180
                if (ppiclf_nid.eq.0)
     >             print*,"Half Cylinder Initialized!"
                ppiclf_xdrange(1,1) = -1.0*rout
                ppiclf_xdrange(2,1) = rout 
                ppiclf_xdrange(1,2) = tan(xangle)*rout
                ppiclf_xdrange(2,2) = rout 
              
              CASE DEFAULT
                call ppiclf_exittr('Invalid Rotational Case!$',0.0d0
     >             ,ppiclf_nid)
              END SELECT

      call ppiclf_solve_InitSolve
      
      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitAngularPlane(i,rin, rout,
     >                                         angle, xangle,
     >                                         dist1, dist2)
!
      implicit none
!
      include "PPICLF"
! Inputs 
!
      integer*4 i
      real*8 rin, rout, angle, xangle
! Local Variables:
!     
      real*8 p1(3), p2(3), p3(3), p4(3), p5(3), p6(3),
     >       v1(3), v2(3), v3(3), v4(3), n1(3), n2(3),
     >       A, B, C, D, E, F, G, H, zt, xp, yp, zp
!
! Outputs
      real*8 dist1, dist2
!
      xp = ppiclf_y(PPICLF_JX,i)
      yp = ppiclf_y(PPICLF_JY,i)
      zp = ppiclf_y(PPICLF_JZ,i)
      
      zt = ppiclf_binb(6) - ppiclf_binb(5) ! bin thickness in z-direction

      !!! upper plane calculation !!! 
      
      ! plane equation
      ! Ax + By + Cz + D = 0

      ! p1, p2, p3 are 3 points in the upper plane
      p1 = (/rin, tan(angle - abs(xangle))*rin, 0.0d0/) 
      p2 = (/rout, tan(angle - abs(xangle))*rout, 0.0d0/) 
      p3 = (/rin, tan(angle - abs(xangle))*rin, zt/) 

      v1 = p2 - p1 ! vector P1P2
      v2 = p3 - p1 ! vector P1P3
      
      ! upper plane normal vector - n1(A,B,C) = v1 x v2
      ! cross product calculation
      A =  v1(2)*v2(3) - v1(3)*v2(2)
      B = -v1(1)*v2(3) + v1(3)*v2(1)
      C =  v1(1)*v2(2) - v1(2)*v2(1)
      n1(1)=A ; n1(2)=B; n1(3)=C
      
      ! values of either p1, p2, or p3 can be used to calculate D
      D = -A*p1(1) - B*p1(2) - C*p1(3)
      
      ! P(xp, yp, zp) arbitrary point
      ! dist = distance between P and upper plane 
      dist1 = abs(A*xp + B*yp + C*zp + D)
      dist1 = dist1/sqrt(A**2 + B**2 + C**2)
      
      !!! lower plane calculation !!! 
      ! plane equation
      ! Ex + Fy + Gz + H = 0

      ! p4, p5, p6 are 3 points in the lower plane
      p4 = (/rin, -tan(angle - abs(xangle))*rin, 0.0d0/)
      p5 = (/rout, -tan(angle - abs(xangle))*rout, 0.0d0/)
      p6 = (/rin, -tan(angle - abs(xangle))*rin, zt/)
      
      v3 = p5 - p4 ! vector P4P5
      v4 = p6 - p4 ! vector P4P6
      
      ! lower plane normal vector - n2(E,F,G) = v3 x v4
      ! cross product calculation
      E =  v3(2)*v4(3) - v3(3)*v4(2)
      F = -v3(1)*v4(3) + v3(3)*v4(1)
      G =  v3(1)*v4(2) - v3(2)*v4(1)
      n2(1)=E ; n2(2)=F; n2(3)=G
      
      ! values of either p4, p5, or p6 can be used to calculate H
      H = -E*p4(1) - F*p4(2) - G*p4(3)

      ! P(xp, yp, zp) arbitrary point
      ! dist = distance between P and lower plane 
      dist2 = abs(E*xp + F*yp + G*zp + H)
      dist2 = dist2/sqrt(E**2 + F**2 + G**2)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InvokeLinearPeriodic(i)
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, j, jj(3), jchk, in_part(PPICLF_LPART)
      ! 02/08/2025 - Thierry Daoud - starts here 
      !                              Included vector ppiclf_y1(PPICLF_LRS*PPICLF_NPART) in periodicity.      
      !                              Solution of array ppiclf_y is stored in ppiclf_y1 for 1st RK stage (iStage==1).
      integer*4 isl   
      ! 02/08/2025 - Thierry Daoud - ends here   
!
      jj(1) = 1 ; jj(2) = 2 ; jj(3) = 3
      ! 02/08/2025 - Thierry Daoud - starts here 
      isl = (i-1) * PPICLF_LRS + 1 ! This is how RemoveParticle does it in original code version.
      ! 02/08/2025 - Thierry Daoud - ends here 


! Case 1 - Linear Periodicity in any of 3 directions ; NO Anuglar Periodicity
      if(((x_per_flag.eq.1).or.(y_per_flag.eq.1).or.(z_per_flag.eq.1))
     >   .and.(ang_per_flag.eq.0)) then

      do j=0,ppiclf_ndim-1
        jchk = jj(j+1)
        ! particle leaving min. periodic face -> move it relative to 
        !                                         max periodic face
        if(ppiclf_y(jchk,i).lt.ppiclf_xdrange(1,j+1)) then
          ppiclf_y(jchk,i) = ppiclf_xdrange(2,j+1) - 
     >                  abs(ppiclf_xdrange(1,j+1) - ppiclf_y(jchk,i))
        ! 02/08/2025 - Thierry Daoud - starts here 
          ppiclf_y1(isl+j)   = ppiclf_xdrange(2,j+1) +
     >            abs(ppiclf_xdrange(1,j+1) - ppiclf_y1(isl+j))
        ! 02/08/2025 - Thierry Daoud - ends here 
          goto 1512
        endif

        ! particle leaving max. periodic face -> move it relative to 
        !                                         min periodic face
        if(ppiclf_y(jchk,i).gt.ppiclf_xdrange(2,j+1)) then
          ppiclf_y(jchk,i) = ppiclf_xdrange(1,j+1) + 
     >                  abs(ppiclf_y(jchk,i) - ppiclf_xdrange(2,j+1))
        ! 02/08/2025 - Thierry Daoud - starts here 
          ppiclf_y1(isl+j)   = ppiclf_xdrange(1,j+1) +
     >             abs(ppiclf_y1(isl+j) - ppiclf_xdrange(2,j+1))
        ! 02/08/2025 - Thierry Daoud - ends here 
          goto 1512
        endif
        
        ! Thierry - I'm not sure what this does but this is how it was implemented
        if (ppiclf_iprop(1,i) .eq. 2) then
             in_part(i) = -1 ! only if periodic check fails it will get here
        endif
 1512 continue
        end do ! j=0, ndim-1
      endif ! Case 1 

! Case 2 - Linear Periodicity in Z-direction only; WITH Anuglar Periodicity
      if((z_per_flag.eq.1).and.(ang_per_flag.eq.1)) then

        ! particle leaving min. or max. x-face -> delete it
        if((ppiclf_y(1,i).lt.ppiclf_xdrange(1,1)).or.
     >  (ppiclf_y(1,i).gt.ppiclf_xdrange(2,1))) then
          call ppiclf_solve_MarkForRemoval(i)
          goto 1515
        endif

        ! particle leaving min. z-periodic face -> move it relative to 
        !                                         max z-periodic face
        if(ppiclf_y(3,i).lt.ppiclf_xdrange(1,3)) then
          ppiclf_y(3,i) = ppiclf_xdrange(2,3) - 
     >                  abs(ppiclf_xdrange(1,3) - ppiclf_y(3,i))
        ! 02/08/2025 - Thierry Daoud - starts here 
          ppiclf_y1(isl+2)   = ppiclf_xdrange(2,3) +
     >            abs(ppiclf_xdrange(1,3) - ppiclf_y1(isl+2))
        ! 02/08/2025 - Thierry Daoud - ends here 
          goto 1515
        endif

        ! particle leaving max. z-periodic face -> move it relative to 
        !                                         min z-periodic face
        if(ppiclf_y(3,i).gt.ppiclf_xdrange(2,3)) then
          ppiclf_y(3,i) = ppiclf_xdrange(1,3) + 
     >                  abs(ppiclf_y(3,i) - ppiclf_xdrange(2,3))
        ! 02/08/2025 - Thierry Daoud - starts here 
          ppiclf_y1(isl+2)   = ppiclf_xdrange(1,3) +
     >             abs(ppiclf_y1(isl+2) - ppiclf_xdrange(2,3))
        ! 02/08/2025 - Thierry Daoud - ends here 
          goto 1515
        endif
        
        if (ppiclf_iprop(1,i) .eq. 2) then
             in_part(i) = -1 ! only if periodic check fails it will get here
        endif
 1515 continue
      endif ! Case 2

      return
      end
!-----------------------------------------------------------------------
       subroutine ppiclf_solve_InvokeAngularPeriodic(i,flag,
     >                                              per_alpha,
     >                                              angle, xangle,
     >                                              register)
!
      implicit none
!
      include "PPICLF"
! :
! Input: 
! 
      ! Thierry -  07/24/24 - modified code begins here
      ! global variables
      integer*4 i, flag
      real*8 rin, rout, per_alpha, angle, xangle
      ! local variables
      real*8 ct, st, ex, ey, ez, local_angle
      real*8 rotmat(3,3) , v(3), x(3)
      integer*4 register
      ! 02/08/2025 - Thierry Daoud - starts here 
      integer*4 isl, j
      ! 02/08/2025 - Thierry Daoud - ends here 
      ! 02/08/2025 - Thierry Daoud - starts here 
      isl = (i-1) * PPICLF_LRS + 1 ! This is how RemoveParticle does it in original code version.
      ! 02/08/2025 - Thierry Daoud - ends here 



        ! Thierry - 07/24/24 - modified code begins here
        ! Implementation of Rotational Periodicity
        ! Just like how Rocflu does it in modflu/RFLU_ModRelatedPatches.F90
        ! this is invoked when particle is leaving x-axis or y-axis
       
!        print*, "!!! Rotational Periodicity Invoked !!!!" 
          
        ! use local angle so the value of angle does not get affected globally
        local_angle = angle
        ! Thierry 
        !    (1) sign convention for theta is +ve when measured CCW
        !           switch angle sign when particle is leaving from 
        !           upper periodic face
        !    (2) 0.5 instead of 1.0 to switch angle for ghost algorithm
        !           since the ghost is being created before the 
        !           particle is leaving domain
        if(per_alpha.gt. 0.5*(xangle+angle)) 
     >    local_angle=-1.0*local_angle
        
        ! Half-cylinder case - particle leaving +ve x-axis 
        !                    - adjust rotation matrix angle accordingly
        if(ang_case .eq. 3) then
          if(per_alpha .lt. xangle) local_angle = 0.0 
        end if

        ! convert from degrees to radians
        ct = cos(local_angle)
        st = sin(local_angle)
        
        SELECT CASE(flag)
          !CASE(1)
          !  ex = 1.0d0
          !  ey = 0.0d0
          !  ez = 0.0d0
          !  print*, "X-Rotational Axis"

          !CASE(2)
          !  ex = 0.0d0
          !  ey = 1.0d0
          !  ez = 0.0d0
          !  print*, "Y-Rotational Axis"

          CASE(1)
            ex = 0.0d0
            ey = 0.0d0
            ez = 1.0d0
!            print*, "Z-Rotational Axis"
          CASE DEFAULT
            call ppiclf_exittr('Invalid Axis of Rotation!$',0.0d0
     >         ,ppiclf_nid)

          END SELECT 
          
          ! Rotation Matrix calculation
          rotmat(1,1) = ct + (1.0d0-ct)*ex*ex
          rotmat(1,2) =      (1.0d0-ct)*ex*ey - st*ez
          rotmat(1,3) =      (1.0d0-ct)*ex*ez + st*ey
          
          rotmat(2,1) =      (1.0d0-ct)*ey*ex + st*ez
          rotmat(2,2) = ct + (1.0d0-ct)*ey*ey
          rotmat(2,3) =      (1.0d0-ct)*ey*ez - st*ex
          
          rotmat(3,1) =      (1.0d0-ct)*ez*ex - st*ey
          rotmat(3,2) =      (1.0d0-ct)*ez*ey + st*ex
          rotmat(3,3) = ct + (1.0d0-ct)*ez*ez

          ! Corrdinates modification
          x(1) = ppiclf_y(PPICLF_JX,i)
          x(2) = ppiclf_y(PPICLF_JY,i)
          x(3) = ppiclf_y(PPICLF_JZ,i)
          
          xrot = MATMUL(rotmat, x)
          
          ! Velocity vector modification

          v(1) = ppiclf_y(PPICLF_JVX,i)
          v(2) = ppiclf_y(PPICLF_JVY,i)
          v(3) = ppiclf_y(PPICLF_JVZ,i)

          vrot = MATMUL(rotmat, v)
          
          ! 08/27/24 - Thierry - we add a register variable to 
          !   choose if we want to register
          !   the angularly modified variables 
          ! register = 1 when called from RemoveParticle -> we want to modify the values
          ! register = 0 when called from AngularCreateGhost -> we don't want to modify values
          
          ! register modified values
          if (register==1) then
            !print*, "Registering values!" 
            ppiclf_y(PPICLF_JX,i) = xrot(1)
            ppiclf_y(PPICLF_JY,i) = xrot(2)
            ppiclf_y(PPICLF_JZ,i) = xrot(3)
            
            ppiclf_y(PPICLF_JVX,i) = vrot(1)
            ppiclf_y(PPICLF_JVY,i) = vrot(2)
            ppiclf_y(PPICLF_JVZ,i) = vrot(3)
            do j=0,ppiclf_ndim-1
              ppiclf_y1(isl+j)   = xrot(j+1) ! coordinates modification
              ppiclf_y1(isl+j+3) = vrot(j+1) ! velocities modification
            end do ! j 
          end if 
       
      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitGaussianFilter(filt,alpha,iwallm)
     > bind(C, name="ppiclc_solve_InitGaussianFilter")
#else
      subroutine ppiclf_solve_InitGaussianFilter(filt,alpha,iwallm)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8    filt
      real*8    alpha
      integer*4 iwallm
! 
! Internal: 
! 
      real*8 rsig
! 
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitFilter$',0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitFilter$',0.0d0
     >                  ,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.0d0
     >                  ,0)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter can only be called once$',0.0d0,0)
      if (iwallm .lt. 0 .or. iwallm .gt. 1)
     >call ppiclf_exittr('0 or 1 must be used to specify filter mirror$'
     >                  ,0.0d0,iwallm)

      ppiclf_filter = filt
      ppiclf_alpha  = alpha 
      ppiclf_iwallm = iwallm

      rsig             = ppiclf_filter/(2.0d0*sqrt(2.0d0*log(2.0d0)))
      ppiclf_d2chk(2)  = rsig*sqrt(-2*log(ppiclf_alpha))

      ! TLJ added 12/21/2024
      if (ppiclf_nid==0) then
         print*,'TLJ recompute d2chk(2) based on Gausian filter = ',
     >     ppiclf_d2chk(2)
      endif

      PPICLF_LSUBBIN = .true.
      if (ppiclf_ngrids .eq. 0) PPICLF_LSUBBIN = .false.

      PPICLF_LFILT      = .true.
      PPICLF_LFILTGAUSS = .true.

      ppiclf_ngrids = 0 ! for now leave sub bin off

      return
      end
!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_InitBoxFilter(filt,iwallm,sngl_elem)
     > bind(C, name="ppiclc_solve_InitBoxFilter")
#else
      subroutine ppiclf_solve_InitBoxFilter(filt,iwallm,sngl_elem)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8    filt
      integer*4 iwallm
      integer*4 sngl_elem
! 
      if (.not.PPICLF_LCOMM)
     >call ppiclf_exittr('InitMPI must be before InitFilter$',0.0d0,0)
      if (.not.PPICLF_LINIT)
     >call ppiclf_exittr('InitParticle must be before InitFilter$',0.0d0
     >                   ,0)
      if (PPICLF_OVERLAP)
     >call ppiclf_exittr('InitFilter must be before InitOverlap$',0.0d0
     >                   ,0)
      if (PPICLF_LFILT)
     >call ppiclf_exittr('InitFilter can only be called once$',0.0d0,0)

c     filt = sqrt(1.5d0*filt**2/log(2.0d0) + 1.0d0)

      ppiclf_filter = filt
      ppiclf_iwallm = iwallm

      ppiclf_d2chk(2)  = filt/2.0d0

      ! TLJ added 12/21/2024
      if (ppiclf_nid==0) then
         print*,'TLJ checking d2chk(2) = ',ppiclf_d2chk(2)
      endif

      PPICLF_LSUBBIN = .true.
      if (ppiclf_ngrids .eq. 0) PPICLF_LSUBBIN = .false.

      PPICLF_LFILT    = .true.
      PPICLF_LFILTBOX = .true.

      ! option to only use the current element (filter width will be 
      ! ignored)
      ! Note that this assumes the element volume is that of
      ! a cuboid... will need to get a better way for general
      ! hexahedral element eventually
      if ( sngl_elem == 1 ) PPICLF_SNGL_ELEM = .true.

      ppiclf_ngrids = 0 ! for now leave sub bin off

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_SetParticleTag(npart)
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      integer*4 npart
! 
! Internal: 
! 
      integer*4 i
!
      do i=ppiclf_npart-npart+1,ppiclf_npart
         ppiclf_iprop(5,i) = ppiclf_nid 
         ppiclf_iprop(6,i) = ppiclf_cycle
         ppiclf_iprop(7,i) = i
      enddo

      return
      end
c----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_WriteVTU(time)
     > bind(C, name="ppiclc_solve_WriteVTU")
#else
      subroutine ppiclf_solve_WriteVTU(time)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      real*8    time
! 
! Internal:
!
      ppiclf_time   = time

      ! already wrote initial conditions
      !if (ppiclf_cycle .ne. 0) then
            call ppiclf_io_WriteParticleVTU('')
            call ppiclf_io_WriteBinVTU('')
      !endif

      if (ppiclf_lsubbin)
     >      call ppiclf_io_WriteSubBinVTU('')

      ! Output diagnostics
      call ppiclf_io_OutputDiagAll

      return
      end
c----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_solve_IntegrateParticle(istep,iostep,dt,time)
     > bind(C, name="ppiclc_solve_IntegrateParticle")
#else
      subroutine ppiclf_solve_IntegrateParticle(istep,iostep,dt,time)
#endif
!
      implicit none
!
      include "PPICLF"
! 
! Input: 
! 
      integer*4 istep
      integer*4 iostep
      real*8    dt
      real*8    time
! 
! Internal:
!
      logical iout
!
      ppiclf_cycle  = istep
      ppiclf_iostep = iostep
      ppiclf_dt     = dt
      ppiclf_time   = time

      ! integerate in time
      if (ppiclf_imethod .eq. 1) 
     >   call ppiclf_solve_IntegrateRK3(iout)
      if (ppiclf_imethod .eq. -1) 
     >   call ppiclf_solve_IntegrateRK3s(iout)
      if (ppiclf_imethod .eq. 2)
     >   call ppiclf_solve_IntegrateRK3s_Rocflu(iout)

!      ! output files
!      if (ppiclf_iostep .gt.0)then
!      if (mod(ppiclf_cycle,ppiclf_iostep) .eq. 0 .and. iout) then
!
!         ! already wrote initial conditions
!         if (ppiclf_cycle .ne. 0) then
!            call ppiclf_io_WriteParticleVTU('')
!            call ppiclf_io_WriteBinVTU('')
!         endif
!
!         if (ppiclf_lsubbin)
!     >      call ppiclf_io_WriteSubBinVTU('')
!      endif
!
!      ! Output diagnostics
!      if (mod(ppiclf_cycle,ppiclf_iostep) .eq. 0 .and. iout) then
!         call ppiclf_io_OutputDiagAll
!      endif
!      endif

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateRK3(iout)
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, ndum, nstage, istage
!
! Output:
!
      logical iout
!
      ! save stage 1 solution
      ndum = PPICLF_NPART*PPICLF_LRS
      do i=1,ndum
         ppiclf_y1(i) = ppiclf_y(i,1)
      enddo

      ! get rk3 coeffs
      call ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      do istage=1,nstage

         ! evaluate ydot
         call ppiclf_solve_SetYdot

         ! rk3 integrate
         do i=1,ndum
            ndum = PPICLF_NPART*PPICLF_LRS
            ppiclf_y(i,1) =  ppiclf_rk3coef(1,istage)*ppiclf_y1   (i)
     >                     + ppiclf_rk3coef(2,istage)*ppiclf_y    (i,1)
     >                     + ppiclf_rk3coef(3,istage)*ppiclf_ydot (i,1)
         enddo
      enddo

      iout = .true.

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateRK3s(iout)
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, ndum, nstage, istage
      integer*4 icalld
      save      icalld
      data      icalld /0/
!
! Output:
!
      logical iout
!
      icalld = icalld + 1


      ! get rk3 coeffs
      call ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      istage = mod(icalld,nstage)
      if (istage .eq. 0) istage = 3
      iout = .false.
      if (istage .eq. nstage) iout = .true.

      ! save stage 1 solution
      if (istage .eq. 1) then
      ndum = PPICLF_NPART*PPICLF_LRS
      do i=1,ndum
         ppiclf_y1(i) = ppiclf_y(i,1)
      enddo
      endif

      ! evaluate ydot
      call ppiclf_solve_SetYdot

      ! rk3 integrate
      ndum = PPICLF_NPART*PPICLF_LRS
      do i=1,ndum
         ppiclf_y(i,1) =  ppiclf_rk3coef(1,istage)*ppiclf_y1   (i)
     >                  + ppiclf_rk3coef(2,istage)*ppiclf_y    (i,1)
     >                  + ppiclf_rk3coef(3,istage)*ppiclf_ydot (i,1)
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_IntegrateRK3s_Rocflu(iout)
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, ndum, nstage, istage
      integer*4 icalld
      integer*4 j
      save      icalld
      data      icalld /0/

!
! Output:
!
      logical iout
!
      icalld = icalld + 1

      ! get rk3 coeffs
      call ppiclf_solve_SetRK3Coeff(ppiclf_dt)

      nstage = 3
      istage = mod(icalld,nstage)
      if (istage .eq. 0) istage = 3
      iout = .false.
      if (istage .eq. nstage) iout = .true.
      
      ! 02/08/2025 - Thierry Daoud - starts here
      if (istage .eq. 1) then
      ndum = 0
      do j=1,PPICLF_NPART
      do i=1,PPICLF_LRS
         ndum = ndum + 1
         ppiclf_y1(ndum) =  ppiclf_y(i,j)
      enddo
      enddo
      endif
      ! 02/08/2025 - Thierry Daoud - ends here

      ! evaluate ydot
      call ppiclf_solve_SetYdot

      ! 02/08/2025 - Thierry Daoud - commented out TLj zeroing below
      
      !Zero out for first stage
      ! TLJ this loop is fine
!      if (istage .eq. 1) then
!        ppiclf_y1 = 0.0d0
!        !ndum = PPICLF_NPART*PPICLF_LRS
!        !do i=1,ndum
!        !  ppiclf_y1(i) = 0.0d0 
!        !enddo
!      endif
      
      ! 02/08/2025 - Thierry Daoud - ends here

      ! TLJ comment Dec 7, 2023
      ! The Rocflu RK3 can be found in equation (7) of:
      ! S. Yu. "Runge-Kutta Methods Combined with Compact
      !   Difference Schemes for the Unsteady Euler Equations".
      !   Center for Modeling of Turbulence and Transition.
      !   Research Briefs, 1991.

      ! TLJ modified loop to prevent -fcheck=all error
      !     must preserve order
      !     ppiclf_y   (PPICLF_LRS, PPICLF_LPART)
      !     ppiclf_ydot(PPICLF_LRS, PPICLF_LPART)
      ndum = 0
      do j=1,PPICLF_NPART
      do i=1,PPICLF_LRS
         ndum = ndum + 1
         ppiclf_y(i,j) =  -ppiclf_rk3coef(1,istage)*ppiclf_y1   (ndum)
     >                   + ppiclf_rk3coef(2,istage)*ppiclf_y    (i,j)
     >                   + ppiclf_rk3coef(3,istage)*ppiclf_ydot (i,j)
      enddo
      enddo
      !ndum = PPICLF_NPART*PPICLF_LRS
      !do i=1,ndum
      !   ppiclf_y(i,1) =  -ppiclf_rk3coef(1,istage)*ppiclf_y1   (i)
      !>                   + ppiclf_rk3coef(2,istage)*ppiclf_y    (i,1)
      !>                   + ppiclf_rk3coef(3,istage)*ppiclf_ydot (i,1)
      !enddo

!Store Current stage RHS for next stage's use
      ! TLJ modified loop to prevent -fcheck=all error
      !     must preserve order
      ndum = 0
      do j=1,PPICLF_NPART
      do i=1,PPICLF_LRS
         ndum = ndum + 1
         ppiclf_y1(ndum) =  ppiclf_ydot(i,j)
      enddo
      enddo
      !ndum = PPICLF_NPART*PPICLF_LRS
      !do i=1,ndum
      !   ppiclf_y1(i) = ppiclf_ydot(i,1)
      !enddo
!
!!WAARNING: Experimental fix to keep particles unsure where to place this
!!          command. Either before or after the storing of the current 
!!          storage
        call ppiclf_solve_RemoveParticle      
!End Experimental fix
      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_SetYdot
!
      implicit none
!
      include "PPICLF"
! 
      call ppiclf_solve_InitSolve
      call ppiclf_user_SetYdot
      call ppiclf_solve_RemoveParticle

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_InitSolve
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 i, j
!
      call ppiclf_comm_CreateBin
      call ppiclf_comm_FindParticle
      call ppiclf_comm_MoveParticle
      if (ppiclf_overlap) 
     >   call ppiclf_comm_MapOverlapMesh
      if ((ppiclf_lintp .and. ppiclf_int_icnt .ne. 0) .or.
     >    (ppiclf_lproj .and. ppiclf_sngl_elem))
     >   call ppiclf_solve_InterpParticleGrid
      call ppiclf_solve_RemoveParticle
      if (ppiclf_lsubsubbin .or. ppiclf_lproj) then
         ! Thierry - standard ghost algorithm when no angular periodicity
         if(ang_per_flag.eq.0) then
           call ppiclf_comm_CreateGhost
         elseif(ang_per_flag.eq.1) then
           call ppiclf_comm_AngularCreateGhost
         endif
         call ppiclf_comm_MoveGhost
      endif

      if (ppiclf_lproj .and. ppiclf_overlap) 
     >   call ppiclf_solve_ProjectParticleGrid
      if (ppiclf_lsubsubbin) 
     >   call ppiclf_solve_SetNeighborBin
      ! Zero 
      do i=1,PPICLF_LPART
      do j=1,PPICLF_LRS
         ppiclf_ydotc(j,i) = 0.0d0
      enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpParticleGrid
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 j
!
      call ppiclf_solve_InitInterp
      do j=1,PPICLF_INT_ICNT
         call ppiclf_solve_InterpField(j)
      enddo
      call ppiclf_solve_FinalizeInterp

      call ppiclf_solve_LocalInterp

      call ppiclf_solve_PostInterp

      PPICLF_INT_ICNT = 0


      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpFieldUser(jp,infld)
!
      implicit none
!
      include "PPICLF"
!
! Input: 
!
      integer*4 jp
      real*8 infld(*)
!
! Internal:
!
      integer*4 n
!
      if (PPICLF_INTERP .eq. 0)
     >call ppiclf_exittr(
     >     'No specified interpolated fields, set PPICLF_LRP_INT$',0.0d0
     >                   ,0)

      PPICLF_INT_ICNT = PPICLF_INT_ICNT + 1

      if (PPICLF_INT_ICNT .gt. PPICLF_LRP_INT)
     >   call ppiclf_exittr('Interpolating too many fields$'
     >                     ,0.0d0,PPICLF_INT_ICNT)
      if (jp .le. 0 .or. jp .gt. PPICLF_LRP)
     >   call ppiclf_exittr('Invalid particle array interp. location$'
     >                     ,0.0d0,jp)

      ! set up interpolation map
      PPICLF_INT_MAP(PPICLF_INT_ICNT) = jp

      ! copy to infld internal storage
      n = PPICLF_NEE*PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      call ppiclf_copy(ppiclf_int_fldu(1,1,1,1,PPICLF_INT_ICNT)
     >                ,infld(1),n)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InitInterp
!
      implicit none
!
      include "PPICLF"
! 
! Internal: 
! 
      integer*4 n, ie, npt_max, np, ndum
      real*8 tol, bb_t
!
      if (.not.ppiclf_overlap)
     >call ppiclf_exittr('Cannot interpolate unless overlap grid$',0.0d0
     >                   ,0)
      if (.not.ppiclf_lintp) 
     >call ppiclf_exittr('To interpolate, set PPICLF_LRP_PRO to ~= 0$'
     >                   ,0.0d0,0)

      n = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_neltb
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,1)
     >                   ,ppiclf_xm1b(1,1,1,1,ie),n)
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,2)
     >                   ,ppiclf_xm1b(1,1,1,2,ie),n)
         call ppiclf_copy(ppiclf_xm1bi(1,1,1,ie,3)
     >                   ,ppiclf_xm1b(1,1,1,3,ie),n)
      enddo

      tol     = 5e-13
      bb_t    = 0.01
      npt_max = 128
      np      = 1
c     ndum    = ppiclf_neltb*n
      ndum    = ppiclf_neltb+2

      ! initiate findpts since mapping can change on next call
      call pfgslib_findpts_setup(ppiclf_fp_hndl
     >                         ,ppiclf_comm_nid
     >                         ,np ! only 1 rank on this comm
     >                         ,ppiclf_ndim
     >                         ,ppiclf_xm1bi(1,1,1,1,1)
     >                         ,ppiclf_xm1bi(1,1,1,1,2)
     >                         ,ppiclf_xm1bi(1,1,1,1,3)
     >                         ,PPICLF_LEX
     >                         ,PPICLF_LEY
     >                         ,PPICLF_LEZ
     >                         ,ppiclf_neltb
     >                         ,2*PPICLF_LEX
     >                         ,2*PPICLF_LEY
     >                         ,2*PPICLF_LEZ
     >                         ,bb_t
     >                         ,ndum
     >                         ,ndum
     >                         ,npt_max
     >                         ,tol)


      ! copy MapOverlapMesh mapping from prior to communicating map
      ppiclf_neltbbb = ppiclf_neltbb
      do ie=1,ppiclf_neltbbb
         call ppiclf_icopy(ppiclf_er_mapc(1,ie),ppiclf_er_maps(1,ie)
     >             ,PPICLF_LRMAX)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_InterpField(j)
!
      implicit none
!
      include "PPICLF"
!
! Input: 
!
      integer*4 jp
!
! Internal:
!
      integer*4 n, ie, iee, j
!
      ! use the map to take original grid and map to fld which will be
      ! sent to mapped processors
      n = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_neltbbb
         iee = ppiclf_er_mapc(1,ie)
         call ppiclf_copy(ppiclf_int_fld (1,1,1,j  ,ie)
     >                   ,ppiclf_int_fldu(1,1,1,iee,j ),n)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_FinalizeInterp
!
      implicit none
!
      include "PPICLF"
!
! Internal: 
!
      real*8 FLD(PPICLF_LEX,PPICLF_LEY,PPICLF_LEZ,PPICLF_LEE)
      integer*4 nkey(2), nl, nii, njj, nxyz, nrr, ix, iy, iz, i, jp, ie
      logical partl
!
      ! send it all
      nl   = 0
      nii  = PPICLF_LRMAX
      njj  = 6
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      nrr  = nxyz*PPICLF_LRP_INT
      nkey(1) = 2
      nkey(2) = 1
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl,ppiclf_neltbbb
     >      ,PPICLF_LEE,ppiclf_er_mapc,nii,partl,nl,ppiclf_int_fld
     >      ,nrr,njj)
      call pfgslib_crystal_tuple_sort    (ppiclf_cr_hndl,ppiclf_neltbbb
     >       ,ppiclf_er_mapc,nii,partl,nl,ppiclf_int_fld,nrr,nkey,2)

      ! find which cell particle is in locally
      ix = 1
      iy = 2
      iz = 1
      if (ppiclf_ndim .eq. 3)
     >iz = 3

      call pfgslib_findpts(PPICLF_FP_HNDL           !   call pfgslib_findpts( ihndl,
     >        , ppiclf_iprop (1 ,1),PPICLF_LIP        !   $             rcode,1,
     >        , ppiclf_iprop (3 ,1),PPICLF_LIP        !   &             proc,1,
     >        , ppiclf_iprop (2 ,1),PPICLF_LIP        !   &             elid,1,
     >        , ppiclf_rprop2(1 ,1),PPICLF_LRP2       !   &             rst,ndim,
     >        , ppiclf_rprop2(4 ,1),PPICLF_LRP2       !   &             dist,1,
     >        , ppiclf_y     (ix,1),PPICLF_LRS        !   &             pts(    1),1,
     >        , ppiclf_y     (iy,1),PPICLF_LRS        !   &             pts(  n+1),1,
     >        , ppiclf_y     (iz,1),PPICLF_LRS ,PPICLF_NPART) !   &             pts(2*n+1),1,n)

      do i=1,PPICLF_INT_ICNT
         jp = PPICLF_INT_MAP(i)

         do ie=1,ppiclf_neltbbb
            call ppiclf_copy(fld(1,1,1,ie)
     >                      ,ppiclf_int_fld(1,1,1,i,ie),nxyz)
         enddo

         ! sam commenting out eval nearest neighbor to use Local Interp instead
         ! leaving findpts call to help with projection, where the element id is
         ! needed. 
         ! interpolate field locally
!         call pfgslib_findpts_eval_local( PPICLF_FP_HNDL
!     >                                  ,ppiclf_rprop (jp,1)
!     >                                  ,PPICLF_LRP
!     >                                  ,ppiclf_iprop (2,1)
!     >                                  ,PPICLF_LIP
!     >                                  ,ppiclf_rprop2(1,1)
!     >                                  ,PPICLF_LRP2
!     >                                  ,PPICLF_NPART
!     >                                  ,fld)

      enddo

      ! free since mapping can change on next call
      call pfgslib_findpts_free(PPICLF_FP_HNDL)

      ! Set interpolated fields to zero again
      ! Sam - commenting out for local routine
      !PPICLF_INT_ICNT = 0

      return
      end
!
!-----------------------------------------------------------------------
!
!     Avery's latest version Sept 26, 2024
!
      SUBROUTINE ppiclf_solve_LocalInterp
      IMPLICIT NONE

      include "PPICLF"

      ! Local Variables
      INTEGER*4 i, j, k, l, ix, iy, iz, ip, ie, iee, nxyz, nnearest, 
     >          inearest(28)
      REAL*8    d2l, d2i, wsum, eps, A(27,4), d2(28), xp(3),  
     >          center(3,28), b(27,1), w(27), centeri(3,ppiclf_neltbbb),
     >          d2i_EleLen(3), MaxPoint(3), MinPoint(3),d2Max_EleLen(3)
      LOGICAL   added, farAway 
      !***************************************************************
 
      eps = 1.0e-12 !machine epsilon
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ !number of points in mesh
                                              !element 
      ! Calculate centroid, max cell lengths
      DO ie = 1,ppiclf_neltbbb !Loop fluid cells on this processor
        ! Initialize as zero for each element
        DO l = 1,3
          centeri(l,ie)   =  0.0 
          MaxPoint(l)     = -1.0E10 
          MinPoint(l)     =  1.0E10 
          d2i_EleLen(l)   =  0.0    
        ENDDO !l
        ! Add all x,y,z mesh points for centroid and find extremes
        DO l = 1,3
          DO k = 1,PPICLF_LEZ
            DO j = 1,PPICLF_LEY
              DO i = 1,PPICLF_LEX
                centeri(l,ie) = centeri(l,ie) + ppiclf_xm1b(i,j,k,l,ie)
                IF (ppiclf_xm1b(i,j,k,l,ie) .GT. MaxPoint(l)) 
     >            MaxPoint(l) = ppiclf_xm1b(i,j,k,l,ie)  
                IF (ppiclf_xm1b(i,j,k,l,ie) .LT. MinPoint(l))
     >            MinPoint(l) = ppiclf_xm1b(i,j,k,l,ie)
              ENDDO !i
            ENDDO !j
          ENDDO !k
          ! Find individual (mesh element length)^2 in all directions
          d2i_EleLen(l) = (MaxPoint(l)-MinPoint(l))**2
          ! Find max (mesh element length)^2 for all mesh elements in
          ! all directions
          IF (d2i_EleLen(l) .GT. d2Max_EleLen(l)) THEN
            d2Max_EleLen(l) = d2i_EleLen(l)
          ENDIF
          ! Divide by number of points in mesh element to find centroid
          centeri(l,ie) = centeri(l,ie) / nxyz
        ENDDO !l
      ENDDO !ie
      DO ip=1,ppiclf_npart !Loop all particles in this bin
        ! particle centers in all directions
        xp(1) = ppiclf_y(PPICLF_JX, ip)
        xp(2) = ppiclf_y(PPICLF_JY, ip)
        xp(3) = ppiclf_y(PPICLF_JZ, ip)
        nnearest = 0 ! number of nearest elements
        DO ie = 1,28
          inearest(ie) = -1 ! index of nearest elements
          d2(ie) = 1E20 ! distance to center of nearest element
        ENDDO !ie
        DO ie = 1,ppiclf_neltbbb
          ! get distance from particle to center
          d2l     = 0.0
          d2i     = 0.0
          farAway = .TRUE.
          DO l=1,3
            d2l  =(centeri(l,ie) - xp(l))**2 
            d2i = d2i + d2l
            IF (d2l < (1.5**2)*d2Max_EleLen(l)) farAway = .FALSE.
          ENDDO !l
          ! skip to next fluid cell if greater than 1.5*max cell
          ! distance in respective x,y,z direction.
          if (farAWAY) CYCLE !ie
          ! Sort closest fluid cell centers
          added = .FALSE.
          DO i=1,27
            j = 27 - i + 1
            IF (d2i .LT. d2(j)) THEN
              d2(j+1) = d2(j)
              inearest(j+1) = inearest(j)
              DO l=1,3
                center(l, j+1) = center(l, j)
              ENDDO
              d2(j) = d2i
              inearest(j) = ie
              DO l=1,3
                center(l, j) = centeri(l,ie)
              ENDDO
              added = .TRUE.
            ELSE ! If not within closest cell list
              EXIT !i
            ENDIF
          ENDDO !i
          IF (added) nnearest = nnearest + 1
        ENDDO ! ie
        nnearest = min(nnearest, 27)
        IF (nnearest .lt. 1) THEN
          PRINT *, 'nnearest', nnearest, ip, ppiclf_npart, xp
          PRINT *, ppiclf_rprop(1:PPICLF_LRP, ip)
          PRINT *, ppiclf_y(1:PPICLF_LRS, ip)
          CALL ppiclf_exittr('Failed to interpolate',0.0d0,nnearest)
        ELSE
            DO i=1,nnearest
              DO j=1,3
                A(i, j) = xp(j) - center(j, i)
              ENDDO !j
              A(i, 4) = 1
            ENDDO !i
            DO i=1,PPICLF_INT_ICNT
              DO k=1,nnearest
                b(k, 1) = 0.0 ! cell averaged properties
                DO iz=1,PPICLF_LEZ
                  DO iy=1,PPICLF_LEY
                    DO ix=1,PPICLF_LEX
                      b(k, 1) = b(k, 1) + ppiclf_int_fld(ix,iy,iz,i,
     >                                                  inearest(k))
                    ENDDO !ix
                  ENDDO !iy
                ENDDO !iz
                b(k, 1) = b(k, 1) / nxyz
              ENDDO ! nnearest
              j = PPICLF_INT_MAP(i)
              ! Inverse Distance Interpolation
              ppiclf_rprop(j, ip) = 0
              wsum = 0
              DO k=1,nnearest
                w(k) = 1.0d0 / (SQRT(d2(k)) + eps)
                ppiclf_rprop(j, ip) = ppiclf_rprop(j, ip) + w(k)*b(k, 1)
                wsum = wsum + w(k)
              ENDDO ! k
              ppiclf_rprop(j, ip) = ppiclf_rprop(j, ip) / wsum
              IF (isnan(ppiclf_rprop(j,ip))) THEN
                PRINT *, ip,ppiclf_nid,xp,nnearest
              ENDIF
            ENDDO ! i
        ENDIF ! nnearest
      ENDDO ! ip
      RETURN
      END

!
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_LocalInterp_old
      implicit none

      include "PPICLF"

      ! internal
      integer*4 i,j,k,l,ix,iy,iz
      integer*4 ip, ie, iee, inearest(28), nnearest, nxyz, neltgg
      real*8 A(27, 4), d2i, d2(28), xp(3), center(3, 28)
      real*8 U(27, 27), SIG(4), Vt(4, 4), b(27, 1)
      real*8 interp(4, 1) ! for SVD
      integer*4 m, n, lda, ldu, ldvt, lwork, info, ierr
      real*8 w(27), wsum, eps
      real*8 work(5*27+4)
      character jobu, jobv
      logical added 
      integer*4 nl, nii, njj, nkey(2), nrr
      logical partl
      ! Avery added
      real*8 centeri(3,ppiclf_neltbbb), MaxEleSize, EleSizei,
     >       MaxPoint(3), MinPoint(3)
      !integer subbin_part, subbin_ele(piclf_neltbbb,1) 

      eps = 1.0e-12
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      
      ! Avery Add.  Calculate centroid, max cell size, cell subbin
      MaxEleSize = 0.0
      do ie=1,ppiclf_neltbbb
        ! calculate centroid and largest cell diagonal length
        EleSizei = 0.0  
        do l=1,3
          centeri(l,ie) =  0.0
          MaxPoint(l)   = -1.0E10 
          MinPoint(l)   =  1.0E10
        enddo !l

        do l=1,3
        do k=1,PPICLF_LEZ
        do j=1,PPICLF_LEY
        do i=1,PPICLF_LEX
          centeri(l,ie) = centeri(l,ie) + ppiclf_xm1b(i, j, k, l, ie)
          if (ppiclf_xm1b(i,j,k,l,ie) .gt. MaxPoint(l)) 
     >        MaxPoint(l) = ppiclf_xm1b(i,j,k,l,ie)  
          if (ppiclf_xm1b(i,j,k,l,ie) .lt. MinPoint(l))
     >        MinPoint(l) = ppiclf_xm1b(i,j,k,l,ie)
        enddo !i
        enddo !j
        enddo !k
        enddo !l
        
        do l=1,3
          EleSizei = EleSizei + (MaxPoint(l)-MinPoint(l))**2
        enddo

        if (EleSizei .gt. MaxEleSize) MaxEleSize = EleSizei

        do l=1,3
          centeri(l,ie) = centeri(l,ie) / nxyz
        enddo !l

      enddo !ie
      ! End Avery Added
      
      ! neighbor search O(Nparticles * Nelements)
      ! this should be switched to a KD tree for Log(Nelements) scaling
  
      ! Avery - I'm not sure that KD tree is best for dynamic point
      ! clouds.  Maybe we have multiple "sub-bins" on each processor.
      ! Seems that this is called an Octree in CS speak.
      ! If we do implement some type of tree, we should probably do the
      ! same thing for particle-particle nearest neighbor search as
      ! well.  The current nearest neighbor does the same cycle method
      ! as below.
  
      do ip=1,ppiclf_npart
        ! particle center
        xp(1) = ppiclf_y(PPICLF_JX, ip)
        xp(2) = ppiclf_y(PPICLF_JY, ip)
        xp(3) = ppiclf_y(PPICLF_JZ, ip)
        
        !start particle in subbin
        !subbin_part = !subbin
        !end particle subbin

        nnearest = 0 ! number of nearby elements

        do ie=1,28
          inearest(ie) = -1 ! index of nearest elements
          d2(ie) = 1E20 ! distance to center of nearest element
        enddo

        do ie=1,ppiclf_neltbbb
                   ! get distance from particle to center
          d2i = 0
          do l=1,3
            d2i = d2i + (centeri(l,ie) - xp(l))**2
          enddo
          ! Avery Added if / cycle
          ! Go to next cell if particle is 1.5*largest grid cell diagonal
          ! direction away from neighboring cell centroid
          if (d2i .gt. 2.25d0*MaxEleSize) cycle !1.5**2 = 2.25

          ! sort
          added = .false.
          do i=1,27
            j = 27 - i + 1

            if (d2i .lt. d2(j)) then
              d2(j+1) = d2(j)
              inearest(j+1) = inearest(j)
              do l=1,3
                center(l, j+1) = center(l, j)
              enddo

              d2(j) = d2i
              inearest(j) = ie
              do l=1,3
                center(l, j) = centeri(l,ie)
              enddo

              added = .true.
            else ! Avery added else/exit
              exit
            endif
          enddo !i

          if (added) nnearest = nnearest + 1
          
        enddo ! ie
        ! Avery added if
        ! Check for at least 16 due to cases when one layer of 9 cells
        ! isn't available because the particle is near the bin boundary.
        !if (nnearest .lt. 16) then
        !    print *, '***WARNING***: Less than 16 interpolated cells'
        !endif
       
        nnearest = min(nnearest, 27)

        if (nnearest .lt. 1) then
          print *, 'nnearest', nnearest, ip, ppiclf_npart, xp
          print *, ppiclf_rprop(1:PPICLF_LRP, ip)
          print *, ppiclf_y(1:PPICLF_LRS, ip)
          print *, ppiclf_y(1:PPICLF_LRS, ip)
          call ppiclf_exittr('Failed to interpolate',0.0d0,nnearest)
        else
            do i=1,nnearest
              do j=1,3
                A(i, j) = xp(j) - center(j, i)
              end do
              A(i, 4) = 1
            enddo
! Avery print
!        write(ppiclf_nid,*) ip, xp, nnearest, ppiclf_nid
            do i=1,PPICLF_INT_ICNT

              do k=1,nnearest
                b(k, 1) = 0.0 ! cell averaged properties
                ! Avery Add
!                if (i ==1) then
!                  write(ppiclf_nid,*) ip,
!     >                     center(1:3,k),indg(k)
!                endif
                do iz=1,PPICLF_LEZ
                do iy=1,PPICLF_LEY
                do ix=1,PPICLF_LEX
                  b(k, 1) = b(k, 1) + ppiclf_int_fld(ix,iy,iz,i,
     >                                inearest(k))
                enddo
                enddo
                enddo

                b(k, 1) = b(k, 1) / nxyz
              enddo ! nnearest

              j = PPICLF_INT_MAP(i)

              ! Nearest neighbor interpolation
              !ppiclf_rprop(j, ip) = b(1, 1) ! nearest neighbor interpolation

              ! "harmonic" interpolation
              ppiclf_rprop(j, ip) = 0
              wsum = 0
              do k=1,nnearest
                w(k) = 1.0d0 / (sqrt(d2(k)) + eps)
                ppiclf_rprop(j, ip) = ppiclf_rprop(j, ip) + w(k)*b(k, 1)
                wsum = wsum + w(k)
              enddo

              ppiclf_rprop(j, ip) = ppiclf_rprop(j, ip) / wsum
            if (isnan(ppiclf_rprop(j,ip))) then
              PRINT *, ip,ppiclf_nid,xp,nnearest
            endif
            enddo

        endif ! nnearest
      enddo ! ip

      ! reset here instead of in finalize
      !PPICLF_INT_ICNT = 0

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_PostInterp
!
      implicit none
!
      include "PPICLF"
!
! Internal: 
!
      integer*4 rstride, istride
      parameter(rstride = 7 + PPICLF_LRP_INT)
      parameter(istride = 3)
      real*8 coord(rstride,PPICLF_LPART)
      integer*4 flag(istride,PPICLF_LPART)
      integer*4 fp_handle, i, j, k, npart
      external ppiclf_iglsum
      integer*4 ppiclf_iglsum
      integer*4 npt_max, np, ndum
      real*8 tol, bb_t
      integer*4 copy_back, jp, nxyz, ie
      real*8 xgrid(PPICLF_LEX, PPICLF_LEY, PPICLF_LEZ, PPICLF_LEE)
     >      ,ygrid(PPICLF_LEX, PPICLF_LEY, PPICLF_LEZ, PPICLF_LEE)
     >      ,zgrid(PPICLF_LEX, PPICLF_LEY, PPICLF_LEZ, PPICLF_LEE)
!
      ! Copy not found particles
      npart = 0
      do i=1,ppiclf_npart
         if (ppiclf_iprop(1,i) .eq. 2) then
            npart = npart + 1
            do j=1,ppiclf_ndim
               coord(j,npart) = ppiclf_y(j,i)
            enddo
         endif
      enddo

      if (ppiclf_iglsum(npart,1) .eq. 0) then
         return
      endif

      ! Copy grid indexing 
      ! TLJ changing loop structure to prevent -fcheck=all error
      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ
      do ie=1,ppiclf_nee
      !do i=1,nxyz
      !   xgrid(i,1,1,ie) = ppiclf_xm1bs(i,1,1,1,ie)
      !   ygrid(i,1,1,ie) = ppiclf_xm1bs(i,1,1,2,ie)
      !   zgrid(i,1,1,ie) = ppiclf_xm1bs(i,1,1,3,ie)
      !enddo
      do k=1,PPICLF_LEZ
      do j=1,PPICLF_LEY
      do i=1,PPICLF_LEX
         xgrid(i,j,k,ie) = ppiclf_xm1bs(i,j,k,1,ie)
         ygrid(i,j,k,ie) = ppiclf_xm1bs(i,j,k,2,ie)
         zgrid(i,j,k,ie) = ppiclf_xm1bs(i,j,k,3,ie)
      enddo
      enddo
      enddo

      enddo

      tol     = 5e-13
      bb_t    = 0.01
      npt_max = 128
      np      = ppiclf_np
c     ndum    = ppiclf_nee*n
      ndum    = ppiclf_nee+2

      ! initiate findpts since mapping can change on next call
      call pfgslib_findpts_setup(fp_handle
     >                         ,ppiclf_comm
     >                         ,np 
     >                         ,ppiclf_ndim
     >                         ,xgrid
     >                         ,ygrid
     >                         ,zgrid
     >                         ,PPICLF_LEX
     >                         ,PPICLF_LEY
     >                         ,PPICLF_LEZ
     >                         ,ppiclf_nee
     >                         ,2*PPICLF_LEX
     >                         ,2*PPICLF_LEY
     >                         ,2*PPICLF_LEZ
     >                         ,bb_t
     >                         ,ndum
     >                         ,ndum
     >                         ,npt_max
     >                         ,tol)

      call pfgslib_findpts(fp_handle           !   call pfgslib_findpts( ihndl,
     >        , flag (1 ,1),istride        !   $             rcode,1,
     >        , flag (3 ,1),istride        !   &             proc,1,
     >        , flag (2 ,1),istride        !   &             elid,1,
     >        , coord(4 ,1),rstride       !   &             rst,ndim,
     >        , coord(7 ,1),rstride       !   &             dist,1,
     >        , coord(1,1) ,rstride        !   &             pts(    1),1,
     >        , coord(2,1) ,rstride        !   &             pts(  n+1),1,
     >        , coord(3,1) ,rstride ,npart) !   &             pts(2*n+1),1,n)

      do i=1,PPICLF_LRP_INT

      ! sam - see finalize interp for note
         ! interpolate field (non-local)
!         call pfgslib_findpts_eval( fp_handle
!     >                                  ,coord (7+i,1)
!     >                                  ,rstride
!     >                                  ,flag (1,1)
!     >                                  ,istride
!     >                                  ,flag (3,1)
!     >                                  ,istride
!     >                                  ,flag (2,1)
!     >                                  ,istride
!     >                                  ,coord(4,1)
!     >                                  ,rstride
!     >                                  ,npart
!     >                                  ,ppiclf_int_fldu(1,1,1,1,i))

      enddo

      ! free since mapping can change on next call
      call pfgslib_findpts_free(fp_handle)

      ! Now copy particles back (assumes same ordering)
      k = 0
      do i=1,ppiclf_npart
         copy_back = 0
         if (ppiclf_iprop(1,i) .eq. 2) then
            k = k + 1
            if (flag(1,k) .lt. 2) then
               copy_back = 1
            endif
         endif

         if (copy_back .eq. 1) then
            ppiclf_iprop(1,i) = flag(1,k)
            do j=1,PPICLF_LRP_INT
               jp = PPICLF_INT_MAP(j)
               ppiclf_rprop(jp,i) = coord(7+j,k)
            enddo
         endif
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_SetRK3Coeff(dt)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      real*8 dt
!


      if (ppiclf_imethod .eq. 2) then
        !BD:Rocflu's rk3 scheme

        !Folowing form
        !rk3(1,:) = Temp storage i.e. Previous Stage RHS
        !rk3(2,:) = Temp storage i.e. Current Stage iteration
        !rk3(3,:) = Temp storage i.e. Current Stage RHS
        ppiclf_rk3ark(1) = 8.0d0/15.0d0
        ppiclf_rk3ark(2) = 5.0d0/12.0d0
        ppiclf_rk3ark(3) = 0.75d0

        ppiclf_rk3coef(1,1) = 0.d00
        ppiclf_rk3coef(2,1) = 1.0d0
        ppiclf_rk3coef(3,1) = dt*8.0d0/15.0d0
        ppiclf_rk3coef(1,2) = dt*17.0d0/60.0d0
        ppiclf_rk3coef(2,2) = 1.0d0
        ppiclf_rk3coef(3,2) = dt*5.0d0/12.0d0
        ppiclf_rk3coef(1,3) = dt*5.0d0/12.0d0
        ppiclf_rk3coef(2,3) = 1.0d0
        ppiclf_rk3coef(3,3) = dt*3.0d0/4.0d0
      else
        !BD:Original Code This follows CMT-nek's rk 3 scheme
        ppiclf_rk3coef(1,1) = 0.d00
        ppiclf_rk3coef(2,1) = 1.0d0 
        ppiclf_rk3coef(3,1) = dt
        ppiclf_rk3coef(1,2) = 3.0d0/4.0d0
        ppiclf_rk3coef(2,2) = 1.0d0/4.0d0 
        ppiclf_rk3coef(3,2) = dt/4.0d0
        ppiclf_rk3coef(1,3) = 1.0d0/3.0d0
        ppiclf_rk3coef(2,3) = 2.0d0/3.0d0
        ppiclf_rk3coef(3,3) = dt*2.0d0/3.0d0
        !BD: Original Code END
      end if

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_MarkForRemoval(i)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 i
!
      ppiclf_iprop(1,i) = 3

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_RemoveParticle
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 in_part(PPICLF_LPART), jj(3), iperiodicx, iperiodicy,
     >          iperiodicz,ndim, i, isl, isr, j, jchk, ic
      ! 08/18/24 - Thierry - added for Angular Periodicity begins here
      real*8 per_alpha
      ! 08/18/24 - Thierry - added for Angular Periodicity ends here
!
      iperiodicx = ppiclf_iperiodic(1)
      iperiodicy = ppiclf_iperiodic(2)
      iperiodicz = ppiclf_iperiodic(3)
      ndim       = ppiclf_ndim

      jj(1) = 1
      jj(2) = 2
      jj(3) = 3

      do i=1,ppiclf_npart
        
         isl = (i -1) * PPICLF_LRS + 1
         in_part(i) = 0
         if (ppiclf_iprop(1,i) .eq. 3) then
            in_part(i) = -1 ! User removed particle
            goto 1513
         endif
!-----------------------------------------------------------------------------------
!!!!!!!!!!!!!!!        Rotational Periodicity Starts Here     !!!!!!!!!!!!!!!!!!!!
            ! currently only supports angular rotation around z-axis
            ! and only check theta component and not radial
            ! applying the radial periodicity is very straightforward, but not needed for now

            if(ang_per_flag .eq. 1) then  ! Angular periodicity
           
            ! particle angle w/ x-axis
            ! per_alpha here is obtained in radians
            ! ang_per_angle & ang_per_xangle are transformed 
            !   to radians in PICL_TEMP_InitFlowSolver.F90
            per_alpha = 
     >         atan2(ppiclf_y(PPICLF_JY,i), ppiclf_y(PPICLF_JX,i))

            ! check if particle leaving through lower face or upper face of wedge
            if ((per_alpha .lt. ang_per_xangle) .or. 
     >          (per_alpha .gt. (ang_per_xangle + ang_per_angle))) then
              call ppiclf_solve_InvokeAngularPeriodic(i, ang_per_flag,
     >                                                per_alpha, 
     >                                                ang_per_angle,
     >                                                ang_per_xangle,
     >                                                1)
            endif ! per_alpha
           endif ! ang_per_flag

        ! Linear Periodicity Invoked
          if((x_per_flag.eq.1) .or. (y_per_flag.eq.1) 
     >                         .or. (z_per_flag.eq.1)) then
!------------------------------------------------------------------------------------------------------
! 10/16/2024 - Thierry - this is now implemented in
! ppiclf_solve_InvokeLinearPeriodic(i). Either delete below or comment
! out
!         do j=0,ndim-1
!            jchk = jj(j+1)
!            ! Thierry - checks if particle is about to leave min. periodic face
!            ! moves it linearly relative to the max. periodic face
!            if (ppiclf_y(jchk,i).lt.ppiclf_xdrange(1,j+1))then
!               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
!     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
!     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
!                   ppiclf_y(jchk,i) = ppiclf_xdrange(2,j+1) - 
!     &                     abs(ppiclf_xdrange(1,j+1) - ppiclf_y(jchk,i))
!!                   ppiclf_y1(isl+j)   = ppiclf_xdrange(2,j+1) +
!!     &                     abs(ppiclf_xdrange(1,j+1) - ppiclf_y1(isl+j))
!                  goto 1512
!                endif
!            endif
!            ! Thierry - checks if particle is about to leave max. periodic face
!            ! moves it relative to the min. periodic face
!            if (ppiclf_y(jchk,i).gt.ppiclf_xdrange(2,j+1))then
!               if (((iperiodicx.eq.0) .and. (j.eq.0)) .or.   ! periodic
!     >             ((iperiodicy.eq.0) .and. (j.eq.1)) .or.     
!     >             ((iperiodicz.eq.0) .and. (j.eq.2)) ) then
!                   ppiclf_y(jchk,i) = ppiclf_xdrange(1,j+1) +
!     &                     abs(ppiclf_y(jchk,i) - ppiclf_xdrange(2,j+1))
!!                   ppiclf_y1(isl+j)   = ppiclf_xdrange(1,j+1) +
!!     &                     abs(ppiclf_y1(isl+j) - ppiclf_xdrange(2,j+1))
!                  goto 1512
!                endif
!            endif
!            if (ppiclf_iprop(1,i) .eq. 2) then
!               in_part(i) = -1 ! only if periodic check fails it will get here
!            endif
! 1512 continue
!         enddo ! j=0,ndim-1
!------------------------------------------------------------------------------------------------------
         call ppiclf_solve_InvokeLinearPeriodic(i)
         endif ! x_per_flag
 1513 continue
      enddo ! i=1,ppiclf_part

      ic = 0
      do i=1,ppiclf_npart
         if (in_part(i).eq.0) then
            ic = ic + 1 
            if (i .ne. ic) then
               isl = (i -1) * PPICLF_LRS + 1
               isr = (ic-1) * PPICLF_LRS + 1
               call ppiclf_copy
     >              (ppiclf_y     (1,ic),ppiclf_y(1,i)     ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_y1    (isr) ,ppiclf_y1(isl)    ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_ydot  (1,ic),ppiclf_ydot(1,i)  ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_ydotc (1,ic),ppiclf_ydotc(1,i) ,PPICLF_LRS)
               call ppiclf_copy
     >              (ppiclf_rprop (1,ic),ppiclf_rprop(1,i) ,PPICLF_LRP)
               call ppiclf_copy
     >              (ppiclf_rprop2(1,ic),ppiclf_rprop2(1,i),PPICLF_LRP2)
               call ppiclf_copy
     >              (ppiclf_rprop3(1,ic),ppiclf_rprop3(1,i),PPICLF_LRP3)
               call ppiclf_icopy
     >              (ppiclf_iprop(1,ic) ,ppiclf_iprop(1,i) ,PPICLF_LIP)
            endif
         endif
      enddo
      ppiclf_npart = ic

      return
      end
!----------------------------------------------------------------------
      subroutine ppiclf_solve_FindWallProject(rx2)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
       real*8 rx2(3)
!
! Internal:
!
      real*8 rnx,rny,rnz,rpx1,rpy1,rpz1,rpx2,rpy2,rpz2,rflip,rd,rdist
      integer*4 j, istride
!
      istride = ppiclf_ndim
      ppiclf_nwall_m = 0
      do j = 1,ppiclf_nwall

         rnx  = ppiclf_wall_n(1,j)
         rny  = ppiclf_wall_n(2,j)
         rnz  = 0.0d0
         rpx1 = rx2(1)
         rpy1 = rx2(2)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(1,j)
         rpy2 = ppiclf_wall_c(2,j)
         rpz2 = 0.0d0
         rpx2 = rpx2 - rpx1
         rpy2 = rpy2 - rpy1

         if (ppiclf_ndim .eq. 3) then
            rnz  = ppiclf_wall_n(3,j)
            rpz1 = rx2(3)
            rpz2 = ppiclf_wall_c(3,j)
            rpz2 = rpz2 - rpz1
         endif
    
         rflip = rnx*rpx2 + rny*rpy2 + rnz*rpz2
         if (rflip .gt. 0.0d0) then
            rnx = -1.0d0*rnx
            rny = -1.0d0*rny
            rnz = -1.0d0*rnz
         endif

         rpx1 = ppiclf_wall_c(1,j)
         rpy1 = ppiclf_wall_c(2,j)
         rpz1 = 0.0d0
         rpx2 = ppiclf_wall_c(istride+1,j)
         rpy2 = ppiclf_wall_c(istride+2,j)
         rpz2 = 0.0d0

         if (ppiclf_ndim .eq. 3) then
            rpz1 = ppiclf_wall_c(3,j)
            rpz2 = ppiclf_wall_c(istride+3,j)
         endif

         rd   = -(rnx*rpx1 + rny*rpy1 + rnz*rpz1)

         rdist = abs(rnx*rx2(1)+rny*rx2(2)
     >              +rnz*rx2(3)+rd)
         rdist = rdist/sqrt(rnx**2 + rny**2 + rnz**2)
         rdist = rdist*2.0d0 ! Mirror

         if (rdist .gt. ppiclf_d2chk(2)) goto 1511

         ppiclf_nwall_m = ppiclf_nwall_m + 1

         ppiclf_xyz_mirror(1,ppiclf_nwall_m) = rx2(1) - rdist*rnx
         ppiclf_xyz_mirror(2,ppiclf_nwall_m) = rx2(2) - rdist*rny
         ppiclf_xyz_mirror(3,ppiclf_nwall_m) = 0.0d0
         if (ppiclf_ndim .eq. 3) then
            ppiclf_xyz_mirror(3,ppiclf_nwall_m) = rx2(3) - rdist*rnz
         endif

 1511 continue
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_ProjectParticleGrid
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 iproj(4,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 ppiclf_jxgp,ppiclf_jygp,ppiclf_jzgp
      logical partl, if3d
      integer*4 nkey(2), nxyz, nxyzdum, i, j, k, idum, ic, ip, iip, jjp,
     >          kkp, ilow, ihigh, jlow, jhigh, klow, khigh, ie, jj, j1,
     >          neltbc, ndum, nl, nii, njj, nrr, nlxyzep, iee, ndumdum
      integer*4 jp
      real*8 pi, d2chk2_sq, rdum, multfci, rsig, rdist2, rexp, rx2(3),
     >       rx22, ry22, rz22, rtmp2, evol

      ! Sam - variables for general hex volume calculation
      real*8 v1(3), v2(3), v3(3), cross(3), centroid(3), voltet
      real*8 face(2,2,3) ! ifacex, ifacey, xyz
      integer*4 face_map(3,2,2,2,3) ! xyz,front/back,ifacex,ifacey,ixyz
      integer*4 ix,ix2,iface,inode,ia,ib
!
      if3d = .false.
      if (ppiclf_ndim .eq. 3) if3d = .true.

      PI=4.D0*DATAN(1.D0)

      nxyz = PPICLF_LEX*PPICLF_LEY*PPICLF_LEZ

      nxyzdum = nxyz*PPICLF_LRP_PRO*PPICLF_LEE
      ! TLJ changed loop structure to prevent -fcheck=all error
      ppiclf_pro_fldb = 0.0d0
      !do i=1,nxyzdum
      !   ppiclf_pro_fldb(i,1,1,1,1) = 0.0d0
      !enddo
      !do jp=1,PPICLF_LEE
      !do ip=1,PPICLF_LRP_PRO
      !do k=1,PPICLF_LEZ
      !do j=1,PPICLF_LEY
      !do i=1,PPICLF_LEX
      !   ppiclf_pro_fldb(i,j,k,ip,jp) = 0.0d0
      !enddo
      !enddo
      !enddo
      !enddo
      !enddo

      d2chk2_sq = ppiclf_d2chk(2)**2

      ! real particles
      ppiclf_jxgp  = 1
      ppiclf_jygp  = 2
      ppiclf_jzgp  = 1
      if (if3d) ppiclf_jzgp  = 3

      rdum = 0.0d0
      if (ppiclf_lfiltgauss) then
         rsig    = ppiclf_filter/(2.0d0*sqrt(2.0d0*log(2.0d0)))
         multfci = 1.0d0/(sqrt(2.0d0*pi)**2 * rsig**2) ! in 2D
         if (if3d) multfci = multfci**(1.5d0) ! in 3D
         rdum   = 1.0d0/(-2.0d0*rsig**2)
      endif

      if (ppiclf_lfiltbox) then
         if (ppiclf_sngl_elem) then
           multfci = 1.0d0
           rdum = multfci
         else
           multfci = 1.0d0/(PI/4.0d0*ppiclf_filter**2)
           if (if3d) multfci = multfci/(1.0d0/1.5d0*ppiclf_filter)
           rdum = multfci
         endif
      endif

      ! real particles
      do ip=1,ppiclf_npart

         rproj(1 ,ip) = rdum
         rproj(2 ,ip) = ppiclf_cp_map(ppiclf_jxgp,ip)
         rproj(3 ,ip) = ppiclf_cp_map(ppiclf_jygp,ip)
         if (if3d) 
     >   rproj(4 ,ip) = ppiclf_cp_map(ppiclf_jzgp,ip)

         idum = PPICLF_LRS+PPICLF_LRP
         ic = 4
         ! TLJ modifed loop to remove out of bounds in first index
         !do j=idum+1,idum+PPICLF_LRP_GP
         do j=idum+1,PPICLF_LRP_GP
            ic = ic + 1
            rproj(ic,ip) = ppiclf_cp_map(j,ip)*multfci
         enddo
                    
         iproj(1,ip)  = ppiclf_iprop(8,ip)
         iproj(2,ip)  = ppiclf_iprop(9,ip)
         if (if3d)
     >   iproj(3,ip)  = ppiclf_iprop(10,ip)
         iproj(4,ip)  = ppiclf_iprop(11,ip)
      enddo

      if (.not. ppiclf_sngl_elem) then
  
        ! ghost particles
        do ip=1,ppiclf_npart_gp
  
           rproj(1 ,ip+ppiclf_npart) = rdum
           rproj(2 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jxgp,ip)
           rproj(3 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jygp,ip)
           if (if3d) 
     >     rproj(4 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jzgp,ip)
  
           idum = PPICLF_LRS+PPICLF_LRP
           ic = 4
           ! TLJ modified loop to remove out of bounds in first index
           !do j=idum+1,idum+PPICLF_LRP_GP
           do j=idum+1,PPICLF_LRP_GP
              ic = ic + 1
              rproj(ic,ip+ppiclf_npart) = ppiclf_rprop_gp(j,ip)*multfci
           enddo
                      
           iproj(1,ip+ppiclf_npart)  = ppiclf_iprop_gp(2,ip)
           iproj(2,ip+ppiclf_npart)  = ppiclf_iprop_gp(3,ip)
           if (if3d)
     >     iproj(3,ip+ppiclf_npart)  = ppiclf_iprop_gp(4,ip)
           iproj(4,ip+ppiclf_npart)  = ppiclf_iprop_gp(5,ip)
        enddo
  
        ndum = ppiclf_npart+ppiclf_npart_gp
  
        do ip=1,ndum
           iip      = iproj(1,ip)
           jjp      = iproj(2,ip)
           if (if3d)
     >     kkp      = iproj(3,ip)
           ndumdum  = iproj(4,ip)
  
           ilow  = iip-1
           ihigh = iip+1
           jlow  = jjp-1
           jhigh = jjp+1
           if (if3d) then
              klow  = kkp-1
              khigh = kkp+1
           endif
  
           ! Find if particle near wall and should mirror itself
           if (ppiclf_iwallm .eq. 1) then
              rx2(1) = rproj(2,ip)
              rx2(2) = rproj(3,ip)
              rx2(3) = rproj(4,ip)
              call ppiclf_solve_FindWallProject(rx2)
           endif
  
           do ie=1,ppiclf_neltb
  
                 if (ppiclf_el_map(1,ie) .gt. ndumdum) exit
                 if (ppiclf_el_map(2,ie) .lt. ndumdum) cycle 
           
                 if (ppiclf_el_map(3,ie) .gt. ihigh) cycle
                 if (ppiclf_el_map(4,ie) .lt. ilow)  cycle
                 if (ppiclf_el_map(5,ie) .gt. jhigh) cycle
                 if (ppiclf_el_map(6,ie) .lt. jlow)  cycle
                 if (if3d) then
                    if (ppiclf_el_map(7,ie) .gt. khigh) cycle
                    if (ppiclf_el_map(8,ie) .lt. klow)  cycle
                 endif
  
           do k=1,PPICLF_LEZ
           do j=1,PPICLF_LEY
           do i=1,PPICLF_LEX
              if (ppiclf_modgp(i,j,k,ie,4).ne.ndumdum) cycle
  
              rdist2  = (ppiclf_xm1b(i,j,k,1,ie) - rproj(2,ip))**2 +
     >                  (ppiclf_xm1b(i,j,k,2,ie) - rproj(3,ip))**2
              if(if3d) rdist2 = rdist2 +
     >                  (ppiclf_xm1b(i,j,k,3,ie) - rproj(4,ip))**2
  
              if (rdist2 .gt. d2chk2_sq) cycle
  
              rexp = 1.0d0  ! for box filter
              if (ppiclf_lfiltgauss)
     >           rexp = exp(rdist2*rproj(1,ip))
  
              ! add wall effects
              if (ppiclf_iwallm .eq. 1) then
                 do jj=1,ppiclf_nwall_m
                    rx22 = (ppiclf_xm1b(i,j,k,1,ie) 
     >                     -ppiclf_xyz_mirror(1,jj))**2
                    ry22 = (ppiclf_xm1b(i,j,k,2,ie)
     >                     -ppiclf_xyz_mirror(2,jj))**2
                    rtmp2 = rx22 + ry22
                    if (if3d) then
                       rz22 = (ppiclf_xm1b(i,j,k,3,ie)
     >                        -ppiclf_xyz_mirror(3,jj))**2
                       rtmp2 = rtmp2 + rz22
                    endif
                    if (ppiclf_lfiltgauss) then
                       rexp = rexp + exp(rtmp2*rproj(1,ip))
                    else
                       rexp = rexp + 1.0d0
                    endif
                 enddo
              endif
  
              
              do jj=1,PPICLF_LRP_PRO
                 j1 = jj+4
                 ppiclf_pro_fldb(i,j,k,jj,ie) = 
     >                           ppiclf_pro_fldb(i,j,k,jj,ie) 
     >                         + rproj(j1,ip)*rexp
              enddo
           enddo
           enddo
           enddo
           enddo
        enddo
      ! sngl elem
      else

!BD: This is where rproj gets stored in fld, this is one of the steps
!that should be tracked

        ! Sam - build map to face coordinates for 3D general hex volume
        ! calculation
        do ix=1,3
          ia = max(3-ix,1)
          ib = min(5-ix,3)
          do iface=1,2
            do j=1,2
            do i=1,2
              face_map(ix,iface,i,j,ix) = iface ! constant
              face_map(ix,iface,i,j,ia) = i ! ix
              face_map(ix,iface,i,j,ib) = j ! iy
            end do
            end do
          end do
        end do

        do ip=1,ppiclf_npart
           !do ie=1,ppiclf_neltb
  
             !if (ie .ne. ppiclf_iprop(2,ip)+1) cycle
             ie = ppiclf_iprop(2, ip) + 1
             if ((ie .lt. 1) .or. (ie .gt. ppiclf_neltb)) cycle

             ! Sam - general hexahedron volume calculation
             if (if3d) then
               ! get centroid of hexahedron
               do ix=1,3
                 centroid(ix) = 0.0
               end do

               do ix=1,3
               do k=1,PPICLF_LEZ
               do j=1,PPICLF_LEY
               do i=1,PPICLF_LEX
                 centroid(ix) = centroid(ix) + ppiclf_xm1b(i,j,k,ix,ie)
               end do
               end do
               end do
               end do

               do ix=1,3
                 centroid(ix) = centroid(ix) / 8.0
               end do


               ! calculate volume based on two contributions from each
               ! face as tetrahedrons
               evol = 0.0
               do ix=1,3
                 do iface=1,2

                   ! get face coordinates
                   do j=1,2
                   do i=1,2
                   do ix2=1,3
                     face(i,j,ix2) = ppiclf_xm1b(
     >                                 face_map(ix,iface,i,j,1),
     >                                 face_map(ix,iface,i,j,2),
     >                                 face_map(ix,iface,i,j,3),
     >                                 ix2,ie)
                   end do
                   end do
                   end do

                   do ix2=1,3
                     v1(ix2) = face(1,2,ix2) - face(2,1,ix2)
                     v2(ix2) = centroid(ix2) - face(2,1,ix2)
                   end do ! ix2

                   ! take cross product
                   cross(1) = v1(2)*v2(3) - v1(3)*v2(2)
                   cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
                   cross(3) = v1(1)*v2(2) - v1(2)*v2(1)

                   ! get contriubtions to volume from each tetrahedron
                   do inode=1,2
                   do ix2=1,3
                     v3(ix2) = face(inode,inode,ix2) - face(2,1,ix2)
                   end do ! ix2

                   ! really 6 times the volume of the tet, but we can
                   ! save an operation by dividing at the end
                   voltet = 0.0
                   do ix2=1,3
                     voltet = voltet + v3(ix2)*cross(ix2)
                   end do ! ix2
                   evol = evol + abs(voltet)
                   end do ! inode
                   
                 end do ! iface
              end do ! ix
               evol = evol / 6.0
             else
               ! Sam - default to naive solution for 2D. ASSUMES
               ! rectangular elements. This will
               ! probably never get used, but if it does throw an error
               ! so the user is absolutely sure of what they're doing.
!               call ppiclf_exittr('Single element projection only
!     >          supported in 3D for general hex elements. Comment and
!     >          ignore this error if your elements are perfect
!     >          rectangles. $',0.0d0,0)

               evol = (ppiclf_xm1b(PPICLF_LEX,1,1,1,ie) 
     >               - ppiclf_xm1b(1,1,1,1,ie))
               evol = evol
     >              * (ppiclf_xm1b(1,PPICLF_LEY,1,2,ie) 
     >               - ppiclf_xm1b(1,1,1,2,ie))
!               if (if3d) evol = evol
!     >              * (ppiclf_xm1b(1,1,PPICLF_LEZ,3,ie) 
!     >               - ppiclf_xm1b(1,1,1,3,ie))
            end if ! if3d

             rexp = 1.0 / evol
           do k=1,PPICLF_LEZ
           do j=1,PPICLF_LEY
           do i=1,PPICLF_LEX
              do jj=1,PPICLF_LRP_PRO
                 j1 = jj+4
                 ppiclf_pro_fldb(i,j,k,jj,ie) = 
     >                           ppiclf_pro_fldb(i,j,k,jj,ie) 
     >                         + rproj(j1,ip)*rexp
              enddo
           enddo
           enddo
           enddo
           enddo
        !enddo ! ppiclf_neltb
      endif ! ppiclf_npart

      ! now send xm1b to the processors in nek that hold xm1

      neltbc = ppiclf_neltb
      ndum = PPICLF_LRMAX*neltbc
      call ppiclf_icopy(ppiclf_er_mapc,ppiclf_er_map,ndum)
      do ie=1,neltbc
         ppiclf_er_mapc(5,ie) = ppiclf_er_mapc(2,ie)
         ppiclf_er_mapc(6,ie) = ppiclf_er_mapc(2,ie)
      enddo
      nl = 0
      nii = PPICLF_LRMAX
      njj = 6
      nrr = nxyz*PPICLF_LRP_PRO
      nkey(1) = 2
      nkey(2) = 1
      call pfgslib_crystal_tuple_transfer(ppiclf_cr_hndl,neltbc,
     >   PPICLF_LEE,ppiclf_er_mapc,nii,partl,nl,ppiclf_pro_fldb,nrr,njj)
      call pfgslib_crystal_tuple_sort    (ppiclf_cr_hndl,neltbc
     >       ,ppiclf_er_mapc,nii,partl,nl,ppiclf_pro_fldb,nrr,nkey,2)

      ! add the fields from the bins to ptw array
      nlxyzep = nxyz*PPICLF_LEE*PPICLF_LRP_PRO
      ! TLJ changed looping to prevent -fcheck=all error
      ppiclf_pro_fld = 0.0d0
      !do i=1,nlxyzep
      !   ppiclf_pro_fld(i,1,1,1,1) = 0.0d0
      !enddo
      !do jp=1,PPICLF_LRP_PRO
      !do ip=1,PPICLF_LEE
      !do k=1,PPICLF_LEZ
      !do j=1,PPICLF_LEY
      !do i=1,PPICLF_LEX
      !   ppiclf_pro_fld(i,j,k,ip,jp) = 0.0d0
      !enddo
      !enddo
      !enddo
      !enddo
      !enddo


      do ie=1,neltbc
         iee = ppiclf_er_mapc(1,ie)
         do ip=1,PPICLF_LRP_PRO
         do k=1,PPICLF_LEZ
         do j=1,PPICLF_LEY
         do i=1,PPICLF_LEX
           ppiclf_pro_fld(i,j,k,iee,ip) = ppiclf_pro_fld(i,j,k,iee,ip) +
     >                                    ppiclf_pro_fldb(i,j,k,ip,ie)

         enddo
         enddo
         enddo
         enddo
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine ppiclf_solve_ProjectParticleSubBin
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      real*8    rproj(1+PPICLF_LRP_GP,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 iproj(4,PPICLF_LPART+PPICLF_LPART_GP)
      integer*4 ppiclf_jxgp,ppiclf_jygp,ppiclf_jzgp, nxyz, nxyzdum,
     >          idum, jdum, kdum, ic, i, j, k, ip, ndum, il, ir, jl, jr,
     >          kl, kr, jj, j1, iip, jjp, kkp
      logical if3d
      real*8 pi, d2chk2_sq, rdum, rsig, multfci, rexp, rdist2
!

      if3d = .false.
      if (ppiclf_ndim .eq. 3) if3d = .true.

      PI=4.D0*DATAN(1.D0)

      nxyz = PPICLF_BX1*PPICLF_BY1*PPICLF_BZ1

      nxyzdum = nxyz*PPICLF_LRP_PRO
      do i=1,nxyzdum
         ppiclf_grid_fld(i,1,1,1) = 0.0d0
      enddo

      d2chk2_sq = ppiclf_d2chk(2)**2

      ! real particles
      ppiclf_jxgp  = 1
      ppiclf_jygp  = 2
      ppiclf_jzgp  = 1
      if (if3d)
     >ppiclf_jzgp  = 3

      rdum = 0.0d0
      if (ppiclf_lfiltgauss) then
         rsig    = ppiclf_filter/(2.0d0*sqrt(2.0d0*log(2.0d0)))
         multfci = 1.0d0/(sqrt(2.0d0*pi)**2 * rsig**2) 
         if (if3d) multfci = multfci**(1.5d0)
         rdum   = 1.0d0/(-2.0d0*rsig**2)
      endif

      if (ppiclf_lfiltbox) then
         multfci = 1.0d0/(PI/4.0d0*ppiclf_filter**2)
         if (if3d) multfci = multfci/(1.0d0/1.5d0*ppiclf_filter)
      endif

      ! real particles
      do ip=1,ppiclf_npart

         rproj(1 ,ip) = rdum
         rproj(2 ,ip) = ppiclf_cp_map(ppiclf_jxgp,ip)
         rproj(3 ,ip) = ppiclf_cp_map(ppiclf_jygp,ip)
         if (if3d)
     >   rproj(4 ,ip) = ppiclf_cp_map(ppiclf_jzgp,ip)

         idum = PPICLF_LRS+PPICLF_LRP
         ic = 4
         ! TLJ modifed loop to remove out of bounds in first index
         !do j=idum+1,idum+PPICLF_LRP_GP
         do j=idum+1,PPICLF_LRP_GP
            ic = ic + 1
            rproj(ic,ip) = ppiclf_cp_map(j,ip)*multfci
         enddo

         iproj(1,ip) = 
     >       floor( (rproj(2,ip) - ppiclf_binx(1,1))/ppiclf_rdx)
         iproj(2,ip) = 
     >       floor( (rproj(3,ip) - ppiclf_biny(1,1))/ppiclf_rdy)
         if (if3d)
     >   iproj(3,ip) = 
     >       floor( (rproj(4,ip) - ppiclf_binz(1,1))/ppiclf_rdz)
      enddo

      ! ghost particles
      do ip=1,ppiclf_npart_gp

         rproj(1 ,ip+ppiclf_npart) = rdum
         rproj(2 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jxgp,ip)
         rproj(3 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jygp,ip)
         if (if3d)
     >   rproj(4 ,ip+ppiclf_npart) = ppiclf_rprop_gp(ppiclf_jzgp,ip)

         idum = PPICLF_LRS+PPICLF_LRP
         ic = 4
         ! TLJ modifed loop to remove out of bounds in first index
         !do j=idum+1,idum+PPICLF_LRP_GP
         do j=idum+1,PPICLF_LRP_GP
            ic = ic + 1
            rproj(ic,ip+ppiclf_npart) = ppiclf_rprop_gp(j,ip)*multfci
         enddo
                    
         iproj(1,ip+ppiclf_npart) = 
     >     floor((rproj(2,ip+ppiclf_npart)-ppiclf_binx(1,1))/ppiclf_rdx)
         iproj(2,ip+ppiclf_npart) = 
     >     floor((rproj(3,ip+ppiclf_npart)-ppiclf_biny(1,1))/ppiclf_rdy)
         if (if3d)
     >   iproj(3,ip+ppiclf_npart) = 
     >     floor((rproj(4,ip+ppiclf_npart)-ppiclf_binz(1,1))/ppiclf_rdz)
      enddo

      ndum = ppiclf_npart+ppiclf_npart_gp

      if (ppiclf_lfiltgauss) then
         idum = floor(ppiclf_filter/2.0d0/ppiclf_rdx
     >    *sqrt(-log(ppiclf_alpha)/log(2.0d0)))+1
         jdum = floor(ppiclf_filter/2.0d0/ppiclf_rdy
     >    *sqrt(-log(ppiclf_alpha)/log(2.0d0)))+1
         kdum = 999999999
         if (if3d)
     >   kdum = floor(ppiclf_filter/2.0d0/ppiclf_rdz
     >    *sqrt(-log(ppiclf_alpha)/log(2.0d0)))+1
      endif

      if (ppiclf_lfiltbox) then
         idum = ppiclf_ngrids/2+1
         jdum = ppiclf_ngrids/2+1
         kdum = 999999999
         if (if3d)
     >   kdum = ppiclf_ngrids/2+1
      endif

      do ip=1,ndum
         iip = iproj(1,ip)
         jjp = iproj(2,ip)
         if (if3d)
     >   kkp = iproj(3,ip)

         il  = max(1     ,iip-idum)
         ir  = min(ppiclf_bx,iip+idum)
         jl  = max(1     ,jjp-jdum)
         jr  = min(ppiclf_by,jjp+jdum)
         kl  = 1
         kr  = 1
         if (if3d) then
         kl  = max(1     ,kkp-kdum)
         kr  = min(ppiclf_bz,kkp+kdum)
         endif

c        do k=kl,kr
c        do j=jl,jr
c        do i=il,ir
         do k=1,ppiclf_bz
         do j=1,ppiclf_by
         do i=1,ppiclf_bx
            rdist2  = (ppiclf_grid_x(i,j,k) - rproj(2,ip))**2 +
     >                (ppiclf_grid_y(i,j,k) - rproj(3,ip))**2
            if(if3d) rdist2 = rdist2 +
     >                (ppiclf_grid_z(i,j,k) - rproj(4,ip))**2

            if (rdist2 .gt. d2chk2_sq) cycle

            rexp = 1.0d0
            if (ppiclf_lfiltgauss)
     >         rexp = exp(rdist2*rproj(1,ip))

            do jj=1,PPICLF_LRP_PRO
               j1 = jj+4
               ppiclf_grid_fld(i,j,k,jj) = 
     >                         ppiclf_grid_fld(i,j,k,jj) 
     >                       + sngl(rproj(j1,ip)*rexp)
            enddo
         enddo
         enddo
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_solve_GetProFldIJKEF(i,j,k,e,m,fld)
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 i,j,k,e,m
!
! Output:
!
      real*8 fld
!
      fld = ppiclf_pro_fld(i,j,k,e,m)

      return
      end
!-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine pfgslib_userExitHandler()
!
      implicit none
!
! This has been added so when linked with another code that links
! to gslib, there are no naming conflicts. See USREXIT=1 flag in
! install script for gslib.
!
      return
      end
c-----------------------------------------------------------------------
      subroutine pfgslib_mxm(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3
      real*8 a(n1,n2),b(n2,n3),c(n1,n3)
!
      call ppiclf_mxmf2(a,n1,b,n2,c,n3)

      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxmf2(A,N1,B,N2,C,N3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3
      real*8 a(n1,n2),b(n2,n3),c(n1,n3)
!
      if (n2.le.8) then
         if (n2.eq.1) then
            call ppiclf_mxf1(a,n1,b,n2,c,n3)
         elseif (n2.eq.2) then
            call ppiclf_mxf2(a,n1,b,n2,c,n3)
         elseif (n2.eq.3) then
            call ppiclf_mxf3(a,n1,b,n2,c,n3)
         elseif (n2.eq.4) then
            call ppiclf_mxf4(a,n1,b,n2,c,n3)
         elseif (n2.eq.5) then
            call ppiclf_mxf5(a,n1,b,n2,c,n3)
         elseif (n2.eq.6) then
            call ppiclf_mxf6(a,n1,b,n2,c,n3)
         elseif (n2.eq.7) then
            call ppiclf_mxf7(a,n1,b,n2,c,n3)
         else
            call ppiclf_mxf8(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.16) then
         if (n2.eq.9) then
            call ppiclf_mxf9(a,n1,b,n2,c,n3)
         elseif (n2.eq.10) then
            call ppiclf_mxf10(a,n1,b,n2,c,n3)
         elseif (n2.eq.11) then
            call ppiclf_mxf11(a,n1,b,n2,c,n3)
         elseif (n2.eq.12) then
            call ppiclf_mxf12(a,n1,b,n2,c,n3)
         elseif (n2.eq.13) then
            call ppiclf_mxf13(a,n1,b,n2,c,n3)
         elseif (n2.eq.14) then
            call ppiclf_mxf14(a,n1,b,n2,c,n3)
         elseif (n2.eq.15) then
            call ppiclf_mxf15(a,n1,b,n2,c,n3)
         else
            call ppiclf_mxf16(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.24) then
         if (n2.eq.17) then
            call ppiclf_mxf17(a,n1,b,n2,c,n3)
         elseif (n2.eq.18) then
            call ppiclf_mxf18(a,n1,b,n2,c,n3)
         elseif (n2.eq.19) then
            call ppiclf_mxf19(a,n1,b,n2,c,n3)
         elseif (n2.eq.20) then
            call ppiclf_mxf20(a,n1,b,n2,c,n3)
         elseif (n2.eq.21) then
            call ppiclf_mxf21(a,n1,b,n2,c,n3)
         elseif (n2.eq.22) then
            call ppiclf_mxf22(a,n1,b,n2,c,n3)
         elseif (n2.eq.23) then
            call ppiclf_mxf23(a,n1,b,n2,c,n3)
         elseif (n2.eq.24) then
            call ppiclf_mxf24(a,n1,b,n2,c,n3)
         endif
      else
         call ppiclf_mxm44_0(a,n1,b,n2,c,n3)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf1(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,1),b(1,n3),c(n1,n3)
!
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf2(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,2),b(2,n3),c(n1,n3)
!
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf3(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,3),b(3,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf4(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,4),b(4,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf5(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,5),b(5,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf6(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,6),b(6,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf7(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,7),b(7,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf8(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,8),b(8,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf9(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,9),b(9,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf10(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,10),b(10,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf11(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,11),b(11,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf12(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,12),b(12,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf13(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,13),b(13,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf14(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,14),b(14,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf15(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,15),b(15,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf16(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,16),b(16,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf17(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,17),b(17,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf18(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,18),b(18,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf19(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,19),b(19,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf20(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,20),b(20,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf21(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,21),b(21,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf22(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,22),b(22,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf23(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,23),b(23,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxf24(a,n1,b,n2,c,n3)
!
      implicit none
!
! Internal:
!
      integer*4 n1, n2, n3, i, j
      real*8 a(n1,24),b(24,n3),c(n1,n3)
c
      n1 = n1
      n2 = n2
      n3 = n3
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
     $             + a(i,17)*b(17,j)
     $             + a(i,18)*b(18,j)
     $             + a(i,19)*b(19,j)
     $             + a(i,20)*b(20,j)
     $             + a(i,21)*b(21,j)
     $             + a(i,22)*b(22,j)
     $             + a(i,23)*b(23,j)
     $             + a(i,24)*b(24,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ppiclf_mxm44_0(a, m, b, k, c, n)
!
      implicit none
!
! Internal:
!
      integer*4 m, k, n, i, j, l, m1, n1, nresid, mresid
      real*8 a(m,k), b(k,n), c(m,n)
      real*8 s11, s12, s13, s14, s21, s22, s23, s24
      real*8 s31, s32, s33, s34, s41, s42, s43, s44
!
      mresid = iand(m,3) 
      nresid = iand(n,3) 
      m1 = m - mresid + 1
      n1 = n - nresid + 1

      do i=1,m-mresid,4
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s41 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s42 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          s43 = 0.0d0
          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0
          s44 = 0.0d0
          do l=1,k
            s11 = s11 + a(i,l)*b(l,j)
            s12 = s12 + a(i,l)*b(l,j+1)
            s13 = s13 + a(i,l)*b(l,j+2)
            s14 = s14 + a(i,l)*b(l,j+3)

            s21 = s21 + a(i+1,l)*b(l,j)
            s22 = s22 + a(i+1,l)*b(l,j+1)
            s23 = s23 + a(i+1,l)*b(l,j+2)
            s24 = s24 + a(i+1,l)*b(l,j+3)

            s31 = s31 + a(i+2,l)*b(l,j)
            s32 = s32 + a(i+2,l)*b(l,j+1)
            s33 = s33 + a(i+2,l)*b(l,j+2)
            s34 = s34 + a(i+2,l)*b(l,j+3)

            s41 = s41 + a(i+3,l)*b(l,j)
            s42 = s42 + a(i+3,l)*b(l,j+1)
            s43 = s43 + a(i+3,l)*b(l,j+2)
            s44 = s44 + a(i+3,l)*b(l,j+3)
          enddo
          c(i,j)     = s11 
          c(i,j+1)   = s12 
          c(i,j+2)   = s13
          c(i,j+3)   = s14

          c(i+1,j)   = s21 
          c(i+2,j)   = s31 
          c(i+3,j)   = s41 

          c(i+1,j+1) = s22
          c(i+2,j+1) = s32
          c(i+3,j+1) = s42

          c(i+1,j+2) = s23
          c(i+2,j+2) = s33
          c(i+3,j+2) = s43

          c(i+1,j+3) = s24
          c(i+2,j+3) = s34
          c(i+3,j+3) = s44
        enddo
* Residual when n is not multiple of 4
        if (nresid .ne. 0) then
          if (nresid .eq. 1) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,n)
              s21 = s21 + a(i+1,l)*b(l,n)
              s31 = s31 + a(i+2,l)*b(l,n)
              s41 = s41 + a(i+3,l)*b(l,n)
            enddo
            c(i,n)     = s11 
            c(i+1,n)   = s21 
            c(i+2,n)   = s31 
            c(i+3,n)   = s41 
          elseif (nresid .eq. 2) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
            enddo
            c(i,j)     = s11 
            c(i,j+1)   = s12

            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 

            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
          else
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0
            s43 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)
              s13 = s13 + a(i,l)*b(l,j+2)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)
              s23 = s23 + a(i+1,l)*b(l,j+2)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)
              s33 = s33 + a(i+2,l)*b(l,j+2)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
              s43 = s43 + a(i+3,l)*b(l,j+2)
            enddo
            c(i,j)     = s11 
            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 
            c(i,j+1)   = s12 
            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
            c(i,j+2)   = s13
            c(i+1,j+2) = s23
            c(i+2,j+2) = s33
            c(i+3,j+2) = s43
          endif
        endif
      enddo

* Residual when m is not multiple of 4
      if (mresid .eq. 0) then
        return
      elseif (mresid .eq. 1) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,j)
            s12 = s12 + a(m,l)*b(l,j+1)
            s13 = s13 + a(m,l)*b(l,j+2)
            s14 = s14 + a(m,l)*b(l,j+3)
          enddo
          c(m,j)     = s11 
          c(m,j+1)   = s12 
          c(m,j+2)   = s13
          c(m,j+3)   = s14
        enddo
* mresid is 1, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n)
          enddo
          c(m,n) = s11
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s12 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-1)
            s12 = s12 + a(m,l)*b(l,n)
          enddo
          c(m,n-1) = s11
          c(m,n) = s12
          return
        else
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-2)
            s12 = s12 + a(m,l)*b(l,n-1)
            s13 = s13 + a(m,l)*b(l,n)
          enddo
          c(m,n-2) = s11
          c(m,n-1) = s12
          c(m,n) = s13
          return
        endif          
      elseif (mresid .eq. 2) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          s21 = 0.0d0
          s22 = 0.0d0
          s23 = 0.0d0
          s24 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,j)
            s12 = s12 + a(m-1,l)*b(l,j+1)
            s13 = s13 + a(m-1,l)*b(l,j+2)
            s14 = s14 + a(m-1,l)*b(l,j+3)

            s21 = s21 + a(m,l)*b(l,j)
            s22 = s22 + a(m,l)*b(l,j+1)
            s23 = s23 + a(m,l)*b(l,j+2)
            s24 = s24 + a(m,l)*b(l,j+3)
          enddo
          c(m-1,j)   = s11 
          c(m-1,j+1) = s12 
          c(m-1,j+2) = s13
          c(m-1,j+3) = s14
          c(m,j)     = s21
          c(m,j+1)   = s22 
          c(m,j+2)   = s23
          c(m,j+3)   = s24
        enddo
* mresid is 2, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n)
          enddo
          c(m-1,n) = s11
          c(m,n) = s21
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-1)
            s12 = s12 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-1)
            s22 = s22 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-1) = s11
          c(m-1,n)   = s12
          c(m,n-1)   = s21
          c(m,n)     = s22
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-2)
            s12 = s12 + a(m-1,l)*b(l,n-1)
            s13 = s13 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-2)
            s22 = s22 + a(m,l)*b(l,n-1)
            s23 = s23 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-2) = s11
          c(m-1,n-1) = s12
          c(m-1,n)   = s13
          c(m,n-2)   = s21
          c(m,n-1)   = s22
          c(m,n)     = s23
          return
        endif
      else
* mresid is 3
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0

          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0

          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0

          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0

          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,j)
            s12 = s12 + a(m-2,l)*b(l,j+1)
            s13 = s13 + a(m-2,l)*b(l,j+2)
            s14 = s14 + a(m-2,l)*b(l,j+3)

            s21 = s21 + a(m-1,l)*b(l,j)
            s22 = s22 + a(m-1,l)*b(l,j+1)
            s23 = s23 + a(m-1,l)*b(l,j+2)
            s24 = s24 + a(m-1,l)*b(l,j+3)

            s31 = s31 + a(m,l)*b(l,j)
            s32 = s32 + a(m,l)*b(l,j+1)
            s33 = s33 + a(m,l)*b(l,j+2)
            s34 = s34 + a(m,l)*b(l,j+3)
          enddo
          c(m-2,j)   = s11 
          c(m-2,j+1) = s12 
          c(m-2,j+2) = s13
          c(m-2,j+3) = s14

          c(m-1,j)   = s21 
          c(m-1,j+1) = s22
          c(m-1,j+2) = s23
          c(m-1,j+3) = s24

          c(m,j)     = s31 
          c(m,j+1)   = s32
          c(m,j+2)   = s33
          c(m,j+3)   = s34
        enddo
* mresid is 3, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n)
          enddo
          c(m-2,n) = s11
          c(m-1,n) = s21
          c(m,n) = s31
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-1)
            s12 = s12 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-1)
            s22 = s22 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-1)
            s32 = s32 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-1) = s11
          c(m-2,n)   = s12
          c(m-1,n-1) = s21
          c(m-1,n)   = s22
          c(m,n-1)   = s31
          c(m,n)     = s32
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-2)
            s12 = s12 + a(m-2,l)*b(l,n-1)
            s13 = s13 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-2)
            s22 = s22 + a(m-1,l)*b(l,n-1)
            s23 = s23 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-2)
            s32 = s32 + a(m,l)*b(l,n-1)
            s33 = s33 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-2) = s11
          c(m-2,n-1) = s12
          c(m-2,n)   = s13
          c(m-1,n-2) = s21
          c(m-1,n-1) = s22
          c(m-1,n)   = s23
          c(m,n-2)   = s31
          c(m,n-1)   = s32
          c(m,n)     = s33
          return
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
