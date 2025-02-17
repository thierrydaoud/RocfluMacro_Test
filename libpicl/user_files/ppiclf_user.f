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
