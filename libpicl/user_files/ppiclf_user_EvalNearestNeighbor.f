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
