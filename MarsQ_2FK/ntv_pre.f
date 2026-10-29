      module NTV_pre
C*****************************************************************************
C   calculate the spectrum of the  magnetic field strength
C       variation for NTV module
C   calculate the collisionality and related frequencies
C
C
C      In this part, we read the data from the MARS-F output
C         
C
C          by Y. Sun    
C          Mar. 01, 2011
C*****************************************************************************
      use NTV_lib
      implicit none
!---------- parameters
      real*8,parameter::npi=3.1415926535897932
      contains

      subroutine MARS_out(nr,n_chi,n_phi,n_prof,keyntv,rho,
     &                    mui,Z,Zi,knc,m_max,n_max,R_axis,nu_prof,
     &                    bc,bs,bc2,bs2)
C****************************************************************
C           load data from MARS-F output
C
C   input:
C     rho(nr) is the radii (s coordinate) will evaluate the NTV torque.
C     n_chi,n_phi,are the number of grids for poloidal and toroidal angle.
C     n_prof, is the plasma profile number 
C     mi=mui*mp, mui=2 for D
C     Z, is the charge number of the particle
C     Zi, is the charge number of ion
C     knc=1.17, is neoclasical coeffiecient for poloidal rotation
C     m_max, is the maximum poloidal mode number, m=m(-m_max:m_max)
C     n_max, is the maximum toroidal mode number, n=n(1:n_max)
C
C   output:
C
C    R_axis is the major radius of the magnetic axis
C   
C    bc and bs are the sin and cos component of the spetrum bmn
C      
C   nu_prof(nr,n_prof)=
C           (rho,q,eps,rho_i,nu,wti,wB0,qwe,qws,qwsT,|dB|_surf_average/B0)
C             1  2  3    4    5   6   7   8  9  10   11
C
C      rho is the radius, q is the safty factor, eps=epsillon~ r/R
C      rho_i=n_i*m_i, nu is the collisionality, wti is the trasit frequency
C      wB0=omega_B0 is the grad B drift of the deeply trapped pariticles with x=1,
C      qwe=q*omega_E is the E \times B drift
C      qws=q*omega_*p, is the diamagnetic frequency
C      qwsT=q*omega_*T,
C
C****************************************************************
      character*(30) f1
      integer M0,NP,NV,i,j,n_chi,n_prof,n_phi,nr,m_max,n_max,keyntv
      integer Z,Zi,KROTWE
      real R0,B0,dum,R_axis,rho(nr),knc
      real,allocatable::R_mars(:,:),Z_mars(:,:)
      real,allocatable::dpsids(:),FF(:),s(:),s_prof(:,:),omega(:)
      real,allocatable::dwds(:),q(:),gg(:,:,:),dsdpsi(:),dBav(:)
      real,allocatable::dBds(:,:),dBdc(:,:),chiH(:,:),dphiH(:,:)
      real,allocatable::Jac(:,:),Beq(:,:),eps(:),tmpB(:),tmpgg(:,:)
      real,allocatable::dm(:,:,:),bth(:),cp(:),pp(:)
      real qdev
      real min_t,max_t,ieps,co,mui
      complex,allocatable::bmns(:,:,:)
      real,allocatable::tmpbmn(:),tmpbmn2(:)
      real,allocatable::tmps(:),tmps2(:)
      real nu_prof(nr,n_prof)
      real bc(2*m_max+1,n_max,nr),bs(2*m_max+1,n_max,nr)
      real bc2(2*m_max+1,n_max,nr),bs2(2*m_max+1,n_max,nr)
      real,allocatable::chi(:),tmpbc(:),tmpr(:),r_prof(:,:)

!------- define chi
      min_t=0.0
      max_t=2.0*npi
      allocate(chi(n_chi)) 
      call xgrid(n_chi,min_t,max_t,chi)

      call get_dimens(NP,R0,B0)

      allocate(tmpbc(nr),tmpr(nr))
      allocate(r_prof(nr,n_prof))

      allocate(R_mars(NP,n_chi))  ! only inside the plasma
      allocate(Z_mars(NP,n_chi))

      call get_RZ(NP,n_chi,chi,R_mars,Z_mars)

C     print*,'R0=',R0,'  B0=',B0
C     co=B0  !old version, incorrect
C     correction following Youwen's e-mail:

      co=1.0

      R_axis=sum(R_mars(1,:))/n_chi*R0
C     print*, 'R_c=',R_axis,'co=',co

!---------------- load plasma profile data
            allocate(s(NP))
            allocate(q(NP))
            allocate(dpsids(NP))
            allocate(FF(NP))
            allocate(omega(NP))
            allocate(dwds(NP))
            allocate(s_prof(NP,n_prof))
       call plasma_prof(KROTWE,NP,n_prof,s_prof,s,dpsids,FF,omega,dwds)
            q=s_prof(:,2)

!---------------- calculate metrics 
            allocate(gg(NP,n_chi,3))
            allocate(dBdc(NP,n_chi))
            allocate(dBds(NP,n_chi))
            allocate(Jac(NP,n_chi))
            allocate(Beq(NP,n_chi))
       call cal_metric(NP,n_chi,s,chi,R_mars,Z_mars,FF,dpsids,q,
     .                  Jac,gg,Beq,dBds,dBdc)
!---------------- div(xi_perp) metric, B^theta and mu0 p'(psi) psi'(s)
!                 for the Hamada-gauge parallel displacement
            allocate(dm(NP,n_chi,7))
            allocate(bth(NP))
            allocate(cp(NP))
            allocate(pp(NP))
       call cal_divmetric(NP,n_chi,s,chi,R_mars,Z_mars,FF,dpsids,
     .                    Jac,Beq,dm)
       call get_ppeq(NP,pp)
       bth=0.0
       qdev=0.0
       do i=2,NP
          bth(i)=n_chi*dpsids(i)/sum(Jac(i,:))
          qdev=max(qdev,abs(FF(i)*sum(Jac(i,:)/R_mars(i,:)**2)
     .         /(n_chi*dpsids(i)*q(i))-1.0))
       end do
       cp=pp*dpsids
       print*,'NTV Hamada chart: max|F<J/R^2>/(psi''q)-1| =',qdev
!---------------- coordinate transform
            allocate(chiH(NP,n_chi))
            allocate(dphiH(NP,n_chi))
       call coord_trans(NP,n_chi,chi,R_mars,q,Jac,
     .                                chiH,dphiH)
!       call coord_map(n_chi,chi,NP,chiH,R_mars)
!       call coord_map(n_chi,chi,NP,chiH,Z_mars)
!---------------- calculate epsillon
       call coord_map(n_chi,chi,NP,chiH,Beq)
              allocate(eps(NP))
              allocate(tmpB(n_chi))         
       do i=1,NP
          tmpB=Beq(i,:)
          call cal_eps(n_chi,chi,tmpB,ieps)
          eps(i)=ieps 
       end do
       s_prof(:,10)=eps
!---------------- calculate collisionalities
            allocate(tmps(NP))
            allocate(tmps2(NP))
        do i=2,n_prof
           tmps=s_prof(:,i)
           call spline1d(tmpr,rho,nr,tmps,s,NP,tmps2)
           r_prof(:,i)=tmpr
        end do
              r_prof(:,1)=rho
        dpsids=dpsids*(R0**2*B0)    
            allocate(dsdpsi(nr))
            call spline1d(dsdpsi,rho,nr,dpsids,s,NP,tmps2) 
          dsdpsi=1.0/dsdpsi
       call cal_nu(KROTWE,nr,rho,dsdpsi,n_prof,r_prof,mui,R_axis,B0,
     .             Z,Zi,knc,nu_prof)
!---------------- displacement
C      print*,'Calculating spectrum...'
       call coord_map(n_chi,chi,NP,chiH,dBds)
       call coord_map(n_chi,chi,NP,chiH,dBdc)

       allocate(tmpgg(NP,n_chi))
       do i=1,3
          tmpgg=gg(:,:,i)
          call coord_map(n_chi,chi,NP,chiH,tmpgg)
          gg(:,:,i)=tmpgg
       end do
       do i=1,7
          tmpgg=dm(:,:,i)
          call coord_map(n_chi,chi,NP,chiH,tmpgg)
          dm(:,:,i)=tmpgg
       end do

       allocate(bmns(2*m_max+1,n_max,NP))
       allocate(dBav(NP))
       call cal_xi_liu(NP,n_chi,chiH,n_phi,dphiH,dBav,
     .              gg,dBds,dBdc,omega,dwds,m_max,n_max,bmns,
     .              s,q,Beq,dm,bth,cp)
         bmns=co*bmns
!           print*,'spec',abs(co*bmns(m_max+2,1,:))
       allocate(tmpbmn(NP))
       allocate(tmpbmn2(NP))

       if (nr.ne.NP-1) then
       call spline1d(tmpbc,rho,nr,dBav,s,NP,tmpbmn2)
       nu_prof(:,11) = tmpbc

       do 20 i=1,(2*m_max+1)
          do 20 j=1,n_max
             tmpbmn=real(bmns(i,j,:))
             call spline1d(tmpbc,rho,nr,tmpbmn,s,NP,tmpbmn2)
             bc(i,j,:)=tmpbc
             tmpbmn=aimag(bmns(i,j,:))
             call spline1d(tmpbc,rho,nr,tmpbmn,s,NP,tmpbmn2)
             bs(i,j,:)=tmpbc
 20    continue
       endif

C      new procedure: assuming rho is half-integer mesh of
C      integer mesh s, avoiding using spline,
C      which can do strange things near rational surfaces
       if (nr.eq.NP-1) then
          do i=1,nr
             bc(:,:,i) =  real(bmns(:,:,i)+bmns(:,:,i+1))*0.5
             bs(:,:,i) = aimag(bmns(:,:,i)+bmns(:,:,i+1))*0.5
             nu_prof(i,11) = (dBav(i)+dBav(i+1))*0.5
          enddo
       endif
99     format(1x,500E15.7)

C     note that nu_prof(i,11)=dBav is never used further in the code
C     therefore no need to store it
      IF (keyntv.EQ.1) THEN
      OPEN(31,FILE='DATATORQNTV.OUT')
      REWIND(31)
      DO i=1,nr
         WRITE(31,130) 
     &                 (bc(j,1,i),j=1,(2*m_max+1)),
     &                 (bs(j,1,i),j=1,(2*m_max+1))
      ENDDO
      CLOSE(31)
 130  FORMAT(900(E15.8,1X))
      ENDIF

      IF (keyntv.EQ.2) THEN
      OPEN(31,FILE='DATATORQNTV_A.IN')
      REWIND(31)
      DO i=1,nr
         READ(31,*) 
     &             (bc(j,1,i),j=1,(2*m_max+1)),
     &             (bs(j,1,i),j=1,(2*m_max+1))
      ENDDO
      CLOSE(31)

      OPEN(31,FILE='DATATORQNTV_B.IN')
      REWIND(31)
      DO i=1,nr
         READ(31,*) 
     &             (bc2(j,1,i),j=1,(2*m_max+1)),
     &             (bs2(j,1,i),j=1,(2*m_max+1))
      ENDDO
      CLOSE(31)
      ENDIF

      deallocate(chi)
      deallocate(tmpbc)
      deallocate(tmpr)
      deallocate(r_prof)
      deallocate(R_mars)  
      deallocate(Z_mars)
      deallocate(s)
      deallocate(q)
      deallocate(dpsids)
      deallocate(dBav)
      deallocate(FF)
      deallocate(omega)
      deallocate(dwds)
      deallocate(s_prof)
      deallocate(gg)
      deallocate(dBdc)
      deallocate(dBds)
      deallocate(Jac)
      deallocate(Beq)
      deallocate(chiH)
      deallocate(dphiH)
      deallocate(eps)
      deallocate(tmpB)         
      deallocate(tmps)
      deallocate(tmps2)
      deallocate(dsdpsi)
      deallocate(tmpgg)
      deallocate(bmns)
      deallocate(tmpbmn)
      deallocate(tmpbmn2)
      deallocate(dm,bth,cp,pp)

      end subroutine MARS_out

      subroutine cal_xi_liu(NP,n_chi,chiH,n_phi,dphiH,dBav,
     .                    gg,dBds,dBdc,omega,dwds,m_max,n_max,bmn,
     .                    s,q,B,dm,bth,cp)
C****************************************************************
C           load the plasma displacement Xi
C           and calculate the spectrum of dB
C
C     dB is |B| in perturbed Hamada coordinates
C       dB = dB_L + lam B.grad B,  dB_L = b.deltaB + xi_perp.grad B,
C     with deltaB the computed perturbed field (B1U,B2U,B3U),
C     xi_perp = X1 (e_s)_perp + X2 J B x grad s/B**2 and the
C     Hamada parallel displacement xi_par = lam B from div xi = 0,
C       B.grad lam = -d,  d = -dB_L/B - kappa.xi_perp,
C       lam_k = i d_k (k+nq)/(B^theta ((k+nq)**2+sigma**2)).
C     Under ideal induction d = div xi_perp and dB_L = -B(d+kappa.xi_perp);
C     the field form avoids radial derivatives of xi, which break the
C     near-cancellation in div xi_perp in resistive layers and at the
C     boundary, where X2U(NRP1) is not defined.
C     kappa.xi_perp = xi_perp.grad B/B + mu0 p'(psi) psi'(s) X1/B**2
C     (force balance). The MHD parallel displacement X3 is not used:
C     it is a fluid velocity component, not the Hamada coordinate gauge.
C****************************************************************
!---------------- load displacement
      integer M1,NP1,n,i,j,NP,n_chi,n_phi,k,m_max,n_max
      real gg(NP,n_chi,3),dBds(NP,n_chi),dBdc(NP,n_chi)
      real s(NP),q(NP),B(NP,n_chi),dm(NP,n_chi,7),bth(NP),cp(NP)
      integer,allocatable::M(:)
      real,allocatable::X1R(:,:),X1I(:,:),X2R(:,:),X2I(:,:)
      real,allocatable::X3R(:,:),X3I(:,:)
      real,allocatable::B1R(:,:),B1I(:,:),B2R(:,:),B2I(:,:),
     &                  B3R(:,:),B3I(:,:)
      real,allocatable::F1R(:),F1I(:),F2R(:),F2I(:)
      real,allocatable::fr(:),fi(:),dphi(:),dqds(:)
      complex,allocatable::dB(:),XI1(:),XI2(:),C1(:),C2(:),C3(:)
      complex,allocatable::divx(:),gB(:),lam(:),divk(:),dBL(:)
      real htheta,chiH(NP,n_chi),dphiH(NP,n_chi)
      real chi(n_chi)
      real omega(NP),dwds(NP),dBav(NP)
      real x,sigma,den
      complex bmn(2*m_max+1,n_max,NP)
      integer kmin

       call get_dimensV(M1,NP1,n)
    
       allocate(M(M1))
       allocate(X1R(NP1,M1))
       allocate(X1I(NP1,M1))
       allocate(X2R(NP1,M1))
       allocate(X2I(NP1,M1))
       allocate(X3R(NP1,M1))
       allocate(X3I(NP1,M1))
       allocate(B1R(NP1,M1),B1I(NP1,M1),B2R(NP1,M1),B2I(NP1,M1),
     &          B3R(NP1,M1),B3I(NP1,M1))

       call get_disp(M1,NP1,M,X1R,X1I,X2R,X2I,X3R,X3I)
       call get_field(M1,NP1,B1R,B1I,B2R,B2I,B3R,B3I)
       allocate(dqds(NP))
       call cald_1d(q,NP,s,dqds)

       bmn = 0.
       ! xgrid() supplies periodic midpoint samples with spacing 2*pi/n_chi.
       ! The previous n_chi-1 trapezoid used the wrong spacing, omitted the
       ! periodic closing cell, and leaked a manufactured single Fourier mode.
       htheta = 2.*npi/n_chi
       kmin = -(n_chi/2)

!-------------- transform to real space and calculate dB and its spectrum
       allocate(F1R(n_chi))
       allocate(F1I(n_chi))
       allocate(F2R(n_chi))
       allocate(F2I(n_chi))
       allocate(XI1(n_chi),XI2(n_chi),C1(n_chi),C2(n_chi),C3(n_chi))
       allocate(divx(n_chi),gB(n_chi),lam(n_chi),divk(n_chi),dBL(n_chi))
       allocate(dB(n_chi))
       allocate(fr(M1))
       allocate(fi(M1))
       allocate(dphi(n_chi))

       do 30 i=2,NP1

C      transform X1,X2 and B1,B2,B3 from Fourier space to real
C      space in chi-angle; note that chi is non-unifrom MARS-F chi,
C      corresponding to uniform theta^H in Hamada coordinate system
       chi=chiH(i,:)
       dphi=dphiH(i,:)

       fr=X1R(i,:)
       fi=X1I(i,:) 
       call trans_1d(fr,fi,M1,M,n,n_chi,chi,dphi,F1R,F1I)
       XI1 = F1R+(0.,1.)*F1I

       fr=X2R(i,:)
       fi=X2I(i,:) 
       call trans_1d(fr,fi,M1,M,n,n_chi,chi,dphi,F2R,F2I)
       XI2 = F2R+(0.,1.)*F2I

       call trans_1d(B1R(i,:),B1I(i,:),M1,M,n,n_chi,chi,dphi,F1R,F1I)
       C1 = F1R+(0.,1.)*F1I
       call trans_1d(B2R(i,:),B2I(i,:),M1,M,n,n_chi,chi,dphi,F1R,F1I)
       C2 = F1R+(0.,1.)*F1I
       call trans_1d(B3R(i,:),B3I(i,:),M1,M,n,n_chi,chi,dphi,F1R,F1I)
       C3 = F1R+(0.,1.)*F1I

C      xi_perp.grad B
       gB = XI1*dBds(i,:) + (XI1*gg(i,:,1) + XI2*gg(i,:,2))*dBdc(i,:)

C      dB_L = b.deltaB + xi_perp.grad B with J deltaB^(s,chi,phi) = C1,C2,C3,
C      b.deltaB = [(psi'/J)(g_schi C1 + g_chichi C2) + F C3]/(J B)
       dBL = B(i,:)*(-gg(i,:,1)*C1/gg(i,:,3) - dm(i,:,6)*C2/dm(i,:,7)
     &     + gg(i,:,2)*C3)/dm(i,:,7) + gB

C      d = -dB_L/B - kappa.xi_perp
       divx = -dBL/B(i,:) - gB/B(i,:) - cp(i)*XI1/B(i,:)**2

C      Hamada gauge B.grad lam = -d on exp(i(k theta + n zeta)),
C      B.grad = B^theta (d_theta + q d_zeta); Park regularisation sigma
C      of one radial cell, sigma = |n q'| ds
       if (i.lt.NP1) then
          sigma = abs(n*dqds(i))*(s(i+1)-s(i-1))*0.5
       else
          sigma = abs(n*dqds(i))*(s(i)-s(i-1))
       endif
       do k=1,n_chi
          divk(k) = 0.
          do j=1,n_chi
             divk(k) = divk(k) +
     &       divx(j)*exp(-(0.,1.)*(kmin+k-1)*(j-.5)*htheta)
          enddo
          x = kmin+k-1+n*q(i)
          den = x**2+sigma**2
          if (den.gt.0.) then
             divk(k) = (0.,1.)*divk(k)/n_chi*x/(bth(i)*den)
          else
             divk(k) = 0.
          endif
       enddo
       do j=1,n_chi
          lam(j) = 0.
          do k=1,n_chi
             lam(j) = lam(j) +
     &       divk(k)*exp((0.,1.)*(kmin+k-1)*(j-.5)*htheta)
          enddo
       enddo

C      dB = dB_L + lam B.grad B
       dB = dBL + lam*gg(i,:,3)*dBdc(i,:)

C      compute surface averaged value of abs(dB/B)
C      to be used to estimate various collisionality regimes
       dBav(i) = sum(abs(dB))/DFLOAT(n_chi)

C      decompose dB in Fourier space in Hamada coordinate system (V,theta^H,phi^H)
C      note that transform to Hamada coordinates does not change the toroidal 
C      mode number n
       do k=-m_max,m_max
          do j=1,n_chi
             bmn(k+m_max+1,1,i) = bmn(k+m_max+1,1,i) +
     &       dB(j)*exp(-(0.,1.)*k*(j-.5)*htheta)
          enddo
       enddo

 30    continue
       bmn = bmn*htheta/2./npi
       bmn(:,:,1)=bmn(:,:,2)

       deallocate(M, X1R, X1I, X2R, X2I, X3R, X3I,
     &            B1R, B1I, B2R, B2I, B3R, B3I,
     &            F1R, F1I, F2R, F2I, XI1, XI2, C1, C2, C3,
     &            divx, gB, lam, divk, dBL, dB, fr, fi, dphi, dqds)

      end subroutine cal_xi_liu


      subroutine cal_xi_sun(NP,n_chi,chiH,n_phi,dphiH,
     .                    gg,dBds,dBdc,omega,dwds,m_max,n_max,bmn)
C****************************************************************
C           load the plasma displacement Xi
C           and calculate the spectrum of dB
C****************************************************************
!---------------- load displacement
      character*(30) f3
      integer M1,NP1,n,i,j,NP,n_chi,n_phi,k,m_max,n_max
      real gg(NP,n_chi,3),dBds(NP,n_chi),dBdc(NP,n_chi)
      integer,allocatable::M(:)
      real,allocatable::X1R(:,:),X1I(:,:),X2R(:,:),X2I(:,:)
      real,allocatable::F3R(:,:),F3I(:,:),V3R(:,:),V3I(:,:)
      real,allocatable::F1R(:,:),F1I(:,:),F2R(:,:),F2I(:,:)
      real,allocatable::fr(:),fi(:),X3(:),phi(:,:),phi0(:)
      real,allocatable::dB(:,:),XI_chi(:)
      real dum,chiH(NP,n_chi),dphiH(NP,n_chi)
      real chi(n_chi)
      real min_t,max_t,omega(NP),dwds(NP)
      complex bmn(2*m_max+1,n_max,NP),tmpbmn(2*m_max+1,n_max)

      call get_dimensV(M1,NP1,n)
    
       allocate(M(M1))
       allocate(X1R(NP1,M1))
       allocate(X1I(NP1,M1))
       allocate(X2R(NP1,M1))
       allocate(X2I(NP1,M1))
       allocate(V3R(NP1,M1))
       allocate(V3I(NP1,M1))

       call get_disp(M1,NP1,M,X1R,X1I,X2R,X2I,V3R,V3I)

!-------------- transform to real space and calculate dB and its spectrum
       allocate(F1R(n_chi,n_phi))
       allocate(F1I(n_chi,n_phi))
       allocate(F2R(n_chi,n_phi))
       allocate(F2I(n_chi,n_phi))
       allocate(F3R(n_chi,n_phi))
       allocate(F3I(n_chi,n_phi))
       allocate(X3(n_chi))
       allocate(XI_chi(n_chi))
       allocate(dB(n_chi,n_phi))
       allocate(fr(M1))
       allocate(fi(M1))
       allocate(phi(n_chi,n_phi))
       allocate(phi0(n_phi))

        min_t=0.0
        max_t=2.0*npi
           call xgrid(n_phi,min_t,max_t,phi0)

       do 30 i=2,NP1
 !        print*,'i=',i,'/',NP1
           chi=chiH(i,:)
          do 20 j=1,n_phi
            phi(:,j)=phi0(j)-dphiH(i,:)
 20       continue 
        !---------transform  
        fr=X1R(i,:)
        fi=X1I(i,:) 
!        print*,M,n
         call trans_2d(fr,fi,M1,M,n,n_chi,n_phi,chi,phi,F1R,F1I)
!         print*,'fr',i,fr
!         print*,'f1r',i,F1R(:,1)
        fr=X2R(i,:)
        fi=X2I(i,:) 
         call trans_2d(fr,fi,M1,M,n,n_chi,n_phi,chi,phi,F2R,F2I)
        fr=V3R(i,:)
        fi=V3I(i,:) 
         call trans_2d(fr,fi,M1,M,n,n_chi,n_phi,chi,phi,F3R,F3I)

         !---------calculate dB 
         do 40 k=1,n_phi
          X3=(F3I(:,k)+gg(i,:,2)*dwds(i)*F1I(:,k))/n/(omega(i)+1.0e-6)
                                                     ! avoid singularity for omega=0
          XI_chi=F1R(:,k)*gg(i,:,1)+F2R(:,k)*gg(i,:,2)+
     .                   X3*gg(i,:,3)
          dB(:,k)=F1R(:,k)*dBds(i,:)+XI_chi*dBdc(i,:)
 40      continue
        !--------- calculate dB spectrum
         call myspec(n_chi,n_phi,dB,m_max,n_max,tmpbmn)
         bmn(:,:,i)=tmpbmn

C        kliuprint
         if (i.eq.186.and.1.eq.0) then
            print*,'kliuprint...'
            print*,X1R(i,31),X1I(i,31)
            print*,X2R(i,31),X2I(i,31)
            print*,V3R(i,31),V3I(i,31)
            print*,F1R(1,1),F2R(1,1),omega(i)  
            print*,dB(1,1),bmn(31,1,i)
         endif

 30    continue
        bmn(:,:,1)=bmn(:,:,2)

      close(16)
99     format(1x,500E15.7)

      deallocate(M, X1R, X1I, X2R, X2I, V3R, V3I, 
     &              F1R, F1I, F2R, F2I, F3R, F3I,
     &           X3, XI_chi, dB, fr, fi, phi, phi0)

      end subroutine cal_xi_sun


      subroutine plasma_prof(KROTWE,NP,n_prof,s_prof,s,ndpsids,FF,
     &                       nomega,dwds)
C****************************************************************
C           load the plasma profile
C****************************************************************
      use GLOBALM
      use RCOMDM

      character*(30) f2
      integer NP,i,n_prof,KROTWE
      real dum,nne0,nTi0,nTe0,nw_phi0
      real s(NP),ndpsids(NP),FF(NP),s_prof(NP,n_prof),nomega(NP),
     &     dwds(NP)
      real,allocatable::w_phi(:),nq(:),Ne(:),Te(:),Ti(:),dNe(:),
     &     dTe(:),dTi(:),nui(:),nue(:)

      allocate(w_phi(NP),nq(NP),Ne(NP),Te(NP),Ti(NP),dNe(NP),
     &         dTe(NP),dTi(NP),nui(NP),nue(NP))

      KROTWE = NPROFWE
!---------------- open files
      do i=1,NP
         s(i) = CS(i)
         nq(i) = Q(i)
         Ne(i) = RHO(i)
         nomega(i) = (ROT(i)+TROTI(i))
         Ti(i) = TEMPI(i)
         Te(i) = TEMPE(i)
         ndpsids(i) = DPSIDS(i)
         FF(i) = T(i)
         s_prof(i,12) = ROTWEI(i)/ZTAUA0 
         nui(i) = GNUI(i)
         nue(i) = GNUE(i)
      end do
      nne0 = ZNE0*1.e-20
      nTi0 = ZTI0*1.e-03
      nTe0 = ZTE0*1.e-03
      nw_phi0 = 1./ZTAUA0
      Ne=nne0*Ne
      Te=nTe0*Te
      Ti=nTi0*Ti
      w_phi=nw_phi0*nomega
      nui = nui * nw_phi0
      nue = nue * nw_phi0
      select case(numodel)
      case (0)
         nui = -1.
         nue = -1.
      end select
C     print*,'ne0,Ti0,Te0,w_phi0',Ne(1),Ti(1),Te(1),w_phi(1)

      call cald_1d(Te,NP,s,dTe)
      call cald_1d(Ti,NP,s,dTi)
      call cald_1d(Ne,NP,s,dNe)
      call cald_1d(nomega,NP,s,dwds)
      s_prof(:,1)=s
      s_prof(:,2)=nq
      s_prof(:,3)=Ne
      s_prof(:,4)=Te
      s_prof(:,5)=Ti
      s_prof(:,6)=w_phi
      s_prof(:,7)=dNe
      s_prof(:,8)=dTe
      s_prof(:,9)=dTi
      s_prof(:,13)=nui
      s_prof(:,14)=nue
      deallocate(w_phi,nq,Ne,Te,Ti,dNe,dTe,dTi,nui,nue)

      end subroutine plasma_prof


      subroutine cal_metric(NP,n_chi,s,chi,R,Z,F,dpsids,q,
     .                       Jac,gg,B,dBds,dBdc)
C************************************************************
C Input
C    R(s,chi),Z(s,chi) is postion of the flux surface
C           for equal arc coordinates (s,chi).
C    dpsids=dpsi_p/ds
C    q, the safty factor
C    F=R*BT
C
C Output:  
C       dB=dBds*Xi^s+dBdchi*Xi^chi 
c     Xi^s=Xi^1
c   Xi^chi=gg(:,:,1)*chi^1+gg(:,:,2)*chi^2+gg(:,:,3)*chi^3
C          gg(:,:,1)=-(dpsids/JB)**2*g_schi
C          gg(:,:,2)=F/B**2
C          gg(:,:,3)=dpsids/J
c
C************************************************************
      integer NP,n_chi,i
      real R(NP,n_chi),Z(NP,n_chi),s(NP),chi(n_chi),F(NP),q(NP)
      real dpsids(NP),dchi,Jac(NP,n_chi)
      real B(NP,n_chi),dBds(NP,n_chi),dBdc(NP,n_chi),gg(NP,n_chi,3)

      real,allocatable::dRds(:,:),dRdchi(:,:),dZds(:,:),
     &                  dZdchi(:,:),g_sc(:,:),g_cc(:,:),
     &                  PJ(:,:),FB(:,:),BT(:,:),BP(:,:)

      allocate(dRds(NP,n_chi),dRdchi(NP,n_chi),dZds(NP,n_chi),
     &         dZdchi(NP,n_chi),g_sc(NP,n_chi),g_cc(NP,n_chi),
     &         PJ(NP,n_chi),FB(NP,n_chi),BT(NP,n_chi),BP(NP,n_chi))

!-------  metrics calculation
      call cald_2d(R,NP,n_chi,s,chi,dRds,dRdchi)
      call cald_2d(Z,NP,n_chi,s,chi,dZds,dZdchi)
       Jac=(dRds*dZdchi-dRdchi*dZds)*R
       g_cc=dRdchi**2+dZdchi**2
       g_sc=dRds*dRdchi+dZds*dZdchi
       do i=1,NP
          PJ(i,:)=dpsids(i)/Jac(i,:)
          BT(i,:)=F(i)/R(i,:)
       end do
          PJ(1,:)=0.0 ! singularity at s=0
        BP=PJ*sqrt(g_cc) 
          BP(1,:)=0.0
        B=sqrt(BT**2+BP**2)
        call cald_2d(B,NP,n_chi,s,chi,dBds,dBdc) ! dB/ds,dB/dchi
       do i=1,NP
          FB(i,:)=F(i)/(B(i,:)**2)
       end do
       gg(:,:,1)=-g_sc*PJ**2/B**2
       gg(:,:,2)=FB             ! F/B**2
       gg(:,:,3)=PJ
       
       deallocate(dRds,dRdchi,dZds,dZdchi,g_sc,g_cc,PJ,FB,BT,BP)

      end subroutine cal_metric

      subroutine cal_divmetric(NP,n_chi,s,chi,R,Z,F,dpsids,Jac,B,dm)
C************************************************************
C Input
C    R(s,chi),Z(s,chi) on the equal arc grid, F=R*BT, dpsids=dpsi_p/ds,
C    Jac, the Jacobian of (s,chi,phi), and B=|B| from cal_metric
C
C Output: coefficients of
C    J div xi_perp = dm1*X1 + dm7*dX1/ds + dm2*dX1/dchi
C                  + dm3*X2 + dm4*dX2/dchi + dX1,X2/dphi terms
C    for xi_perp = X1 (e_s - B_s B/B**2) + X2 J B x grad s/B**2:
C          J xi^s   = J X1                       (dm7=J)
C          J xi^chi = dm2*X1 + dm4*X2,  dm2=-g_schi psi'**2/(J B**2),
C                                       dm4=J F/B**2
C          J xi^phi = dm5*X1 + dm6*X2,  dm5=-g_schi psi' F/(R B)**2,
C                                       dm6=-g_chichi psi'/B**2
C          dm1 = dJ/ds + d(dm2)/dchi,  dm3 = d(dm4)/dchi
C    chi derivatives are periodic centred differences
C************************************************************
      integer NP,n_chi,i,j,jp,jm
      real R(NP,n_chi),Z(NP,n_chi),s(NP),chi(n_chi),F(NP)
      real dpsids(NP),Jac(NP,n_chi),B(NP,n_chi),dm(NP,n_chi,7)
      real hchi

      real,allocatable::dRds(:,:),dRdchi(:,:),dZds(:,:),
     &                  dZdchi(:,:),dJds(:,:),dJdchi(:,:),
     &                  g_sc(:,:),g_cc(:,:)

      allocate(dRds(NP,n_chi),dRdchi(NP,n_chi),dZds(NP,n_chi),
     &         dZdchi(NP,n_chi),dJds(NP,n_chi),dJdchi(NP,n_chi),
     &         g_sc(NP,n_chi),g_cc(NP,n_chi))

      call cald_2d(R,NP,n_chi,s,chi,dRds,dRdchi)
      call cald_2d(Z,NP,n_chi,s,chi,dZds,dZdchi)
      call cald_2d(Jac,NP,n_chi,s,chi,dJds,dJdchi)
       g_cc=dRdchi**2+dZdchi**2
       g_sc=dRds*dRdchi+dZds*dZdchi
       dm=0.0
       do i=2,NP
          dm(i,:,2)=-g_sc(i,:)*dpsids(i)**2/(Jac(i,:)*B(i,:)**2)
          dm(i,:,4)=Jac(i,:)*F(i)/B(i,:)**2
          dm(i,:,5)=-g_sc(i,:)*dpsids(i)*F(i)/(R(i,:)*B(i,:))**2
          dm(i,:,6)=-g_cc(i,:)*dpsids(i)/B(i,:)**2
          dm(i,:,7)=Jac(i,:)
       end do
       hchi=chi(2)-chi(1)
       do j=1,n_chi
          jp=mod(j,n_chi)+1
          jm=mod(j+n_chi-2,n_chi)+1
          dm(:,j,1)=dJds(:,j)+(dm(:,jp,2)-dm(:,jm,2))/(2.0*hchi)
          dm(:,j,3)=(dm(:,jp,4)-dm(:,jm,4))/(2.0*hchi)
       end do
       dm(1,:,:)=0.0 ! singularity at s=0

      deallocate(dRds,dRdchi,dZds,dZdchi,dJds,dJdchi,g_sc,g_cc)

      end subroutine cal_divmetric

      subroutine cal_eps(n_chi,chi,B,eps)
C************************************************************
C     calculate epsillon
C************************************************************
      integer n_chi,i
      real chi(n_chi),B(n_chi),eps
      eps=-sum(B*cos(chi))/sum(cos(chi)**2)
      end subroutine cal_eps

      subroutine coord_trans(NP,n_chi,chi,R,q,Jac,chiH,dphiH)
C************************************************************
C Input
C    R(s,chi),Z(s,chi) is postion of the flux surface
C           for equal arc coordinates (s,chi).
C    q, the safty factor
C    Jac is the Jacobin of (s,chi,phi)
C
C Output:  
C     chiH=chi(theta^H)
C     phiH=phi(theta^H,zeta^H)
C       zeta^H=phiH+dphiH
C    dphiH=dphi(theta^H)
c 
C************************************************************
      integer NP,n_chi,i
      real R(NP,n_chi),chi(n_chi),q(NP)
      real dchi,Jac(NP,n_chi),chiH(NP,n_chi),dphiH(NP,n_chi)
      real,dimension(:),allocatable::eps,VS,theHc,theH,chi_H,tmp,
     &                  tmp_chiH,fac,phic,dphi,tmpdphiH,tmpchi

      allocate(eps(NP),VS(NP),theHc(n_chi+1),theH(n_chi+1),
     &         chi_H(n_chi+1),tmp(n_chi+1),tmp_chiH(n_chi),
     &         fac(n_chi+1),phic(n_chi+1),dphi(n_chi+1),
     &         tmpdphiH(n_chi),tmpchi(n_chi+2))

!----------- mid points for theta^H calculation
         chi_H(1)=0.0
         chi_H(n_chi+1)=2.0*npi
       do i=2,n_chi
         chi_H(i)=(chi(i)+chi(i-1))/2.0
       end do
       dchi=chi(2)-chi(1)
!----------- coordinate trans
       do i=2,NP
        ! -----calculate chi(theta^H)
          VS(i)=sum(Jac(i,:))*dchi
          theHc(1)=0.0
          theHc(2:(n_chi+1))=2*npi*Jac(i,:)/VS(i)
          theH=cumsum(n_chi+1,theHc)*dchi   ! calculate theta^H
        !    print*,'theH',i,theH(1),theH(n_chi+1)
          call spline1d(tmp_chiH,chi,n_chi,chi_H,theH,n_chi+1,tmp)
          chiH(i,:)=tmp_chiH  ! tranform to correspoinding chi for equally spaced theta^H
        ! -----calculate dphi(theta^H)  
          fac(2:(n_chi+1))=1.0-1.0/R(i,:)**2*
     .                      VS(i)/sum(Jac(i,:)/R(i,:)**2*dchi)
          fac(1)=0.0
          fac=q(i)*fac*theHc
           dphi=cumsum(n_chi+1,fac)*dchi    !  zeta=phi+dphi
           call spline1d(tmpdphiH,chi,n_chi,dphi,theH,n_chi+1,tmp)
           dphiH(i,:)=tmpdphiH ! tranform to correspoinding chi for equally spaced theta^H
       end do
       chiH(1,:)=chi
       dphiH(1,:)=dphiH(2,:)    ! singularity at s=0

       deallocate(eps,VS,theHc,theH,chi_H,tmp,tmp_chiH,fac,phic,
     &            dphi,tmpdphiH,tmpchi)

      end subroutine coord_trans


      subroutine coord_map(n_chi,chi,NP,chiH,B)
C*********************************************************************
C        eq mapping
C         B(NP,chi) to B(NP,theta^H)
C*********************************************************************
      integer n_chi,NP,i
      real B(NP,n_chi),chi(n_chi),chiH(NP,n_chi)
      real,dimension(:),allocatable::tmpB,tmpchi,tmp_chiH,tmp2,NB

      allocate(tmpB(n_chi+2),tmpchi(n_chi+2),tmp_chiH(n_chi),
     &         tmp2(n_chi+2),NB(n_chi))

!----------- B(chi) for interpolating, to include the point in -h<chi<h
         tmpchi(1)=chi(n_chi)-2.0*npi
         tmpchi(2:(n_chi+1))=chi
         tmpchi(n_chi+2)=chi(1)+2.0*npi

        do i=1,NP
         tmpB(1)=B(i,n_chi)
         tmpB(2:(n_chi+1))=B(i,:)
         tmpB(n_chi+2)=B(i,1)
         tmp_chiH=chiH(i,:)
         call spline1d(NB,tmp_chiH,n_chi,tmpB,tmpchi,n_chi+2,tmp2)
         B(i,:)=NB
        end do

        deallocate(tmpB,tmpchi,tmp_chiH,tmp2,NB)

      end subroutine coord_map

      subroutine cal_nu(KROTWE,nr,r,drdpsi,n_prof,r_prof,mui,R0,B0,
     &                  Z,Zi,knc,nu_prof)
C*********************************************************************
C    calculate the collisionality and some frequencies
C        for NTV torque calcualtion.
C
C    Input:
C        if Z=-1 (Z<0), then it always for electrons.
C
C        mi=mui*mp
C        R0,major radius at magnetic axis
C        B0,magnetic field strength at the magnetic axis
C        knc=1.17 is the neoclassical coefficent for poloidal flow 
C                  without momentum input.
C        
C
C  radial profiles:
C   drdpsi=dr/dpsi_p
C            drdpsi=q/(rB_0) for r=r_t
C   r_prof(nr,n_prof)=
C           (rho,q,Ne,Te,Ti,w_phi,dNe/dr,dTe/dr,dTi/dr,eps,dB,we,nui,nue)
C             1  2  3  4  5  6      7      8       9    10 11 12  13  14
C
C   output:
C
C   nu_prof(nr,n_prof)=
C           (rho,q,eps,rho_i,nu,wti,wB0,qwe,qws,qwsT,dB,we,nui,nue)
C             1  2  3    4    5   6   7   8  9  10   11 12  13  14
C
C      rho is the radius, q is the safty factor, eps=epsillon~ r/R
C      rho_i=n_i*m_i, nu is the collisionality, wti is the trasit frequency
C      wB0=omega_B0 is the grad B drift of the deeply trapped pariticles with x=1,
C      qwe=q*omega_E is the E \times B drift
C      qws=q*omega_*p, is the diamagnetic frequency
C      qwsT=q*omega_*T,
C
C         Ne=Ne(e+20)
C         T=T(keV)
C         W=W(rad/s),  omega
C
C         by Y. Sun
C*********************************************************************
      implicit none
      integer nr,Z,n_prof
      integer Zi,KROTWE
      real knc
      real r(nr),drdpsi(nr),r_prof(nr,n_prof),nu_prof(nr,n_prof)
      real R0,B0,mui,wg,me

      real,dimension(:),allocatable::eps,q,Ne,Te,Ti,w_phi,Pe,Pi,
     &                  dNe,dTe,dTi,dPe,dPi,qws,qwsT,qwsi,qwsTi,
     &                  qwthe,qwe,nu,wti,wt,wte,wB0,rho_i

      allocate(eps(nr),q(nr),Ne(nr),Te(nr),Ti(nr),w_phi(nr),Pe(nr),
     &         Pi(nr),dNe(nr),dTe(nr),dTi(nr),dPe(nr),dPi(nr),
     &         qws(nr),qwsT(nr),qwsi(nr),qwsTi(nr),qwthe(nr),
     &         qwe(nr),nu(nr),wti(nr),wt(nr),wte(nr),wB0(nr),
     &         rho_i(nr))

!----------
!         Zi=1
!         knc=1.17
!---------- get profiles
          r=r_prof(:,1)
          q=r_prof(:,2)
         Ne=r_prof(:,3)
         Te=r_prof(:,4)
         Ti=r_prof(:,5)
      w_phi=r_prof(:,6)
        dNe=r_prof(:,7)
        dTe=r_prof(:,8)
        dTi=r_prof(:,9)       
        eps=r_prof(:,10)
!---------- diamagnetic frequencies
      Pe=Ne*Te        ! Pe/e0, e0=1.6e-19, Te(keV)
      Pi=Ne*Ti
        dPe=Ne*dTe+Te*dNe
        dPi=Ne*dTi+Ti*dNe
      rho_i=1.67e-7*mui*Ne

!      call cald_1d(Te,nr,r,dTe)
!      call cald_1d(Ti,nr,r,dTi)
!      call cald_1d(Pe,nr,r,dPe)
!      call cald_1d(Pi,nr,r,dPi)

      dTe=1.0e+3*dTe*drdpsi
      dTi=1.0e+3*dTi*drdpsi
      dPe=1.0e+3*dPe*drdpsi        
      dPi=1.0e+3*dPi*drdpsi

        qwsi=dPi/Ne/Zi
         qwsTi=dTi/Zi
      if (Z.GT.0)then        ! diamagnetic rotation
         qws=qwsi
         qwsT=qwsTi
      else
          qws=dPe/Ne/Z
         qwsT=dTe/Z
      end if
      
      qwthe=knc*qwsTi        ! poloidal rotation
      if (KROTWE.eq.4) then
         qwe=q*r_prof(:,12)
      else
         qwe=w_phi-qwthe+qwsi   ! E * B rotation
      endif
!---------- collisionalities
    !  Ne=Ne/1.0e20
    !  Te=Te/1.0e3
    !  Ti=Ti/1.0e3
    !  print*,R0
      wti=4.4e5*sqrt(Ti)/sqrt(mui)/R0/q
      wte=1.89e7*sqrt(Te)/sqrt(mui)/R0/q

      if (Z.GT.0)then     ! collisionality and grad B drift
         if (r_prof(1,13).lt.0) then
        !sqrt(2)*pi*1e+20*e^4/(4*pi*epsilon0)^2/sqrt(mp)/(e*1e+3)^1.5=2.85e+2
            nu=2.85e2*Ne*18.7/sqrt(mui)/Ti**1.5
         else
            nu=r_prof(:,13)
         endif
          wg=9.6e+7*Z/mui*B0
        !  print*,'wg',wg
         wB0=q**3*wti**2/2.0/eps/wg
        ! print*,'wB0',wB0
      else
         if (r_prof(1,14).lt.0) then
            nu=1.22e4*Ne*18.7/Te**1.5
         else
            nu=r_prof(:,14)
         endif
          wg=1.78e+11*Z*B0
         wB0=q**3*wte**2/2.0/eps/wg
      end if    
C   nu_prof(nr,10)=
C           (rho,q,eps,rho_i,nu,wti,wB0,qwe,qws,qwsT)
C             1  2  3    4    5   6   7   8  9  10
      nu_prof(:,1)=r
      nu_prof(:,2)=q
      nu_prof(:,3)=eps
      nu_prof(:,4)=rho_i
      nu_prof(:,5)=nu
      nu_prof(:,6)=wti
      nu_prof(:,7)=wB0
      nu_prof(:,8)=qwe
      nu_prof(:,9)=qws
      nu_prof(:,10)=qwsT
      nu_prof(:,12)=r_prof(:,12)
      nu_prof(:,13)=r_prof(:,13)
      nu_prof(:,14)=r_prof(:,14)

      deallocate(eps,q,Ne,Te,Ti,w_phi,Pe,Pi, dNe,dTe,dTi,
     &           dPe,dPi,qws,qwsT,qwsi,qwsTi,qwthe,qwe,nu,
     &           wti,wt,wte,wB0,rho_i)

      end subroutine cal_nu

      subroutine myspec(nt,np,db,m_max,n_max,bmn)
C************************************************************
C        calculate the spectrum of bmn(-m_max:m_max,1:n_max) 
C        from  db(the(nt),phi(np))
C          by Y. Sun    
C************************************************************
      use NTV_lib
      integer nt,np,m_max,n_max,i
      real    db(nt,np)
      complex bmn(2*m_max+1,n_max)
      integer,parameter::ndim=2
      integer nn(2)
      
      complex,allocatable::fk(:,:)

      allocate(fk(nt,np))

      if ((m_max.GT.(nt/2-1)).or.(n_max.GT.(np/2-1))) then
         print*,'ERROR: m_max or n_max is set too large.'
         print*,'m_max=',m_max,',   n_max=',n_max
          stop
      end if
      nn(1)=nt
      nn(2)=np
      call fftn(2,nn,db,fk)
      fk=2.0*fk/dble(nt*np)  ! the factor 2 comes from the m/n part

      do i=1,n_max
         bmn((m_max+1):(2*m_max+1),i)=fk(1:m_max,np-i+1)  ! bmn(0:m_max,-n)
         bmn(1:m_max,i)=fk((nt-m_max+1):nt,np-i+1)        ! bmn(-m_max:-1,-n)
      end do

      deallocate(fk)

      end subroutine myspec


      subroutine get_dimens(NP,R0,B0)
      use GLOBALM
      use DIMENSIM
      use TORQUEM
      integer NP
      real*8  R0,B0

      NP = NRP1
      R0 = R0EXP
      B0 = B0EXP

      end subroutine get_dimens


      subroutine get_RZ(NP,n_chi,chi,R_mars,Z_mars)
      use DIMENSIM
      use RCOMDM

      integer NP,n_chi,i,j
      real chi(n_chi)
      real R_mars(NP,n_chi),Z_mars(NP,n_chi)

      real hchi
      real,allocatable::rchi(:),tmpc(:),tmpc2(:),tmpr(:)

      allocate(rchi(NCHI+1),tmpc(NCHI+1),tmpc2(NCHI+1))
      allocate(tmpr(n_chi))

      hchi = 2*npi/NCHI
      rchi(1) = 0.0
      do j=1,NCHI
         rchi(j+1) = rchi(j) + hchi
      enddo

      do i=1,NP
         do j=1,NCHI
            tmpc(j) = REQ(i,j)
         enddo
         tmpc(NCHI+1) = tmpc(1)
         call spline1d(tmpr,chi,n_chi,tmpc,rchi,NCHI+1,tmpc2)
         do j=1,n_chi
            R_mars(i,j) = tmpr(j)
         enddo
 
         do j=1,NCHI
            tmpc(j) = ZEQ(i,j)
         enddo
         tmpc(NCHI+1) = tmpc(1)
         call spline1d(tmpr,chi,n_chi,tmpc,rchi,NCHI+1,tmpc2)
         do j=1,n_chi
            Z_mars(i,j) = tmpr(j)
         enddo
      enddo

      deallocate(rchi,tmpc,tmpc2,tmpr)

      end subroutine get_RZ


      subroutine get_dimensV(MM1,NP1,n)
      use DIMENSIM
      use GLOBALM

      integer MM1,NP1,n
   
      MM1  = MSMAX
      NP1 = NRP1
      n   = RNTOR

      end subroutine get_dimensV

      subroutine get_ppeq(NP,pp)
C     mu0 dp/dpsi on the integer mesh (MARS normalisation)
      use RCOMDM

      integer NP
      real    pp(NP)

      pp = PPEQ(1:NP)

      end subroutine get_ppeq

      subroutine get_disp(M10,NP1,M,X1R,X1I,X2R,X2I,X3R,X3I)
      use DIMENSIM
      use GLOBALM

      integer M10,NP1,M(M10)
      real    X1R(NP1,M10),X1I(NP1,M10),X2R(NP1,M10),X2I(NP1,M10),
     &        X3R(NP1,M10),X3I(NP1,M10)

      real,dimension(:,:),allocatable::tmp

      integer i,j,jx,jy

      do j=1,M10
         M(j) = RM(j,2)
      enddo
   
      do i=1,NP1
         do j=1,M10
            jx = (j-1)*NXCOMP
            jy = (j-1)*NYCOMP
            X1R(i,j) = real(X1U(i,j))
            X1I(i,j) = imag(X1U(i,j))
            X2R(i,j) = real(X2U(i,j))
            X2I(i,j) = imag(X2U(i,j))
            X3R(i,j) = real(X3U(i,j))
            X3I(i,j) = imag(X3U(i,j))
         enddo
      enddo

C     recompute X2 and X3 from half-integer to integer mesh
C     using linear interpolation
      allocate(tmp(NP1,M10))

      tmp = X2R
      do i=1,NP1-1
      X2R(i+1,:)=(CSH(i+1)*tmp(i,:)+CSH(i)*tmp(i+1,:))/(CSH(i)+CSH(i+1))
      enddo

      tmp = X2I
      do i=1,NP1-1
      X2I(i+1,:)=(CSH(i+1)*tmp(i,:)+CSH(i)*tmp(i+1,:))/(CSH(i)+CSH(i+1))
      enddo

      tmp = X3R
      do i=1,NP1-1
      X3R(i+1,:)=(CSH(i+1)*tmp(i,:)+CSH(i)*tmp(i+1,:))/(CSH(i)+CSH(i+1))
      enddo

      tmp = X3I
      do i=1,NP1-1
      X3I(i+1,:)=(CSH(i+1)*tmp(i,:)+CSH(i)*tmp(i+1,:))/(CSH(i)+CSH(i+1))
      enddo

      deallocate(tmp)
      end subroutine get_disp


      subroutine get_field(M10,NP1,B1R,B1I,B2R,B2I,B3R,B3I)
C     perturbed field J deltaB^(s,chi,phi) on the integer mesh; B2U and
C     B3U live on the half mesh and are centred as X2U, X3U in get_disp
      use DIMENSIM
      use GLOBALM

      integer M10,NP1
      real    B1R(NP1,M10),B1I(NP1,M10),B2R(NP1,M10),B2I(NP1,M10),
     &        B3R(NP1,M10),B3I(NP1,M10)

      complex*16,dimension(:,:),allocatable::tmp
      integer i

      B1R = real(B1U(1:NP1,1:M10))
      B1I = imag(B1U(1:NP1,1:M10))
      allocate(tmp(NP1,M10))
      tmp = B2U(1:NP1,1:M10)
      do i=1,NP1-1
      tmp(i+1,:)=(CSH(i+1)*B2U(i,1:M10)+CSH(i)*B2U(i+1,1:M10))
     &           /(CSH(i)+CSH(i+1))
      enddo
      B2R = real(tmp)
      B2I = imag(tmp)
      tmp = B3U(1:NP1,1:M10)
      do i=1,NP1-1
      tmp(i+1,:)=(CSH(i+1)*B3U(i,1:M10)+CSH(i)*B3U(i+1,1:M10))
     &           /(CSH(i)+CSH(i+1))
      enddo
      B3R = real(tmp)
      B3I = imag(tmp)
      deallocate(tmp)

      end subroutine get_field

      subroutine get_alpha(NRHO,NTHETA,X3,alpha1,alpha2)
      use GLOBALM
      integer NRHO,NTHETA
      complex alpha1,X3(NRHO,NTHETA)
      real    alpha2

      alpha1 = CALPHA1
      alpha2 = CALPHA2
    
      end subroutine get_alpha
  
      subroutine set_alpha(NRHO,NTHETA,X3)
      use GLOBALM
      integer NRHO,NTHETA
      complex X3(NRHO,NTHETA)

      end subroutine set_alpha
  
      end
      
