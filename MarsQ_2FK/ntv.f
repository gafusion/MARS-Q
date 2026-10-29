      Module NTV
C***********************************************************************
C      calculate the NTV torque density
C      according to Shaing's Formula.
C
C        based on the version in MATLAB
C
C  important subroutines:
C        gntv1d(), calculate the radial profile 
C                      of the NTV torque density
C
C         gntv(), calculate the NTV torque density
C    Ik_shaing(), calculate the pitch angle integral
C                 according to Shaing's formula
C       dbspec(), calculate related bounce averaged perturbation
C
C
C        by Y. Sun
C        Mar. 01, 2011
C***********************************************************************
      use NTV_lib
      implicit none

      contains

      subroutine gntv1d
C**********************************************************************
C      calculate the NTV torque density profile
C
C     need to call:
C           xgrid, generate the grids for energy space x.
C           f_nud, calculate the coeffcients in the deflection frequency.
C           gntv, calculate the NTV torque density.
C
C   Input:
C
C   nx,nk,nr are the grid numbers in the energy, kappa and real spaces.
C   R_axis is the major radius of the magnetic axis.
C   Z is the charge number of the particles.
C   nrom_Z=Z_j/Z_i.
C
C    radial profiles:
C   nu_prof(nr,10)=
C           (rho,q,eps,rho_i,nu,wti,wB0,qwe,qws,qwsT)
C             1  2  3    4    5   6   7   8  9  10
C      rho is the radius, q is the safty factor, eps=epsillon~ r/R
C      rho_i=n_i*m_i, nu is the collisionality, wti is the trasit frequency
C      wB0=omega_B0 is the grad B drift of the deeply trapped pariticles with x=1,
C      qwe=q*omega_E is the E \times B drift
C      qws=q*omega_*p, is the diamagnetic frequency
C      qwsT=q*omega_*T,
C   nu_prof(:,15) is the magnetic shear s=(r/q)dq/dr.
C      
C    bc(n_m,n_n,nr),bs(n_m,n_n,nr) is the spectrum of the magnetic
C                      field strength. 
C
C  in(n_in) is the toroidal mode numbers included in the NTV calculation
C
C
C    ic is a control parameter. 
C          if ic=1, Shaing's formula will be used.                       
C   
C   Output:
C                T_ntv_non(nr): non-resonant part of NTV torque 
C                T_ntv_res(nr): resonant part of NTV torque 
C
C             by Y. Sun
C**********************************************************************
      use common_ntv
      implicit none
      
      integer i,j,n,kprint
      real q,nu,eps,qwe,wB0,shear
      real rho_i,wti,qws,qwsT
      real x(nx),f_nu(nx),nud(nx)
      real t_bc(n_m,2),t_bs(n_m,2)
      real k0,k_ind(nk),delk(nk),k2(nk),KK1(nk),EE1(nk),KK2(nk),EE2(nk)
      complex*16 T_tmp1,T_tmp2,T_tmp3,I_kn(nx),lambda(6)

C     generate x grid
      call xgrid(nx,0.0,15.0,x)

C     calculate the coefficient 
      call f_nud(x,nx,Z,f_nu)

C     calculate (KK1,EE1) for later usage
      k0=-5.0                   
      call xgrid_f(nk,0.0,k0,k_ind)
      delk=10.0**k_ind
      k2=1.0-delk
      k2(1)=0.0                
      call ellipke(k2,nk,KK1,EE1)

C     calculate (KK2,EE2) for later usage
      call xgrid(nk,0.0,1.0,k2)
      call ellipke(k2,nk,KK2,EE2)

C     calculate the NTV torque density
      do 11 i=1,nr
         !-------get profiles
            q=nu_prof(i,2)
          eps=nu_prof(i,3)
        rho_i=nu_prof(i,4)
 	   nu=nu_prof(i,5)
          wti=nu_prof(i,6)
          wB0=nu_prof(i,7)
          qwe=nu_prof(i,8)
          qws=nu_prof(i,9)
         qwsT=nu_prof(i,10)
        shear=nu_prof(i,15)
         !-------
         nud=0.5*nu/eps*f_nu
	 T_ntv_non(i)=(0.0,0.0)
	 T_ntv_res(i)=(0.0,0.0)

C        printing...
         if (i.eq.50) then
            kprint=1
         else
            kprint=0
         endif         

      do 11 j=1,n_in
         n=nin(j)
         t_bc(:,1)=bc(:,j,i)
         t_bs(:,1)=bs(:,j,i)
         t_bc(:,2)=bc2(:,j,i)
         t_bs(:,2)=bs2(:,j,i)
         call gntv(n_m,mm,t_bc,t_bs,q,n,nk,nx,shear,
     &             x,nud,qwe,wB0,qws,qwsT,wti,rho_i,R_axis,
     &             eps,norm_Z,ic,I_kn,lambda,T_tmp1,T_tmp2,T_tmp3,
     &             KK1,EE1,KK2,EE2,kprint,keyntv)        
         T_ntv_non(i)=T_ntv_non(i)+T_tmp1    
         T_ntv_res(i)=T_ntv_res(i)+T_tmp2    
         rot_diff(i,j)=real(T_tmp3)
11    continue
C      print*,'... end T_ntv calculation.'
      end subroutine gntv1d

      subroutine gntv(n_m,mm,bc,bs,q,n,nk,nx,shear,x,nud,
     &                qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z,
     &                ic,I_kn,lambda,T_ntv1,T_ntv2,T_ntv3,
     &                KK1,EE1,KK2,EE2,kprint,keyntv)
C***************************************************************************
C      calculate the NTV torque density
C
C      need to call:
C           Ik_shaing, Pitch angle integral from Shaing''s connected formula.
C           intE, energy integral to calculate lambda.
C     
C      Input:
C
C     n, is the toroidal mode number
C     mm(n_m) is the poloidal mode number
C     bc(n_m),bs(n_m) are the n^th harmonic of the spectrum 
C                     of the magnetic field strength. 
C
C     nk, is the number of the grids in kappa space
C     shear=(r/q)dq/dr sets the precession of the resonant pitch.
C
C      nx is the number of the grids in energy space.
C      x=v^2/v_t^2, is the normalized energy
C      q is the safty factor, eps=epsillon~ r/R
C      rho_i=n_i*m_i, nu is the collisionality, wti is the trasit frequency
C      wB0=omega_B0 is the grad B drift of the deeply trapped pariticles with x=1,
C      qwe=q*omega_E is the E \times B drift
C      qws=q*omega_*p, is the diamagnetic frequency
C      qwsT=q*omega_*T,
C      nud(x) is the deflection frequency
C      nrom_Z=Z_j/Z_i.
C
C      ic is a control parameter. 
C          if ic=1, Shaing''s formula will be used
C
C       Output:
C          I_kn, pitch angle integral
C        lambda, energy integral
C         T_ntv1, non-resonant part of NTV torque density
C         T_ntv2, resonant part of NTV torque density
C
C        by Y. Sun
C***************************************************************************
      implicit none
      real*8,parameter::pi=3.1415926535897932
      integer n_m,mm(n_m),n,nk,nx,ic,nt,kprint,k,keyntv
      real bc(n_m,2),bs(n_m,2),q,shear
      real x(nx),nud(nx)
      real qwe,wB0,qws,qwsT,wti,rho_i,R0,eps,norm_Z
      logical resx(nx)
      real min_t,max_t
      real,allocatable::theta(:),bnc(:,:),bns(:,:)
      real T0,w1,w2
      real KK1(nk),EE1(nk),KK2(nk),EE2(nk)
      complex*16 T_ntv1,T_ntv2,T_ntv3,I_kn(nx),lambda(6)

C     generate grids for theta and spectrum transform
C     form bn(m) to bn(theta)
      nt=2*nk
      allocate(theta(nt))
      allocate(bnc(nt,2))
      allocate(bns(nt,2))
      min_t=-pi
      max_t=pi
      call xgrid_f(nt,min_t,max_t,theta)
      call spec_trans(n_m,mm,bc,bs
     &     ,q,n,nt,theta,bnc,bns)
!-------calculate pitch angle integral
      if (ic.EQ.1)then          !  Shaing's formula
        call Ik_shaing(n,nt,theta,bnc,bns,x,nx,nud,qwe,wB0,shear,
     &                 nk,resx,I_kn,
     &                 KK1,EE1,KK2,EE2,kprint,keyntv)

      elseif (ic.EQ.2)then
                  !     numerical
      end if
!--------calculate energy integeral
      call intE(nx,x,resx,I_kn,lambda)
!--------calculate torque density
      w1=qwe-qws
      w2=-qwsT 
      T0=-0.5/sqrt(2.0)/pi**(1.5)*rho_i*R0**2
     .     *sqrt(eps)*abs(norm_Z)*q**2*wti**2
      T_ntv1=T0*(lambda(1)*w1+lambda(2)*w2)
      T_ntv2=T0*(lambda(4)*w1+lambda(5)*w2)
      T_ntv3=w1+w2*(lambda(2)+lambda(5))/(lambda(1)+lambda(4))

      deallocate(theta,bnc,bns)

      end subroutine gntv

      subroutine Ik_shaing(n,nt,theta,bnc,bns,x,nx,nud,qwe,wB0,shear,
     &                     nk,resx,I_kn,
     &                     KK1,EE1,KK2,EE2,kprint,keyntv)
C***************************************************
C      calculate the Pitch angle integeral 
C      according to Shaing''s Formula.
C        based on the version in MATLAB
C       
C      input:
C
C
C      output: I_kn(x), Pitch angle integeral 
C
C        by Y. Sun
C        Mar. 01, 2011
C***************************************************
      use NTV_lib
      implicit none
      real*8,parameter::pi=3.1415926535897932
      integer,parameter::maxroot=4
      integer i,n,nx,nk,nt,keyntv,r,nroot(nx)
      real qwe,wB0,shear,xmin,fac_refine
      real,dimension(1:nx)::x,nud
     &     ,fac_1nu,fac_nu,nusd,delk_i,fac_sbp
      real,dimension(1:nk)::tmp
      real kap(nx,maxroot),KX(nx,maxroot),dGdk(nx,maxroot)
      real bncb(nx,2),bnsb(nx,2)
      logical resx(nx)
      real k0,k_ind(nk),k2(nk),K(nk),E(nk),F(nk),delk(nk)
      real theta(nt),min_t,max_t,phi(nt),bnc(nt,2),bns(nt,2)
      real KK1(nk),EE1(nk),KK2(nk),EE2(nk)
      complex*16,dimension(1:nx)::k_1nu,k_nu,k_sbp,k_non,k_res,I_kn
      complex*16 int_eta,int_bdy(nx)

      integer kprint
C------------------------prepare
!------generate grids for kappa (k2,delk,K,F)for bounce average.
      k0=-5.0                   !  boundary for kappa^2, 1-10**k0 
      call xgrid_f(nk,0.0,k0,k_ind)
      delk=10.0**k_ind
      k2=1.0-delk
      k2(1)=0.0                ! avoid numerical error
C     call ellipke(k2,nk,K,E)
      K = KK1
      E = EE1
      F=2.0*(E/K-1.0+k2)
!------generate grids for phi 
      min_t=-pi/2.0
      max_t=pi/2.0
      call xgrid(nt,min_t,max_t,phi)
!-------calculate x_min
      xmin=qwe/max(abs(wB0),1.0e-15) !  avoid sigularity 
      xmin=sign(xmin,wB0*qwe)
C------------------------I_kn calculation
C---- 1/nu regime
      call bavgb_eta(n,nt,bnc,bns,theta,phi,nk,int_eta,KK2,EE2,
     &               kprint,keyntv)
      fac_1nu=1.0/nud
      k_1nu=fac_1nu*int_eta
C----nu-sqrt(nu) regime
      fac_refine=0.0  
      nusd=nud/n/abs(-qwe+fac_refine*wB0*x)  
      ! refinement in Shaing's 2010 NF2 paper
      fac_nu=nusd**2/nud
      nusd=nusd/(1.0+8.0*n*nusd) 
      ! It is used to avoid problem in the big nusd case, in which I_nu is not useful.  
      !  The 8*n fator comes from the difference in the diffinition of nu_{*d}.
      !  It keeps it equivalent to the formula in Shaing's 2010 NF paper. 
      delk_i=sqrt(8.0*nusd/log(16.0/sqrt(8.0*nusd))) 
      call bavgb_bdy(n,nt,bnc,bns,theta,phi,nk
     .               ,k2,K,F,nx,delk_i,int_bdy,keyntv)
      k_nu=fac_nu*int_bdy                
C---- superbanana-plateau regime
C     Resonance x*wB0*G(kappa0) = -qwe.  The shear part of G is the q'
C     term of dJ/dpsi (Martitsch 2016, Albert 2016); for |s| > 1/2 G is
C     not monotone, so a pitch may have two resonances and x < |xmin|
C     may resonate.
      call kappa0_roots(nx,x,xmin,shear,nk,k2,K,E,maxroot,
     &                  nroot,kap,KX,dGdk)
      resx=nroot.gt.0
      k_sbp=(0.,0.)
      do r=1,maxval(nroot)
         call bavgb(nt,bnc,bns,theta,phi,nx,kap(:,r),KX(:,r),
     &              bncb,bnsb)
         fac_sbp=abs(4.0*pi*KX(:,r)*n/(wB0*x*dGdk(:,r)))
         where (nroot.lt.r) fac_sbp=0.0
         if (keyntv.ne.2) k_sbp=k_sbp+fac_sbp*(bncb(:,1)**2+
     &                                         bnsb(:,1)**2)
         if (keyntv.eq.2) k_sbp=k_sbp+fac_sbp*(bncb(:,1)*bncb(:,2)+
     &                                bnsb(:,1)*bnsb(:,2) +
     &                       (0.,1.)*(bnsb(:,1)*bncb(:,2)-
     &                                bncb(:,1)*bnsb(:,2))) 
      end do
C---- combined formula
      k_non=k_nu/(1.0+k_nu/k_1nu)
      k_res=k_sbp/(1.0+k_sbp/k_1nu)    
      do i=1,nx
         if(.not.resx(i))then
            I_kn(i)=k_non(i)
         else
            I_kn(i)=k_res(i)
         end if
      end do

      end subroutine Ik_shaing

      subroutine check_k0(nx,k02,delk1,cc)
C***************************************************
C    check the kappa_0^2 value near k^2=1
C            avoid problem in the interpolating
C***************************************************
      implicit none
      integer nx,i,cc,j
      real k02(nx),delk1
      j=0
      do i=1,nx
         if (k02(i).LT.delk1)then
            k02(i)=delk1
            j=j+1
         end if
      end do
         if(j.GT.0)then
           if (cc.EQ.1)then
c$$$            print*,'Warning, some of the 1-k_0^2 < min(dk^2)=',delk1
c$$$            print*,'it is forced to become',delk1
           else
            print*,'Warning, some of the dk_bdy^2 < min(dk^2)=',delk1
            print*,'it is changed to',delk1
           end if
         end if
      
      end subroutine check_k0

      subroutine intE(nx,x,resx,k,lambda)
C***************************************************
C    energy integral 
C***************************************************
      implicit none
      integer nx,i,j
      real x(nx),dx
      logical resx(nx)
      complex*16 k(nx),lambda(6),f
      dx=x(2)-x(1)
      lambda = (0.,0.)
      do i=1,3
           do j=1,nx
              f=0.5*k(j)*(x(j)-2.5)**(i-1.0)*x(j)**(2.5)*exp(-x(j))
              if (.not.resx(j)) then
                 lambda(i)=lambda(i) + f*dx
              else
                 lambda(i+3)=lambda(i+3) + f*dx
              endif
           enddo
      enddo
      end subroutine intE

      subroutine kappa0_roots(nx,x,xmin,shear,nk,k2,K,E,maxroot,
     &                        nroot,kap,KX,dGdk)
C***************************************************
C      trapped pitches kappa0^2 in [0,1) with G(kappa0^2) = -xmin/x,
C      G = 2E/K-1 + 4 shear (E/K+kappa^2-1), dGdk = dG/dkappa^2.
C      Unused slots hold a valid dummy pitch.
C***************************************************
      implicit none
      integer nx,nk,maxroot,nroot(nx),i,j,it
      real x(nx),xmin,shear,k2(nk),K(nk),E(nk)
      real kap(nx,maxroot),KX(nx,maxroot),dGdk(nx,maxroot)
      real m(nk+1),g(nk+1),target,ua,ub,fa,fm,dedm,m1(1),k1(1),e1(1)
      m(1:nk)=k2
      g(1:nk)=2.0*E/K-1.0+4.0*shear*(E/K+k2-1.0)
C     analytic limit kappa^2 -> 1, E/K -> 0
      m(nk+1)=1.0
      g(nk+1)=-1.0
      nroot=0
      kap=0.5
      m1=0.5
      call ellipke(m1,1,k1,e1)
      KX=k1(1)
      dGdk=1.0
      do i=1,nx
         target=-xmin/x(i)
         do j=1,nk
            if (((g(j).gt.target).eqv.(g(j+1).gt.target)).or.
     &          nroot(i).eq.maxroot) cycle
C           bisection in log(1-kappa^2), down to 1-kappa^2 = 1e-16
            ua=log(1.0-m(j))
            ub=log(max(1.0-m(j+1),1.0e-16))
            fa=g(j)-target
            do it=1,60
               m1=1.0-exp(0.5*(ua+ub))
               call ellipke(m1,1,k1,e1)
               fm=2.0*e1(1)/k1(1)-1.0+4.0*shear*(e1(1)/k1(1)+m1(1)-1.0)
     &            -target
               if ((fm.gt.0.0).eqv.(fa.gt.0.0)) then
                  ua=0.5*(ua+ub)
                  fa=fm
               else
                  ub=0.5*(ua+ub)
               end if
            end do
            m1=1.0-exp(0.5*(ua+ub))
            call ellipke(m1,1,k1,e1)
            nroot(i)=nroot(i)+1
            kap(i,nroot(i))=m1(1)
            KX(i,nroot(i))=k1(1)
            dedm=-((e1(1)/k1(1)-1.0)**2/m1(1)+2.0*e1(1)/k1(1)-1.0)
     &           /(2.0*(1.0-m1(1)))
            dGdk(i,nroot(i))=(2.0+4.0*shear)*dedm+4.0*shear
         end do
      end do
      end subroutine kappa0_roots


      subroutine f_nud(x,nx,z,f)
C***************************************************
C      nu_d=f*nu_0, deflection frequency
C       z>0  for ion
C       z<0  for electron
C***************************************************
      use NTV_lib
      implicit none
      real*8,parameter::pi=3.1415926535897932
      integer nx,i,z
      real x(nx),f(nx),tmp    
      if(x(1).lt.1.0e-15) stop 'energy x is set too small'
      do i=1,nx
         f(i)=(1-0.5/x(i))*erf(sqrt(x(i)))+exp(-x(i))/sqrt(pi*x(i))
      end do
      if(z>0)then
         f=x**(-1.5)*f
      else
         f=x**(-1.5)*(1+f)
      endif
      end subroutine f_nud

      subroutine bavgb_bdy(n,ng,bnc,bns,theta,phi
     &                     ,nk,k2,K,F,ndk,delk,int_bdy,keyntv)
C***************************************************
C    bounce averaged bn*(1-exp(-(1+i)*y))
C      for the NTV calculation in nu-sqrt(nu) regime
C***************************************************
      implicit none
      integer ng,nk,i,j,n,ndk,l,keyntv
      real bnc(ng,2),bns(ng,2),theta(ng),phi(ng)
      real delk(ndk),K(nk),bncb(nk,2),bnsb(nk,2)
      real bncbs(nk,ndk,2),bnsbs(nk,ndk,2)
      real y(nk),y1(nk),y2(nk),k2(nk),F(nk),c0(nk)
      real dbncb(nk,ndk,2),dbnsb(nk,ndk,2)
      complex*16 int_bdy(ndk),dbnb2(nk,ndk),tmp2(nk),temp_bdy

      call bavgb(ng,bnc,bns,theta,phi,nk,k2,K,bncb,bnsb)
    
      do 10 l=1,2
      do 10 j=1,ndk
         y=(1.0-k2)/delk(j)
         y1=1.0-exp(-y)*cos(y)
         y2=exp(-y)*sin(y)
         bncbs(:,j,l)=bncb(:,l)*y1-bnsb(:,l)*y2
         bnsbs(:,j,l)=bnsb(:,l)*y1+bncb(:,l)*y2
 10   continue

      do l=1,2
         call cald(bncbs(:,:,l),nk,ndk,k2,dbncb(:,:,l))
         call cald(bnsbs(:,:,l),nk,ndk,k2,dbnsb(:,:,l))
      enddo
      
      if (keyntv.ne.2) dbnb2=dbncb(:,:,1)**2+dbnsb(:,:,1)**2
      if (keyntv.eq.2) dbnb2=dbncb(:,:,1)*dbncb(:,:,2)+
     &                       dbnsb(:,:,1)*dbnsb(:,:,2) +
     &              (0.,1.)*(dbnsb(:,:,1)*dbncb(:,:,2)-
     &                       dbncb(:,:,1)*dbnsb(:,:,2))
      
      c0=4.0*K*F
      do i=1,ndk
         tmp2=n**2*c0*dbnb2(:,i)
         call calint(tmp2,k2,nk,temp_bdy)
         int_bdy(i)=temp_bdy
      end do
      end subroutine bavgb_bdy

      subroutine bavgb_eta(n,ng,bnc,bns,theta,phi,nk,int_eta,
     &                     KK2,EE2,kprint,keyntv)
C***********************************************************
C    bounce averaged bn*eta for NTV in the 1/nu regime
C***********************************************************
      use NTV_lib
      implicit none
      integer ng,nk,i,n,j,kprint,l,keyntv
      real bnc(ng,2),bns(ng,2),theta(ng),phi(ng)
      real k2(nk),K(nk),E(nk),F(nk),dk,ik2,iK
      real tmp(ng),the(ng),bnc0(ng,2),bns0(ng,2),bncb(2),bnsb(2),eta(ng)
      real tmp_f(ng,2),nthe(ng)
      real tmp2c(ng,2),tmp2s(ng,2)
      real KK2(nk),EE2(nk)
      real,parameter::yp1=-1.0e30,ypn=-1.0e-30
      complex*16 int_eta,bnb2(nk)

      do l=1,2
      call spline(theta,bnc(:,l),ng,yp1,ypn,tmp2c(:,l))
      call spline(theta,bns(:,l),ng,yp1,ypn,tmp2s(:,l))
      enddo
      call xgrid(nk,0.0,1.0,k2)

C     call ellipke(k2,nk,K,E)
      K=KK2
      E=EE2

      F=2.0*(E/K-1.0+k2)
      
      do 10 i=1,nk
         ik2=k2(i)
         iK=K(i)
         the=2.0*asin(sqrt(ik2)*sin(phi))
         eta=abs(ik2-sin(the/2.0)**2)
         nthe=the
         do 20 l=1,2
         do 20 j=1,ng
            call splint(theta,bnc(:,l),tmp2c(:,l),ng,the(j),bnc0(j,l))
            call splint(theta,bns(:,l),tmp2s(:,l),ng,the(j),bns0(j,l))
 20      continue    
         tmp_f(:,1)=eta*bnc0(:,1)
         tmp_f(:,2)=eta*bnc0(:,2)
         call bavg(ik2,iK,ng,phi,tmp_f,bncb)
         tmp_f(:,1)=eta*bns0(:,1)
         tmp_f(:,2)=eta*bns0(:,2)
         call bavg(ik2,iK,ng,phi,tmp_f,bnsb)

         if (keyntv.ne.2) bnb2(i)=bncb(1)**2+bnsb(1)**2
         if (keyntv.eq.2) bnb2(i)=bncb(1)*bncb(2)+bnsb(1)*bnsb(2) + 
     &                   (0.,1.)*(bnsb(1)*bncb(2)-bncb(1)*bnsb(2))
 10   continue    

      dk=k2(2)-k2(1)
      int_eta=n**2*sum(16.0*K*bnb2/F*dk)

      end subroutine bavgb_eta

      subroutine bavgb(ng,bnc,bns,theta,phi,nk,k2,K,bncb,bnsb)
C***************************************************
C    bounce averaged bn
C      for the ntv calculation in the sb-p regime
C***************************************************
      implicit none
      integer ng,nk,i,j,l
      real bnc(ng,2),bns(ng,2),theta(ng),phi(ng),k2(nk),K(nk)
      real tmp(ng),the(ng),bnc0(ng,2),bns0(ng,2),tmpc(2),tmps(2),ik2,iK
      real bncb(nk,2),bnsb(nk,2)
      real tmp2c(ng,2),tmp2s(ng,2)
      real,parameter::yp1=-1.0e30,ypn=-1.0e-30

      do l=1,2
      call spline(theta,bnc(:,l),ng,yp1,ypn,tmp2c(:,l))
      call spline(theta,bns(:,l),ng,yp1,ypn,tmp2s(:,l))
      enddo

      do 10 i=1,nk
         ik2=k2(i)
         iK=K(i)
         the=2.0*asin(sqrt(ik2)*sin(phi))
         do 20 l=1,2
         do 20 j=1,ng
          call splint(theta,bnc(:,l),tmp2c(:,l),ng,the(j),bnc0(j,l))      
          call splint(theta,bns(:,l),tmp2s(:,l),ng,the(j),bns0(j,l))      
 20      continue    
         call bavg(ik2,iK,ng,phi,bnc0,tmpc)
         call bavg(ik2,iK,ng,phi,bns0,tmps)
         do l=1,2
            bncb(i,l)=tmpc(l)
	    bnsb(i,l)=tmps(l)
         enddo
 10 	continue
      
      end subroutine bavgb

      subroutine bavg(ik2,iK,ng,phi,f,bf)
C***************************************************
C    bonce averge of f
C       output: bf
C***************************************************
      implicit none
      integer ng
      real ik2,iK,phi(ng),dphi,f(ng,2),bf(2)
      dphi=phi(2)-phi(1)
      bf(1)=sum(0.5/iK*f(:,1)/sqrt(abs(1.0-ik2*sin(phi)**2))*dphi)
      bf(2)=sum(0.5/iK*f(:,2)/sqrt(abs(1.0-ik2*sin(phi)**2))*dphi)
      end subroutine bavg

      subroutine spec_trans(n_m,mm,bc,bs,q,n,ng,theta,bnc,bns)
C***************************************************
C    transform the spectrum bmn(m)
C              to bn(theta)
C***************************************************
      implicit none
      integer n_m,mm(n_m),ng,i,n,k
      real q,bc(n_m,2),bs(n_m,2),theta(ng),bnc(ng,2),bns(ng,2)
      do 10 k=1,2
      do 10 i=1,ng
         bnc(i,k)=sum(bc(:,k)*cos((mm-n*q)*theta(i))-
     &                bs(:,k)*sin((mm-n*q)*theta(i)))
         bns(i,k)=sum(bs(:,k)*cos((mm-n*q)*theta(i))+
     &                bc(:,k)*sin((mm-n*q)*theta(i)))
 10   continue
      end subroutine spec_trans


      end




