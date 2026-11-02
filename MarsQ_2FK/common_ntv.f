      module common_ntv
C****************************************************************
C define common variants of the input for the NTV module
C         prepare data for NTV calculation
C
C
C      cal_nu(), calculate the related collisionality and frequencies
C
C        by Y. Sun
C        Mar. 01, 2011
C****************************************************************
C NR      =      65,  ! No. of the radial grids for NTV
C NX      =     101,  ! No. of the grids in energy space
C NK      =     201,  ! No. of the grids in pitch angle space
C N_CHI   =     128,  ! number of poloidal grids, 2^n
C N_PHI   =      64,  ! number of toroidal grids, 2^n
C Z       =       1,  ! charge number of particles considered. e_j=Z*e0. Z=-1 always for electrons.
C Zi      =       1,  ! charge number of the main ions. e_i=Zi*e0.
C MUI     =     2.0,  ! ion mass, mi=mui*mp
C IC      =       1,  ! ic=1 use Shaings formula
C N_MAX   =       1,  ! maximum toroidal number in the bmn spectrum. nin=n          
C M_MAX   =      25,  ! maximum poloidal number in the bmn spectrum. mm=mm(-M_max:M_max)
C N_IN    =       1,  ! number of the toroidal modes included in the calculation.
C KNC     =    1.17,  !  neoclassical coefficient for poloidal rotation
C IRMP    =    21.0   ! current(kA*turns) of the rmp field 
C NIN     =       1   ! the toroidal mode number inclued in the calculation. NIN=NIN(N_IN)

      use NTV_lib
      implicit none
!---------- parameters
      real*8,parameter::pi=3.1415926535897932
      integer,parameter::n_prof=15
!----------- evaluated in basic_input
      integer,parameter::nx=101
      integer,parameter::nk=201
      integer,parameter::Zi=1
      integer,parameter::ic=1
      integer,parameter::n_max=1
      integer,parameter::n_in=1
      real,parameter::mui=2.0
      integer,save :: Z,nr,m_max,n_chi,keyntv
      integer,save :: n_m,n_n,nin(1)
      real,save :: norm_Z
      integer,allocatable,save:: mm(:),nn(:)
      integer,parameter::n_phi=64
!----------- evaluated in MARS_output
      real,parameter::knc=1.17
C     real,parameter::knc=0.0
      real,allocatable,save::bc(:,:,:),bs(:,:,:)
      real,allocatable,save::bc2(:,:,:),bs2(:,:,:)
      real,allocatable,save::rho(:),nu_prof(:,:)
      real,save :: R_axis
!----------- evaluated in gntv1d and output to MARS-Q
      real,allocatable,save::rot_diff(:,:)
      complex*16,allocatable,save::T_ntv_non(:),T_ntv_res(:)

      contains

      subroutine basic_input(Z_mars,nr_mars,nc_mars,n_mars,
     &                       M_mars,keyntv_mars,s_mars)
C****************************************************************
C           load basic input about dimensions
C****************************************************************
      integer i,Z_mars,nr_mars,nc_mars,n_mars,M_mars,keyntv_mars
      REAL*8,DIMENSION(nr_mars+1)::s_mars

!------- load dimension data
      n_chi = nc_mars + 1
      Z  = Z_mars
      nr = nr_mars
      nin(1) = n_mars
      m_max  = M_mars + 5
      keyntv = keyntv_mars
      if (.not.allocated(rho)) allocate(rho(nr))
      do i=1,nr
         rho(i) = (s_mars(i+1)+s_mars(i))*0.5
      enddo

      norm_Z=1.0*Z/Zi
      n_m=2*m_max+1
      n_n=n_max
      if (.not.allocated(mm)) allocate(mm(n_m))
      if (.not.allocated(nn)) allocate(nn(n_n))
        do i=1,n_m
           mm(i)=i-1-m_max
        end do
        do i=1,n_n
           nn(i)=i
        end do
      if (.not.allocated(bc)) allocate(bc(n_m,n_n,nr))
      if (.not.allocated(bs)) allocate(bs(n_m,n_n,nr))
      if (.not.allocated(bc2)) allocate(bc2(n_m,n_n,nr))
      if (.not.allocated(bs2)) allocate(bs2(n_m,n_n,nr))
      if (.not.allocated(nu_prof)) allocate(nu_prof(nr,n_prof))
      if (.not.allocated(T_ntv_non)) allocate(T_ntv_non(nr),
     &                                        T_ntv_res(nr),
     &                                        rot_diff(nr,n_n))

      bc        = 0.0
      bs        = 0.0
      bc2       = 0.0
      bs2       = 0.0
      nu_prof   = 0.0
      T_ntv_non = (0.0,0.0)
      T_ntv_res = (0.0,0.0)
      rot_diff  = 0.0

      end subroutine basic_input

      subroutine ntv_predata
C*******************************************************************
C    load the plasma data and the bn spectrum for NTV calculation
C*******************************************************************
      use NTV_pre

      call MARS_out(nr,n_chi,n_phi,n_prof,keyntv,rho,
     &              mui,Z,Zi,knc,m_max,n_max,R_axis,nu_prof,
     &              bc,bs,bc2,bs2)    ! load the data from the MARS output
  
      end subroutine ntv_predata


      end
