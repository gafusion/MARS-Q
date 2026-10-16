module manufactured_force
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    implicit none
    complex(dp),allocatable :: rhs_x(:,:),rhs_y(:,:)
    integer :: calls=0
end module

program lifecycle
    use manufactured_force
    implicit none
    integer,parameter::md=2,my=1,nd=2
    complex(dp)::a(md,md,nd+3),b(md,md,nd+3),c(md,md,nd+2)
    complex(dp)::d(my,my,nd+2),e(md,my,nd+2),f(my,md,nd+2)
    complex(dp)::g(my,md,nd+2),h(md,my,nd+3)
    complex(dp)::dx(md,nd+1),dy(my,nd+1),x(md,nd+1),y(my,nd+1)
    complex(dp)::dix(md,nd+1),diy(my,nd),r(md,nd+1),ry(my,nd)
    complex(dp)::bt1(md,md),bty1(md,md),cte1(md,md),cty1(md,my),dt1(md),dty1(my)
    complex(dp)::sc1(md,2*md+1),scy1(my,2*my+1)
    complex(dp)::bt2(md,md),bty2(md,md),cte2(md,md),cty2(md,my),dt2(md),dty2(my)
    complex(dp)::sc2(md,2*md+1),scy2(my,2*my+1)
    complex(dp)::xold(md,nd+1),yold(my,nd),xpivot(md),ypivot(my)
    complex(dp)::tmp1(md),tmp2(md),tmpy1(my),tmpy2(my),al0,alam,alnorm
    complex(dp)::work(10000),phase,asave(size(a)),bsave(size(b)),hsave(size(h))
    integer::neg,negeig,nit,nxdev,nydev,nxpneg,nypneg,nev,m,n,mm_y
    common /auxint/ neg,negeig,nit,nxdev,nydev,nxpneg,nypneg,nev,m,n,mm_y
    integer::k,i,nxc,nyc,ncase,nitmax,noncon
    real(dp)::epspam,epsdet
    character(64)::backend,argument
    real(dp)::family_eps
    call get_command_argument(1,argument)
    if(len_trim(argument)>0)then
        read(argument,*)family_eps
        call run_family(family_eps)
        stop 0
    end if
    allocate(rhs_x(md,nd+1),rhs_y(my,nd))
    m=md;mm_y=my;n=nd
    nxc=1;nyc=1;ncase=2;nitmax=1
    epspam=1e-7_dp;epsdet=1e-199_dp
    alam=1;al0=cmplx(0.2_dp,0.3_dp,dp);alnorm=cmplx(7.0_dp,-9.0_dp,dp)
    a=0;b=0;c=0;d=0;e=0;f=0;g=0;h=0
    ! Poison unused guards: real callers pass A/H starting at radial index2.
    a(:,:,1)=999;h(:,:,1)=-777
    do k=1,nd+1
        do i=1,md
            b(i,i,k)=4
            if(k<=nd)then
                a(i,i,k+1)=cmplx(-0.1_dp,0.2_dp,dp)
                c(i,i,k)=cmplx(0.2_dp,-0.1_dp,dp)
            end if
        end do
        if(k>nd)cycle
        d(1,1,k)=3
        e(:,1,k)=[cmplx(0.2_dp,0.1_dp,dp),cmplx(0.3_dp,0.0_dp,dp)]
        h(:,1,k+1)=[cmplx(0.2_dp,0.0_dp,dp),cmplx(0.0_dp,0.3_dp,dp)]
        f(1,:,k)=[cmplx(0.4_dp,0.0_dp,dp),cmplx(0.0_dp,0.2_dp,dp)]
        g(1,:,k)=[cmplx(-0.1_dp,0.0_dp,dp),cmplx(0.2_dp,0.0_dp,dp)]
    end do
    phase=exp(cmplx(0.0_dp,0.37_dp,dp))
    do k=1,nd+1
        x(:,k)=phase*[cmplx(k,0.4_dp,dp),cmplx(-2.0_dp*k,0.1_dp,dp)]
        if(k<=nd)y(1,k)=phase*cmplx(-k,0.2_dp*k,dp)
    end do
    xold=x;yold=y(:,1:nd)
    call muamat(md,my,nd,nxc,nyc,ncase,nitmax,epspam,epsdet,al0,alam,alnorm,noncon, &
        a(1,1,2),b,c,d,e,f,g,h(1,1,2),dx,dy,x,y,dix,diy,r,ry, &
        bt1,bty1,cte1,cty1,dt1,dty1,sc1,scy1,bt2,bty2,cte2,cty2,dt2,dty2,sc2,scy2, &
        xold,yold,xpivot,ypivot,tmp1,tmp2,tmpy1,tmpy2)
    rhs_x=r;rhs_y=ry
    asave=reshape(a,[size(a)]);bsave=reshape(b,[size(b)]);hsave=reshape(h,[size(h)])
    x=0;y=0;dx=huge(1.0_dp);dy=-huge(1.0_dp);work=0
    nit=-99;noncon=-99
    call calpam(md,my,nd,nxc,nyc,ncase,nitmax,epspam,epsdet,al0,alam,alnorm,noncon, &
        a(1,1,2),b,c,d,e,f,g,h(1,1,2),dx,dy,x,y,work,size(work))
    if(calls/=1)error stop 'Native UBXRHS did not form forcing exactly once'
    if(nit/=1.or.nitmax/=1)error stop 'Forced-response iteration lifecycle changed'
    if(alam/=cmplx(1.0_dp,0.0_dp,dp))error stop 'ALAM changed'
    if(alnorm/=cmplx(7.0_dp,-9.0_dp,dp))error stop 'ALNORM changed'
    ! xold/yold here are independent caller arrays, not CALPAM workspace views.
    if(any(abs(x-xold)>1e-12_dp*abs(xold)))error stop 'X mapping or complex phase changed'
    if(any(abs(y(:,1:nd)-yold)>1e-12_dp*abs(yold)))error stop 'Y mapping or complex phase changed'
    call get_environment_variable('MARS_PAMS_BACKEND',backend)
    if(trim(backend)=='UMFPACK_IPC')then
        if(noncon/=0)error stop 'Sparse success must report NONCON=0'
        if(any(reshape(a,[size(a)])/=asave).or.any(reshape(b,[size(b)])/=bsave).or. &
           any(reshape(h,[size(h)])/=hsave))error stop 'Original operator was mutated'
    else
        if(noncon/=3)error stop 'Default NCASE2 lifecycle changed'
    end if
    write(*,*)'PASS native CALPAM/PAMERA/UBXRHS lifecycle ',trim(backend)
end program

subroutine run_family(eps)
    use manufactured_force
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(dp),intent(in)::eps
    integer,parameter::md=1,my=1,nd=1,wide=selected_real_kind(30)
    complex(dp)::a(md,md,nd+3),b(md,md,nd+3),c(md,md,nd+2)
    complex(dp)::d(my,my,nd+2),e(md,my,nd+2),f(my,md,nd+2)
    complex(dp)::g(my,md,nd+2),h(md,my,nd+3)
    complex(dp)::dx(md,nd+1),dy(my,nd+1),x(md,nd+1),y(my,nd+1),work(10000)
    complex(dp)::al0,alam,alnorm,full(3,3),copy(3,3),rhs(3,1),exact(3),got(3)
    complex(wide)::aq(3,3),bq(3),xq(3),residual(3)
    real(dp)::forward_error,backward_error
    integer::noncon,info,ipiv(3)
    integer::neg,negeig,nit,nxdev,nydev,nxpneg,nypneg,nev,m,n,mm_y
    common /auxint/ neg,negeig,nit,nxdev,nydev,nxpneg,nypneg,nev,m,n,mm_y
    allocate(rhs_x(md,nd+1),rhs_y(my,nd))
    ! Known well-conditioned original system in order (X1,Y1,X2).
    full(1,:)=[(1.0_dp,0.0_dp),(1.0_dp,0.0_dp),(0.0_dp,0.0_dp)]
    full(2,:)=[(1.0_dp,0.0_dp),cmplx(eps,0.0_dp,dp),(1.0_dp,0.0_dp)]
    full(3,:)=[(0.0_dp,0.0_dp),(1.0_dp,0.0_dp),(1.0_dp,0.0_dp)]
    exact=[(1.0_dp,0.0_dp),(2.0_dp,0.0_dp),(3.0_dp,0.0_dp)]
    rhs(:,1)=matmul(full,exact)
    rhs_x(1,:)=[rhs(1,1),rhs(3,1)]
    rhs_y(1,1)=rhs(2,1)
    aq=cmplx(real(full,dp),aimag(full),wide)
    bq=cmplx(real(rhs(:,1),dp),aimag(rhs(:,1)),wide)
    copy=full
    call zgetrf(3,3,copy,3,ipiv,info)
    if(info/=0)error stop 'Independent global LU factorization failed'
    call zgetrs('N',3,1,copy,3,ipiv,rhs,3,info)
    if(info/=0.or.maxval(abs(rhs(:,1)-exact))>1e-14_dp) &
        error stop 'Independent global LU solution failed'
    write(*,'(a,es25.17)')'FAMILY_GLOBAL_LU_MAXERR ',maxval(abs(rhs(:,1)-exact))
    a=0;b=0;c=0;d=0;e=0;f=0;g=0;h=0
    a(:,:,1)=999;h(:,:,1)=-777
    b(1,1,1:2)=1;d(1,1,1)=eps
    e(1,1,1)=1;f(1,1,1)=1;g(1,1,1)=1;h(1,1,2)=1
    x=0;y=0;dx=huge(1.0_dp);dy=-huge(1.0_dp);work=0
    alam=1;al0=cmplx(0.2_dp,0.3_dp,dp);alnorm=cmplx(7.0_dp,-9.0_dp,dp)
    call calpam(md,my,nd,1,1,2,1,1e-7_dp,1e-199_dp,al0,alam,alnorm,noncon, &
        a(1,1,2),b,c,d,e,f,g,h(1,1,2),dx,dy,x,y,work,size(work))
    got=[x(1,1),y(1,1),x(1,2)]
    xq=cmplx(real(got,dp),aimag(got),wide)
    residual=bq-matmul(aq,xq)
    forward_error=maxval(abs(got-exact))
    backward_error=real(maxval(abs(residual)/(abs(bq)+matmul(abs(aq),abs(xq)))),dp)
    write(*,'(a,3es25.17)')'FAMILY_ORIGINAL_ERRORS ',eps,forward_error,backward_error
    if(.not.ieee_is_finite(forward_error).or..not.ieee_is_finite(backward_error)) &
        error stop 'Native PAMS family nonfinite original-system error'
    if(forward_error>1e-11_dp.or.backward_error>1e-12_dp) &
        error stop 'Native PAMS family original-system accuracy failure'
    if(calls/=1.or.nit/=1)error stop 'Family forcing/iteration lifecycle changed'
    if(alam/=cmplx(1.0_dp,0.0_dp,dp).or.alnorm/=cmplx(7.0_dp,-9.0_dp,dp)) &
        error stop 'Family forcing normalization changed'
    write(*,*)'PASS_NATIVE_CALLER_FAMILY'
end subroutine
