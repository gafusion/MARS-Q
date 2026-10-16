program pams_global_csc_cases
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer :: i
    real(dp) :: eps
    character(40)::argument
    call get_command_argument(1,argument)
    if(trim(argument)=='fatal-singular')then
        call run_case(1,1,1,5,0.0_dp)
        error stop 'Singular wrapper unexpectedly returned'
    end if
    do i=0,5
        eps=0.0_dp
        if(i>0)eps=10.0_dp**(-6-2*i)
        call run_case(1,1,1,1,eps)
    end do
    call run_case(2,3,3,2,0.0_dp)
    call run_case(3,2,1,2,0.0_dp)
    call run_case(2,1,1,3,0.0_dp)
    call run_case(1,1,1,4,0.0_dp)
    call index_overflow_case()
    write(*,*)'PASS native PAMS global CSC cases'
end program

subroutine run_case(md,my,nd,variant,eps)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use pams_global_csc, only: pams_csc_pack,pams_csc_solve,pams_csc_errors
    use fortsparse_csc, only: csc_z_t,csc_matvec
    implicit none
    integer,intent(in)::md,my,nd,variant
    real(dp),intent(in)::eps
    complex(dp)::a(md,md,nd+2),b(md,md,nd+3),c(md,md,nd+2)
    complex(dp)::d(my,my,nd+2),e(md,my,nd+2),f(my,md,nd+2)
    complex(dp)::g(my,md,nd+2),h(md,my,nd+2)
    complex(dp)::dx(md,nd+1),dy(my,nd+1),x(md,nd+1),y(my,nd+1)
    complex(dp)::dix(md,nd+1),diy(my,nd),r(md,nd+1),ry(my,nd)
    complex(dp)::bt1(md,md),bty1(md,md),cte1(md,md),cty1(md,my),dt1(md),dty1(my)
    complex(dp)::sc1(md,2*md+1),scy1(my,2*my+1)
    complex(dp)::bt2(md,md),bty2(md,md),cte2(md,md),cty2(md,my),dt2(md),dty2(my)
    complex(dp)::sc2(md,2*md+1),scy2(my,2*my+1)
    complex(dp)::xold(md,nd+1),yold(my,nd),xpivot(md),ypivot(my)
    complex(dp)::tmp1(md),tmp2(md),tmpy1(my),tmpy2(my),al0,alam,alnorm
    integer::nxc,nyc,ncase,nitmax,noncon
    real(dp)::epspam,epsdet
    integer::neg,negeig,nit,nxdev,nydev,nxpneg,nypneg,nev,m,n,mm_y
    common /auxint/ neg,negeig,nit,nxdev,nydev,nxpneg,nypneg,nev,m,n,mm_y
    type(csc_z_t)::matrix
    complex(dp),allocatable::vector(:),action(:)
    real(dp)::errors(nd*(md+my)+md),scale
    integer::i,j,k,offset,ierr,pos
    m=md
    mm_y=my
    n=nd
    a=0;b=0;c=0;d=0;e=0;f=0;g=0;h=0
    if(variant==1.or.variant==4.or.variant==5)then
        b(1,1,1:nd+1)=1
        d(1,1,1)=eps
        e(1,1,1)=1;f(1,1,1)=1;g(1,1,1)=1;h(1,1,1)=1
        if(variant==4.or.variant==5)then
            b=0;d=0;e=0;f=0;g=0;h=0
        end if
    else if(variant==2)then
        ! Distinct complex entries exercise every block and both radial boundaries.
        do k=1,nd+1
            do j=1,md
                do i=1,md
                    b(i,j,k)=cmplx(0.01_dp*i/(j+k),0.02_dp*j/(i+k),dp)
                    if(i==j)b(i,j,k)=b(i,j,k)+3
                    if(k<=nd)then
                        a(i,j,k)=cmplx(0.03_dp/(i+j+k),-0.01_dp,dp)
                        c(i,j,k)=cmplx(-0.02_dp,0.01_dp/(i+j+k),dp)
                    end if
                end do
            end do
        end do
        do k=1,nd
            do j=1,my
                do i=1,my
                    d(i,j,k)=cmplx(0.02_dp*i/(j+k),-0.01_dp*j/(i+k),dp)
                    if(i==j)d(i,j,k)=d(i,j,k)+2
                end do
                do i=1,md
                    e(i,j,k)=cmplx(0.03_dp/(i+j+k),0.02_dp,dp)
                    h(i,j,k)=cmplx(-0.01_dp,0.03_dp/(i+j+k),dp)
                    f(j,i,k)=cmplx(0.02_dp/(i+j+k),-0.04_dp,dp)
                    g(j,i,k)=cmplx(-0.04_dp,0.02_dp/(i+j+k),dp)
                end do
            end do
        end do
        ! Apply independently varying row scales; exact powers of two avoid noise.
        do k=1,nd+1
            do i=1,md
                scale=2.0_dp**(40*(i-k))
                b(i,:,k)=b(i,:,k)*scale
                if(k>1)then
                    a(i,:,k-1)=a(i,:,k-1)*scale
                    h(i,:,k-1)=h(i,:,k-1)*scale
                end if
                if(k<=nd)then
                    c(i,:,k)=c(i,:,k)*scale
                    e(i,:,k)=e(i,:,k)*scale
                end if
            end do
            if(k>nd)cycle
            do i=1,my
                scale=2.0_dp**(40*(k-i))
                d(i,:,k)=d(i,:,k)*scale
                f(i,:,k)=f(i,:,k)*scale
                g(i,:,k)=g(i,:,k)*scale
            end do
        end do
    else
        ! Two isolated tiny X rows coexist with ordinary nonzero response rows.
        do k=1,nd+1
            b(1,1,k)=1
            b(2,2,k)=1e-200_dp
        end do
        d(1,1,1)=1
    end if
    do k=1,nd+1
        do i=1,md
            x(i,k)=cmplx(i+k,0.1_dp*i,dp)
        end do
        if(k<=nd)then
            do i=1,my
                y(i,k)=cmplx(-i-k,0.2_dp*k,dp)
            end do
        end if
    end do
    if(variant==3)x(2,:)=cmplx(1e-100_dp,0.0_dp,dp)
    xold=x;yold=y(:,1:nd)
    ! Independent oracle: unmodified native MUAMAT on the original blocks.
    call muamat(md,my,nd,nxc,nyc,ncase,nitmax,epspam,epsdet,al0,alam,alnorm,noncon, &
        a,b,c,d,e,f,g,h,dx,dy,x,y,dix,diy,r,ry, &
        bt1,bty1,cte1,cty1,dt1,dty1,sc1,scy1, &
        bt2,bty2,cte2,cty2,dt2,dty2,sc2,scy2, &
        xold,yold,xpivot,ypivot,tmp1,tmp2,tmpy1,tmpy2)
    dx=r;dy(:,1:nd)=ry
    if(variant==5)call pams_global_solve(md,my,nd,a,b,c,d,e,f,g,h,dx,dy(:,1:nd),x,y(:,1:nd))
    call pams_csc_pack(md,my,nd,a,b,c,d,e,f,g,h,matrix,ierr)
    if(ierr/=0)error stop 'CSC assembly failed'
    allocate(vector(matrix%nrow))
    do k=1,nd+1
        offset=(k-1)*(md+my)
        vector(offset+1:offset+md)=x(:,k)
        if(k<=nd)vector(offset+md+1:offset+md+my)=y(:,k)
    end do
    action=csc_matvec(matrix,vector)
    do k=1,nd+1
        offset=(k-1)*(md+my)
        if(any(abs(action(offset+1:offset+md)-r(:,k))>1e-14_dp*max(abs(r(:,k)),tiny(1.0_dp)))) &
            error stop 'CSC X action differs from native MUAMAT'
        if(k>nd)cycle
        if(any(abs(action(offset+md+1:offset+md+my)-ry(:,k))>1e-14_dp*max(abs(ry(:,k)),tiny(1.0_dp)))) &
            error stop 'CSC Y action differs from native MUAMAT'
    end do
    call pams_csc_solve(md,my,nd,a,b,c,d,e,f,g,h,dx,dy(:,1:nd),x,y(:,1:nd),errors,ierr)
    if(variant==4)then
        if(ierr/=1)error stop 'True singular matrix did not return typed singular status'
        write(*,*)'PASS singular status',ierr
        return
    end if
    if(ierr/=0)then
        write(*,*)'FAILED variant,dimensions,eps,status',variant,md,my,nd,eps,ierr
        error stop 1
    end if
    if(any(abs(x-xold)>1e-11_dp*abs(xold)))error stop 'X manufactured solution mismatch'
    if(any(abs(y(:,1:nd)-yold)>1e-11_dp*abs(yold)))error stop 'Y manufactured solution mismatch'
    if(variant==3)then
        x(2,1)=1.05_dp*x(2,1)
        call pams_csc_errors(md,my,nd,a,b,c,d,e,f,g,h,dx,dy(:,1:nd),x,y(:,1:nd),errors)
        if(errors(2)<0.02_dp)error stop 'Tiny row error was hidden by a denominator floor'
        ! Also verify the explicitly defined all-zero row/zero RHS convention.
        b(2,2,1)=0;dx(2,1)=0
        call pams_csc_errors(md,my,nd,a,b,c,d,e,f,g,h,dx,dy(:,1:nd),x,y(:,1:nd),errors)
        if(errors(2)/=0)error stop 'Zero row convention failed'
    end if
    write(*,*)'PASS variant,dimensions,eps',variant,md,my,nd,eps
end subroutine

subroutine index_overflow_case()
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use pams_global_csc, only: pams_csc_pack
    use fortsparse_csc, only: csc_z_t
    implicit none
    complex(dp)::zero(1,1,1)
    type(csc_z_t)::matrix
    integer::md,ierr
    zero=0
    md=huge(0)
    call pams_csc_pack(md,1,1,zero,zero,zero,zero,zero,zero,zero,zero,matrix,ierr)
    if(ierr/=7)error stop 'Index overflow was not rejected before allocation'
    if(allocated(matrix%val).or.allocated(matrix%col_ptr))error stop 'Overflow allocated CSC storage'
    write(*,*)'PASS public index overflow before allocation'
end subroutine
