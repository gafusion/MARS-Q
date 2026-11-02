! Exercise the production KIA_PASS energy integral with the passing magnetic
! drift (KPASSDRIFT=1) on the production KENERGYMESH tangent mesh. Each line
! prints the case parameters and the native result; the Python oracle
! integrates the same physical integrand independently.
program passing_drift
    use rcomdm
    use dimensim
    use globalm
    use kineticm
    use anisotropicm
    implicit none
    integer :: i, j, k, mode, flag
    character(len=16) :: argument
    complex*16 :: value
    real*8, parameter :: drifts(4) = [0.d0, 0.8d0, -0.8d0, 3.d0]
    real*8, parameter :: rotations(2) = [0.37d0, -0.37d0]
    real*8, parameter :: widths(3) = [2.d-1, 1.d-3, 0.d0]

    call get_command_argument(1, argument)
    read(argument, *) nepk
    nepk2 = 2*nepk-1
    allocate(epk(nepk), zepko(nepk2), zepkn(nepk2), zepk(nepk2))
    allocate(zresu(nepk2,3,5,2), rlm(1), slam0(1,2))
    allocate(nueff(3,2,2), epsalpha(3,2,2), omegae0(3,2))
    allocate(depsalphadpsi(3,2,2), cs(3), dpsids(3), t(3))
    allocate(especies_tem(3,2,2), alphaa1(3,2,2), ispecies_ek(2))
    call kenergymesh()
    zresu = 0.d0
    do j = 1, 3, 2
        ! Floor F_E above underflow: KIA_PASS forms F_psi/F_E, which is 0/0 in a
        ! zero tail (gafusion/MARS-Q#38); the floored tail adds nothing visible.
        zresu(:,j,2,1) = -exp(-min(zepkn, 700.d0))
        zresu(:,j,1,1) = 0.4d0*zresu(:,j,2,1)
    end do
    rlm = 1.d0
    slam0 = -1.d0
    cs = [0.d0, 0.5d0, 1.d0]
    dpsids = 1.d0
    t = 1.d0
    depsalphadpsi = 0.d0
    especies_tem = 1.d0
    alphaa1 = 1.d0
    especies_z(1:2) = [1.d0, -1.d0]
    especies_m(1:2) = [2.d0, 0.00054463d0]
    ispecies_f0(1:2) = 0
    ispecies_ek = 0
    pspecies_fowp(1:2) = 0.d0
    ifowpsi0 = 1
    rntor = -3.d0
    b0k = 1.d0
    rqk = 0.7d0
    omegab = 0.8d0
    omegaci0 = 10.d0
    epsalpha = 1.d0
    call set_test_omega()

    do flag = 0, 1
        kpassdrift = flag
        do mode = 0, 1
            inutype = mode
            do i = 1, size(drifts)
                drift = drifts(i)
                do j = 1, size(rotations)
                    omegae0 = rotations(j)
                    do k = 1, size(widths)
                        if (mode == 1 .and. widths(k) < 1.d-2) cycle
                        nueff = widths(k)
                        call kia_pass(2,1,1,0,1,0.1d0,0,value)
                        write(*,'(I2,1X,I2,1X,5(ES25.16E3,1X))') flag, mode, &
                            drift, rotations(j), widths(k), dble(value), aimag(value)
                    end do
                end do
            end do
        end do
    end do
end program passing_drift
