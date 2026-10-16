module pams_global_csc
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortsparse_csc, only: csc_z_t
    use fortsparse_umfpack_ipc, only: umfpack_ipc_backend_t
    use fortsparse_status, only: fortsparse_status_t
    use fortsparse_backend, only: sparse_factor_stats_t
    implicit none
    private
    integer, parameter :: wide = selected_real_kind(18, 400)
    real(dp), parameter, public :: pams_csc_tolerance = 1.0e-11_dp
    public :: pams_csc_pack, pams_csc_solve, pams_csc_errors
contains
    subroutine pams_csc_pack(md, my, nd, a, b, c, d, e, f, g, h, matrix, ierr, report_unit)
        integer, intent(in) :: md, my, nd
        complex(dp), intent(in) :: a(md, md, *), b(md, md, *), c(md, md, *)
        complex(dp), intent(in) :: d(my, my, *), e(md, my, *), f(my, md, *)
        complex(dp), intent(in) :: g(my, md, *), h(md, my, *)
        type(csc_z_t), intent(out) :: matrix
        integer, intent(out) :: ierr
        integer, intent(in), optional :: report_unit
        integer(int64), allocatable :: counts(:)
        integer(int64) :: size64, total
        integer :: n, col, k, j, stat
        integer(int64) :: p
        logical :: fill
        ierr = 0
        if (md < 1 .or. my < 1 .or. nd < 0) then
            ierr = 2
            return
        end if
        size64 = int(nd, int64)*(int(md, int64) + my) + md
        if (size64 > huge(n) - 1_int64 .or. &
            int(md, int64) + my > huge(n)) then
            if (present(report_unit)) write(report_unit, *) &
                'PAMS counted original dimension exceeds CSC32:', size64
            ierr = 7
            return
        end if
        n = int(size64)
        allocate (counts(n), matrix%col_ptr(n + 1), stat=stat)
        if (stat /= 0) then
            ierr = 6
            return
        end if
        fill = .false.
        !$omp parallel do private(k,j,p)
        do col = 1, n
            k = (col - 1)/(md + my) + 1
            j = mod(col - 1, md + my) + 1
            p = 0
            call column(k, j, p)
            counts(col) = p
        end do
        !$omp end parallel do
        total = sum(counts)
        if (present(report_unit)) write(report_unit, *) &
            'PAMS counted original matrix n,nnz:', size64, total
        if (total > int(huge(n), int64) - 1_int64) then
            ierr = 7
            return
        end if
        total = 1_int64
        do col = 1, n
            matrix%col_ptr(col) = int(total)
            total = total + counts(col)
        end do
        matrix%col_ptr(n + 1) = int(total)
        matrix%nrow = n
        matrix%ncol = n
        matrix%nnz = int(total - 1)
        allocate (matrix%row_idx(matrix%nnz), matrix%val(matrix%nnz), stat=stat)
        if (stat /= 0) then
            ierr = 6
            return
        end if
        fill = .true.
        !$omp parallel do private(k,j,p)
        do col = 1, n
            k = (col - 1)/(md + my) + 1
            j = mod(col - 1, md + my) + 1
            p = matrix%col_ptr(col) - 1
            call column(k, j, p)
        end do
        !$omp end parallel do
    contains
        subroutine column(k, j, p)
            integer, intent(in) :: k, j
            integer(int64), intent(inout) :: p
            if (j <= md) then
                if (k > 1) then
                    call block(c(:, j, k - 1), (k - 2)*(md + my), p)
                    call block(g(:, j, k - 1), (k - 2)*(md + my) + md, p)
                end if
                call block(b(:, j, k), (k - 1)*(md + my), p)
                if (k <= nd) then
                    call block(f(:, j, k), (k - 1)*(md + my) + md, p)
                    call block(a(:, j, k), k*(md + my), p)
                end if
            else
                call block(e(:, j - md, k), (k - 1)*(md + my), p)
                call block(d(:, j - md, k), (k - 1)*(md + my) + md, p)
                call block(h(:, j - md, k), k*(md + my), p)
            end if
        end subroutine
        subroutine block(values, offset, p)
            complex(dp), intent(in) :: values(:)
            integer, intent(in) :: offset
            integer(int64), intent(inout) :: p
            integer :: row
            do row = 1, size(values)
                ! No numerical threshold: retain every nonzero coefficient.
                if (values(row) == cmplx(0.0_dp, 0.0_dp, dp)) cycle
                p = p + 1
                if (fill) then
                    matrix%row_idx(p) = offset + row
                    matrix%val(p) = values(row)
                end if
            end do
        end subroutine
    end subroutine pams_csc_pack

    subroutine pams_csc_errors(md, my, nd, a, b, c, d, e, f, g, h, bx, by, x, y, errors)
        integer, intent(in) :: md, my, nd
        complex(dp), intent(in) :: a(md, md, *), b(md, md, *), c(md, md, *)
        complex(dp), intent(in) :: d(my, my, *), e(md, my, *), f(my, md, *)
        complex(dp), intent(in) :: g(my, md, *), h(md, my, *)
        complex(dp), intent(in) :: bx(md, nd + 1), by(my, nd), x(md, nd + 1), y(my, nd)
        real(dp), intent(out) :: errors(:)
        complex(wide) :: value, rhs
        real(wide) :: denominator
        integer :: k, row, pos
        !$omp parallel do private(row,pos,value,rhs,denominator)
        do k = 1, nd + 1
            do row = 1, md
                rhs = bx(row, k)
                value = 0
                denominator = abs(rhs)
                if (k > 1) then
                    call accumulate(a(row, :, k - 1), x(:, k - 1), value, denominator)
                    call accumulate(h(row, :, k - 1), y(:, k - 1), value, denominator)
                end if
                call accumulate(b(row, :, k), x(:, k), value, denominator)
                if (k <= nd) then
                    call accumulate(c(row, :, k), x(:, k + 1), value, denominator)
                    call accumulate(e(row, :, k), y(:, k), value, denominator)
                end if
                pos = (k - 1)*(md + my) + row
                call finish(value, rhs, denominator, pos)
            end do
            if (k > nd) cycle
            do row = 1, my
                rhs = by(row, k)
                value = 0
                denominator = abs(rhs)
                call accumulate(f(row, :, k), x(:, k), value, denominator)
                call accumulate(g(row, :, k), x(:, k + 1), value, denominator)
                call accumulate(d(row, :, k), y(:, k), value, denominator)
                pos = (k - 1)*(md + my) + md + row
                call finish(value, rhs, denominator, pos)
            end do
        end do
        !$omp end parallel do
    contains
        subroutine accumulate(coefficients, vector, value, denominator)
            complex(dp), intent(in) :: coefficients(:), vector(:)
            complex(wide), intent(inout) :: value
            real(wide), intent(inout) :: denominator
            complex(wide) :: coefficient, element
            integer :: j
            do j = 1, size(vector)
                coefficient = coefficients(j)
                element = vector(j)
                value = value + coefficient*element
                denominator = denominator + abs(coefficient)*abs(element)
            end do
        end subroutine
        subroutine finish(value, rhs, denominator, pos)
            complex(wide), intent(in) :: value, rhs
            real(wide), intent(in) :: denominator
            integer, intent(in) :: pos
            if (denominator == 0) then
                if (abs(rhs - value) == 0) then
                    errors(pos) = 0
                else
                    errors(pos) = huge(1.0_dp)
                end if
            else
                errors(pos) = real(abs(rhs - value)/denominator, dp)
            end if
        end subroutine
    end subroutine pams_csc_errors

    subroutine pams_csc_solve(md,my,nd,a,b,c,d,e,f,g,h,bx,by,x,y,errors,ierr, &
                              memory_budget_bytes, report_unit)
        integer, intent(in) :: md, my, nd
        complex(dp), intent(in) :: a(md, md, *), b(md, md, *), c(md, md, *)
        complex(dp), intent(in) :: d(my, my, *), e(md, my, *), f(my, md, *)
        complex(dp), intent(in) :: g(my, md, *), h(md, my, *)
        complex(dp), intent(in) :: bx(md, nd + 1), by(my, nd)
        complex(dp), intent(out) :: x(md, nd + 1), y(my, nd)
        real(dp), intent(out) :: errors(:)
        integer, intent(out) :: ierr
        integer(int64), intent(in), optional :: memory_budget_bytes
        integer, intent(in), optional :: report_unit
        type(csc_z_t) :: matrix
        type(umfpack_ipc_backend_t) :: solver
        type(fortsparse_status_t) :: status, stats_status
        type(sparse_factor_stats_t) :: stats
        complex(dp), allocatable :: rhs(:), solution(:)
        integer :: k, offset, n, stat
        call pams_csc_pack(md, my, nd, a, b, c, d, e, f, g, h, matrix, ierr, report_unit)
        if (ierr /= 0) return
        if (any(.not. ieee_is_finite(real(matrix%val, dp))) .or. &
            any(.not. ieee_is_finite(aimag(matrix%val)))) then
            ierr = 2
            return
        end if
        n = matrix%nrow
        if (size(errors) /= n) then
            ierr = 2
            return
        end if
        allocate (rhs(n), solution(n), stat=stat)
        if (stat /= 0) then
            ierr = 6
            return
        end if
        do k = 1, nd + 1
            offset = (k - 1)*(md + my)
            rhs(offset + 1:offset + md) = bx(:, k)
            if (k <= nd) rhs(offset + md + 1:offset + md + my) = by(:, k)
        end do
        if (any(.not. ieee_is_finite(real(rhs, dp))) .or. &
            any(.not. ieee_is_finite(aimag(rhs)))) then
            ierr = 2
            return
        end if
        solver%options%pivot_tolerance = 1.0_dp
        if (present(memory_budget_bytes)) then
            if (memory_budget_bytes < 0_int64) then
                ierr = 2
                return
            end if
            solver%options%memory_budget_bytes = memory_budget_bytes
        end if
        if (present(report_unit)) then
            write (report_unit, *) 'PAMS UMFPACK requested budget bytes:', &
                solver%options%memory_budget_bytes
            write (report_unit, *) 'PAMS global matrix n,nnz:', matrix%nrow, matrix%nnz
            write (report_unit, *) 'PAMS caller CSC bytes:', &
                4_int64*(int(matrix%ncol, int64) + 1) + 20_int64*matrix%nnz
        end if
        call solver%factor_complex(matrix, .true., status)
        call solver%factor_stats(stats, stats_status)
        if (present(report_unit)) then
            write (report_unit, *) 'PAMS UMFPACK status:', status%code, trim(status%msg)
            write (report_unit, *) 'PAMS UMFPACK statistics status:', stats_status%code
          write (report_unit, *) 'PAMS UMFPACK symbolic,estimated,numeric,IPC bytes:', &
                stats%symbolic_peak_bytes, stats%estimated_peak_bytes, &
                stats%numeric_peak_bytes, stats%ipc_bytes
            write (report_unit, *) 'PAMS UMFPACK effective pivot tolerance:', &
                stats%effective_pivot_tolerance
           write (report_unit, *) 'PAMS UMFPACK numeric started:', stats%numeric_started
        end if
        if (status%code == 0) call solver%solve_complex(rhs, solution, status)
        ierr = status%code
        call solver%free()
        if (ierr /= 0) return
        if (any(.not. ieee_is_finite(real(solution, dp))) .or. &
            any(.not. ieee_is_finite(aimag(solution)))) then
            ierr = 5
            return
        end if
        do k = 1, nd + 1
            offset = (k - 1)*(md + my)
            x(:, k) = solution(offset + 1:offset + md)
            if (k <= nd) y(:, k) = solution(offset + md + 1:offset + md + my)
        end do
        call pams_csc_errors(md, my, nd, a, b, c, d, e, f, g, h, bx, by, x, y, errors)
        if (any(.not. ieee_is_finite(errors))) then
            ierr = 5
        else if (maxval(errors) > pams_csc_tolerance) then
            ierr = 100
        end if
    end subroutine pams_csc_solve
end module pams_global_csc

subroutine pams_global_requested(requested)
    implicit none
    logical, intent(out) :: requested
    character(64) :: backend
    integer :: stat
    call get_environment_variable('MARS_PAMS_BACKEND', backend, status=stat)
    if (stat == 1) backend = 'PAMS'
    if (stat /= 0 .and. stat /= 1) error stop 'Invalid MARS_PAMS_BACKEND value'
    select case (trim(backend))
    case ('', 'PAMS')
        requested = .false.
    case ('UMFPACK_IPC')
        requested = .true.
    case default
        error stop 'Unknown MARS_PAMS_BACKEND'
    end select
end subroutine

subroutine pams_global_solve(md, my, nd, a, b, c, d, e, f, g, h, bx, by, x, y)
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, output_unit
    use pams_global_csc, only: pams_csc_solve
    implicit none
    integer, intent(in)::md, my, nd
    complex(dp), intent(in)::a(md, md, *), b(md, md, *), c(md, md, *), d(my, my, *)
    complex(dp), intent(in)::e(md, my, *), f(my, md, *), g(my, md, *), h(md, my, *)
    complex(dp), intent(in)::bx(md, nd + 1), by(my, nd)
    complex(dp), intent(out)::x(md, nd + 1), y(my, nd)
    real(dp), allocatable::errors(:)
    integer::ierr, unit, row, stat
    integer(int64)::size64, budget_bytes
    character(4096)::path
    character(64)::budget_text
    if (md < 1 .or. my < 1 .or. nd < 0) error stop 'Invalid PAMS CSC dimensions'
    size64 = int(nd, int64)*(int(md, int64) + my) + md
  if (size64 < 1 .or. size64 > huge(row) - 1_int64) error stop 'PAMS CSC index overflow'
    allocate (errors(int(size64)), stat=stat)
    if (stat /= 0) error stop 'PAMS CSC error-vector allocation failed'
 call get_environment_variable('MARS_PAMS_UMFPACK_BUDGET_BYTES',budget_text,status=stat)
    if (stat /= 0) &
        error stop 'Set a positive MARS_PAMS_UMFPACK_BUDGET_BYTES ' // &
        'before sparse solve'
    read (budget_text, *, iostat=stat) budget_bytes
    if (stat /= 0) error stop 'Invalid MARS_PAMS_UMFPACK_BUDGET_BYTES'
    if (budget_bytes <= 0) error stop 'UMFPACK factor memory budget must be positive'
   call pams_csc_solve(md, my, nd, a, b, c, d, e, f, g, h, bx, by, x, y, errors, ierr, &
                        budget_bytes, output_unit)
    ! Only success and the finite row-quality failure own a valid certificate.
    if (ierr == 0 .or. ierr == 100) then
        write (*, *) 'PAMS original componentwise backward error,row:', &
            maxval(errors), maxloc(errors)
        call get_environment_variable('MARS_PAMS_ROW_ERRORS', path, status=stat)
        if (stat == 0 .and. len_trim(path) > 0) then
            open (newunit=unit, file=trim(path), status='replace', action='write')
            do row = 1, size(errors)
                write (unit, '(i0,1x,es25.17)') row, errors(row)
            end do
            close (unit)
        end if
    end if
    if (ierr /= 0) then
        write (*, *) 'PAMS global sparse failure:', ierr
        flush (output_unit)
        error stop 1
    end if
end subroutine
