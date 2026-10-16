subroutine pams_global_requested(requested)
    implicit none
    logical, intent(out)::requested
    character(64)::backend
    integer::stat
    call get_environment_variable('MARS_PAMS_BACKEND', backend, status=stat)
    if (stat == 1) backend = 'PAMS'
    if (stat /= 0 .and. stat /= 1) error stop 'Invalid MARS_PAMS_BACKEND value'
    if (trim(backend) /= '' .and. trim(backend) /= 'PAMS') then
        error stop 'Build with PAMS_CSC=1 to enable MARS_PAMS_BACKEND=UMFPACK_IPC'
    end if
    requested = .false.
end subroutine

subroutine pams_global_solve()
    implicit none
    error stop 'Global PAMS backend not compiled'
end subroutine
