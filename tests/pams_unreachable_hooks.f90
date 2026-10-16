subroutine get_inormsol(inorm)
implicit none
integer,intent(out)::inorm
inorm=0
end subroutine
subroutine get_adapcrit()
error stop 'Unexpected unrelated hook: get_adapcrit'
end subroutine
subroutine torqjxb()
error stop 'Unexpected unrelated hook: torqjxb'
end subroutine
subroutine torqntv()
error stop 'Unexpected unrelated hook: torqntv'
end subroutine
subroutine torqrey()
error stop 'Unexpected unrelated hook: torqrey'
end subroutine
subroutine torqergo()
error stop 'Unexpected unrelated hook: torqergo'
end subroutine
subroutine calcdens()
error stop 'Unexpected unrelated hook: calcdens'
end subroutine
subroutine calcdispnorm()
error stop 'Unexpected unrelated hook: calcdispnorm'
end subroutine
subroutine feedevol()
error stop 'Unexpected unrelated hook: feedevol'
end subroutine
subroutine solvemoment()
error stop 'Unexpected unrelated hook: solvemoment'
end subroutine
subroutine solvedntran()
error stop 'Unexpected unrelated hook: solvedntran'
end subroutine
