"""Native MAPPIN axis-F extrapolation on nonuniform radial coordinates."""
from support import NativeTest, NATIVE, routine


class AxisField(NativeTest):
    def test_signed_polynomial_axis_and_field_units(self):
        source = routine("MAPPIN")
        start = source.index("T0    = FCCCC0(TMF(1)")
        assignment = "         " + "\n".join(source[start:].splitlines()[:2])
        native = """      DOUBLE PRECISION FUNCTION AXIS_F(TMF,CSM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TMF(4),CSM(4)
      PARAMETER (RC0P=0D0)
      INCLUDE 'CUCCCC.inc'
""" + assignment + "\n      AXIS_F=T0\n      END\n"
        oracle = """program polynomial_axis
implicit none
real(8),external :: axis_f
real(8) :: x(4),y(4),expected,actual,scale
integer :: grid,profile,units,failures
failures=0
 do grid=1,2
  if(grid==1)then
   x=[.11d0,.24d0,.39d0,.62d0]
  else
   x=[.0078273d0,.023482d0,.039137d0,.054791d0]
  endif
  do units=1,2
   scale=real(units*2-1,8)
   do profile=1,4
    expected=scale
    select case(profile)
    case(1)
     y=scale
    case(2)
     y=scale*(1+.3d0*x**2)
    case(3)
     y=scale*(1+.3d0*x**2-.1d0*x**3)
    case(4)
     y=-scale*(1+.3d0*x**2)
     expected=-scale
    end select
    actual=axis_f(y,x)
    print *,grid,units,profile,expected,actual
    if(abs(actual-expected)>2.d-13*scale)failures=failures+1
   enddo
  enddo
 enddo
 if(failures/=0)error stop 1
end program
"""
        exe = self.compile({"CUCCCC.inc": (NATIVE / "CUCCCC.inc").read_text(),
                            "axis.f": native, "oracle.f90": oracle})
        self.execute(exe)
