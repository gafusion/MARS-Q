#!/usr/bin/env python3
"""Run native periodic-spline routines against exact integrated cubics.

The nonuniform oracle integrates the continuous periodic piecewise-linear
second derivative [0, -3, 0, 5], starting with first derivative 15/8.
It uses four distinct knots, without oversized storage or a duplicate endpoint.
"""

import os
import argparse
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest


DRIVER = """program periodic_endpoint
  implicit none
  real :: x(4), y(4), m(4), work(4,3), result(11), q(11)
  real :: a(11), b(11), c(11), d(11), expected(11), err, tol
  real :: epnon0, epsmch, ezmag, rc0p, rc1p, rc1m12, rc1m13
  real :: rc1m14, rc2p, rc2pi
  integer :: index(11), case_id
  common /comdec/ epnon0, epsmch, ezmag, rc0p, rc1p, rc1m12, &
    rc1m13, rc1m14, rc2p, rc2pi
  rc1m14 = 1.e-10
  tol = 5.e-14
  do case_id = 1, 2
    if (case_id == 1) then
      x = [.25, 1.25, 2.25, 3.25]
      y = [0., 1., 0., -1.]
      q = [.75, 1.75, 2.75, 3.75, -.25, 0., 4., x]
      expected = [11./16., 11./16., -11./16., -11./16., &
        -11./16., -47./128., -47./128., y]
    else
      x = [.25, 1.25, 2.75, 3.25]
      y = [0., 11./8., -5./16., -25./24.]
      q = [.75, 2., 3., 3.75, -.25, 0., 4., x]
      expected = [7./8., 61./64., -145./192., -5./6., &
        -5./6., -175./384., -175./384., y]
    end if
    call splchc(x,4,4.,q,11,index,a,b,c,d)
    if (any(index < 1) .or. any(index > 4)) stop 1
    if (any(index(4:7) /= 4)) stop 2
    call splchi(x,y,4,4.,m,work,result,11,index,a,b,c,d)
    err = maxval(abs(result-expected))
    if (err > tol) then
      print *, 'evaluation error', case_id, err
      stop 3
    end if
    if (case_id == 1) then
      err = maxval(abs(m-[0.,-3.,0.,3.]))
    else
      err = maxval(abs(m-[0.,-3.,0.,5.]))
    end if
    if (err > tol) stop 4
    if (abs(a(6)-1./4.) > tol .or. abs(b(6)-3./4.) > tol) stop 5
    if (abs(c(6)+5./128.) > tol .or. abs(d(6)+7./128.) > tol) stop 6
    if (abs(a(4)-.5) > tol .or. abs(b(4)-.5) > tol) stop 7
    if (abs(c(4)+1./16.) > tol .or. abs(d(4)+1./16.) > tol) stop 8
    print *, 'PASS periodic spline exact oracle', case_id, err
  end do
end program
"""


def run_oracle(source_dir=None, legacy_coefficients=False):
    source_dir = source_dir or Path(__file__).resolve().parents[1] / "CheaseMerge"
    source = (source_dir / "chease.f").read_text()
    routines = []
    names = ["SPLCHC", "SPLCHI", "SPLCY", "TRICYC", "ISRCHFGE"]
    if legacy_coefficients:
        # Execute the old native GIJREA block with exactly KCHI-sized arrays.
        start = source.index("         SUBROUTINE GIJREA(KPSI,PS)")
        start = source.index("         DO I=1,INCHI\n           ICHISO", start)
        end = source.index("         ENDDO", start) + len("         ENDDO")
        routines.append("""      SUBROUTINE SPLCHC(CHIISO,IGCHISO,PERIOD,ZCHI,INCHI,
     + ICHIISO,ZA,ZB,ZC,ZD)
      DIMENSION CHIISO(IGCHISO),ZCHI(INCHI),ICHIISO(INCHI),
     + ZA(INCHI),ZB(INCHI),ZC(INCHI),ZD(INCHI)
""" + source[start:end] + "\n      END\n")
        names.remove("SPLCHC")
    for name in names:
        start = re.search(
            rf"^ +(?:SUBROUTINE|FUNCTION) {name}\(", source, re.MULTILINE
        )
        if start is None:
            raise RuntimeError(f"native routine {name} missing")
        end = re.search(r"^ +END *$", source[start.start():], re.MULTILINE)
        routines.append(source[start.start():start.start() + end.end()] + "\n")
    with tempfile.TemporaryDirectory(prefix="chease-periodic-") as directory:
        path = Path(directory)
        (path / "native.f").write_text("\n".join(routines))
        (path / "oracle.f90").write_text(DRIVER)
        compiler = os.environ.get("FC", "gfortran")
        subprocess.run(
            [compiler, "-O0", "-g", "-fcheck=all", "-fbacktrace",
             "-fdefault-real-8", "-fdefault-double-8", "-std=legacy",
             "-I" + str(source_dir), str(path / "native.f"),
             str(path / "oracle.f90"), "-o", str(path / "oracle")],
            check=True,
        )
        subprocess.run([str(path / "oracle")], check=True)


class TestCheasePeriodicEndpoint(unittest.TestCase):
    def test_native_exact_periodic_splines(self):
        run_oracle()


if __name__ == "__main__":
    if "--source-dir" in sys.argv or "--legacy-coefficients" in sys.argv:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--source-dir", type=Path)
        parser.add_argument("--legacy-coefficients", action="store_true")
        args = parser.parse_args()
        run_oracle(args.source_dir, args.legacy_coefficients)
    else:
        unittest.main()
