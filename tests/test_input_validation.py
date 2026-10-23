"""Input guards exercise the native rule and interpolation implementations."""
from pathlib import Path
import re
import subprocess
from support import NativeTest, NATIVE, routine


class InputValidation(NativeTest):
    def test_quadrature_capacity_counts_and_axis(self):
        init = routine("INITIA")
        prefix = init.split("CALL RESETI(MPLA1", 1)[0]
        match = re.search(r"^.*(?:CALL CHEASE_GAUSS_CHECK|IF *\(NSGAUS *\.EQ\. *1\))", prefix, re.M)
        guard = prefix[match.start():] if match else ""
        start = re.search(r"^.*CALL GAUSS\(NSGAUS", init, re.M).start()
        block = re.split(r"(?m)^ *RETURN", init[start:], maxsplit=1)[0]
        declarations = "      INCLUDE 'DECLAR.inc'\n      INCLUDE 'COMDIM.inc'\n"
        assignments = ""
        program = """      program tensor_rule
DECLARATIONS
      integer :: nsgaus,ntgaus,nwgaus,j6,j7,j8,j9,kpoint,k,ierr
      real(kind(1.d0)),allocatable :: zracs(:),zract(:),zwgts(:),zwgtt(:),cw(:),dzeta(:,:)
      real(kind(1.d0)) :: value,expected
      integer,parameter :: rkind=kind(1.d0)
      character(32) :: arg
      call get_command_argument(1,arg)
      read(arg,*) nsgaus
      call get_command_argument(2,arg)
      read(arg,*) ntgaus
      nwgaus=nsgaus*ntgaus
ASSIGNMENTS
      allocate(zracs(npsgs+1),zwgts(npsgs+1),zract(nptgs+1),zwgtt(nptgs+1))
      allocate(cw(2*npgaus),dzeta(2*npgaus,2))
      zracs=-17
      zract=-17
      zwgts=-17
      zwgtt=-17
      cw=-17
      dzeta=-17
GUARD
BLOCK
      if (any(cw(npgaus+1:)/=-17)) stop 5
      if (any(zracs(npsgs+1:)/=-17).or.any(zract(nptgs+1:)/=-17)) stop 5
      if (nwgaus/=nsgaus*ntgaus) stop 6
      do k=0,1
         value=sum(cw(:nwgaus)*dzeta(:nwgaus,1)**k*dzeta(:nwgaus,2)**k)
         expected=1.d0/(k+1)**2
         if (.not.(abs(value-expected)<2.d-13)) stop 7
      enddo
      print *, 'PASS: capacity, effective count, exact tensor moments'
      end
""".replace("DECLARATIONS", declarations).replace("ASSIGNMENTS", assignments).replace("GUARD", guard).replace("BLOCK", block)
        files = {}
        for name in ("DECLAR", "COMDIM"):
            files[name+".inc"] = (NATIVE/(name+".inc")).read_text()
        files["gauss.f"] = routine("GAUSS")
        if "CHEASE_GAUSS_CHECK" in guard:
            files["check.f"] = routine("CHEASE_GAUSS_CHECK")
        files["oracle.f"] = program
        exe = self.compile(files)
        for radial, angular in ((2, 1), (4, 4)):
            with self.subTest(radial=radial, angular=angular):
                self.execute(exe, radial, angular)
        rejected = [(1, 4), (8, 4), (4, 8), (0, 4), (11, 4)]
        for radial, angular in rejected:
            with self.subTest(rejected=(radial, angular)):
                result = subprocess.run([str(exe), str(radial), str(angular)], cwd=self.work,
                                        capture_output=True, text=True, timeout=300)
                self.assertEqual(result.returncode, 1, result.stdout+result.stderr)
                self.assertRegex(result.stderr, "undefined at the magnetic axis|Invalid Gaussian order/capacity")
