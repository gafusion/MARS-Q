"""All four smoothed Hermite components follow the node permutation."""
from pathlib import Path
from support import NativeTest, NATIVE, routine


class Smoothing(NativeTest):
    def test_four_jet_components(self):
        files = {"oracle.f90": Path(__file__).with_name("smoothing.f90").read_text(),
                 "wrapper.f": Path(__file__).with_name("smoothing_wrapper.f").read_text()}
        for name in ("COMDIM", "DECLAR"):
            files[name + ".inc"] = (NATIVE / (name + ".inc")).read_text()
        for name in ("COMCON", "COMINT", "COMPHY"):
            files[name + ".inc"] = "C Unused by these routines\n"
        for name, block in {"COMBLA": "COMMON /FIXC/ CPSI(NP4NST),CPSICL(NP4NST)",
                            "COMESH": "COMMON /FIXM/ CSIG(NSP1),CT(NTP1)",
                            "COMNUM": "COMMON /FIXN/ NS,NT,NS1,NT1",
                            "COMSOL": "COMMON /FIXP/ NUPDWN(NSNT)"}.items():
            files[name + ".inc"] = "         " + block + "\n"
        files["native.f"] = "\n".join(routine(name) for name in ("SMOOTH", "MSPLINE", "MSPLCY", "TRIDAGM", "TRICYCM", "ISAMIN")).replace("SUBROUTINE SMOOTH", "SUBROUTINE NATIVE_SMOOTH")
        self.execute(self.compile(files))
