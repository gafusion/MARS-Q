"""MARS-Q must not read uninitialised heap at the plasma edge (element NRP1 of the 1D profiles).

Builds CHEASE and MARS-Q from source, computes the EXAMPLE equilibrium, and runs MARS-Q on it
with every large heap block pre-filled with 0xFF bytes (a double NaN), emulating heap garbage.
A correct run is independent of the malloc'd contents: every exported field and
profile must be finite. Opt in with MARS_HALF_MESH_SYSTEM_TEST=1 on a cluster.
MARS_EXE reuses a qualified native executable; MARS_TEST_ARTIFACTS retains outputs.
MARS_HALF_MESH_NEGATIVE_EXE optionally reproduces the pre-fix numerical failure
on the same generated equilibrium, without rebuilding that qualified executable.
"""
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(os.environ.get("MARSQ_ROOT") or Path(__file__).resolve().parents[1])
MPIF90 = os.environ.get("MPIF90", "mpif90")
FFLAGS = ("-O2 -fdefault-real-8 -fdefault-double-8 -fopenmp -ffixed-line-length-none "
          "-ffree-line-length-none -fallow-argument-mismatch -std=legacy -w")

SHIM = r"""
#define _GNU_SOURCE
#include <dlfcn.h>
#include <string.h>
#include <stddef.h>
void *malloc(size_t n) {
    static void *(*real)(size_t);
    if (!real) real = dlsym(RTLD_NEXT, "malloc");
    void *p = real(n);
    if (p && n >= 64) memset(p, 0xFF, n);
    return p;
}
"""

RUN_IN = """&BASIC
 NCASE = 2, NITMAX = 1, NSWEEP = 1, GAMMA = 1.66667, P0 = 1e-06, RNTOR = -1,
 M1 = -6, M2 = 6, NPROFN = 4, NPROFIE = 4, NPROFT = 3, NPROFR = 0, ETA = 0.0,
 ROTE = 1.0e-04, NV = 160, NWALL = 0, TAUW = 10000, IWALL = 161, IWO = 1, IVISC = 1,
 PVISC = 1.5, FRACPTH = 1, TALPHA1 = (0.0, 0.0),
&END
&FEEDBACK
 INCFEED = 8, NCOIL = 1, IFEED = 8, FCCHI = 0.0, FWCHI = 0.2, SCCHI = 0.0, SWCHI = 0.05,
 BTCHI = 0.05, ISENS = 8, NSENS = 0, KKF = -2, KREADECA = 0, FEEDI = (1.0e-3,0.0),
&END
&KINETIC
 INCKIN = 0, ALPHAP = 0.5, ALPHAD = 1, OMEGACI0 = 214.7434, IPERTURB = 2, KFASTRUN = 0,
 KENORM = 2, NSPECIES = 2, IFOWPSI0 = 1, ISPECIES_F0 = 0, 0, 3, KNBI = 0, DZETA0 = 1,
 S0TYPE4 = 0.5, R0TYPE4 = 1.5, HHTYPE4C = 0.5, ESPECIES_M = 2, 0.00054463, 2,
 ESPECIES_Z = 1, -1, 1, PSPECIES_AP = 0, 0, 1, PSPECIES_AT = 0, 0, 0, PSPECIES_NP = 0, 0, 0,
 PSPECIES_NTB = 0, 0, 0, PSPECIES_NTD = 0, 0, 0, PSPECIES_FOWP = 0, 0, 1, PSPECIES_FOWT = 0, 0, 0,
 NPROFK = 0, INUTYPE = 0, NPROFUI = 0, NUEFFIA = 0, NPROFUE = 0, NUEFFEA = 0,
 ATAU = (1.0E-04, 0.0), ALTAU = (0.5, 0.0),
&END
&QLIN
 TALPHA2 = 1, TALPHA3 = 1, TALPHA6 = 1, TALPHA8 = 0.9, TDELTAMIN = 0.04, TDELTAMAX = 0.08,
 TDELTALIM = 100, TDELTALOW = 2, CTJXB = 1, CTNTV = 1, CTREY = 1, ITSATURAT = 2, NPROFVD = 3,
 TCHIM0 = 8.4571e-07, NPROFVP = 0, TVPINCH0 = 0, CTEDGE = 0.995, KSOLSAVE = 0, KSOLREAD = 0,
&END
&NUMERIC
 PTRAPI = 1, PTRAPH = 0.9, NRES = 7, NFIT = 3, NTORQ = 5, IPDIVB = 1, NCONVB1 = 1,
 NCONVCS = 0, IOMPNUM = 1, NVACJ = 0, V2XKEY = 2,
&END
&OUTOPT
 JSOUT = 50, ORMIN = 0, ORMAX = 2, OZMIN = -2, OZMAX = 2, NORR = 0, NOZZ = 0,
&END
"""


def build(src, dst, target, flags):
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("*.o", "*.mod", "*.x", "tests"))
    result = subprocess.run(["make", "-j4", "--no-print-directory", "PROG=" + target,
                             "F95=" + MPIF90, "F95FLAGS=" + flags],
                            cwd=dst, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    return dst / target


@unittest.skipUnless(os.environ.get("MARS_HALF_MESH_SYSTEM_TEST") == "1",
                     "opt-in native CHEASE/fluid system regression")
class HalfmeshUninitialised(unittest.TestCase):
    def test_edge_profiles_independent_of_heap_contents(self):
        tmp = tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR"))
        work = Path(tmp.name)
        artifacts = os.environ.get("MARS_TEST_ARTIFACTS")
        if artifacts:
            destination = Path(artifacts)
            def retain():
                destination.mkdir(parents=True, exist_ok=True)
                chease_log = work / "equi" / "log_chease"
                if chease_log.exists():
                    shutil.copy2(chease_log, destination / "log_chease")
                for run in work.glob("run_*"):
                    retained = destination / run.name
                    retained.mkdir(exist_ok=True)
                    for name in ("RUN.IN", "log", "XPLASMA.OUT", "BPLASMA.OUT",
                                 "JPLASMA.OUT", "PROFEQ.OUT"):
                        source = run / name
                        if source.exists():
                            shutil.copy2(source, retained / name)
                tmp.cleanup()
            self.addCleanup(retain)
        else:
            self.addCleanup(tmp.cleanup)
        # The current owning GNU target requires initialized legacy CHEASE locals.
        chease = build(ROOT / "CheaseMerge", work / "chease", "chease.x",
                       "-O2 -mcmodel=medium -finit-local-zero -fdefault-real-8 -fdefault-double-8 "
                       "-ffixed-line-length-none -ffree-line-length-none "
                       "-fallow-argument-mismatch -std=legacy -w")
        executable = os.environ.get("MARS_EXE")
        marsq = Path(executable).resolve() if executable else build(
            ROOT / "MarsQ_2FK", work / "marsq", "marsq.x", FFLAGS)
        self.assertTrue(marsq.is_file(), str(marsq))
        # equilibrium: the EXAMPLE CHEASE namelist, taken from the shipped run script
        equi = work / "equi"
        equi.mkdir()
        shutil.copy(ROOT / "EXAMPLE/Data/EXPEQ_EXAMPLE", equi / "EXPEQ")
        script = (ROOT / "EXAMPLE/Equi/Rchease_EXAMPLE_base").read_text()
        datain = re.search(r"<<'EOF'\n(.*?)\nEOF", script, re.S).group(1)
        (equi / "datain").write_text(datain + "\n")
        with open(equi / "datain") as fin, open(equi / "log_chease", "w") as fout:
            subprocess.run([str(chease)], cwd=equi, stdin=fin, stdout=fout, check=True)
        # MARS-Q run with all heap blocks pre-filled with NaN bytes
        (work / "shim.c").write_text(SHIM)
        subprocess.run([os.environ.get("CC", "gcc"), "-shared", "-fPIC", "-O1", "-o",
                        str(work / "nanmalloc.so"), str(work / "shim.c"), "-ldl"], check=True)
        variants = [("positive", marsq, True)]
        negative = os.environ.get("MARS_HALF_MESH_NEGATIVE_EXE")
        if negative:
            variants.insert(0, ("negative", Path(negative).resolve(), False))
        receipts = {}
        for tag, executable, expected_finite in variants:
            self.assertTrue(executable.is_file(), str(executable))
            run = work / ("run_" + tag)
            run.mkdir()
            for name in ("OUTRMAR", "OUTVMAR"):
                shutil.copy(equi / name, run / name)
            (run / "RUN.IN").write_text(RUN_IN)
            (run / "CURHARMO.IN").write_text("1\n3 1.0 0.0\n")
            prof = "11 1\n" + "".join("%.6f 1.000000\n" % (0.1 * i) for i in range(11))
            for name in ("PROFDEN", "PROFROT", "PROFTE", "PROFTI"):
                (run / (name + ".IN")).write_text(prof)
            env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", LD_PRELOAD=str(work / "nanmalloc.so"))
            env.pop("MALLOC_PERTURB_", None)
            with open(run / "log", "w") as log:
                result = subprocess.run(["bash", "-c", 'ulimit -s unlimited; exec "$1"',
                                         "mars-halfmesh-test", str(executable)],
                                        cwd=run, env=env, stdout=log,
                                        stderr=subprocess.STDOUT, timeout=800)
            checked = {}
            for name in ("BPLASMA.OUT", "XPLASMA.OUT", "JPLASMA.OUT", "PROFEQ.OUT"):
                path = run / name
                text = path.read_text() if path.exists() else ""
                rows = text.splitlines() if name == "PROFEQ.OUT" else text.splitlines()[1:]
                try:
                    values = [float(token.replace("D", "E").replace("d", "e"))
                              for row in rows for token in row.split()]
                except ValueError:
                    values = []
                nonfinite = bool(re.search(r"nan|infinity|(?<![a-z])inf(?![a-z])", text, re.I))
                checked[name] = {"rows": len(rows), "values": len(values),
                                 "nonfinite_text": nonfinite,
                                 "finite": bool(values) and all(map(math.isfinite, values))}
            log_text = (run / "log").read_text()
            finite = result.returncode == 0 and all(c["finite"] for c in checked.values())
            numerical_failure = any(c["nonfinite_text"] for c in checked.values()) or bool(
                re.search(r"ZGETRF ERROR|STOP FMIND|floating.point exception|NaN", log_text, re.I))
            receipts[tag] = {"mars_exe": str(executable), "returncode": result.returncode,
                             "finite": finite, "numerical_failure": numerical_failure,
                             "expected_finite": expected_finite, "checked": checked}
            if artifacts:
                destination.mkdir(parents=True, exist_ok=True)
                (destination / "heap_poison_receipt.json").write_text(json.dumps({
                    "source": str(ROOT), "poison": "malloc >=64 bytes filled 0xFF",
                    "variants": receipts,
                    "oracle": "native linear-response fields/profiles must be finite independently of heap"
                }, indent=2) + "\n")
            if expected_finite:
                self.assertTrue(finite, str(receipts[tag]) + "\n" + log_text[-2000:])
            else:
                self.assertFalse(finite, "negative control did not reproduce the defect")
                self.assertTrue(numerical_failure, "negative failed without a numerical witness: "
                                + log_text[-2000:])


if __name__ == "__main__":
    unittest.main()
