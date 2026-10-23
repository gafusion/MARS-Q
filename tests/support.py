"""Helpers for native CHEASE regression tests in MARS-Q (standard library only)."""
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROOT = Path(os.environ.get("MARSQ_ROOT") or Path(__file__).resolve().parents[1])
NATIVE = ROOT / "CheaseMerge"
EXE = Path(os.environ.get("CHEASE_EXECUTABLE") or NATIVE / "chease.x")


def routine(name):
    """Extract one production routine from the fixed-form CHEASE source."""
    source = (NATIVE / "chease.f").read_text()
    match = re.search(r"^ +(?:SUBROUTINE|(?:[A-Z*0-9]+ +)?FUNCTION) +" + name + r"\b",
                      source, re.M | re.I)
    if match is None:
        raise ValueError("Native routine not found: " + name)
    end = source.find("\nC*DECK", match.end())
    return source[match.start():end if end >= 0 else len(source)]


class NativeTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR", ROOT))
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)

    def compile(self, files):
        """Compile test fixtures together with unmodified production routines."""
        paths = []
        for name, content in files.items():
            path = self.work / name
            path.write_text(content)
            if path.suffix.lower() in (".f", ".f90"):
                paths.append(str(path))
        flags = ["-O0", "-fcheck=all", "-ffree-line-length-none", "-ffixed-line-length-none",
                 "-std=legacy", "-fdefault-real-8", "-fdefault-double-8",
                 "-fallow-argument-mismatch", "-mcmodel=medium"]
        result = subprocess.run([os.environ.get("FC", "gfortran"), *flags, *paths, "-o", "oracle"],
                                cwd=self.work, capture_output=True, text=True, timeout=300)
        self.assertEqual(result.returncode, 0, result.stderr)
        return self.work / "oracle"

    def execute(self, exe, *args):
        result = subprocess.run([str(exe), *map(str, args)], cwd=self.work,
                                capture_output=True, text=True, timeout=300)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result
