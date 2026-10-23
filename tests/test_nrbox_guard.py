"""NRBOX/NZBOX above NPISO must be rejected before the EQDSK box overruns EQDSPSI."""
import resource
import shutil
import subprocess
import unittest

from support import EXE, ROOT, NativeTest


def unlimited_stack():
    resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))


class NrboxGuard(NativeTest):
    def run_deck(self, nrbox):
        if not EXE.exists():
            self.skipTest("set CHEASE_EXECUTABLE to a Make-built chease.x")
        deck = (ROOT / "EXAMPLE/Equi/Rchease_EXAMPLE_base").read_text()
        deck = deck.split("<<'EOF'\n", 1)[1].split("\nEOF", 1)[0]
        deck = deck.replace(" NEQDSK=0,", " NEQDSK=0, NIDEAL=6, NRBOX=%d, NZBOX=%d," % (nrbox, nrbox), 1)
        shutil.copyfile(ROOT / "EXAMPLE/Data/EXPEQ_EXAMPLE", self.work / "EXPEQ")
        (self.work / "input").write_text(deck + "\n")
        with (self.work / "input").open() as stream:
            return subprocess.run([str(EXE)], stdin=stream, cwd=self.work, capture_output=True,
                                  text=True, timeout=300, preexec_fn=unlimited_stack)

    def test_box_within_npiso_still_runs(self):
        result = self.run_deck(1025)
        self.assertEqual(result.returncode, 0, result.stdout[-1500:])
        self.assertTrue((self.work / "EQDSK.OUT").stat().st_size > 0)

    def test_box_above_npiso_is_rejected(self):
        result = self.run_deck(2049)
        self.assertGreaterEqual(result.returncode, 0, "terminated by signal %d" % -result.returncode)
        self.assertIn("OUT OF RANGE", result.stdout)
        self.assertFalse((self.work / "EQDSK.OUT").exists())


if __name__ == "__main__":
    unittest.main()
