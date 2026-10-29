"""Coil-inductance solve (INCFEED=4) on a grid that reaches the magnetic axis.

FEEDM replaces the half-cell B2/B3 rows by J1=0 and div b=0 for the inductance
problem. On the axis half-cell BOUNDC already sets J1=0, so replacing the B2 row
there leaves B2 and B3 with one equation and PAMS stops on a singular first
block. Runs CHEASE on EXAMPLE/Data/EXPEQ_EXAMPLE at reduced resolution and MARS
on the EXAMPLE/RmpLin vacuum deck (MARS_EXE); skipped without both executables.
"""
import os
from pathlib import Path
import re
import subprocess
import unittest

from support import EXE as CHEASE, NativeTest, ROOT


def setting(text, key, value):
    if re.search(r'(?i)\b' + key + r'\s*=', text):
        return re.sub(r'(?i)\b' + key + r'\s*=\s*[^,\n]+', lambda m: key + '=' + str(value), text)
    return text.replace('&END', f' {key}={value},\n &END')


def heredoc(path):
    return path.read_text().split("<<'EOF'\n", 1)[1].split('\nEOF', 1)[0]


def mars():
    value = os.environ.get('MARS_EXE')
    return Path(value).expanduser().resolve() if value else None


class AxisInductanceRows(NativeTest):
    def setUp(self):
        super().setUp()
        if mars() is None:
            raise unittest.SkipTest('set MARS_EXE for the full-program axis test')
        probe = subprocess.run([str(CHEASE)], input='', cwd=self.work, capture_output=True,
                               text=True, timeout=60) if CHEASE.exists() else None
        if probe is None or probe.returncode == 127:
            raise unittest.SkipTest('set CHEASE_EXECUTABLE to a runnable native build')

    def test_vacuum_inductance_solve_on_axis_grid(self):
        deck = heredoc(ROOT / 'EXAMPLE/Equi/Rchease_EXAMPLE_base')
        for key, value in dict(NS=24, NT=48, NPSI=48, NCHI=64, NV=40, NPLOT=0).items():
            deck = setting(deck, key, value)
        (self.work / 'EXPEQ').write_text((ROOT / 'EXAMPLE/Data/EXPEQ_EXAMPLE').read_text())
        env = {**os.environ, 'OMP_NUM_THREADS': '1'}
        result = subprocess.run([str(CHEASE)], input=deck + '\n', cwd=self.work,
                                capture_output=True, text=True, timeout=300, env=env)
        self.assertEqual(result.returncode, 0, result.stdout[-2000:] + result.stderr)

        run = heredoc(ROOT / 'EXAMPLE/RmpLin/Rmars_EXAMPLE_RmpLin_vac')
        self.assertRegex(run, r'INCFEED\s*=\s*4\b')
        run = re.sub(r'&KINETIC.*?&END', '&KINETIC\n INCKIN = 0,\n&END', run, flags=re.S)
        for key, value in dict(M1=-4, M2=4, NV=40, IWALL=10).items():
            run = setting(run, key, value)
        (self.work / 'RUN.IN').write_text(run + '\n')
        for source, target in (('PROFDEN_EXAMPLE.IN', 'PROFDEN.IN'),
                               ('PROFROT_EXAMPLE.IN', 'PROFROT.IN'),
                               ('PANGLE_EXAMPLE_PEST', 'PANGLE_PEST.IN'),
                               ('PANGLE_EXAMPLE_GEOM', 'PANGLE_EQAC.IN')):
            (self.work / target).write_text((ROOT / 'EXAMPLE/Data' / source).read_text())
        result = subprocess.run([str(mars())], cwd=self.work, capture_output=True,
                                text=True, timeout=600, env=env)
        output = result.stdout + result.stderr
        self.assertFalse('ZGETRF ERROR' in output, output[-1500:])
        self.assertEqual(result.returncode, 0, output[-1500:])
        self.assertIn('AFTER PAMS', output)


if __name__ == '__main__':
    unittest.main()
