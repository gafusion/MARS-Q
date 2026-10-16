#!/usr/bin/env python3
"""Compile one focused native PAMS/CSC qualification with an existing fortsparse build."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fortsparse-prefix',type=Path,required=True)
    parser.add_argument('--build-parent',type=Path,required=True)
    parser.add_argument('--link-library',action='append',default=[])
    parser.add_argument('--baseline-revision',help='Compare identical caller fixture against this original git revision')
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    source=Path(__file__).resolve().parents[1]
    prefix=args.fortsparse_prefix.resolve()
    args.build_parent.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1',
             FORTSPARSE_UMFPACK_HELPER=str(prefix/'bin/fortsparse_umfpack_helper'),
             MARS_PAMS_UMFPACK_BUDGET_BYTES=str(64*1024**2))
    flags=['-O2','-fopenmp','-fcheck=bounds','-fdefault-real-8',
           '-fdefault-double-8','-ffixed-line-length-none','-ffree-line-length-none',
           '-fallow-argument-mismatch','-ffunction-sections','-fdata-sections']
    with tempfile.TemporaryDirectory(prefix='pams-csc-',dir=args.build_parent) as root:
        env['TMPDIR']=root
        env['MARS_PAMS_ROW_ERRORS']=str(Path(root)/'row-errors.txt')
        def run(command):
            subprocess.run(command,cwd=root,env=env,check=True,timeout=60)
        run(['gfortran',*flags,'-I'+str(prefix/'include'),'-c',
             str(source/'MarsQ_2FK/pams_global_csc.f90'),
             str(source/'MarsQ_2FK/pams.f'),str(source/'MarsQ_2FK/smlpck.f')])
        run(['gfortran',*flags,'-I'+str(prefix/'include'),
             str(source/'tests/pams_global_csc_cases.f90'),
             'pams_global_csc.o','pams.o','smlpck.o',
             str(prefix/'lib/libfortsparse.a'),'-Wl,--gc-sections',
             '-lblas','-ldl','-lrt',*['-l'+x for x in args.link_library],'-o','cases'])
        run([str(Path(root)/'cases')])
        singular=subprocess.run([str(Path(root)/'cases'),'fatal-singular'],
                                cwd=root,env=env,text=True,capture_output=True,timeout=10)
        if singular.returncode==0 or 'PAMS global sparse failure:' not in singular.stdout:
            raise AssertionError('Singular production wrapper must fail with nonzero exit')
        print('PASS singular production wrapper nonzero exit')
        run(['gfortran',*flags,'-I'+str(prefix/'include'),
             str(source/'tests/pams_global_lifecycle.f90'),
             str(source/'tests/pams_manufactured_force.f'),
             str(source/'tests/pams_unreachable_hooks.f90'),
             'pams_global_csc.o','pams.o','smlpck.o',
             str(prefix/'lib/libfortsparse.a'),'-Wl,--gc-sections',
             '-lblas','-ldl','-lrt',*['-l'+x for x in args.link_library],'-o','lifecycle'])
        for backend in ['PAMS','UMFPACK_IPC']:
            env['MARS_PAMS_BACKEND']=backend
            run([str(Path(root)/'lifecycle')])
            if backend=='UMFPACK_IPC':
                rows=[line.split() for line in (Path(root)/'row-errors.txt').read_text().splitlines()]
                if len(rows)!=8 or any(int(row[0])!=i+1 or float(row[1])>1e-11 for i,row in enumerate(rows)):
                    raise AssertionError('Original per-row certificate is incomplete')
        env['MARS_PAMS_UMFPACK_BUDGET_BYTES']='1'
        limited=subprocess.run([str(Path(root)/'lifecycle')],cwd=root,env=env,
                               text=True,capture_output=True,timeout=10)
        if limited.returncode==0 or 'PAMS global sparse failure:           8' not in limited.stdout or \
           'PAMS UMFPACK numeric started: F' not in limited.stdout:
            raise AssertionError('Budget must reject before numerical factorization')
        del env['MARS_PAMS_UMFPACK_BUDGET_BYTES']
        missing=subprocess.run([str(Path(root)/'lifecycle')],cwd=root,env=env,
                               text=True,capture_output=True,timeout=10)
        if missing.returncode==0 or 'Set a positive MARS_PAMS_UMFPACK_BUDGET_BYTES' not in missing.stderr:
            raise AssertionError('Production wrapper requires an explicit positive budget')
        print('PASS complete row certificate, pre-numeric budget and missing-budget failure')
        env['MARS_PAMS_UMFPACK_BUDGET_BYTES']=str(64*1024**2)
        run(['gfortran',*flags,str(source/'tests/pams_global_lifecycle.f90'),
             str(source/'tests/pams_manufactured_force.f'),
             str(source/'tests/pams_unreachable_hooks.f90'),
             str(source/'MarsQ_2FK/pams_global_stub.f90'),
             'pams.o','smlpck.o','-Wl,--gc-sections','-lblas','-o','default-lifecycle'])
        env['MARS_PAMS_BACKEND']='PAMS'
        run([str(Path(root)/'default-lifecycle')])
        env['MARS_PAMS_BACKEND']='UMFPACK_IPC'
        unavailable=subprocess.run([str(Path(root)/'default-lifecycle')],
                                   cwd=root,env=env,text=True,capture_output=True,timeout=10)
        if unavailable.returncode==0 or 'Build with PAMS_CSC=1' not in unavailable.stderr:
            raise AssertionError('Uncompiled backend selection must fail explicitly')
        print('PASS default dependency-free build and unavailable-backend failure')
        records=[]
        variants=[('after','lifecycle')]
        baseline_sha=None
        if args.baseline_revision:
            before=subprocess.check_output(['git','-C',str(source),'show',
                                            args.baseline_revision+':MarsQ_2FK/pams.f'])
            baseline_sha=hashlib.sha256(before).hexdigest()
            (Path(root)/'pams_before.f').write_bytes(before)
            run(['gfortran',*flags,'-c','pams_before.f','-o','pams_before.o'])
            run(['gfortran',*flags,'-I'+str(prefix/'include'),
                 str(source/'tests/pams_global_lifecycle.f90'),
                 str(source/'tests/pams_manufactured_force.f'),
                 str(source/'tests/pams_unreachable_hooks.f90'),
                 'pams_global_csc.o','pams_before.o','smlpck.o',
                 str(prefix/'lib/libfortsparse.a'),'-Wl,--gc-sections',
                 '-lblas','-ldl','-lrt',*['-l'+x for x in args.link_library],'-o','before-lifecycle'])
            variants.insert(0,('before','before-lifecycle'))
        for tag,executable in variants:
            for eps in ['0','1e-12','1e-16']:
                result=subprocess.run([str(Path(root)/executable),eps],cwd=root,env=env,
                                      text=True,capture_output=True,timeout=10)
                completed='PASS_NATIVE_CALLER_FAMILY' in result.stdout
                passed=result.returncode==0 and completed
                if (tag=='after')!=passed:
                    raise AssertionError(f'Unexpected native caller result {tag} eps={eps}: '+result.stdout+result.stderr)
                if tag=='before' and 'STOP FMIND' not in result.stderr and \
                   'Native PAMS family original-system accuracy failure' not in result.stderr:
                    raise AssertionError('Baseline must fail numerically, not through missing adapter features')
                oracle=[float(line.split()[1]) for line in result.stdout.splitlines()
                        if line.startswith('FAMILY_GLOBAL_LU_MAXERR ')]
                if len(oracle)!=1 or oracle[0]>1e-14:
                    raise AssertionError('Independent global LU oracle did not pass')
                errors=[list(map(float,line.split()[1:])) for line in result.stdout.splitlines()
                        if line.startswith('FAMILY_ORIGINAL_ERRORS ')]
                records.append(dict(variant=tag,eps=eps,passed=passed,returncode=result.returncode,
                                    completed=completed,global_lu_maxerr=oracle[0],errors=errors,
                                    stdout=result.stdout,stderr=result.stderr))
                print('PASS native caller before/after gate',tag,eps,errors or 'native FMIND STOP')
        if args.report:
            args.report.write_text(json.dumps(dict(baseline_revision=args.baseline_revision,
                baseline_pams_sha256=baseline_sha,native_caller_family=records),indent=2)+'\n')


if __name__=='__main__':
    main()
