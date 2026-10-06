"""One explicit D13 build; only constant native metadata is invoked.

Run under an outer 60-second supervisor. Intermediate receipts and compiler
streams survive rejection; an existing build directory is never overwritten.
"""
import ctypes as C
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import sysconfig
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
QCB = ROOT.parent
REFERENCE = QCB/'specialization_optimization_20260930'
CONFIG = REFERENCE/'build/flint_header_config'
SOURCE = REFERENCE/'sources/flint-3.6.0/src'
SCIENTIFIC = ('block.c', 'block_interop.py', 'build.py', 'block_acquisition.py', 'diagnose.py')
DOCUMENTS = ('README.md', 'PROTOCOL.md', 'API_CONTRACT.md', 'REVIEW.md')
LIBRARY_SHA = '871a4132fd1e9f3638391b2208e07088f8e3e72a10e41d45f58b150a60c2a1a9'
HEADERS_SHA = 'a94c463fafd1060ed192935311949a54cd341ca53abe3e0edc19f7aa2528fcfb'
PREMISES = {
    ROOT/'BLOCK_ELIMINATION_DESIGN_DRAFT.md': '8b07f2f1a2e7db77dcdc316d209de7d7d57b9c0b5cf922c716b050bdd57c8123',
    ROOT/'BLOCK_ELIMINATION_DESIGN_REVIEW.md': '163e7e6506b102a94ecf77d797cba824c6829ca886561c2461196ed621be5609',
    ROOT/'BLOCK_ELIMINATION_API_NOTE.md': '01088dc06238a2b907aa032ce167cee23950bb0e230f994e82f58ae8d042d513',
    ROOT/'BLOCK_ELIMINATION_MATRIX_API_CONTRACT.md': 'eb5193db7a3f659b7411406d218b71def6f5a95e156e91b403bef389b818a2cc',
    ROOT/'BLOCK_ELIMINATION_SCHEDULE_DRAFT.md': 'c0a206c050b37cbf9352bc3934b6caf07b180543a6411496161c7b91644937bc',
    ROOT/'BLOCK_ELIMINATION_SCHEDULE_REVIEW.md': 'b1b388fbee0d8429c7c7c1cbd953afd1c133d9a60a79a443d50a7a3956a2c057',
}
D11 = HERE.parent/'terminal_readout'
D12 = HERE.parent/'native_second'
EVIDENCE = {
    D12/'AUDIT.md': 'e1a90bee49168a8d69ec46a8bd5431101f6f9de95501de5d661dce83b908fdf1',
    D12/'RESULTS_ZH.md': 'e9a3cd1efd0e5b176f4b3cf010614f2f0dda1e2f24dfa33bd7eee39182005f37',
    D12/'audit/d12_v1_receipt.json': '6611172b8c727e5cc3256c98d535369f2b84e7f4ce471939bc09d0441b72b236',
    D12/'runs/d12_v1/manifest.json': '821ec23a9493e0be21fb61e905895864448327ddcad1bc51a44eefc71fd7da42',
}
PRIMARY = {
    'fmpz_mod_mat.h': '9dda40e3470de26a23d07ca7ed5b9cda40332ed7ada3a254e3ee2d76b3d329a2',
    'fmpz_mod_mat/init.c': 'fbf49afa9cc8bbe22b4b088403b3c0df1ac5fa15ef4a94bc121ae2a079f5ba52',
    'fmpz_mod_mat/clear.c': 'a48a22a8e5afdacf7547da878d0309737e0ae13ecdd8012dde36f12de6490432',
    'fmpz_mod_mat/set_get.c': '639ec8b5441ca785e455e78cff45eed5e5ec5fa2e73af9b3cd7c892c76c9c31c',
    'fmpz_mod_mat/mul.c': 'c7deab8b82a2aaa514a9f1c84c4c1b33e8d1f54fbdf047dd865b42e5fd042300',
    'fmpz_mat/init.c': 'a056ee15924fc4849924fef4ded7d85772077a51c3dcd1fed11f6aaac2009bcb',
    'fmpz_mod.h': '8a940231631182b3522052540b5214ec45db58e0be3d426d10230e2ee71b8f61',
    'fmpz_mod/ctx.c': '73a8f5c7090304e997a0059f5b60fd249247d67860e78f1af2eb68f729d4670a',
    'fmpz_mod/inv.c': '7dd703caeeb623b5d99b62b67bb01cc9d9157aa6b48bb09b74d1e4dada9a1fb3',
    'fmpz_mod/is_invertible.c': '53d67fe983efda92e00412f5fb44afa1a7d347963d25376cf74175d74463a0de',
    'fmpz_mod/set_fmpz.c': 'cb03a6095d221668bcb4d79d6a544ef9938b5eb0802470fb1fe54f608780a3a5',
    'fmpz_mod_poly/set_get_coeff.c': 'f4deef9602378df8110638485195bf51c6f8febdb49b3f54017283afb6f19424',
    'fmpz_mod/sub.c': 'ec2fd39807e29e8f2b26468d75d0ef5b2fef7f1e0ac0bdf2ab88dce2ae00d110',
    'fmpz_mod/mul.c': '45c248560d5480ed7023ffc3a5942dd69ac3a286765504260d69afc013bf55e1',
}
WORDS = dict(char_bits=8, flint_bits=64, long_bytes=8, slong_bytes=8, ulong_bytes=8, size_t_bytes=8)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = perf_counter()
    resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
    build = HERE/'build'
    if build.exists():
        raise RuntimeError('refusing to overwrite the registered D13 build')
    build.mkdir()
    receipt_path = build/'receipt.json'
    receipt = dict(status='STARTED', passed=False, mathematical_operation_executed=False,
                   complete_build_budget_seconds=60, memory_limit_bytes=8*1024**3)

    def persist():
        receipt['build_to_receipt_seconds'] = perf_counter()-started
        receipt['clock_scope'] = 'Entry through the sample before this receipt write; final write/print/exit is covered by the outer supervisor.'
        receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')

    def remaining():
        value = 60-(perf_counter()-started)
        if value <= 0:
            raise TimeoutError('D13 complete build deadline exceeded')
        return value

    persist()
    try:
        gate = (HERE/'REVIEW.md').read_text()
        review = (HERE/'CODE_REVIEW.md').read_text()
        if 'code_status: CODE_PERMITTED' not in gate.splitlines():
            raise RuntimeError('D13 registered code gate absent')
        if 'build_status: BUILD_PERMITTED' not in review.splitlines():
            raise RuntimeError('D13 independent static code review does not permit build')
        scientific = {name: digest(HERE/name) for name in SCIENTIFIC}
        documents = {name: digest(HERE/name) for name in DOCUMENTS}
        if any(value not in review for value in (*scientific.values(), *documents.values())):
            raise RuntimeError('D13 source review does not bind every source and contract')
        if any(documents[name] not in gate for name in DOCUMENTS if name != 'REVIEW.md'):
            raise RuntimeError('D13 prerequisite review does not bind the formal contracts')
        for path, expected in (PREMISES | EVIDENCE).items():
            if digest(path) != expected or expected not in gate:
                raise RuntimeError('D13 changed or unreviewed premise: '+str(path))
        inherited = json.loads((D12/'runs/d12_v1/manifest.json').read_text())
        audit = json.loads((D12/'audit/d12_v1_receipt.json').read_text())
        if (audit['status'] != 'PASS_SCOPED_READ_ONLY_ARTIFACT_AUDIT' or audit['workers'] != 96
                or audit['outputs'] != 384 or audit['errors']):
            raise RuntimeError('D12 evidence is incomplete or failed')
        for path, expected in inherited['sources'].items():
            if digest(QCB/path) != expected:
                raise RuntimeError('changed inherited science: '+path)
            remaining()
        for path, expected in inherited['dependencies'].items():
            if digest(Path(path)) != expected:
                raise RuntimeError('changed inherited installed dependency: '+path)
            remaining()
        for name, expected in PRIMARY.items():
            if digest(SOURCE/name) != expected:
                raise RuntimeError('changed D13 primary source: '+name)
        repo = Path(subprocess.check_output(['git', 'rev-parse', '--show-toplevel'],
                    cwd=HERE, text=True, timeout=remaining()).strip())
        bound = [HERE/name for name in (*SCIENTIFIC, *DOCUMENTS, 'CODE_REVIEW.md')]
        bound += list(PREMISES)+list(EVIDENCE)
        # D12's compact audit records counts, not a new tracked-path list.
        # Reuse the frozen D11 disposition and add every D12 owned source.
        old_audit = json.loads((D11/'audit/d11_v1_receipt_v2.json').read_text())
        bound += [QCB/name for name in old_audit['sections']['provenance']['source_dispositions']['tracked']]
        bound += [D12/name for name in ('second.c', 'second_interop.py', 'native_acquisition.py',
            'build.py', 'diagnose.py', 'README.md', 'PROTOCOL.md', 'CHECK_SCHEDULE.md',
            'API_CONTRACT.md', 'REVIEW.md', 'CODE_REVIEW.md', 'BUILD_REVIEW.md', 'build/receipt.json')]
        bound += [ROOT/name for name in ('SECOND_SUPPORT_CUT_PROPOSAL.md',
            'SECOND_SUPPORT_CUT_REVIEW.md', 'SECOND_NATIVE_API_DESIGN.md')]
        paths = list(dict.fromkeys(str(p.relative_to(repo)) for p in bound))
        subprocess.run(['git', 'ls-files', '--error-unmatch', '--', *paths], cwd=repo,
                       stdout=subprocess.DEVNULL, check=True, timeout=remaining())
        subprocess.run(['git', 'diff', '--exit-code', 'HEAD', '--', *paths], cwd=repo,
                       stdout=subprocess.DEVNULL, check=True, timeout=remaining())
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo,
                    text=True, timeout=remaining()).strip()
        headers_path = REFERENCE/'build/headers_receipt.json'
        if digest(headers_path) != HEADERS_SHA:
            raise RuntimeError('changed matching FLINT header receipt')
        headers = json.loads(headers_path.read_text())['headers_sha256']
        if len(headers) != 177:
            raise RuntimeError('unreviewed header inventory')
        for path, expected in headers.items():
            if digest(REFERENCE/path) != expected:
                raise RuntimeError('changed matching header: '+path)
        import flint
        if (flint.__version__, flint.__FLINT_VERSION__, flint.ctx.threads) != ('0.9.0', '3.6.0', 1):
            raise RuntimeError('unreviewed build environment')
        if any(C.sizeof(t) != 8 for t in (C.c_long, C.c_ulong, C.c_size_t, C.c_void_p)):
            raise RuntimeError('unreviewed Python ABI')
        libs = Path(sysconfig.get_path('platlib'))/'python_flint.libs'
        library = libs/'libflint-6839011d.so.24.0.0'
        if digest(library) != LIBRARY_SHA:
            raise RuntimeError('changed bundled FLINT library')
        compiler = Path(sys.prefix)/'bin/gcc'
        output = build/'libqcb_block.so'
        command = [str(compiler), '-O3', '-std=c11', '-fPIC', '-shared', '-Wall', '-Wextra',
            '-Werror=implicit-function-declaration', '-D_POSIX_C_SOURCE=200809L',
            '-I'+str(CONFIG/'src'), '-I'+str(CONFIG), '-I'+str(SOURCE),
            '-I'+str(Path(sys.prefix)/'include'), str(HERE/'block.c'), str(library),
            '-Wl,-rpath,'+str(libs), '-Wl,-rpath-link,'+str(libs), '-o', str(output)]
        receipt.update(command=command, commit=commit, scientific_sources_sha256=scientific,
            sources_sha256=scientific, documents_sha256=documents,
            source_sha256=scientific['block.c'], build_script_sha256=scientific['build.py'],
            code_review_sha256=digest(HERE/'CODE_REVIEW.md'),
            premises_sha256={str(p): v for p, v in PREMISES.items()},
            evidence_sha256={str(p): v for p, v in EVIDENCE.items()},
            primary_sources_sha256={str(SOURCE/p): v for p, v in PRIMARY.items()},
            inherited_sources_sha256=inherited['sources'], inherited_dependencies_sha256=inherited['dependencies'],
            headers_receipt=str(headers_path), headers_receipt_sha256=HEADERS_SHA, headers_checked=len(headers),
            compiler_sha256=digest(compiler.resolve()), compiler_version=None,
            linked_library=str(library), library_sha256=LIBRARY_SHA,
            dependencies_sha256={str(p): digest(p) for p in sorted(libs.glob('*.so.*'))},
            status='COMPILING')
        persist()
        t = perf_counter()
        try:
            completed = subprocess.run(command, capture_output=True, text=True, timeout=remaining())
            stdout, stderr, returncode = completed.stdout, completed.stderr, completed.returncode
        except subprocess.TimeoutExpired as error:
            stdout, stderr, returncode = error.stdout or b'', error.stderr or b'', None
            raise
        finally:
            # The timeout branch retains partial streams without a rerun.
            for name, value in (('compile.stdout', locals().get('stdout', '')),
                                ('compile.stderr', locals().get('stderr', ''))):
                if isinstance(value, bytes):
                    value = value.decode(errors='replace')
                (build/name).write_text(value)
            receipt.update(compile_seconds=perf_counter()-t, returncode=locals().get('returncode'))
        if returncode:
            raise RuntimeError('D13 compilation failed')
        receipt['status'] = 'COMPILED'
        receipt['artifact_sha256'] = digest(output)
        persist()
        receipt['compiler_version'] = subprocess.check_output([str(compiler), '--version'], text=True, timeout=remaining())
        remaining()
        # Only these constant ABI/version getters are called, never solve/apply.
        lib = C.CDLL(str(output), mode=os.RTLD_NOW | os.RTLD_LOCAL)
        for name in ('header_version', 'runtime_version'):
            getter = getattr(lib, 'qcb_block_'+name)
            getter.argtypes, getter.restype = [], C.c_char_p
        lib.qcb_block_abi.argtypes, lib.qcb_block_abi.restype = [], C.c_int
        lib.qcb_block_word_size.argtypes, lib.qcb_block_word_size.restype = [C.c_char_p], C.c_long
        versions = (lib.qcb_block_header_version().decode(), lib.qcb_block_runtime_version().decode())
        words = {k: int(lib.qcb_block_word_size(k.encode('ascii'))) for k in WORDS}
        linkage = subprocess.check_output(['ldd', str(output)], text=True, timeout=remaining())
        receipt.update(header_version=versions[0], runtime_version=versions[1],
            bridge_abi=int(lib.qcb_block_abi()), word_sizes=words, ldd=linkage)
        if (versions != ('3.6.0', '3.6.0') or receipt['bridge_abi'] != 2
                or words != WORDS or str(library) not in linkage or 'not found' in linkage
                or any(path not in linkage for path in receipt['dependencies_sha256'])):
            raise RuntimeError('D13 metadata/linkage mismatch')
        if (any(digest(HERE/name) != value for name, value in scientific.items())
                or any(digest(HERE/name) != value for name, value in documents.items())
                or digest(output) != receipt['artifact_sha256']):
            raise RuntimeError('D13 build identity changed during compilation')
        final_pairs = [(HERE/'CODE_REVIEW.md', receipt['code_review_sha256']),
                       (headers_path, HEADERS_SHA)]
        final_pairs += list((PREMISES | EVIDENCE).items())
        final_pairs += [(SOURCE/p, v) for p, v in PRIMARY.items()]
        final_pairs += [(REFERENCE/p, v) for p, v in headers.items()]
        final_pairs += [(QCB/p, v) for p, v in inherited['sources'].items()]
        final_pairs += [(Path(p), v) for p, v in
                        (inherited['dependencies'] | receipt['dependencies_sha256']).items()]
        for path, expected in final_pairs:
            if digest(path) != expected:
                raise RuntimeError('D13 source/dependency changed during build: '+str(path))
            remaining()
        if subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo,
                text=True, timeout=remaining()).strip() != commit:
            raise RuntimeError('D13 build commit changed during compilation')
        remaining()
        receipt.update(status='BUILD_ACCEPTABLE', passed=True)
        persist()
        remaining()
    except Exception as error:
        receipt.update(status='CENSORED_TIME' if isinstance(error, (TimeoutError, subprocess.TimeoutExpired)) else 'ERROR',
            passed=False, error=type(error).__name__+': '+str(error))
        persist()
        raise
    print(json.dumps(dict(status=receipt['status'], passed=receipt['passed'], receipt=str(receipt_path),
        mathematical_operation_executed=False, build_to_receipt_seconds=receipt['build_to_receipt_seconds']),
        allow_nan=False), flush=True)


if __name__ == '__main__':
    main()

