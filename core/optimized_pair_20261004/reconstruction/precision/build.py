"""Explicit D5 compilation and identity checks; no mathematical query."""
import ctypes as C
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import sysconfig
from time import perf_counter

HERE = Path(__file__).resolve().parent
OWNER = HERE.parent.parent
REFERENCE = OWNER.parent/'specialization_optimization_20260930'
CONFIG = REFERENCE/'build/flint_header_config'
SOURCE = REFERENCE/'sources/flint-3.6.0/src'
SOURCES = ('guard.c', 'guard_interop.py', 'build.py')
SCIENTIFIC = SOURCES+('diagnose.py',)
LIBRARY_SHA = '871a4132fd1e9f3638391b2208e07088f8e3e72a10e41d45f58b150a60c2a1a9'
HEADERS_SHA = 'a94c463fafd1060ed192935311949a54cd341ca53abe3e0edc19f7aa2528fcfb'
PREMISES = {
    str(OWNER/'PRECISION_REDUCTION_PROPOSAL.md'):
        '540e1b66b57e50139f36cf96b4f668dccef75e89d658de67deb94b02bca7565d',
    str(OWNER/'PRECISION_REDUCTION_REVIEW.md'):
        'f88bd8175fb7f73b71fb487d4127724b7e85b87d5d3538583be41c6d16f32e88',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    build_started = perf_counter()

    def remaining():
        value = 60-(perf_counter()-build_started)
        if value <= 0:
            raise TimeoutError('D5 complete build deadline exceeded')
        return value

    protocol, api = HERE/'PROTOCOL.md', HERE.parent/'first_transform/API_REVIEW.md'
    gate, code_review = HERE/'REVIEW.md', HERE/'CODE_REVIEW.md'
    protocol_sha, api_sha = digest(protocol), digest(api)
    gate_text, review = gate.read_text(), code_review.read_text()
    if ('code_status: CODE_PERMITTED' not in gate_text.splitlines()
            or protocol_sha not in gate_text or api_sha not in gate_text):
        raise RuntimeError('D5 premise gate missing or changed')
    if 'build_status: BUILD_PERMITTED' not in review.splitlines():
        raise RuntimeError('D5 static code review has not permitted build')
    sources = {name: digest(HERE/name) for name in SOURCES}
    scientific = {name: digest(HERE/name) for name in SCIENTIFIC}
    if any(value not in review for value in scientific.values()):
        raise RuntimeError('D5 static review does not bind all scientific sources')
    if protocol_sha not in review or api_sha not in review:
        raise RuntimeError('D5 static review does not bind prerequisites')
    premises = {path: digest(Path(path)) for path in PREMISES}
    if (premises != PREMISES or any(value not in gate_text or value not in review
                                    for value in premises.values())):
        raise RuntimeError('D5 precision proof identity is changed or unreviewed')
    repository = Path(subprocess.check_output(
        ['git', 'rev-parse', '--show-toplevel'], cwd=HERE, text=True, timeout=remaining()).strip())
    bound = ([HERE/name for name in SCIENTIFIC]
             + [protocol, api, HERE/'README.md', gate, code_review]
             + [Path(path) for path in premises])
    relative = [str(path.relative_to(repository)) for path in bound]
    tracked = subprocess.run(['git', 'ls-files', '--error-unmatch', '--']+relative,
                             cwd=repository, capture_output=True, text=True, timeout=remaining())
    if tracked.returncode:
        raise RuntimeError('D5 build sources/reviews are not all tracked')
    unchanged = subprocess.run(['git', 'diff', '--quiet', 'HEAD', '--']+relative,
                               cwd=repository, capture_output=True, text=True, timeout=remaining())
    if unchanged.returncode:
        raise RuntimeError('D5 build sources/reviews differ from committed HEAD')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repository, text=True, timeout=remaining()).strip()
    header_path = REFERENCE/'build/headers_receipt.json'
    if digest(header_path) != HEADERS_SHA:
        raise RuntimeError('changed matching FLINT header receipt')
    headers = json.loads(header_path.read_text())
    if len(headers['headers_sha256']) != 177:
        raise RuntimeError('unreviewed FLINT header inventory')
    for name, expected in headers['headers_sha256'].items():
        if digest(REFERENCE/name) != expected:
            raise RuntimeError('changed matching FLINT header: '+name)
    import flint
    if (flint.__version__, flint.__FLINT_VERSION__) != ('0.9.0', '3.6.0'):
        raise RuntimeError('unreviewed D5 Python/native version')
    libs = Path(sysconfig.get_path('platlib'))/'python_flint.libs'
    library = libs/'libflint-6839011d.so.24.0.0'
    if digest(library) != LIBRARY_SHA:
        raise RuntimeError('changed bundled FLINT library')
    build = HERE/'build'
    if build.exists():
        raise RuntimeError('refusing to overwrite an existing D5 build directory')
    build.mkdir()
    output = build/'libqcb_guard.so'
    compiler = Path(sys.prefix)/'bin/gcc'
    command = [str(compiler), '-O3', '-std=c11', '-fPIC', '-shared', '-Wall', '-Wextra',
               '-Werror=implicit-function-declaration', '-D_POSIX_C_SOURCE=200809L',
               '-I'+str(CONFIG/'src'), '-I'+str(CONFIG), '-I'+str(SOURCE),
               '-I'+str(Path(sys.prefix)/'include'), str(HERE/'guard.c'), str(library),
               '-Wl,-rpath,'+str(libs), '-Wl,-rpath-link,'+str(libs), '-o', str(output)]
    started = perf_counter()
    timed_out = False
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=remaining())
        stdout, stderr, returncode = completed.stdout, completed.stderr, completed.returncode
    except (subprocess.TimeoutExpired, TimeoutError) as error:
        timed_out = True
        stdout, stderr = getattr(error, 'stdout', '') or '', getattr(error, 'stderr', '') or ''
        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors='replace')
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors='replace')
        returncode = None
    elapsed = perf_counter()-started
    (build/'compile.stdout').write_text(stdout)
    (build/'compile.stderr').write_text(stderr)
    receipt = dict(command=command, compile_seconds=elapsed, returncode=returncode,
                   status='CENSORED_TIME' if timed_out else 'COMPILED' if returncode == 0 else 'ERROR',
                   compiler_sha256=digest(compiler.resolve()),
                   compiler_version=None,
                   sources_sha256=sources, scientific_sources_sha256=scientific, commit=commit,
                   source_sha256=sources['guard.c'], build_script_sha256=sources['build.py'],
                   headers_receipt=str(header_path), headers_receipt_sha256=HEADERS_SHA,
                   headers_checked=len(headers['headers_sha256']),
                   dependencies_sha256={p.name: digest(p) for p in sorted(libs.glob('*.so.*'))},
                   library_sha256=LIBRARY_SHA, linked_library=str(library),
                   gate_sha256=digest(gate), code_review_sha256=digest(code_review),
                   protocol_sha256=protocol_sha, api_sha256=api_sha,
                   api_path=str(api), premises_sha256=premises, passed=False,
                   mathematical_operation_executed=False)
    receipt_path = build/'receipt.json'
    receipt['build_seconds'] = perf_counter()-build_started
    receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
    if timed_out or returncode:
        raise RuntimeError('D5 compilation failed or timed out; receipt retained')
    try:
        receipt['compiler_version'] = subprocess.check_output(
            [str(compiler), '--version'], text=True, timeout=remaining())
        remaining()
        lib = C.CDLL(str(output), mode=os.RTLD_NOW | os.RTLD_LOCAL)
        for name in ('header_version', 'runtime_version'):
            getter = getattr(lib, 'qcb_first_'+name)
            getter.argtypes = []
            getter.restype = C.c_char_p
        lib.qcb_first_abi.argtypes = []
        lib.qcb_first_abi.restype = C.c_int
        lib.qcb_first_slong_size.argtypes = []
        lib.qcb_first_slong_size.restype = C.c_size_t
        versions = (lib.qcb_first_header_version().decode(), lib.qcb_first_runtime_version().decode())
        linkage = subprocess.check_output(['ldd', str(output)], text=True, timeout=remaining())
        receipt.update(header_version=versions[0], runtime_version=versions[1],
                       bridge_abi=lib.qcb_first_abi(), slong_size=lib.qcb_first_slong_size(),
                       artifact_sha256=digest(output), ldd=linkage)
        receipt['passed'] = (versions == ('3.6.0', '3.6.0') and receipt['bridge_abi'] == 2
                             and receipt['slong_size'] == C.sizeof(C.c_long)
                             and str(library) in linkage)
        remaining()
        receipt['status'] = 'BUILD_ACCEPTABLE' if receipt['passed'] else 'ERROR'
    except Exception as error:
        receipt['passed'] = False
        receipt['status'] = 'CENSORED_TIME' if isinstance(error, (TimeoutError, subprocess.TimeoutExpired)) else 'ERROR'
        receipt['identity_error'] = type(error).__name__+': '+str(error)
        raise
    finally:
        receipt['build_seconds'] = perf_counter()-build_started
        receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
    if not receipt['passed']:
        raise RuntimeError('D5 identity check failed; receipt retained')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
