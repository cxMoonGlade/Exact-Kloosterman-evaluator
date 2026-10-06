"""Explicit D2 build and ABI checks only; no candidate mathematical operation."""
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
SOURCES = ('compose.c', 'interop.py', 'build.py')
SCIENTIFIC = SOURCES+('binary_batch.py', 'reference.py', 'diagnose.py')
LIBRARY_SHA = '871a4132fd1e9f3638391b2208e07088f8e3e72a10e41d45f58b150a60c2a1a9'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    protocol = HERE/'PROTOCOL.md'
    gate = HERE/'REVIEW.md'
    code_review = HERE/'CODE_REVIEW.md'
    protocol_sha = digest(protocol)
    gate_text = gate.read_text()
    if 'code_status: CODE_PERMITTED' not in gate_text.splitlines() or protocol_sha not in gate_text:
        raise RuntimeError('D2 premise gate missing or changed')
    review = code_review.read_text()
    if not any(marker in review.splitlines() for marker in (
            'build_status: BUILD_PERMITTED', 'execution_status: EXECUTION_PERMITTED')):
        raise RuntimeError('D2 code review has not permitted build/execution')
    sources = {name: digest(HERE/name) for name in SOURCES}
    if any(value not in review for value in sources.values()):
        raise RuntimeError('D2 code review does not bind all bridge sources')
    scientific = {name: digest(HERE/name) for name in SCIENTIFIC}
    repository = Path(subprocess.check_output(
        ['git', 'rev-parse', '--show-toplevel'], cwd=HERE, text=True).strip())
    bound = [HERE/name for name in SCIENTIFIC]+[protocol, gate, code_review]
    relative = [str(path.relative_to(repository)) for path in bound]
    tracked = subprocess.run(['git', 'ls-files', '--error-unmatch', '--']+relative,
                             cwd=repository, capture_output=True, text=True)
    if tracked.returncode:
        raise RuntimeError('D2 build sources/reviews are not all tracked')
    unchanged = subprocess.run(['git', 'diff', '--quiet', 'HEAD', '--']+relative,
                               cwd=repository, capture_output=True, text=True)
    if unchanged.returncode:
        raise RuntimeError('D2 build sources/reviews differ from committed HEAD')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repository, text=True).strip()
    header_path = REFERENCE/'build/headers_receipt.json'
    headers = json.loads(header_path.read_text())
    if len(headers['headers_sha256']) != 177:
        raise RuntimeError('unreviewed FLINT header inventory')
    for name, expected in headers['headers_sha256'].items():
        if digest(REFERENCE/name) != expected:
            raise RuntimeError('changed matching header: '+name)
    import flint
    if (flint.__version__, flint.__FLINT_VERSION__) != ('0.9.0', '3.6.0'):
        raise RuntimeError('unreviewed Python/native version')
    libs = Path(sysconfig.get_path('platlib'))/'python_flint.libs'
    library = libs/'libflint-6839011d.so.24.0.0'
    if digest(library) != LIBRARY_SHA:
        raise RuntimeError('changed bundled FLINT')
    build = HERE/'build'
    if build.exists():
        raise RuntimeError('refusing to overwrite an existing batch build directory')
    build.mkdir()
    output = build/'libqcb_batch.so'
    compiler = Path(sys.prefix)/'bin/gcc'
    command = [str(compiler), '-O3', '-std=c11', '-fPIC', '-shared', '-Wall', '-Wextra',
               '-Werror=implicit-function-declaration', '-D_POSIX_C_SOURCE=200809L',
               '-I'+str(CONFIG/'src'), '-I'+str(CONFIG), '-I'+str(SOURCE),
               '-I'+str(Path(sys.prefix)/'include'), str(HERE/'compose.c'), str(library),
               '-Wl,-rpath,'+str(libs), '-Wl,-rpath-link,'+str(libs), '-o', str(output)]
    started = perf_counter()
    completed = subprocess.run(command, capture_output=True, text=True)
    elapsed = perf_counter()-started
    (build/'compile.stdout').write_text(completed.stdout)
    (build/'compile.stderr').write_text(completed.stderr)
    receipt = dict(command=command, compile_seconds=elapsed, returncode=completed.returncode,
                   compiler_sha256=digest(compiler.resolve()),
                   compiler_version=subprocess.check_output([str(compiler), '--version'], text=True),
                   sources_sha256=sources, source_sha256=sources['compose.c'],
                   scientific_sources_sha256=scientific, commit=commit,
                   build_script_sha256=sources['build.py'],
                   headers_receipt=str(header_path), headers_receipt_sha256=digest(header_path),
                   headers_checked=len(headers['headers_sha256']),
                   dependencies_sha256={p.name: digest(p) for p in sorted(libs.glob('*.so.*'))},
                   library_sha256=LIBRARY_SHA, linked_library=str(library),
                   gate_sha256=digest(gate), code_review_sha256=digest(code_review),
                   protocol_sha256=protocol_sha, passed=False,
                   mathematical_operation_executed=False)
    if completed.returncode:
        (build/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
        raise RuntimeError('batch bridge compilation failed; build receipt retained')
    lib = C.CDLL(str(output), mode=os.RTLD_NOW | os.RTLD_LOCAL)
    for name in ('header_version', 'runtime_version'):
        getattr(lib, 'qcb_batch_'+name).restype = C.c_char_p
    lib.qcb_batch_slong_size.restype = C.c_size_t
    versions = (lib.qcb_batch_header_version().decode(), lib.qcb_batch_runtime_version().decode())
    linkage = subprocess.check_output(['ldd', str(output)], text=True)
    receipt.update(header_version=versions[0], runtime_version=versions[1],
                   bridge_abi=lib.qcb_batch_abi(), slong_size=lib.qcb_batch_slong_size(),
                   artifact_sha256=digest(output), ldd=linkage)
    receipt['passed'] = (versions == ('3.6.0', '3.6.0') and receipt['bridge_abi'] == 1
                         and receipt['slong_size'] == C.sizeof(C.c_long)
                         and str(library) in linkage)
    (build/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    if not receipt['passed']:
        raise RuntimeError('batch bridge identity check failed; receipt retained')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
