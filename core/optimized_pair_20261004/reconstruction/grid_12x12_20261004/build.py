"""Grid12x12 identity inventory and one explicit, metadata-only native build.

Import is standard-library only and performs no I/O or scientific operation.
An outer supervisor supplies the 60-second wall and 8 GiB process limits.
Each source identity has its own immutable attempt directory and static review.
"""
import argparse
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
PREVIOUS = HERE.parent/'rank_extension'
SCIENTIFIC = ('grid_resident.c', 'grid_interop.py', 'grid_query.py',
              'grid_parameter.py', 'grid_parameter_query.py', 'build.py',
              'validation.py', 'diagnose.py')
DOCUMENTS = ('README.md', 'PROTOCOL.md', 'REVIEW.md')
LIBRARY_SHA = '871a4132fd1e9f3638391b2208e07088f8e3e72a10e41d45f58b150a60c2a1a9'
HEADERS_SHA = 'a94c463fafd1060ed192935311949a54cd341ca53abe3e0edc19f7aa2528fcfb'
# Frozen D13 build/receipt.json compiler_sha256, retained by D15's manifest.
COMPILER_SHA = '1b866d8450249db5231194834ff440ad20b1f5e9316c94cd6b59c2da575545e0'
PREMISES = {
    HERE/'CAPACITY_REVIEW.md': '08961c9f098bea59e28e295bdbedfc15748b331c7167adb1a974af1f6879ace5',
    HERE/'UPSTREAM_REVIEW.md': '8ba98aa522b9f5209741486da126b6a4fb7d1ce4d3d31ef417f715b2651691f7',
    HERE/'HARNESS_REVIEW.md': 'bf0e7586b253fc04dcbf99e0992d304fd7c89c15b249cc875999dbe348182f7f',
}
EVIDENCE = {
    PREVIOUS/'runs/rank_v1/manifest.json': '1e2ad1b4d238901401f557816ef9e0ed228968900e487a37008543d5b6f53802',
    PREVIOUS/'audit/rank_v1_receipt_v2.json': '6f1a82996b11d8bfc25958fe3fc8181b4e3e3d97403a89649df44828a01c7eb3',
    PREVIOUS/'AUDIT.md': 'a31b54709f636b2d242ed70793886b8554ff11daca34162cd9ecad93b3739939',
    PREVIOUS/'runs/rank_v1/summary.json': '2e4ff1b307d5bb2271b2896fc5b3837bd406f8752b4b26c4709989ff9c0636c4',
}

# Direct public API identities are unchanged from D15. Its manifest also
# retains the mathematical premises and transitive D3/D6/D11 source identities.
PRIMARY = {
    SOURCE/'fmpz/size.c': '834b301b39e1ea99ed9bb232f50e5984e794192405534cde3dd47a57f50d3015',
    SOURCE/'fmpz_poly_mat/mul.c': 'dbeb8ee0afbcb5da491962ea82e9fe2b62cd6b8efe1fc947098827163f37fcb7',
    SOURCE/'fmpz_poly_mat/mul_KS.c': '13f6e01c067c3058a2f67f507b33306bc5cf6dd6b00f6e9f81cfe65ccf511f9c',
    SOURCE/'fmpz_mod_poly/rem.c': 'e67a1693a930924144f00329fec37e7a87e4fa46d58a99e4a6b0e3722847e131',
    SOURCE/'fmpz_mod_poly/set.c': '94b6df700dab237ceb67e151408fbcf60ad552b5263b47fc24bed777c3dd6f7f',
    SOURCE/'fmpz_mod_poly/get.c': '582fc6a4cb0c29f3b682caee7f62522699c6943d31c756ab8760b199213a5fa9',
    SOURCE/'fmpz_mod_poly/compose_mod_brent_kung_vec_preinv.c': '4e5ef82286ef326acdbf4144aa4cc3a6b6316a32e1a32f690fc7f9b505580d36',
    SOURCE.parent/'doc/source/fmpz_mod_poly.rst': '77a3ebe036d657ea166bc15b2194fc7b646b184025f476f0ee971a3615ef6f4a',
    SOURCE/'fmpz_mod_poly/divrem_f.c': '2394f4685fd108b9ed45d519697bbeb483645f34f774202cd84807d5e7db548f',
    SOURCE/'fmpz_poly_mat/init.c': '4161f7e8464934d60abba0187f03749340f5755ca01411386a04b7d7305e9205',
    SOURCE/'fmpz_poly_mat/clear.c': '5432fef57dbc241eca7d8b954f9495eec4712f24f05dfa6f31babc82629cd70e',
    SOURCE/'fmpz_poly/bit_pack.c': '81263115f1847ae1150b1bbebcd13880b095c126c1f7e5e00bf604ae969821b5',
    SOURCE/'fmpz_poly/bit_unpack.c': '13d247f21e53c0a6aa99f4017a2e7d8861623f7007b2d5f320dc47f98e60f790',
    SOURCE/'fmpz_poly/fit_length.c': '5912298d57ff9cb98e98f1d5a6e92c718aa669bbf652fa1f804320a28c7b4309',
    SOURCE/'fmpz_mod_poly/powers_mod_naive.c': '19aceca7f92a729d586699d75b3bc975a735ca136d9f3251fda92957384b0873',
    SOURCE/'fmpz_mod_poly/mulmod_preinv.c': '0499430ef30f36967134c60e17f3e6dd3bb7695ac68efa9b32a2e805c06603ef',
    SOURCE/'fmpz_mat/init.c': 'a056ee15924fc4849924fef4ded7d85772077a51c3dcd1fed11f6aaac2009bcb',
    Path('environment/include/gmp.h'): '748ca6ea93d9f4e21d2888c36fcc62984794d2dd23916412c961c3cd62565676',
}


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def review_path(attempt):
    if type(attempt) is not int or not 1 <= attempt <= 3:
        raise ValueError('attempt must be 1, 2 or 3')
    return HERE/('CODE_REVIEW.md' if attempt == 1 else f'CODE_REVIEW_ATTEMPT_{attempt}.md')


def selected_build():
    """Read the exact, unique selection; never infer it from directory times."""
    review = (HERE/'BUILD_REVIEW.md').read_text()
    lines = review.splitlines()
    selected = [line for line in lines if line.startswith('selected_attempt:')]
    if len(selected) != 1 or selected[0] not in ('selected_attempt: 1', 'selected_attempt: 2', 'selected_attempt: 3'):
        raise RuntimeError('Grid12x12 requires exactly one explicit selected_attempt')
    if 'execution_status: EXECUTION_PERMITTED' not in lines:
        raise RuntimeError('Grid12x12 selected build lacks the independent execution gate')
    attempt = int(selected[0][-1])
    directory = HERE/'build'/f'attempt_{attempt}'
    receipt_path, artifact = directory/'receipt.json', directory/'libqcb_grid_resident.so'
    if any(digest(p) not in review for p in (receipt_path, artifact, review_path(attempt))):
        raise RuntimeError('Grid12x12 build review does not bind its selected identities')
    return attempt, receipt_path, artifact


def inherited_evidence():
    for path, expected in EVIDENCE.items():
        if digest(path) != expected:
            raise RuntimeError('changed rank-extension evidence: '+str(path))
    manifest = json.loads((PREVIOUS/'runs/rank_v1/manifest.json').read_text())
    audit = json.loads((PREVIOUS/'audit/rank_v1_receipt_v2.json').read_text())
    inventory = manifest['inventory']
    if (audit['status'] != 'PASS_SCOPED_READ_ONLY_ARTIFACT_AUDIT' or
            audit['errors'] or audit['sections']['cold']['workers'] != 36 or
            audit['sections']['cold']['outputs'] != 144 or len(inventory['sources']) != 328 or
            len(inventory['dependencies']) != 19 or len(inventory['tracked']) != 269 or
            len(inventory['hash_only']) != 59):
        raise RuntimeError('rank-extension prerequisite evidence is incomplete')
    tracked, hash_only = set(inventory['tracked']), set(inventory['hash_only'])
    if (tracked & hash_only or tracked | hash_only != set(inventory['sources']) or
            len(tracked) != 269 or len(hash_only) != 59):
        raise RuntimeError('rank-extension manifest source disposition is invalid')
    # The immutable manifest supplies classifications; the read-only audit
    # supplies successful evidence counts and binds the exact raw summary.
    for name in ('manifest.json', 'summary.json'):
        key = 'runs/rank_v1/'+name
        if audit['artifact_sha256'][key] != EVIDENCE[PREVIOUS/key]:
            raise RuntimeError('rank-extension audit does not bind '+key)
    return inventory, audit


def source_inventory(include_build=True):
    """Paid stdlib inventory for the parent harness; no scientific imports.

    QCB-relative source names retain the inherited tracked/hash-only distinction.
    The installed GMP header is an absolute dependency, never a repository file.
    """
    inherited, audit = inherited_evidence()
    sources, dependencies = dict(inherited['sources']), dict(inherited['dependencies'])
    hash_only = set(inherited['hash_only'])
    tracked = set(inherited['tracked'])
    owned = [HERE/name for name in (*SCIENTIFIC, *DOCUMENTS, 'CODE_REVIEW.md')]
    owned += list(PREMISES)+list(EVIDENCE)
    if include_build:
        attempt, receipt, artifact = selected_build()
        owned += [HERE/'BUILD_REVIEW.md', review_path(attempt), receipt]
        sources[str(artifact.relative_to(QCB))] = digest(artifact)
        hash_only.add(str(artifact.relative_to(QCB)))
    for path in owned:
        name = str(path.relative_to(QCB))
        sources[name] = digest(path)
        tracked.add(name)
    for path, expected in (PREMISES | PRIMARY).items():
        if digest(path) != expected:
            raise RuntimeError('changed Grid12x12 premise or primary source: '+str(path))
        if path.is_relative_to(QCB):
            name = str(path.relative_to(QCB))
            sources[name] = expected
            if path in PRIMARY and name not in tracked:
                hash_only.add(name)
        else:
            dependencies[str(path)] = expected
    for name, expected in sources.items():
        if digest(QCB/name) != expected:
            raise RuntimeError('changed inherited/source bytes: '+name)
    for name, expected in dependencies.items():
        if digest(Path(name)) != expected:
            raise RuntimeError('changed installed dependency: '+name)
    if tracked & hash_only or tracked | hash_only != set(sources):
        raise RuntimeError('Grid12x12 source disposition is incomplete or conflicting')
    return dict(sources=sources, dependencies=dependencies,
                tracked=sorted(tracked), hash_only=sorted(hash_only))


def main(argv=None):
    started = perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--attempt', required=True, type=int, choices=(1, 2, 3))
    args = parser.parse_args(argv)
    resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
    directory = HERE/'build'/f'attempt_{args.attempt}'
    directory.mkdir(parents=True, exist_ok=False)
    receipt_path = directory/'receipt.json'
    receipt = dict(status='STARTED', passed=False, attempt=args.attempt,
                   mathematical_operation_executed=False,
                   complete_build_budget_seconds=60, aggregate_budget_seconds=180,
                   memory_limit_bytes=8*1024**3)

    def persist():
        elapsed = perf_counter()-started
        receipt.update(build_to_receipt_seconds=elapsed,
            termination_overrun_seconds=max(0.0, elapsed-60),
            clock_scope='Entry through the sample before this write; the outer supervisor includes final receipt write/print/exit.')
        receipt_path.write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')

    def remaining():
        value = 60-(perf_counter()-started)
        if value <= 0:
            raise TimeoutError('Grid12x12 complete build deadline exceeded')
        return value

    persist()
    try:
        code_review = review_path(args.attempt)
        gate, review = (HERE/'REVIEW.md').read_text(), code_review.read_text()
        if 'code_status: CODE_PERMITTED' not in gate.splitlines() or 'build_status: BUILD_PERMITTED' not in review.splitlines():
            raise RuntimeError('Grid12x12 exact code/build gate absent')
        scientific = {name: digest(HERE/name) for name in SCIENTIFIC}
        documents = {name: digest(HERE/name) for name in DOCUMENTS}
        if any(value not in review for value in (*scientific.values(), *documents.values())):
            raise RuntimeError('Grid12x12 static review does not bind every source and contract')
        if any(documents[name] not in gate for name in DOCUMENTS if name != 'REVIEW.md'):
            raise RuntimeError('Grid12x12 preregistration review does not bind the contracts')
        previous = []
        for k in range(1, args.attempt):
            path = HERE/'build'/f'attempt_{k}'/'receipt.json'
            old = json.loads(path.read_text())
            if old.get('attempt') != k or old.get('status') not in ('BUILD_ACCEPTABLE', 'ERROR', 'CENSORED_TIME'):
                raise RuntimeError('previous build attempt has no terminal receipt')
            if old.get('scientific_sources_sha256') == scientific:
                raise RuntimeError('cannot repeat a previously attempted source identity')
            if digest(path) not in review:
                raise RuntimeError('incremental review does not bind the previous attempt')
            previous.append(dict(attempt=k, receipt=str(path), sha256=digest(path),
                                 seconds=old['build_to_receipt_seconds']))
        prior_seconds = sum(row['seconds'] for row in previous)
        if prior_seconds+60 > 180:
            raise RuntimeError('Grid12x12 aggregate build budget has no full attempt remaining')
        receipt.update(previous_attempts=previous, prior_recorded_build_seconds=prior_seconds)
        for path, expected in (PREMISES | EVIDENCE | PRIMARY).items():
            if digest(path) != expected:
                raise RuntimeError('changed Grid12x12 premise/primary: '+str(path))
            if path in PREMISES or path in EVIDENCE:
                if expected not in gate:
                    raise RuntimeError('Grid12x12 gate omits a premise/evidence identity')
            remaining()
        inventory = source_inventory(include_build=False)
        inherited, _ = inherited_evidence()
        repo = Path(subprocess.check_output(['git', 'rev-parse', '--show-toplevel'], cwd=HERE,
                    text=True, timeout=remaining()).strip())
        bound = [QCB/name for name in inventory['tracked']]+[code_review]
        paths = list(dict.fromkeys(str(p.relative_to(repo)) for p in bound))
        subprocess.run(['git', 'ls-files', '--error-unmatch', '--', *paths], cwd=repo,
                       stdout=subprocess.DEVNULL, check=True, timeout=remaining())
        subprocess.run(['git', 'diff', '--exit-code', 'HEAD', '--', *paths], cwd=repo,
                       stdout=subprocess.DEVNULL, check=True, timeout=remaining())
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo,
                    text=True, timeout=remaining()).strip()
        headers_path = REFERENCE/'build/headers_receipt.json'
        if digest(headers_path) != HEADERS_SHA:
            raise RuntimeError('changed matching header receipt')
        headers = json.loads(headers_path.read_text())['headers_sha256']
        if len(headers) != 177 or any(digest(REFERENCE/p) != v for p, v in headers.items()):
            raise RuntimeError('changed matching FLINT headers')
        import flint
        if (flint.__version__, flint.__FLINT_VERSION__, flint.ctx.threads) != ('0.9.0', '3.6.0', 1):
            raise RuntimeError('unreviewed build environment')
        if (any(C.sizeof(t) != 8 for t in (C.c_long, C.c_ulong, C.c_size_t, C.c_void_p)) or C.sizeof(C.c_int) != 4):
            raise RuntimeError('unreviewed Python ABI')
        libs = Path(sysconfig.get_path('platlib'))/'python_flint.libs'
        library = libs/'libflint-6839011d.so.24.0.0'
        if digest(library) != LIBRARY_SHA:
            raise RuntimeError('changed bundled FLINT library')
        compiler = Path(sys.prefix)/'bin/gcc'
        if digest(compiler.resolve()) != COMPILER_SHA:
            raise RuntimeError('changed D13-pinned compiler')
        output = directory/'libqcb_grid_resident.so'
        command = [str(compiler), '-O3', '-std=c11', '-fPIC', '-shared', '-Wall', '-Wextra',
            '-Werror=implicit-function-declaration', '-D_POSIX_C_SOURCE=200809L',
            '-I'+str(CONFIG/'src'), '-I'+str(CONFIG), '-I'+str(SOURCE),
            '-I'+str(Path(sys.prefix)/'include'), str(HERE/'grid_resident.c'), str(library),
            '-Wl,-rpath,'+str(libs), '-Wl,-rpath-link,'+str(libs), '-o', str(output)]
        receipt.update(status='COMPILING', command=command, commit=commit,
            scientific_sources_sha256=scientific, documents_sha256=documents,
            code_review=str(code_review), code_review_sha256=digest(code_review),
            premises_sha256={str(p): v for p, v in PREMISES.items()},
            evidence_sha256={str(p): v for p, v in EVIDENCE.items()},
            primary_sources_sha256={str(p): v for p, v in PRIMARY.items()},
            inherited_sources_sha256=inherited['sources'], inherited_dependencies_sha256=inherited['dependencies'],
            headers_receipt=str(headers_path), headers_receipt_sha256=HEADERS_SHA, headers_checked=len(headers),
            compiler_sha256=COMPILER_SHA, compiler_version=None,
            linked_library=str(library), library_sha256=LIBRARY_SHA,
            dependencies_sha256={str(p): digest(p) for p in sorted(libs.glob('*.so.*'))})
        persist()
        t = perf_counter()
        try:
            completed = subprocess.run(command, capture_output=True, text=True, timeout=remaining())
            stdout, stderr, returncode = completed.stdout, completed.stderr, completed.returncode
        except subprocess.TimeoutExpired as error:
            stdout, stderr, returncode = error.stdout or b'', error.stderr or b'', None
            raise
        finally:
            for name, value in (('compile.stdout', locals().get('stdout', '')),
                                ('compile.stderr', locals().get('stderr', ''))):
                if isinstance(value, bytes):
                    value = value.decode(errors='replace')
                (directory/name).write_text(value)
            receipt.update(compile_seconds=perf_counter()-t, returncode=locals().get('returncode'))
        if returncode:
            raise RuntimeError('Grid12x12 compilation failed')
        receipt.update(status='COMPILED', artifact_sha256=digest(output))
        persist()
        receipt['compiler_version'] = subprocess.check_output([str(compiler), '--version'], text=True, timeout=remaining())
        remaining()
        lib = C.CDLL(str(output), mode=os.RTLD_LOCAL | os.RTLD_NOW)
        for name in ('header_version', 'runtime_version'):
            getter = getattr(lib, 'qcb_resident_'+name)
            getter.argtypes, getter.restype = [], C.c_char_p
        lib.qcb_resident_abi.argtypes, lib.qcb_resident_abi.restype = [], C.c_int
        lib.qcb_resident_slong_size.argtypes, lib.qcb_resident_slong_size.restype = [], C.c_size_t
        versions = (lib.qcb_resident_header_version().decode(), lib.qcb_resident_runtime_version().decode())
        linkage = subprocess.check_output(['ldd', str(output)], text=True, timeout=remaining())
        receipt.update(header_version=versions[0], runtime_version=versions[1],
            bridge_abi=int(lib.qcb_resident_abi()), slong_size=int(lib.qcb_resident_slong_size()), ldd=linkage)
        if (versions != ('3.6.0', '3.6.0') or receipt['bridge_abi'] != 1 or receipt['slong_size'] != 8 or
                str(library) not in linkage or 'not found' in linkage or
                any(path not in linkage for path in receipt['dependencies_sha256'])):
            raise RuntimeError('Grid12x12 constant metadata/linkage mismatch')
        final = [(HERE/n, v) for n, v in (scientific | documents).items()]
        final += [(code_review, receipt['code_review_sha256']), (output, receipt['artifact_sha256']), (headers_path, HEADERS_SHA)]
        final += [(compiler.resolve(), COMPILER_SHA)]
        final += list((PREMISES | EVIDENCE | PRIMARY).items())
        final += [(REFERENCE/p, v) for p, v in headers.items()]
        final += [(QCB/p, v) for p, v in inventory['sources'].items()]
        final += [(Path(p), v) for p, v in (inventory['dependencies'] | receipt['dependencies_sha256']).items()]
        for path, expected in final:
            if digest(path) != expected:
                raise RuntimeError('Grid12x12 identity changed during build: '+str(path))
            remaining()
        if subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo,
                text=True, timeout=remaining()).strip() != commit:
            raise RuntimeError('Grid12x12 source commit changed during build')
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
