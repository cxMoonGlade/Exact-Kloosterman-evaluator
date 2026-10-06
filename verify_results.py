#!/usr/bin/env python3
"""Verify the focused export and recompute statistics from saved measurements."""
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_csv(name):
    with (ROOT / 'results' / name).open(newline='') as stream:
        return list(csv.DictReader(stream))


def check_number(actual, expected):
    require(math.isclose(float(actual), float(expected), rel_tol=1e-12, abs_tol=1e-12),
            f'Numeric mismatch: {actual} versus {expected}')


def main():
    manifest = json.loads((ROOT / 'EXPORT_MANIFEST.json').read_text())
    for name, record in manifest['files'].items():
        path = (ROOT / name).resolve()
        require(path.is_relative_to(ROOT), 'Manifest path escapes repository')
        require(path.stat().st_size == record['bytes'], f'Size changed: {name}')
        require(hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256'],
                f'Hash changed: {name}')
    sources = json.loads((ROOT / 'core/SOURCE_MANIFEST.json').read_text())
    for name, record in sources['files'].items():
        require(hashlib.sha256((ROOT / 'core' / name).read_bytes()).hexdigest()
                == record['distributed_sha256'], f'Core source changed: {name}')

    extraction = json.loads((ROOT / 'results/EXTRACTION_CHECK.json').read_text())
    times = read_csv('timing_samples.csv')
    outputs = read_csv('exact_outputs.csv')
    workers, groups = {}, defaultdict(dict)
    for row in times:
        key = (row['dataset'], row['attempt_id'])
        require(key not in workers, 'Duplicate worker identity')
        workers[key] = row
        wall = float(row['wall_seconds'])
        require(math.isfinite(wall) and wall > 0, 'Invalid process wall time')
        group = (row['dataset'], row['N'], row['L'], row['method'])
        repeat = int(row['repeat'])
        require(repeat not in groups[group], 'Duplicate repetition')
        groups[group][repeat] = wall
    require(dict(Counter(row['dataset'] for row in times)) == extraction['dataset_workers'],
            'Dataset worker count changed')
    require(len(outputs) == extraction['exact_output_rows'], 'Output count changed')

    provenance = {(r['sample_id'], r['source_sha256']) for r in read_csv('output_sources.csv')}
    query_rows, exact_values, repeated_values = defaultdict(dict), {}, {}
    for row in outputs:
        key = (row['dataset'], row['attempt_id'])
        require(key in workers, 'Output has no timing record')
        worker = workers[key]
        require(all(row[k] == worker[k] for k in ('N', 'L', 'method', 'repeat', 'sample_id')),
                'Output and timing identities differ')
        query = int(row['query_index'])
        require(query not in query_rows[key], 'Duplicate query index')
        query_rows[key][query] = row
        n, rank, field, parameter, value = [int(row[k]) for k in ('N', 'L', 'f', 'a', 'value')]
        require(n >= 2 and rank >= 2 and 0 < parameter < (1 << n), 'Invalid input labels')
        require(field.bit_length() == n + 1, 'Wrong field-polynomial degree')
        signature = (row['dataset'], n, rank, field, parameter)
        require(exact_values.setdefault(signature, value) == value,
                'Methods or repetitions disagree on an exact integer')
        repeated = (row['sample_id'], query)
        identity = (n, rank, field, parameter, value)
        require(repeated_values.setdefault(repeated, identity) == identity,
                'An overlapping original sample changed between datasets')
        require((row['sample_id'], row['source_sha256']) in provenance,
                'Output source hash absent from provenance')
    require(set(query_rows) == set(workers), 'Worker missing exact outputs')
    require(all(set(rows) == {0, 1, 2, 3} for rows in query_rows.values()),
            'Worker does not contain exactly four ordered queries')

    summaries = read_csv('timing_summary.csv')
    require(len(summaries) == len(groups), 'Timing summary cell count differs')
    seen = set()
    for row in summaries:
        key = (row['dataset'], row['N'], row['L'], row['method'])
        require(key in groups and key not in seen, 'Unexpected timing summary cell')
        seen.add(key)
        repeats = groups[key]
        require(set(repeats) == {0, 1, 2}, 'Incomplete three-process sample')
        for r in range(3):
            check_number(row[f'repeat_{r}'], repeats[r])
        for column, function in [('median', statistics.median), ('minimum', min), ('maximum', max)]:
            check_number(row[column], function(repeats.values()))

    pari = [r for r in times if r['dataset'] == 'aligned_baseline' and r['method'] == 'pari_cyclic']
    require(len(pari) == 75, 'PARI worker count differs')
    pari_table = read_csv('pari_25_cells.csv')
    require(len(pari_table) == 25, 'PARI grid is incomplete')
    for row in pari_table:
        key = ('aligned_baseline', row['N'], row['L'], 'pari_cyclic')
        samples = list(groups[key].values())
        require(sorted(json.loads(row['samples'])) == sorted(samples), 'PARI samples differ')
        check_number(row['median'], statistics.median(samples))
        check_number(row['minimum'], min(samples))
        check_number(row['maximum'], max(samples))

    for row in read_csv('six_case_times_ms.csv'):
        for method in ('ah_block_second', 'ordinary_int'):
            samples = groups[('six_case', row['N'], row['L'], method)].values()
            check_number(row[method + '_median_ms'], 1000 * statistics.median(samples))

    print(json.dumps({'status': 'PASS', 'scientific_execution': False,
        'verified_release_files': manifest['file_count'], 'core_source_files': len(sources['files']),
        'dataset_workers': extraction['dataset_workers'], 'timing_cells_recomputed': len(groups),
        'exact_output_rows': len(outputs), 'unique_original_samples': len({r['sample_id'] for r in times}),
        'PARI_workers': len(pari), 'PARI_exact_outputs': 4 * len(pari),
        'scope': 'Distributed measurements and hashes verified; mathematical kernels not rerun.'},
        indent=2), flush=True)


if __name__ == '__main__':
    main()
