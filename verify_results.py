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


def verify_direct_pari():
    """Check archived low-rank measurements without importing their workers."""
    data = ROOT / 'results/pari_direct_20260930'
    code = ROOT / 'core/pari_direct_20260930'
    moduli = {2: 7, 4: 19, 6: 67, 8: 283, 12: 4179}
    cases = {(n, rank, parameter) for n in moduli for rank in (2, 3)
             for parameter in (1, 2, 3)}
    small_cases = {case for case in cases if case[0] in (2, 4, 6)}

    def read_json(name):
        return json.loads((data / name).read_text())

    def read_jsonl(name):
        return [json.loads(line) for line in (data / name).read_text().splitlines()]

    def identity(row):
        return row['N'], row['L'], row['a']

    def check_measurements(rows, algorithms):
        expected = {(case, algorithm, repeat) for case in cases
                    for algorithm in algorithms for repeat in range(3)}
        seen, values = set(), {}
        for row in rows:
            case = identity(row)
            key = case, row['algorithm'], row['repetition']
            require(key in expected and key not in seen,
                    'Unexpected or duplicate low-rank measurement')
            seen.add(key)
            require(row['status'] == 'PASS', 'Low-rank worker did not pass')
            require(row['modulus'] == moduli[case[0]], 'Low-rank modulus differs')
            require(type(row['value']) is int, 'Low-rank output is not an exact integer')
            require(values.setdefault(case, row['value']) == row['value'],
                    'Low-rank methods or repetitions disagree')
            for metric in ('compute_seconds', 'process_seconds', 'max_rss_kib'):
                require(math.isfinite(row[metric]) and row[metric] > 0,
                        f'Invalid low-rank metric: {metric}')
            require(row['matches_direct'] is (True if case in small_cases else None),
                    'Saved direct-reference flag disagrees with validation coverage')
        require(seen == expected, 'Low-rank measurement grid is incomplete')
        return values

    direct = read_jsonl('direct_raw.jsonl')
    historical = read_jsonl('raw.jsonl')
    direct_values = check_measurements(direct, ('pari_direct',))
    historical_values = check_measurements(historical, ('ah_unit_pivot', 'pari_binary'))
    require(direct_values == historical_values, 'Direct PARI and historical outputs differ')
    require(all(row['matches_other_runs'] is True for row in historical),
            'Historical exact-comparison flag failed')

    validation = read_json('validation.json')
    require(validation['status'] == 'PASS', 'Saved low-rank validation did not pass')
    saved_references = {}
    for row in validation['direct_checks']:
        case = identity(row)
        require(case in small_cases and case not in saved_references,
                'Unexpected or duplicate saved character-sum reference')
        require(row['direct_sum'] == row['ah_result'] == direct_values[case],
                'Saved independent character sum disagrees with measurements')
        saved_references[case] = row['direct_sum']
    require(set(saved_references) == small_cases,
            'Saved independent character-sum coverage differs')
    for n in moduli:
        for parameter in (1, 2, 3):
            require(direct_values[(n, 3, parameter)] ==
                    direct_values[(n, 2, parameter)] ** 2 - (1 << n),
                    'Saved low-rank outputs violate the rank-three identity')

    for row in direct:
        require(row['matches_original_algorithms'] is True,
                'Direct PARI exact-comparison flag failed')
        require(set(row['stages']) == {'field_parameter_seconds', 'curve_seconds',
                                     'point_count_recovery_seconds'},
                'Direct PARI stage names differ')
        require(all(math.isfinite(value) and value >= 0
                    for value in row['stages'].values()), 'Invalid Direct PARI stage time')
        check_number(sum(row['stages'].values()), row['compute_seconds'])
        startup = row['entry_import_runtime_seconds']
        require(math.isfinite(startup) and startup >= 0, 'Invalid Direct PARI startup time')
        require(row['process_seconds'] >= startup + row['compute_seconds'],
                'Direct PARI total process time is shorter than recorded stages')

    summary = read_json('direct_summary.json')
    require(summary['status'] == 'PASS' and summary['workers'] == 90
            and summary['case_count'] == 30, 'Direct PARI summary counts differ')
    groups = defaultdict(list)
    for row in direct:
        groups[identity(row)].append(row)
    seen = set()
    for group in summary['groups']:
        case = identity(group)
        require(case in cases and case not in seen, 'Unexpected Direct PARI summary group')
        seen.add(case)
        require(group['algorithm'] == 'pari_direct' and group['passed_repetitions'] == 3,
                'Direct PARI summary identity differs')
        require(group['value'] == direct_values[case], 'Direct PARI summary output differs')
        for metric in ('compute_seconds', 'process_seconds', 'max_rss_kib'):
            values = [row[metric] for row in groups[case]]
            for statistic, function in [('median', statistics.median), ('min', min), ('max', max)]:
                check_number(group[metric][statistic], function(values))
    require(seen == cases, 'Direct PARI summary grid is incomplete')

    # The unexported amendment is a historical protocol, not an exported artifact.
    # The failed identity refers to the old worker under its original filename.
    identity_names = {'AMENDMENT_DIRECT_PARI.md', 'direct_worker.py', 'direct_run.py',
                      'raw.jsonl', 'environment.json'}
    for filename, worker_name in [('direct_identity.json', 'direct_worker.py'),
                                  ('direct_failed_identity.json', 'direct_worker_failed.py')]:
        recorded = read_json(filename)
        require(set(recorded) == identity_names, 'Historical Direct PARI identity keys differ')
        paths = {'direct_worker.py': code / worker_name,
                 'direct_run.py': code / 'direct_run.py',
                 'raw.jsonl': data / 'raw.jsonl',
                 'environment.json': data / 'environment.json'}
        for name, path in paths.items():
            require(hashlib.sha256(path.read_bytes()).hexdigest() == recorded[name],
                    f'Historical Direct PARI identity differs: {filename}: {name}')
    failed = read_jsonl('direct_failed_attempt.jsonl')
    require(len(failed) == 1 and failed[0]['status'] == 'ERROR',
            'Expected one retained failed Direct PARI metadata attempt')
    require(identity(failed[0]) == (2, 2, 1) and failed[0]['repetition'] == 0,
            'Retained Direct PARI failure identity differs')

    return {'workers': len(direct), 'cases': len(cases),
            'historical_comparison_workers': len(historical),
            'saved_direct_character_sum_checks': len(saved_references),
            'summary_cells_recomputed': len(seen),
            'retained_failed_attempts': len(failed),
            'historical_identity_hashes_checked': 8,
            'unexported_identity_artifacts_not_checked': ['AMENDMENT_DIRECT_PARI.md']}


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
        'direct_PARI_archive': verify_direct_pari(),
        'scope': 'Distributed measurements and hashes verified; mathematical kernels not rerun.'},
        indent=2), flush=True)


if __name__ == '__main__':
    main()
