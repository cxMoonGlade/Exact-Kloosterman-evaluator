#!/usr/bin/env python3
"""Export compact verified results from retained local snapshots; stdlib only.

Does not run mathematical kernels, benchmarks, or commands from the snapshots.
The snapshots are input evidence and need not be published with this export.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_csv(path, rows):
    require(bool(rows), 'Empty output: ' + path.name)
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class Snapshot:
    def __init__(self, path, label):
        self.path = Path(path).resolve()
        self.label = label
        self.manifest_bytes = (self.path / 'MANIFEST.json').read_bytes()
        self.inventory = json.loads(self.manifest_bytes)['files']
        self.checked = {}

    def data(self, name):
        path = (self.path / name).resolve()
        require(path.is_relative_to(self.path), 'Source escapes snapshot')
        data = path.read_bytes()
        item = self.inventory[name]
        require(digest(data) == item['sha256'], 'Source hash mismatch: ' + name)
        self.checked[name] = item
        return data

    def read(self, name):
        return json.loads(self.data(name))

    def events(self, name):
        events = []
        for line in self.data(name).decode('utf-8').splitlines():
            try:
                value = json.loads(line)
            except ValueError:
                continue
            if isinstance(value, dict):
                events.append(value)
        return events

    def receipts(self, directory):
        return [(str(p.relative_to(self.path)), self.read(str(p.relative_to(self.path))))
                for p in sorted((self.path / directory).glob('*.json'))
                if not p.name.endswith(('.launch.json', '.progress.json'))]


def export(artifact, extension, pari_results, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    main = Snapshot(artifact, 'artifact')
    extra = Snapshot(extension, 'artifact_extension')
    require(all(not output.resolve().is_relative_to(s.path) for s in (main, extra)),
            'Output must be outside the retained snapshots')
    records, clocks, sources = [], [], {}
    six = 'repository/optimized_pair_20261004/reconstruction/rank_extension/runs/rank_v1'
    base = 'repository/baseline_alignment_20261006'
    even = 'repository/even_rank_extension_20261006'

    def add(snapshot, dataset, receipt_name, row, raw_name, attempt=None):
        require(row['status'] == 'COMPLETE', 'Incomplete measurement')
        require(row.get('exit', row.get('exit_code')) == 0, 'Failed process')
        events = snapshot.events(raw_name)
        require(not any(e.get('event') == 'error' for e in events), 'Raw error')
        payload = row['result']
        require([e for e in events if e.get('event') == 'result'] == [payload],
                'Receipt/raw output mismatch')
        n, rank = row['N'], row['L']
        method = row.get('backend', row.get('method'))
        attempt = attempt or row.get('attempt_id') or Path(raw_name).stem
        parameters = payload['parameters']
        values = payload.get('outputs')
        if values is None:
            values = [q['value'] for q in payload['queries']]
        require(len(parameters) == len(values) == 4, 'Expected four queries')
        require(all(type(x) is int for x in values), 'Noninteger output')
        if 'queries' in payload:
            require([q['a'] for q in payload['queries']] == parameters, 'Query order')
            require([q['value'] for q in payload['queries']] == values, 'Query values')
        raw_info = snapshot.inventory[raw_name]
        original = raw_info['original_sha256']
        sample_id = original[:16]
        if sample_id in sources:
            require(sources[sample_id]['original_sha256'] == original, 'Sample ID collision')
        else:
            sources[sample_id] = dict(sample_id=sample_id, original_sha256=original,
                source_sha256=raw_info['sha256'], source_snapshot=snapshot.label,
                source_file=raw_name)
        wall = row['wall_seconds']
        require(type(wall) in (int, float) and math.isfinite(wall) and wall > 0, 'Invalid clock')
        clocks.append(dict(dataset=dataset, attempt_id=attempt, N=n, L=rank,
            method=method, repeat=row['repeat'], wall_seconds=wall, sample_id=sample_id,
            receipt_sha256=snapshot.inventory[receipt_name]['sha256']))
        for index, (parameter, value) in enumerate(zip(parameters, values)):
            records.append(dict(dataset=dataset, attempt_id=attempt, N=n, L=rank,
                method=method, repeat=row['repeat'], f=payload['f'], query_index=index,
                a=parameter, value=value, sample_id=sample_id,
                source_sha256=raw_info['sha256']))

    for row in main.events(six + '/results.jsonl'):
        add(main, 'six_case', six + '/results.jsonl', row, six + '/' + row['raw_file'])
    refs = main.read(base + '/reference/references.json')
    for cell in refs['cells']:
        for method, details in cell['methods'].items():
            for ref in details['repetitions']:
                name = base + '/reference/' + ref['frozen_receipt_path']
                row = main.read(name)
                add(main, 'aligned_reference', name, row,
                    base + '/reference/' + ref['frozen_raw_path'])
    for name, row in main.receipts(base + '/runs/v1/attempts'):
        if row['status'] == 'COMPLETE':
            require(not row.get('execution_epoch') or row.get('interval_complete'), 'Incomplete interval')
            add(main, 'aligned_baseline', name, row, name[:-5] + '.jsonl')
    odd = 'repository/odd_rank_benchmark_20261006/runs/v1/attempts'
    for name, row in main.receipts(odd):
        add(main, 'odd_even', name, row, name[:-5] + '.stdout')
    for name, row in main.receipts(even + '/existing_serial_42/raw'):
        add(main, 'even_serial', name, row, name[:-5] + '.jsonl')
    for name, row in main.receipts(even + '/runs/n36_v1/attempts'):
        add(main, 'even_n36', name, row, name[:-5] + '.stdout')
    extension_dir = even + '/figure1_extension_evidence/matrix_raw'
    for name, row in extra.receipts(extension_dir):
        add(extra, 'figure1_extension', name, row, name[:-5] + '.jsonl')

    counts = collections.Counter(row['dataset'] for row in clocks)
    expected = dict(six_case=36, aligned_reference=150, aligned_baseline=375,
                    odd_even=240, even_serial=252, even_n36=42, figure1_extension=120)
    require(dict(counts) == expected, 'Dataset worker counts differ from completed scope')
    require(len(records) == 4 * sum(expected.values()), 'Output count')
    by_fixture = {}
    for row in records:
        # Never compare across datasets or different field representations.
        key = (row['dataset'], row['N'], row['L'], row['f'], row['a'])
        require(by_fixture.setdefault(key, row['value']) == row['value'], 'Within-dataset integer mismatch')
    identities = [(r['dataset'], r['N'], r['L'], r['method'], r['repeat']) for r in clocks]
    require(len(identities) == len(set(identities)), 'Duplicate repeat identity')
    groups = collections.defaultdict(list)
    for row in clocks:
        groups[row['dataset'], row['N'], row['L'], row['method']].append(row)
    summaries = []
    for (dataset, n, rank, method), rows in sorted(groups.items()):
        rows.sort(key=lambda r: r['repeat'])
        require([r['repeat'] for r in rows] == [0, 1, 2], 'Missing repeat')
        times = [r['wall_seconds'] for r in rows]
        summaries.append(dict(dataset=dataset, N=n, L=rank, method=method,
            repeat_0=times[0], repeat_1=times[1], repeat_2=times[2],
            median=statistics.median(times), minimum=min(times), maximum=max(times)))

    records.sort(key=lambda r: (r['dataset'], r['N'], r['L'], r['method'], r['repeat'], r['query_index']))
    clocks.sort(key=lambda r: (r['dataset'], r['N'], r['L'], r['method'], r['repeat']))
    write_csv(output / 'exact_outputs.csv', records)
    (output / 'exact_outputs.csv.gz').write_bytes(gzip.compress((output / 'exact_outputs.csv').read_bytes(), mtime=0))
    write_csv(output / 'timing_samples.csv', clocks)
    write_csv(output / 'timing_summary.csv', summaries)
    write_csv(output / 'output_sources.csv', sorted(sources.values(), key=lambda r: r['sample_id']))

    copied = {}
    for folder, label in [(main.path / 'recomputed', 'artifact/recomputed'),
                          (extra.path / 'reproduced_figure1_extension', 'artifact_extension/reproduced_figure1_extension'),
                          (Path(pari_results), 'results')]:
        for path in sorted(folder.iterdir()):
            if path.is_file() and (path.suffix == '.csv' or path.name.endswith('CHECK.json')):
                destination = output / path.name
                shutil.copyfile(path, destination)
                copied[path.name] = dict(source=label + '/' + path.name,
                    source_sha256=digest(path.read_bytes()))
    figure_base = 'repository/paper_low_math/report_kloosterman/figures/'
    (output / 'figures').mkdir(exist_ok=True)
    for stem in ('aligned_times', 'odd_even_times', 'even_scaling_times'):
        for suffix in ('.svg', '.csv'):
            name = figure_base + stem + suffix
            relative = 'figures/' + stem + suffix
            data = main.data(name)
            (output / relative).write_bytes(data)
            copied[relative] = dict(source='artifact/' + name,
                source_sha256=digest(data), original_sha256=main.inventory[name]['original_sha256'])
    pari = [r for r in clocks if r['method'] == 'pari_cyclic']
    require(len(pari) == 75, 'PARI worker scope')
    shared = collections.defaultdict(set)
    for row in clocks:
        shared[row['sample_id']].add(row['dataset'])
    result = dict(status='PASS', scientific_execution=False,
        dataset_workers=dict(counts), exact_output_rows=len(records),
        unique_original_raw_streams=len(sources), pari_workers=len(pari), pari_exact_outputs=len(pari)*4,
        overlapping_source_streams=sum(len(v)>1 for v in shared.values()),
        overlap_policy='Repeated sample_id denotes the same original raw stream reused by multiple datasets; dataset totals must not be summed as independent repetitions.',
        precision_policy='Signed integer outputs copied exactly; decimal CSV fields must be parsed as arbitrary-precision integers.',
        consistency_scope='Compared only within (dataset, N, L, f, a); no cross-field or cross-dataset oracle claim.',
        source_manifest_sha256={s.label:digest(s.manifest_bytes) for s in (main, extra)},
        checked_source_files={s.label:len(s.checked) for s in (main, extra)},
        copied_result_sources=copied)
    (output / 'EXTRACTION_CHECK.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact', type=Path, required=True)
    parser.add_argument('--extension', type=Path, required=True)
    parser.add_argument('--pari-results', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(export(args.artifact, args.extension, args.pari_results, args.output), indent=2))
