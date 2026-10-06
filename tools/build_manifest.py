#!/usr/bin/env python3
"""Record the explicit code-and-results release contents."""
import hashlib
import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    files = {}
    for path in sorted(root.rglob('*')):
        relative = path.relative_to(root)
        if (not path.is_file() or '.git' in relative.parts or
                '__pycache__' in relative.parts or relative.name == 'EXPORT_MANIFEST.json'):
            continue
        files[relative.as_posix()] = {'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    record = {'scope': 'Core mathematical source and completed experiment results only.',
              'file_count': len(files), 'total_bytes': sum(v['bytes'] for v in files.values()),
              'files': files}
    (root / 'EXPORT_MANIFEST.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({k: v for k, v in record.items() if k != 'files'}), flush=True)


if __name__ == '__main__':
    main()
