import gzip
import hashlib
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def project_path(p):
    p = Path(p)
    return p if p.is_absolute() else PROJECT_ROOT / p


def load_yaml(path):
    with open(project_path(path)) as f:
        return yaml.safe_load(f)


def load_config(path="config/analysis.yaml"):
    return load_yaml(path)


def open_maybe_gz(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def file_digest(path, algorithm="sha256", chunk=1 << 20):
    h = hashlib.new(algorithm)
    with open(project_path(path), "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def read_fasta(path):
    """Yield (identifier, description, sequence) tuples."""
    name, desc, chunks = None, "", []
    with open_maybe_gz(project_path(path)) as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(chunks)
                header = line[1:].strip()
                parts = header.split(None, 1)
                name = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, desc, "".join(chunks)


def write_fasta(records, path, width=80):
    path = project_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")
    return str(path)
