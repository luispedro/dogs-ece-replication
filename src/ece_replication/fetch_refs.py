"""Download and verify analysis inputs and reference databases.

Every file is recorded in ``data/manifest.tsv`` with its URL, the checksum
published by the source (when there is one), a SHA-256 computed locally, the
database version and the retrieval date.  Once a SHA-256 is recorded it is
treated as pinned: later runs verify against it and fail on mismatch.

Usage (always through pixi)::

    pixi run fetch-refs            # download whatever is missing, then verify
    pixi run python -m ece_replication.fetch_refs --verify-only
"""
import argparse
import csv
import datetime
import gzip
import io
import multiprocessing
import os
import shutil
import ssl
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.request
from pathlib import Path

from .util import file_digest, load_config, project_path

ZENODO = "https://zenodo.org/api/records/16356977/files/{}/content"
PFAM_RELEASE = "38.2"
PFAM = f"https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{PFAM_RELEASE}/{{}}"
PHROG_BASE = "https://phrogs.lmge.uca.fr/downloads_from_website/{}"
UNIPROT = "https://rest.uniprot.org/uniprotkb/{}.fasta"
EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=nuccore&rettype=fasta&retmode=text&id={}")

MANIFEST_FIELDS = ["name", "kind", "path", "url", "source_md5", "sha256",
                   "version", "retrieved", "note"]

# name -> (kind, relative path, url, source md5, version, note)
STATIC_ITEMS = {
    "SHD1_EC.fna.gz": ("input", "data/raw/SHD1_EC.fna.gz", ZENODO.format("SHD1_EC.fna.gz"),
                       "d9869388c1779dfc35a05e2d76d7e895", "zenodo:16356977", "185 circular ECEs"),
    "SHD1_EC_props.tsv.gz": ("input", "data/raw/SHD1_EC_props.tsv.gz", ZENODO.format("SHD1_EC_props.tsv.gz"),
                             "aab2c418539272cbdee60f652b51512f", "zenodo:16356977", "ECE metadata"),
    "Pfam-A.hmm.gz": ("reference", "data/refs/Pfam-A.hmm.gz", PFAM.format("Pfam-A.hmm.gz"),
                      "7ab3c4e215d0daaea3004e37c4e24f8a", f"Pfam {PFAM_RELEASE}", "profile HMMs"),
    "Pfam-A.hmm.dat.gz": ("reference", "data/refs/Pfam-A.hmm.dat.gz", PFAM.format("Pfam-A.hmm.dat.gz"),
                          "aab364689372b1846e8e880976e7f750", f"Pfam {PFAM_RELEASE}", "family metadata"),
    "phrog_annot_v4.tsv": ("reference", "data/refs/phrog_annot_v4.tsv", PHROG_BASE.format("phrog_annot_v4.tsv"),
                           "", "PHROGs v4", "PHROG functional annotation; server TLS certificate is invalid, "
                           "integrity relies on the pinned SHA-256"),
}


def _ssl_context(insecure):
    if not insecure:
        return None
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def download(url, dest, insecure=False, retries=3):
    dest = project_path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ece-replication/0.1"})
            with urllib.request.urlopen(req, context=_ssl_context(insecure), timeout=600) as r, open(tmp, "wb") as out:
                shutil.copyfileobj(r, out, 1 << 20)
                headers = dict(r.headers)
            tmp.rename(dest)
            return headers
        except Exception as e:  # pragma: no cover - network dependent
            print(f"download failed ({attempt + 1}/{retries}) {url}: {e}", file=sys.stderr)
            time.sleep(5)
    raise RuntimeError(f"could not download {url}")


def read_manifest(path=None):
    path = project_path(path or load_config()["paths"]["manifest"])
    if not path.exists():
        return {}
    with open(path) as f:
        return {r["name"]: r for r in csv.DictReader(f, delimiter="\t")}


def write_manifest(entries, path=None):
    path = project_path(path or load_config()["paths"]["manifest"])
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for name in sorted(entries, key=lambda n: (entries[n]["kind"], n)):
            w.writerow({k: entries[name].get(k, "") for k in MANIFEST_FIELDS})


def _record(entries, name, kind, path, url, source_md5, version, note):
    p = project_path(path)
    if source_md5:
        md5 = file_digest(p, "md5")
        if md5 != source_md5:
            raise ValueError(f"{name}: md5 {md5} != published {source_md5}")
    sha = file_digest(p)
    old = entries.get(name)
    if old and old.get("sha256") and old["sha256"] != sha:
        raise ValueError(f"{name}: sha256 changed ({old['sha256']} -> {sha}); refusing to overwrite pin")
    entries[name] = dict(name=name, kind=kind, path=str(path), url=url, source_md5=source_md5,
                         sha256=sha, version=version,
                         retrieved=(old or {}).get("retrieved") or datetime.date.today().isoformat(),
                         note=note)


def fetch_static(entries):
    for name, (kind, path, url, md5, version, note) in STATIC_ITEMS.items():
        if not project_path(path).exists():
            print(f"downloading {name}")
            download(url, path, insecure=url.startswith("https://phrogs"))
        _record(entries, name, kind, path, url, md5, version, note)


def _hmmbuild_one(args):
    msa, out = args
    name = Path(msa).stem
    subprocess.run(["hmmbuild", "--amino", "--informat", "afa", "-n", name, "--cpu", "1", "-o", os.devnull,
                    out, msa], check=True)
    return out


def fetch_phrog_hmms(entries, refs_dir, threads=16):
    """Build HMMER3 profiles from the PHROG v4 multiple sequence alignments.

    PHROGs distributes profiles only in HH-suite/MMseqs2 formats, so the
    HMMER profiles are rebuilt from the published MSAs (hmmbuild defaults).
    """
    out = Path(refs_dir) / "phrogs_v4.hmm.gz"
    url = PHROG_BASE.format("MSA_phrogs.tar.gz")
    entries.pop("HMM_phrog.tar.gz", None)
    if not project_path(out).exists():
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            tarball = Path(tmp) / "MSA_phrogs.tar.gz"
            print("downloading PHROG MSAs")
            download(url, tarball, insecure=True)
            tar_sha = file_digest(tarball)
            old = entries.get("MSA_phrogs.tar.gz")
            if old and old.get("sha256") and old["sha256"] != tar_sha:
                raise ValueError("MSA_phrogs.tar.gz checksum changed")
            entries["MSA_phrogs.tar.gz"] = dict(
                name="MSA_phrogs.tar.gz", kind="reference-source", path="(temporary; not kept)", url=url,
                source_md5="", sha256=tar_sha, version="PHROGs v4",
                retrieved=datetime.date.today().isoformat(),
                note="PHROG MSAs (aligned FASTA); input to hmmbuild. Server TLS certificate is invalid")
            with tarfile.open(tarball) as tar:
                tar.extractall(tmp, filter="data")
            msas = sorted(Path(tmp).rglob("*.fma"), key=lambda p: int(p.stem.split("_")[1]))
            hmm_dir = Path(tmp) / "hmm"
            hmm_dir.mkdir()
            jobs = [(str(m), str(hmm_dir / (m.stem + ".hmm"))) for m in msas]
            with multiprocessing.Pool(threads) as pool:
                built = pool.map(_hmmbuild_one, jobs, chunksize=64)
            with gzip.open(project_path(out), "wt") as dst:
                for h in built:
                    dst.write(Path(h).read_text())
    hmmbuild_version = subprocess.run(["hmmbuild", "-h"], capture_output=True, text=True).stdout.splitlines()[1]
    _record(entries, "phrogs_v4.hmm.gz", "reference-derived", out, url, "", "PHROGs v4",
            f"HMMER3 profiles built from PHROG MSAs ({hmmbuild_version.strip('# ').strip()}), ordered by PHROG id")


def fetch_reference_proteins(entries, cfg):
    table = project_path(cfg["references"]["reference_proteins_table"])
    out = project_path(cfg["references"]["reference_proteins_fasta"])
    release = entries.get("reference_proteins.faa", {}).get("version", "")
    if not out.exists():
        rows = list(csv.DictReader(open(table), delimiter="\t"))
        buf = io.StringIO()
        for r in rows:
            req = urllib.request.Request(UNIPROT.format(r["accession"]))
            with urllib.request.urlopen(req, timeout=120) as resp:
                text = resp.read().decode()
                release = resp.headers.get("X-UniProt-Release", release)
            seq = "".join(text.split("\n", 1)[1].split())
            if not seq:
                raise ValueError(f"empty sequence for {r['accession']}")
            buf.write(f">{r['accession']}\n{seq}\n")
            time.sleep(0.2)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(buf.getvalue())
    _record(entries, "reference_proteins.faa", "reference-derived", cfg["references"]["reference_proteins_fasta"],
            UNIPROT.format("<accession>"), "", f"UniProt {release}".strip(),
            f"curated characterized proteins listed in {cfg['references']['reference_proteins_table']} "
            f"(sha256 of table {file_digest(table)[:16]})")


def fetch_validation_genomes(entries, cfg):
    table = project_path(cfg["references"]["validation_table"])
    out = project_path(cfg["references"]["validation_fasta"])
    if not out.exists():
        rows = list(csv.DictReader(open(table), delimiter="\t"))
        buf = io.StringIO()
        for r in rows:
            with urllib.request.urlopen(EFETCH.format(r["accession"]), timeout=120) as resp:
                text = resp.read().decode()
            seq = "".join(text.split("\n", 1)[1].split()).upper()
            if not seq:
                raise ValueError(f"empty sequence for {r['accession']}")
            buf.write(f">{r['element']} {r['accession']}\n{seq}\n")
            time.sleep(0.4)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(buf.getvalue())
    _record(entries, "validation_elements.fna", "reference-derived", cfg["references"]["validation_fasta"],
            EFETCH.format("<accession>"), "", "NCBI nuccore (versioned accessions)",
            f"characterized replicons listed in {cfg['references']['validation_table']}")


def verify(manifest=None):
    """Return {name: sha256} after checking every pinned file. Raises on mismatch/missing."""
    entries = manifest if manifest is not None else read_manifest()
    if not entries:
        raise FileNotFoundError("manifest is empty; run `pixi run fetch-refs` first")
    result = {}
    for name, e in sorted(entries.items()):
        if e["path"].startswith("("):
            result[name] = e["sha256"]
            continue
        p = project_path(e["path"])
        if not p.exists():
            raise FileNotFoundError(f"{name}: {p} missing; run `pixi run fetch-refs`")
        sha = file_digest(p)
        if sha != e["sha256"]:
            raise ValueError(f"{name}: sha256 mismatch ({sha} != pinned {e['sha256']})")
        result[name] = sha
    return result


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify-only", action="store_true")
    args = ap.parse_args(argv)
    cfg = load_config()
    entries = read_manifest()
    if not args.verify_only:
        fetch_static(entries)
        fetch_phrog_hmms(entries, cfg["paths"]["refs_dir"])
        fetch_reference_proteins(entries, cfg)
        fetch_validation_genomes(entries, cfg)
        write_manifest(entries)
    for name, sha in verify(entries).items():
        print(f"OK\t{name}\t{sha[:16]}")


if __name__ == "__main__":
    main()
