#!/usr/bin/env python3
"""Fetch only the public DTU training archives needed for pseudo references.

The DTU training-data pages publish the pseudo reference inside large tarballs.
This helper uses checked HTTP byte ranges so interrupted downloads can resume and
then extracts only the two small pseudo files.  It does not download or install
the separately licensed NetMHCpan/NetMHCIIpan executable packages.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import tarfile
import urllib.request
from pathlib import Path
from typing import Mapping


ARCHIVES: Mapping[str, Mapping[str, object]] = {
    "netmhcpan-4.2": {
        "url": (
            "https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/"
            "training_data/NetMHCpan_train.tar.gz"
        ),
        "size": 269075958,
        "member": "NetMHCpan_train/pseudoseqs",
        "output": "NetMHCpan-4.2-pseudoseqs.dat",
        "upstream_updated": "2026-07-06",
    },
    "netmhciipan-4.3": {
        "url": (
            "https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/"
            "training_data/NetMHCIIpan_train.tar.gz"
        ),
        "size": 139440263,
        "member": "NetMHCIIpan_train/pseudosequence.2023.dat",
        "output": "NetMHCIIpan-4.3-pseudosequence.2023.dat",
        "upstream_updated": "2026-07-06",
    },
}


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch_range(url: str, start: int, end: int, path: Path) -> None:
    expected = end - start + 1
    if path.is_file() and path.stat().st_size == expected:
        return
    request = urllib.request.Request(
        url,
        headers={
            "Range": f"bytes={start}-{end}",
            "User-Agent": "MHCtransform-reference-builder/0.4",
        },
    )
    temporary = path.with_suffix(path.suffix + ".tmp")
    with urllib.request.urlopen(request, timeout=300) as response:
        content_range = response.headers.get("Content-Range", "")
        if response.status != 206 or not content_range.startswith(
            f"bytes {start}-{end}/"
        ):
            raise RuntimeError(
                f"server did not honor byte range {start}-{end}: "
                f"status={response.status}, Content-Range={content_range!r}"
            )
        with temporary.open("wb") as handle:
            shutil.copyfileobj(response, handle, length=1024 * 1024)
    if temporary.stat().st_size != expected:
        raise RuntimeError(
            f"short range {start}-{end}: {temporary.stat().st_size} != {expected}"
        )
    os.replace(temporary, path)


def fetch_archive(
    url: str,
    size: int,
    destination: Path,
    *,
    jobs: int,
    chunk_size: int,
) -> None:
    if destination.is_file() and destination.stat().st_size == size:
        return
    parts_dir = destination.with_suffix(destination.suffix + ".parts")
    parts_dir.mkdir(parents=True, exist_ok=True)
    ranges = []
    for index, start in enumerate(range(0, size, chunk_size)):
        end = min(start + chunk_size, size) - 1
        ranges.append((index, start, end, parts_dir / f"{index:05d}.part"))
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = [
            executor.submit(fetch_range, url, start, end, path)
            for _, start, end, path in ranges
        ]
        for future in concurrent.futures.as_completed(futures):
            future.result()
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("wb") as output:
        for _, _, _, path in ranges:
            with path.open("rb") as handle:
                shutil.copyfileobj(handle, output, length=1024 * 1024)
    if temporary.stat().st_size != size:
        raise RuntimeError(f"assembled archive size mismatch: {temporary.stat().st_size}")
    os.replace(temporary, destination)


def extract_member(archive: Path, member: str, destination: Path) -> None:
    with tarfile.open(archive, "r:gz") as handle:
        extracted = handle.extractfile(member)
        if extracted is None:
            raise RuntimeError(f"archive member not found: {member}")
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        with temporary.open("wb") as output:
            shutil.copyfileobj(extracted, output)
        os.replace(temporary, destination)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=24)
    parser.add_argument("--chunk-mib", type=int, default=2)
    args = parser.parse_args()
    if args.jobs < 1 or args.chunk_mib < 1:
        parser.error("--jobs and --chunk-mib must be positive")
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    report = {}
    for key, spec in ARCHIVES.items():
        archive = args.cache_dir / f"{key}-training.tar.gz"
        fetch_archive(
            str(spec["url"]),
            int(spec["size"]),
            archive,
            jobs=args.jobs,
            chunk_size=args.chunk_mib * 1024 * 1024,
        )
        output = args.cache_dir / str(spec["output"])
        extract_member(archive, str(spec["member"]), output)
        report[key] = {
            "archive": str(archive),
            "archive_sha256": sha256_path(archive),
            "pseudo_file": str(output),
            "pseudo_sha256": sha256_path(output),
            "upstream_updated": spec["upstream_updated"],
            "url": spec["url"],
        }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
