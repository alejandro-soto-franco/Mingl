"""Snakemake script: download one manuscript dataset and verify its checksum.

Shared by the intestine and melanoma download rules. Uses an exclusive flock
around the shared download directory so concurrent HickeyLab-fork workflows
on this machine cannot race the same file (see README/workflow spec).

Snakemake variables:
    params.url, params.sha256, params.description
    output.raw

If `params.sha256` is falsy (no checksum recorded in config yet), the
downloaded file's SHA-256 is printed so it can be copied into config.yaml,
and the rule fails deliberately rather than silently trusting an
unverified file.
"""

import fcntl
import hashlib
import sys
import urllib.request
from pathlib import Path

snakemake = globals()["snakemake"]

url = snakemake.params.url
expected_sha256 = snakemake.params.sha256
description = snakemake.params.description
out_path = Path(snakemake.output.raw)
out_path.parent.mkdir(parents=True, exist_ok=True)

lock_path = out_path.parent / f".{out_path.name}.lock"
with open(lock_path, "w") as lock_fh:
    fcntl.flock(lock_fh, fcntl.LOCK_EX)

    if out_path.exists():
        print(f"{out_path} already present, reusing it instead of re-downloading.", file=sys.stderr)
    else:
        print(f"Downloading {description}\n  {url}\n  -> {out_path}", file=sys.stderr)
        try:
            urllib.request.urlretrieve(url, out_path)
        except Exception as exc:
            raise RuntimeError(
                f"Could not download {description} from {url}: {exc}. "
                "This dataset's programmatic download was blocked (HTTP 403 from Dryad's "
                "bot protection) when this workflow was last checked from this machine; see "
                "README.md's 'Data availability' section for a manual-download workaround."
            ) from exc

    digest = hashlib.sha256(out_path.read_bytes()).hexdigest()

if not expected_sha256:
    out_path.unlink()
    raise ValueError(
        f"No sha256 recorded in config.yaml for {description}. Downloaded file's SHA-256 is:\n"
        f"  {digest}\n"
        "Copy this into config.yaml's `sha256` field for this dataset and re-run."
    )

if digest != expected_sha256:
    out_path.unlink()
    raise ValueError(f"SHA-256 mismatch for {description}: expected {expected_sha256}, got {digest}")
