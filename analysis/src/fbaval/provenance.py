"""Record which code version produced a result.

Every analysis writes ``provenance.json`` into its output folder with the commit of this repository, whether the
working tree had uncommitted changes, and the run parameters, so that each figure and table can be traced back to
the exact code that produced it.
"""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Optional

from fbaval.paths import REPO

DEFAULT_REPOS = {"Fvenenatum-GEM": REPO}
OUTPUT_DIRS = (":(exclude)results", ":(exclude)figures")


def _git(path: Path, *args: str) -> Optional[str]:
    try:
        out = subprocess.run(["git", "-C", str(path), *args],
                             capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    return out.stdout.strip()


def repo_state(path: Path) -> dict:
    """Commit hash and dirty flag of the git repository at *path* (commit is None if it is not one).

    The output folders ``results/`` and ``figures/`` are left out, so that writing results does not mark the code
    that produced them as changed.
    """
    path = Path(path)
    commit = _git(path, "rev-parse", "HEAD") if path.is_dir() else None
    if commit is None:
        return {"commit": None, "dirty": None}
    dirty = bool(_git(path, "status", "--porcelain", "--", ".", *OUTPUT_DIRS))
    return {"commit": commit, "dirty": dirty}


def write_provenance(out_dir: Path, repos: Optional[Mapping[str, Path]] = None,
                     params: Optional[dict] = None) -> Path:
    """Write ``provenance.json`` (repository commit + run parameters) into *out_dir*."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "repos": {name: repo_state(p) for name, p in (repos or DEFAULT_REPOS).items()},
        "params": params or {},
    }
    path = out_dir / "provenance.json"
    path.write_text(json.dumps(data, indent=2))
    return path
