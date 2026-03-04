"""Sync feedstock recipes from fork repos into vendored directories.

Usage:
    python sync_feedstock.py <dep-name> [--branch <branch>]
    python sync_feedstock.py --all
"""

import argparse
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

DEPS_DIR = Path(__file__).resolve().parent

FEEDSTOCKS = {
    "libmed": ("git@github.com:Krande/libmed-feedstock.git", "feat/win-fortran"),
    "medcoupling": ("git@github.com:Krande/medcoupling-feedstock.git", "main"),
    "mumps": ("git@github.com:Krande/mumps-feedstock.git", "main"),
    "scotch": ("git@github.com:Krande/scotch-feedstock.git", "feat/rattler-build"),
}

# Directories to copy from each feedstock repo
SYNC_DIRS = ["recipe", ".ci_support"]


# Keys to strip from .ci_support yaml files (conflict with rattler-build -c flags)
_STRIP_KEYS = re.compile(r"^(channel_sources|channel_targets):.*\n(?:- .*\n)*", re.MULTILINE)


def _strip_channel_keys(ci_support_dir: Path):
    """Remove channel_sources/channel_targets from all yaml files in a directory."""
    if not ci_support_dir.exists():
        return
    for yaml_file in ci_support_dir.glob("*.yaml"):
        text = yaml_file.read_text()
        cleaned = _STRIP_KEYS.sub("", text)
        if cleaned != text:
            yaml_file.write_text(cleaned)


def sync_feedstock(name: str, branch: str | None = None):
    """Clone a feedstock repo (shallow) and copy recipe/ and .ci_support/ into deps/<name>/."""
    if name not in FEEDSTOCKS:
        raise ValueError(f"Unknown feedstock: {name!r}. Known: {list(FEEDSTOCKS)}")

    url, default_branch = FEEDSTOCKS[name]
    branch = branch or default_branch
    dest = DEPS_DIR / name

    with tempfile.TemporaryDirectory() as tmp:
        print(f"Cloning {name} ({url} @ {branch}) ...")
        subprocess.run(
            ["git", "clone", "--depth", "1", "--branch", branch, url, tmp],
            check=True,
        )

        for dirname in SYNC_DIRS:
            src = Path(tmp) / dirname
            if not src.exists():
                print(f"  warning: {dirname}/ not found in {name}, skipping")
                continue
            target = dest / dirname
            if target.exists():
                shutil.rmtree(target)
            shutil.copytree(src, target)
            print(f"  synced {dirname}/")

    # Strip channel_sources/channel_targets from .ci_support configs
    # (they conflict with rattler-build -c flags)
    _strip_channel_keys(dest / ".ci_support")

    print(f"Done: {name} -> {dest.relative_to(DEPS_DIR)}")


def init_local_channel():
    """Create bare repodata.json so rattler-build/pixi can solve before packages are built."""
    output_dir = DEPS_DIR / "output"
    for subdir in ("win-64", "noarch"):
        d = output_dir / subdir
        d.mkdir(parents=True, exist_ok=True)
        repodata = d / "repodata.json"
        if not repodata.exists():
            repodata.write_text(json.dumps({"info": {}, "packages": {}, "packages.conda": {}}))
            print(f"Initialized {repodata}")


def main():
    parser = argparse.ArgumentParser(description="Sync feedstock recipes from forks")
    parser.add_argument("name", nargs="?", help="Feedstock name to sync")
    parser.add_argument("--branch", help="Override the default branch")
    parser.add_argument("--all", action="store_true", help="Sync all feedstocks")
    args = parser.parse_args()

    if not args.all and not args.name:
        parser.error("Provide a feedstock name or --all")

    if args.all:
        for name in FEEDSTOCKS:
            sync_feedstock(name, args.branch)
    else:
        sync_feedstock(args.name, args.branch)

    init_local_channel()
    print("\nSync complete.")


if __name__ == "__main__":
    main()
