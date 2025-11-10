#!/usr/bin/env python3
"""
Caching module for def generation to avoid expensive dumpbin calls.

This module provides caching functionality for object file symbol extraction.
It caches both file metadata (hash or mtime) and the extracted symbols from dumpbin.
"""

import json
import hashlib
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional


class DumpbinCache:
    """Cache for dumpbin symbol extraction results."""

    def __init__(self, cache_file: Path):
        """Initialize the cache.

        Args:
            cache_file: Path to the JSON cache file
        """
        self.cache_file = cache_file
        self.cache: Dict[str, Dict] = {}
        self.modified = False
        self._load_cache()

    def _load_cache(self):
        """Load the cache from disk."""
        if self.cache_file.exists():
            try:
                with open(self.cache_file, 'r', encoding='utf-8') as f:
                    self.cache = json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                print(f"Warning: Failed to load cache from {self.cache_file}: {e}")
                self.cache = {}

    def save(self):
        """Save the cache to disk if it was modified."""
        if self.modified:
            try:
                # Ensure parent directory exists
                self.cache_file.parent.mkdir(parents=True, exist_ok=True)
                with open(self.cache_file, 'w', encoding='utf-8') as f:
                    json.dump(self.cache, f, indent=2)
            except IOError as e:
                print(f"Warning: Failed to save cache to {self.cache_file}: {e}")

    def _get_file_hash(self, file_path: Path) -> Optional[str]:
        """Compute SHA256 hash of a file.

        Args:
            file_path: Path to the file

        Returns:
            Hex string of the hash, or None if file doesn't exist
        """
        try:
            sha256 = hashlib.sha256()
            with open(file_path, 'rb') as f:
                # Read in chunks to handle large files
                for chunk in iter(lambda: f.read(8192), b''):
                    sha256.update(chunk)
            return sha256.hexdigest()
        except IOError:
            return None

    def _get_file_mtime(self, file_path: Path) -> Optional[float]:
        """Get file modification time.

        Args:
            file_path: Path to the file

        Returns:
            Modification time as float, or None if file doesn't exist
        """
        try:
            return file_path.stat().st_mtime
        except OSError:
            return None

    def get_cached_symbols(self, obj_file: Path, use_hash: bool = False) -> Optional[List[str]]:
        """Get cached symbols for an object file if cache is valid.

        Args:
            obj_file: Path to the object file
            use_hash: If True, use file hash for validation; if False, use mtime

        Returns:
            List of symbols if cache is valid, None otherwise
        """
        obj_key = str(obj_file.resolve())

        if obj_key not in self.cache:
            return None

        cached_entry = self.cache[obj_key]

        # Validate cache based on method
        if use_hash:
            current_hash = self._get_file_hash(obj_file)
            if current_hash is None or cached_entry.get('hash') != current_hash:
                return None
        else:
            current_mtime = self._get_file_mtime(obj_file)
            if current_mtime is None or cached_entry.get('mtime') != current_mtime:
                return None

        return cached_entry.get('symbols')

    def cache_symbols(self, obj_file: Path, symbols: List[str], use_hash: bool = False):
        """Cache symbols for an object file.

        Args:
            obj_file: Path to the object file
            symbols: List of symbols extracted from the file
            use_hash: If True, use file hash for validation; if False, use mtime
        """
        obj_key = str(obj_file.resolve())

        entry = {
            'symbols': symbols
        }

        if use_hash:
            file_hash = self._get_file_hash(obj_file)
            if file_hash is not None:
                entry['hash'] = file_hash
        else:
            mtime = self._get_file_mtime(obj_file)
            if mtime is not None:
                entry['mtime'] = mtime

        self.cache[obj_key] = entry
        self.modified = True

    def extract_symbols_with_cache(self, obj_files: List[Path],
                                   extract_fn,
                                   use_hash: bool = False,
                                   show_progress: bool = True) -> List[str]:
        """Extract symbols from object files using cache when possible.

        Args:
            obj_files: List of object file paths
            extract_fn: Function that takes a list of object files and returns symbols
                       This function should call dumpbin and return raw dumpbin output lines
            use_hash: If True, use file hash for cache validation; if False, use mtime
            show_progress: If True, show progress information

        Returns:
            List of all symbols extracted
        """
        all_symbols = []
        cache_hits = 0
        cache_misses = 0
        files_to_process = []

        # Check cache for each file
        for obj_file in obj_files:
            cached = self.get_cached_symbols(obj_file, use_hash)
            if cached is not None:
                all_symbols.extend(cached)
                cache_hits += 1
            else:
                files_to_process.append(obj_file)
                cache_misses += 1

        if show_progress:
            print(f"Cache: {cache_hits} hits, {cache_misses} misses ({len(files_to_process)} files to process)")

        # Process uncached files
        if files_to_process:
            if show_progress:
                print(f"Running dumpbin on {len(files_to_process)} object files...")

            # Process files and cache results individually
            for obj_file in files_to_process:
                symbols = extract_fn([obj_file])
                all_symbols.extend(symbols)
                self.cache_symbols(obj_file, symbols, use_hash)

        return all_symbols


def run_dumpbin_for_file(obj_file: Path) -> Tuple[bool, List[str]]:
    """Run dumpbin on a single object file and return raw output lines.

    Args:
        obj_file: Path to the object file

    Returns:
        Tuple of (success, lines) where lines are the raw output lines
    """
    try:
        result = subprocess.run(
            ["dumpbin", "/SYMBOLS", str(obj_file)],
            capture_output=True,
            text=True,
            check=True,
        )
        return True, result.stdout.splitlines()
    except FileNotFoundError:
        print("Error: dumpbin.exe not found. Make sure MSVC tools are in PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError:
        return False, []

