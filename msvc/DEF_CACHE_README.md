# Def Generation Caching

## Overview

The def generation scripts now include caching support to significantly speed up the symbol extraction process by avoiding expensive `dumpbin` calls on object files that haven't changed.

## How It Works

Each def generation script (`def_gen_c.py`, `def_gen_cpp.py`, `def_gen_fc.py`, `def_gen_bibfor_ext.py`, `def_gen_astergc.py`) now uses a JSON cache file that stores:

- **File metadata**: Either the modification time (`mtime`) or SHA256 hash of each object file
- **Extracted symbols**: The raw symbols extracted from `dumpbin` for that file (before filtering)

When a script runs:
1. It checks the cache for each object file
2. If the file hasn't changed (based on mtime/hash), it uses the cached symbols
3. If the file has changed or isn't in the cache, it runs `dumpbin` and updates the cache
4. The cache is automatically saved after processing

## Performance Impact

For a typical incremental build where only a few files have changed:
- **Without cache**: Runs `dumpbin` on all ~9800+ object files (very slow)
- **With cache**: Runs `dumpbin` only on changed files (much faster)

Example from bibc (39 files):
```
First run:  Cache: 0 hits, 39 misses (all files processed)
Second run: Cache: 39 hits, 0 misses (instant, no dumpbin calls)
```

## Cache Files

Cache files are stored in the same directory as the output `.def` file:

- `.bibc_defgen_cache.json` - for bibc library
- `.bibcxx_defgen_cache.json` - for bibcxx library  
- `.bibfor_defgen_cache.json` - for bibfor library
- `.bibfor_ext_defgen_cache.json` - for bibfor_ext library
- `.astergc_defgen_cache.json` - for asterGC library

These files are automatically ignored by git (see `msvc/.gitignore`).

## Command Line Options

All def generation scripts support the following options:

### `--cache <path>`
Specify a custom cache file location. If not provided, uses the default location next to the output file.

```bash
python msvc/def_gen_c.py --cache my_custom_cache.json
```

### `--no-cache`
Disable caching entirely (always run dumpbin on all files).

```bash
python msvc/def_gen_c.py --no-cache
```

### `--use-hash`
Use file SHA256 hash instead of modification time for cache validation. This is slower but more reliable if you're concerned about modification time issues (e.g., files restored from backup, clock skew).

```bash
python msvc/def_gen_c.py --use-hash
```

## When to Clear the Cache

The cache should be automatically invalidated when files change. However, you may want to manually delete the cache files if:

1. You change the symbol filtering logic in the def generation scripts
2. You suspect the cache is corrupted
3. You want to force a complete rebuild of the def file

To clear all caches:
```bash
del msvc\.*.json
```

## Integration with Build System

The caching is transparent to the build system. The WAF tasks continue to work as before, but now benefit from the cache automatically.

## Technical Details

The caching implementation is in `msvc/def_gen_cache.py` and provides:

- **`DumpbinCache`** class: Manages loading, validating, and saving the cache
- **`run_dumpbin_for_file()`**: Wrapper for running dumpbin on a single file
- **`extract_symbols_with_cache()`**: Orchestrates the cache-aware symbol extraction

The cache stores symbols as a list of tuples for the C library (to distinguish DATA symbols) and as simple lists for other libraries.

