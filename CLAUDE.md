# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`TagPorter` (binary: `tagporter`) is a Rust CLI tool that transfers auxiliary BAM file tags from unaligned BAM files to aligned BAM files. Originally built for methylation tags (MM, ML) but generalized to support any BAM auxiliary tags.

## Build & Run Commands

```bash
# Build release binary
cargo build --release

# Run the tool
./target/release/tagporter \
  --aln aligned.bam \
  --source-bams unaligned1.bam unaligned2.bam \
  --output output.bam \
  --tags MM,ML \
  --threads 8 \
  --tmp-folder /dev/shm   # optional: use ramdisk for temp Sled DB

# Run tests (always use --release — C libs are only pre-built for release)
cargo test --release

# Run a single test by name
cargo test --release test_mm_ml_values_correct

# Run benchmarks
cargo bench

# Check compilation without building
cargo check
```

Progress output goes to stderr via `env_logger`. Control verbosity with `RUST_LOG=warn/info/debug` or `--quiet`.

## Releasing

[.github/workflows/release.yml](.github/workflows/release.yml) builds `tagporter-{linux-x86_64,linux-aarch64,macos-aarch64,macos-x86_64}` and publishes them with a `SHA256SUMS` file.

- Bump `version` in `Cargo.toml`, commit, then `git tag vX.Y.Z && git push origin vX.Y.Z`. The workflow fails if the tag and `Cargo.toml` version differ.
- A suffix (`vX.Y.Z-rc1`) publishes a pre-release; the suffix is ignored for the version check.
- "Run workflow" in the Actions tab (`workflow_dispatch`) builds artifacts without releasing.
- Linux builds run in `manylinux_2_28` containers and fail if the binary needs a glibc newer than 2.28. macOS builds target macOS 11 (`MACOSX_DEPLOYMENT_TARGET`).

## Architecture

The tool uses a **three-phase pipeline** with an embedded [Sled](https://github.com/spacejam/sled) key-value database as intermediate storage (avoids loading large datasets into RAM):

1. **`collect_tags()`** — Reads N unaligned BAMs, extracts requested auxiliary tags, serializes them with `postcard`, and batch-writes to the Sled database. The key is the read name followed by a one-byte mate suffix (`record_key()`), so the two mates of a pair — which share a QNAME — do not collide. Secondary/supplementary source records are skipped so their clipped tags cannot overwrite the primary's. Batches of 10,000 records are written atomically.

2. **`write_new_bam()`** — Reads the aligned BAM sequentially, looks up each record's mate-aware key in Sled, deserializes matched tags, pushes them onto the record, and writes to the output BAM. Only primary alignments are tagged; secondary/supplementary alignments are written through unchanged (counted in `TransferStats::reads_secondary_supplementary`). CRAM I/O honors `--reference` via `set_reference()` on the reader/writer. Optionally builds a BAI/CSI index afterward (skipped for SAM output, which cannot be indexed).

3. **Cleanup** — The binary owns the Sled DB as a uniquely-named `tempfile::TempDir` (created next to `--aln`, or in `--tmp-folder`), so it is removed on every exit path — success or error — via RAII. The library functions never delete a caller-supplied directory. Unique names let concurrent runs share a directory without colliding.

### Key Data Structures

- **`TagData`** — `HashMap<[u8; 2], SerializableAux>` for flexible multi-tag storage per read.
- **`SerializableAux`** — Enum mirroring all 15+ `rust_htslib::bam::record::Aux` variants, enabling `serde`/`postcard` serialization. Conversion in/out via `convert_aux_to_serializable()` and `push_aux()`.

### Performance Design Decisions

- **Direct access fast path** (in `collect_tags()`, [src/lib.rs](src/lib.rs)): For ≤3 requested tags, uses `record.aux(tag_bytes)` direct lookups instead of iterating all auxiliary fields.
- **Sled `HighThroughput` mode** with configurable cache (`--cache-mb`, default 4096 MB) and manual flushes.
- BAM I/O threading delegated to htslib thread pools (`--threads` flag), not Rayon.

### Source Files

- [src/lib.rs](src/lib.rs) — **Library crate**: all business logic (`collect_tags`, `write_new_bam`, public types). Import as `tagporter::` in tests.
- [src/main.rs](src/main.rs) — **Binary entry point**: CLI definition (`Cli` struct, clap) and `main()` only.
- [tests/integration_test.rs](tests/integration_test.rs) — Integration tests against BAM fixtures in `tests/data/`.
- [benches/my_benchmark.rs](benches/my_benchmark.rs) — Criterion benchmark scaffold (placeholder).

## CLI Surface

```
Required:
  --aln <ALN>               Aligned BAM/CRAM to add tags to
  --source-bams <BAMS>...   Source BAMs containing the tags to transfer (alias: --methyl-bams)
  --output <OUTPUT>         Output path; format auto-detected from extension (.bam/.cram/.sam)

Optional:
  --tags <TAGS>             Comma-separated 2-char tag names [default: MM,ML]
  --threads <N>             htslib thread count [default: 1]
  --no-index                Skip index building
  --index-format bai|csi    Index type — use csi for contigs >512 Mb [default: bai]
  --on-conflict skip|overwrite|error  Tag conflict behavior [default: skip]
  --reference <FASTA>       Reference FASTA (required for CRAM input/output)
  --cache-mb <MB>           Sled cache in MB [default: 4096]
  --tmp-folder <DIR>        Temp dir for Sled DB (e.g. /dev/shm)
  --quiet                   Suppress informational output
```

## Dependencies

| Crate | Purpose |
|-------|---------|
| `rust-htslib` | BAM/SAM/CRAM file I/O (htslib bindings) |
| `sled` | Embedded key-value database for intermediate tag storage |
| `clap` (derive) | CLI argument parsing |
| `postcard` + `serde` | Binary serialization of tag data for Sled |
| `tempfile` | Self-cleaning temp directory for the Sled DB |
| `anyhow` | Error propagation with context |
| `log` + `env_logger` | Structured logging to stderr |
| `criterion` (dev) | Benchmarking framework |
