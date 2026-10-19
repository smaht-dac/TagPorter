# TagPorter

A high-performance Rust tool for transferring auxiliary tags from source BAM files to aligned BAM files. Originally designed for methylation tags (MM, ML) but generalized to support any BAM auxiliary tags.

## Features

- **Fast**: Optimized for performance with batched database operations and direct tag access
- **Flexible**: Copy any combination of BAM auxiliary tags (MM, ML, RG, BC, UB, etc.)
- **Memory Efficient**: Uses Sled embedded database for intermediate storage, handling datasets larger than available RAM

## Installation

### Option 1: Download Pre-built Binaries (Recommended)

Pre-built binaries are available for Linux and macOS on the [Releases](https://github.com/smaht-dac/TagPorter/releases) page. They have no runtime dependencies beyond the operating system:

- **Linux** (x86_64, ARM64): glibc 2.28 or newer — e.g. RHEL/Rocky/Alma 8+, Ubuntu 20.04+, Debian 10+
- **macOS** (Apple Silicon, Intel): macOS 11 or newer

```bash
# Download the latest release for your platform
# Linux (x86_64)
wget https://github.com/smaht-dac/TagPorter/releases/latest/download/tagporter-linux-x86_64
chmod +x tagporter-linux-x86_64
sudo mv tagporter-linux-x86_64 /usr/local/bin/tagporter

# Linux (ARM64)
wget https://github.com/smaht-dac/TagPorter/releases/latest/download/tagporter-linux-aarch64
chmod +x tagporter-linux-aarch64
sudo mv tagporter-linux-aarch64 /usr/local/bin/tagporter

# macOS (Apple Silicon)
wget https://github.com/smaht-dac/TagPorter/releases/latest/download/tagporter-macos-aarch64
chmod +x tagporter-macos-aarch64
sudo mv tagporter-macos-aarch64 /usr/local/bin/tagporter

# macOS (Intel)
wget https://github.com/smaht-dac/TagPorter/releases/latest/download/tagporter-macos-x86_64
chmod +x tagporter-macos-x86_64
sudo mv tagporter-macos-x86_64 /usr/local/bin/tagporter
```

Each release also includes a `SHA256SUMS` file for verifying downloads (`sha256sum -c SHA256SUMS --ignore-missing`).

Verify installation:
```bash
tagporter --version
```

### Option 2: Building from Source

#### Prerequisites

- Rust 1.85 or higher (the crate uses the 2024 edition)
- Cargo (comes with Rust)

#### Build Instructions

```bash
git clone https://github.com/smaht-dac/TagPorter.git
cd TagPorter
cargo build --release
# Binary: target/release/tagporter
```

#### Running the Tests

Tests are in `tests/integration_test.rs` and use pre-built BAM fixtures from `tests/data/`. Because the C libraries (htslib, OpenSSL) are only compiled once for the release profile, always run tests with `--release`:

```bash
cargo test
# or
cargo test --release
```

To run a single test by name:

```bash
cargo test test_mm_ml_values_correct
# or
cargo test --release test_mm_ml_values_correct
```

## Usage

### Basic Command

```bash
tagporter \
  --aln aligned.bam \
  --source-bams unmapped1.bam unmapped2.bam \
  --output output.bam
```

### Command-Line Options

| Option | Description | Default | Required |
|--------|-------------|---------|----------|
| `--aln` | Aligned BAM/CRAM to add tags to | — | Yes |
| `--source-bams` | Space-separated source BAM files containing the tags to transfer | — | Yes |
| `--output` | Output path (format auto-detected from extension: `.bam`, `.cram`, `.sam`) | — | Yes |
| `--tags` | Comma-separated list of 2-character tag names to copy | `MM,ML` | No |
| `--threads` | Number of threads for BAM I/O | `1` | No |
| `--no-index` | Skip building an index for the output file | `false` | No |
| `--index-format` | Index format for BAM output: `bai` or `csi` (use `csi` for contigs >512 Mb). CRAM output always gets a `.crai` index; SAM output is not indexed | `bai` | No |
| `--on-conflict` | Behavior when a tag already exists on the aligned record: `skip`, `overwrite`, or `error` | `skip` | No |
| `--reference` | Reference FASTA file (required when reading or writing CRAM) | — | No |
| `--cache-mb` | Sled database cache capacity in MB | `4096` | No |
| `--tmp-folder` | Directory for the temporary Sled database (use `/dev/shm` on Linux for best performance) | Next to `--aln` | No |
| `--quiet` | Suppress informational output (warnings and errors still shown) | `false` | No |

### Examples

#### Copy methylation tags (MM/ML) — default behavior

```bash
tagporter \
  --aln aligned.bam \
  --source-bams unmapped.bam \
  --output output.bam
```

#### Copy custom tags

```bash
# Copy read group, barcode, and UMI tags
tagporter \
  --aln aligned.bam \
  --source-bams unmapped.bam \
  --output output.bam \
  --tags RG,BC,UB
```

#### Multiple source files with parallelized I/O

```bash
tagporter \
  --aln aligned.bam \
  --source-bams unaligned1.bam unaligned2.bam unaligned3.bam \
  --output output.bam \
  --threads 8
```

#### Large genome — CSI index (required for contigs >512 Mb)

```bash
tagporter \
  --aln aligned.bam \
  --source-bams unmapped.bam \
  --output output.bam \
  --index-format csi
```

#### CRAM output

```bash
tagporter \
  --aln aligned.bam \
  --source-bams unmapped.bam \
  --output output.cram \
  --reference genome.fa
```

#### Use a ramdisk for the temporary database (Linux)

```bash
tagporter \
  --aln aligned.bam \
  --source-bams unmapped.bam \
  --output output.bam \
  --tmp-folder /dev/shm
```

#### Silent mode for pipeline use

```bash
tagporter \
  --aln aligned.bam \
  --source-bams unmapped.bam \
  --output output.bam \
  --quiet
```

## How It Works

1. **Tag Collection Phase**: Reads source BAM files and extracts specified auxiliary tags
   - Stores tags in a Sled embedded key-value database (key: read name, value: serialized tags)
   - Uses batched writes of 10,000 records for optimal throughput
   - Progress updates logged every 100,000 reads

2. **Tag Transfer Phase**: Reads aligned BAM file and adds tags from the database
   - Looks up each read name in the database
   - Adds matching tags to the aligned record (conflict behavior controlled by `--on-conflict`)
   - Writes modified record to output file
   - Reports a transfer summary (reads tagged, unmatched, conflicts)

3. **Indexing Phase** (optional): Creates an index for the output file
   - BAM output gets a BAI index by default; use `--index-format csi` for genomes with contigs larger than 512 Mb
   - CRAM output always gets a `.crai` index (`--index-format` is ignored)
   - SAM output cannot be indexed and is skipped with a warning
   - Skip entirely with `--no-index`

## Verbosity

Progress output goes to stderr and is controlled via the `RUST_LOG` environment variable or `--quiet`:

```bash
RUST_LOG=warn tagporter ...   # suppress info, show only warnings/errors
RUST_LOG=debug tagporter ...  # show all debug output
tagporter --quiet ...         # equivalent to RUST_LOG=warn
```

## Troubleshooting

- **Out of disk space**: Ensure sufficient space for the temporary Sled database; use `--tmp-folder /dev/shm` on Linux to store it in RAM
- **Corrupted BAM files**: Verify input files with `samtools quickcheck`
- **Missing read names**: Tags are only transferred for reads present in both files; check the transfer summary for the unmatched count
- **Tag name mismatch**: Tag names are case-sensitive and must be exactly 2 ASCII characters (e.g., `MM` not `mm`)
- **CSI index needed**: If downstream tools report index errors on large genomes, rebuild with `--index-format csi`

