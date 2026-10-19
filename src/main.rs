use anyhow::{anyhow, bail, Context, Result};
use tagporter::{collect_tags, write_new_bam, ConflictMode, IndexFormat};
use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Aligned BAM/CRAM to add tags to
    #[arg(long)]
    aln: PathBuf,

    /// Source BAM files containing the tags to transfer
    #[arg(long, required = true, value_parser, num_args = 1.., value_delimiter = ' ')]
    source_bams: Vec<PathBuf>,

    /// Comma-separated list of 2-character tag names to copy (e.g., MM,ML,RG)
    #[arg(long, value_delimiter = ',', default_value = "MM,ML")]
    tags: Vec<String>,

    /// Output BAM/CRAM/SAM path (format auto-detected from extension)
    #[arg(long)]
    output: PathBuf,

    /// Number of threads for BAM I/O
    #[arg(long, default_value_t = 1)]
    threads: usize,

    /// Skip building an index for the output file
    #[arg(long, default_value_t = false)]
    no_index: bool,

    /// Index format to build
    #[arg(long, value_enum, default_value = "bai")]
    index_format: IndexFormat,

    /// Behavior when a tag already exists on an aligned record
    #[arg(long, value_enum, default_value = "skip")]
    on_conflict: ConflictMode,

    /// Reference FASTA file (required when reading or writing CRAM)
    #[arg(long)]
    reference: Option<PathBuf>,

    /// Sled database cache capacity in MB
    #[arg(long, default_value_t = 4096)]
    cache_mb: usize,

    /// Directory for the temporary Sled database (e.g., /dev/shm on Linux for best performance)
    #[arg(long)]
    tmp_folder: Option<PathBuf>,

    /// Suppress all informational output (warnings and errors are still shown)
    #[arg(long, default_value_t = false)]
    quiet: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    if cli.quiet {
        log::set_max_level(log::LevelFilter::Warn);
    }

    log::debug!("{:#?}", &cli);

    for tag in &cli.tags {
        if tag.len() != 2 || !tag.is_ascii() {
            bail!(
                "Invalid tag name '{}': must be exactly 2 ASCII characters",
                tag
            );
        }
    }

    // Place the temp Sled DB next to --aln by default (same, usually large,
    // filesystem as the data), or in --tmp-folder if given. A uniquely-named
    // TempDir avoids collisions between concurrent runs and is removed on drop —
    // i.e. on every return path, including errors — so nothing is leaked.
    let tmp_parent: PathBuf = if let Some(tmp_folder) = cli.tmp_folder.clone() {
        tmp_folder
    } else {
        let parent = cli
            .aln
            .parent()
            .ok_or_else(|| anyhow!("Failed to get parent directory of '{}'", cli.aln.display()))?;
        if parent.as_os_str().is_empty() {
            PathBuf::from(".")
        } else {
            parent.to_path_buf()
        }
    };

    let tmp_dir = tempfile::Builder::new()
        .prefix("_tagporter_")
        .tempdir_in(&tmp_parent)
        .with_context(|| {
            format!("creating temporary directory in '{}'", tmp_parent.display())
        })?;
    let sled_db_path = tmp_dir.path();
    log::debug!("Sled DB path: {}", sled_db_path.display());

    log::info!(
        "Collecting tags {:?} from {} source file(s)",
        &cli.tags,
        cli.source_bams.len()
    );
    let t0 = std::time::Instant::now();
    collect_tags(
        &cli.source_bams,
        sled_db_path,
        cli.threads,
        &cli.tags,
        cli.cache_mb,
        cli.reference.as_deref(),
    )
    .context("collecting tags from source BAMs")?;
    let t_collect = t0.elapsed();

    log::info!("Transferring tags to aligned file...");
    let t1 = std::time::Instant::now();
    let stats = write_new_bam(
        sled_db_path,
        cli.threads,
        cli.aln,
        &cli.output,
        cli.no_index,
        cli.index_format,
        cli.on_conflict,
        cli.cache_mb,
        cli.reference.as_deref(),
    )
    .context("writing output file")?;
    let t_write = t1.elapsed();

    log::info!(
        "Phase timing: collect={:.1}s  write={:.1}s  total={:.1}s",
        t_collect.as_secs_f64(),
        t_write.as_secs_f64(),
        t0.elapsed().as_secs_f64()
    );
    log::info!(
        "Transfer complete: {} aligned reads processed",
        stats.total_aligned
    );
    if stats.total_aligned > 0 {
        let tagged_pct = stats.reads_tagged as f64 / stats.total_aligned as f64 * 100.0;
        log::info!(
            "  {} reads received tags ({:.1}%)",
            stats.reads_tagged,
            tagged_pct
        );
        let unmatched_pct = stats.reads_unmatched as f64 / stats.total_aligned as f64 * 100.0;
        log::info!(
            "  {} reads had no matching source record or tags were absent ({:.1}%)",
            stats.reads_unmatched,
            unmatched_pct
        );
    }
    if stats.reads_matched_no_new_tags > 0 {
        log::info!(
            "  {} matched read(s) received no new tag (all requested tags already present)",
            stats.reads_matched_no_new_tags
        );
    }
    if stats.tags_skipped_conflict > 0 {
        log::info!(
            "  {} tag(s) skipped due to conflicts (--on-conflict=skip)",
            stats.tags_skipped_conflict
        );
    }
    if stats.reads_secondary_supplementary > 0 {
        log::info!(
            "  {} secondary/supplementary alignment(s) written without tags",
            stats.reads_secondary_supplementary
        );
    }

    // Explicit cleanup so a removal failure is surfaced on the success path.
    // On any earlier error return, `tmp_dir`'s Drop still removes it best-effort.
    tmp_dir
        .close()
        .context("removing temporary Sled database directory")?;

    Ok(())
}
