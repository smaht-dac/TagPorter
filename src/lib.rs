use anyhow::{bail, Context, Result};
use rust_htslib::bam::index;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};
use serde::{Deserialize, Serialize};
use sled;
use std::collections::HashMap;
use std::collections::HashSet;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};

// ── internal serialization types ─────────────────────────────────────────────

#[derive(Serialize, Deserialize, Debug)]
struct TagData {
    tags: HashMap<[u8; 2], SerializableAux>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
enum SerializableAux {
    String(String),
    Char(u8),
    I8(i8),
    U8(u8),
    I16(i16),
    U16(u16),
    I32(i32),
    U32(u32),
    Float(f32),
    Double(f64),
    ArrayU8(Vec<u8>),
    ArrayU16(Vec<u16>),
    ArrayU32(Vec<u32>),
    ArrayI8(Vec<i8>),
    ArrayI16(Vec<i16>),
    ArrayI32(Vec<i32>),
    ArrayFloat(Vec<f32>),
    HexByteArray(Vec<u8>),
}

// ── public API types ──────────────────────────────────────────────────────────

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum IndexFormat {
    Bai,
    Csi,
}

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum ConflictMode {
    /// Keep the existing tag on the aligned record; skip the source value
    Skip,
    /// Replace the existing tag with the source value
    Overwrite,
    /// Abort with an error on the first conflict
    Error,
}

/// Counts reported after a completed tag-transfer run. For a completed run the
/// four outcome buckets partition `total_aligned` exactly:
/// `reads_tagged + reads_unmatched + reads_matched_no_new_tags
/// + reads_secondary_supplementary == total_aligned`.
#[derive(Debug)]
pub struct TransferStats {
    pub total_aligned: u64,
    /// Primary reads that received at least one new tag.
    pub reads_tagged: u64,
    /// Primary reads with no matching source record in the DB.
    pub reads_unmatched: u64,
    /// Individual tags skipped because they already existed (`--on-conflict=skip`).
    pub tags_skipped_conflict: u64,
    /// Primary reads matched in the DB but which received no new tag (every
    /// candidate tag already existed and was skipped).
    pub reads_matched_no_new_tags: u64,
    /// Secondary/supplementary alignments written unchanged (never tag-eligible).
    pub reads_secondary_supplementary: u64,
}

// ── public functions ──────────────────────────────────────────────────────────

/// Read all source BAMs and store the requested tags in a Sled database keyed
/// by read name.
pub fn collect_tags(
    source_bams: &[PathBuf],
    db_path: &Path,
    threads: usize,
    tags_to_collect: &[String],
    cache_mb: usize,
    reference: Option<&Path>,
) -> Result<()> {
    let db = sled::Config::default()
        .path(db_path)
        .cache_capacity((cache_mb * 1024 * 1024) as u64)
        .flush_every_ms(None)
        .mode(sled::Mode::HighThroughput)
        .open()
        .with_context(|| format!("opening Sled DB at '{}'", db_path.display()))?;

    let tag_bytes: Vec<[u8; 2]> = tags_to_collect
        .iter()
        .filter(|t| t.len() == 2)
        .map(|t| {
            let b = t.as_bytes();
            [b[0], b[1]]
        })
        .collect();

    let tags_set: HashSet<[u8; 2]> = tag_bytes.iter().copied().collect();

    let batch_size = 10_000;
    let use_direct_access = tag_bytes.len() <= 3;

    for path in source_bams {
        let mut bam_reader = bam::Reader::from_path(path)
            .with_context(|| format!("opening source BAM '{}'", path.display()))?;
        bam_reader.set_threads(threads)?;
        if let Some(reference) = reference {
            bam_reader
                .set_reference(reference)
                .with_context(|| format!("setting reference for '{}'", path.display()))?;
        }

        let mut batch_data: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(batch_size);
        let mut total_reads: u64 = 0;

        // Allocate once outside the loop and reuse via .clear() to avoid O(n) allocations.
        let mut found_tags: HashMap<[u8; 2], SerializableAux> =
            HashMap::with_capacity(tag_bytes.len());

        for result in bam_reader.records() {
            let record = result
                .with_context(|| format!("reading record from '{}'", path.display()))?;
            total_reads += 1;

            // Only primary alignments hold the canonical full-length read; skip
            // secondary/supplementary so their clipped tags can't overwrite it.
            if record.is_secondary() || record.is_supplementary() {
                continue;
            }

            found_tags.clear();

            if use_direct_access {
                for tag_arr in &tag_bytes {
                    if let Ok(value) = record.aux(tag_arr) {
                        found_tags.insert(*tag_arr, convert_aux_to_serializable(value));
                    }
                }
            } else {
                for aux in record.aux_iter() {
                    if found_tags.len() == tag_bytes.len() {
                        break;
                    }
                    if let Ok((tag, value)) = aux {
                        if tag.len() == 2 {
                            let tag_arr: [u8; 2] = [tag[0], tag[1]];
                            if tags_set.contains(&tag_arr) {
                                found_tags.insert(tag_arr, convert_aux_to_serializable(value));
                            }
                        }
                    }
                }
            }

            if !found_tags.is_empty() {
                let tag_data = TagData {
                    tags: found_tags.clone(),
                };
                let encoded =
                    postcard::to_allocvec(&tag_data).context("serializing tag data")?;
                batch_data.push((record_key(&record), encoded));

                if batch_data.len() >= batch_size {
                    let mut batch = sled::Batch::default();
                    for (key, value) in batch_data.drain(..) {
                        batch.insert(key, value);
                    }
                    db.apply_batch(batch).context("writing batch to Sled DB")?;
                }
            }

            if total_reads % 100_000 == 0 {
                log::info!("Collected {} reads from {}", total_reads, path.display());
            }
        }

        if !batch_data.is_empty() {
            let mut batch = sled::Batch::default();
            for (key, value) in batch_data {
                batch.insert(key, value);
            }
            db.apply_batch(batch).context("writing final batch to Sled DB")?;
        }

        log::info!(
            "Finished collecting from {}: {} reads processed",
            path.display(),
            total_reads
        );
    }

    db.flush().context("flushing Sled DB")?;
    Ok(())
}

/// Transfer tags from the Sled database to an aligned BAM/CRAM/SAM file,
/// writing results to `output`. Returns transfer statistics.
pub fn write_new_bam(
    db_path: &Path,
    threads: usize,
    aligned_bam: PathBuf,
    output: &PathBuf,
    no_index: bool,
    index_format: IndexFormat,
    on_conflict: ConflictMode,
    cache_mb: usize,
    reference: Option<&Path>,
) -> Result<TransferStats> {
    let db = sled::Config::default()
        .path(db_path)
        .cache_capacity((cache_mb * 1024 * 1024) as u64)
        .flush_every_ms(None)
        .mode(sled::Mode::HighThroughput)
        .open()
        .with_context(|| format!("opening Sled DB at '{}'", db_path.display()))?;

    let mut reader = bam::Reader::from_path(&aligned_bam)
        .with_context(|| format!("opening aligned file '{}'", aligned_bam.display()))?;
    reader.set_threads(threads)?;
    if let Some(reference) = reference {
        reader.set_reference(reference).with_context(|| {
            format!("setting reference for '{}'", aligned_bam.display())
        })?;
    }
    let header = bam::Header::from_template(reader.header());

    let out_fmt = match output.extension().and_then(OsStr::to_str) {
        Some("cram") => bam::Format::Cram,
        Some("sam") => bam::Format::Sam,
        _ => bam::Format::Bam,
    };

    let mut writer = bam::Writer::from_path(output, &header, out_fmt)
        .with_context(|| format!("creating output file '{}'", output.display()))?;
    writer.set_threads(threads)?;
    if let Some(reference) = reference {
        writer.set_reference(reference).with_context(|| {
            format!("setting reference for output '{}'", output.display())
        })?;
    }

    let mut stats = TransferStats {
        total_aligned: 0,
        reads_tagged: 0,
        reads_unmatched: 0,
        tags_skipped_conflict: 0,
        reads_matched_no_new_tags: 0,
        reads_secondary_supplementary: 0,
    };

    for result in reader.records() {
        let mut record = result.context("reading aligned record")?;
        stats.total_aligned += 1;

        // Tags are transferred only to primary alignments; secondary/supplementary
        // records share the read name but carry clipped/reverse sequence. Write
        // them through unchanged.
        if record.is_secondary() || record.is_supplementary() {
            stats.reads_secondary_supplementary += 1;
            writer
                .write(&record)
                .context("writing record to output file")?;
            continue;
        }

        if let Some(bytes) = db
            .get(record_key(&record))
            .context("looking up read name in Sled DB")?
        {
            let data: TagData =
                postcard::from_bytes(&bytes).context("deserializing tag data from Sled DB")?;
            let mut got_tag = false;

            for (tag_array, tag_value) in data.tags {
                let already_exists = record.aux(&tag_array).is_ok();
                if already_exists {
                    match on_conflict {
                        ConflictMode::Skip => {
                            stats.tags_skipped_conflict += 1;
                            continue;
                        }
                        ConflictMode::Overwrite => {
                            let _ = record.remove_aux(&tag_array);
                        }
                        ConflictMode::Error => {
                            bail!(
                                "Tag '{}' already exists on read '{}' \
                                 (use --on-conflict=skip or --on-conflict=overwrite)",
                                std::str::from_utf8(&tag_array).unwrap_or("??"),
                                String::from_utf8_lossy(record.qname())
                            );
                        }
                    }
                }

                let push_result = match tag_value {
                    SerializableAux::String(s) => record.push_aux(&tag_array, Aux::String(&s)),
                    SerializableAux::Char(c) => record.push_aux(&tag_array, Aux::Char(c)),
                    SerializableAux::I8(v) => record.push_aux(&tag_array, Aux::I8(v)),
                    SerializableAux::U8(v) => record.push_aux(&tag_array, Aux::U8(v)),
                    SerializableAux::I16(v) => record.push_aux(&tag_array, Aux::I16(v)),
                    SerializableAux::U16(v) => record.push_aux(&tag_array, Aux::U16(v)),
                    SerializableAux::I32(v) => record.push_aux(&tag_array, Aux::I32(v)),
                    SerializableAux::U32(v) => record.push_aux(&tag_array, Aux::U32(v)),
                    SerializableAux::Float(v) => record.push_aux(&tag_array, Aux::Float(v)),
                    SerializableAux::Double(v) => record.push_aux(&tag_array, Aux::Double(v)),
                    SerializableAux::ArrayU8(v) => {
                        record.push_aux(&tag_array, Aux::ArrayU8(v.as_slice().into()))
                    }
                    SerializableAux::ArrayU16(v) => {
                        record.push_aux(&tag_array, Aux::ArrayU16(v.as_slice().into()))
                    }
                    SerializableAux::ArrayU32(v) => {
                        record.push_aux(&tag_array, Aux::ArrayU32(v.as_slice().into()))
                    }
                    SerializableAux::ArrayI8(v) => {
                        record.push_aux(&tag_array, Aux::ArrayI8(v.as_slice().into()))
                    }
                    SerializableAux::ArrayI16(v) => {
                        record.push_aux(&tag_array, Aux::ArrayI16(v.as_slice().into()))
                    }
                    SerializableAux::ArrayI32(v) => {
                        record.push_aux(&tag_array, Aux::ArrayI32(v.as_slice().into()))
                    }
                    SerializableAux::ArrayFloat(v) => {
                        record.push_aux(&tag_array, Aux::ArrayFloat(v.as_slice().into()))
                    }
                    SerializableAux::HexByteArray(v) => {
                        let hex_string =
                            v.iter().map(|b| format!("{:02X}", b)).collect::<String>();
                        record.push_aux(&tag_array, Aux::HexByteArray(&hex_string))
                    }
                };
                push_result.with_context(|| {
                    format!(
                        "pushing tag '{}'",
                        std::str::from_utf8(&tag_array).unwrap_or("??")
                    )
                })?;
                got_tag = true;
            }

            if got_tag {
                stats.reads_tagged += 1;
            } else {
                // Matched in the DB but every candidate tag already existed and
                // was skipped, so nothing new was applied.
                stats.reads_matched_no_new_tags += 1;
            }
        } else {
            stats.reads_unmatched += 1;
        }

        if stats.total_aligned % 100_000 == 0 {
            log::info!("Processed {} aligned records", stats.total_aligned);
        }

        writer
            .write(&record)
            .context("writing record to output file")?;
    }
    drop(writer);

    if !no_index {
        if matches!(out_fmt, bam::Format::Sam) {
            log::warn!("Skipping index build: SAM output cannot be indexed (use BAM/CRAM output, or pass --no-index to silence this)");
        } else {
            log::info!("Building index...");
            let idx_type = match index_format {
                IndexFormat::Bai => index::Type::Bai,
                IndexFormat::Csi => index::Type::Csi(14),
            };
            index::build(output, None, idx_type, threads as u32).context("building index")?;
        }
    }

    // The temporary Sled DB directory is owned and removed by the caller.
    Ok(stats)
}

// ── private helpers ───────────────────────────────────────────────────────────

/// One reserved byte identifying which mate of a template a record belongs to.
/// The two mates of a pair share a QNAME, so keying on the name alone collides
/// them. Valid SAM QNAMEs never contain bytes 0x00–0x03, so appending this suffix
/// keeps first/second-mate keys distinct and unambiguous.
fn mate_suffix(record: &bam::record::Record) -> u8 {
    match (record.is_first_in_template(), record.is_last_in_template()) {
        (false, false) => 0, // unpaired / single-end
        (true, false) => 1,  // first in pair
        (false, true) => 2,  // second in pair
        (true, true) => 3,   // both bits set (non-standard)
    }
}

/// Key-value store key for a record: its QNAME followed by the mate suffix.
fn record_key(record: &bam::record::Record) -> Vec<u8> {
    let qname = record.qname();
    let mut key = Vec::with_capacity(qname.len() + 1);
    key.extend_from_slice(qname);
    key.push(mate_suffix(record));
    key
}

fn convert_aux_to_serializable(value: Aux) -> SerializableAux {
    match value {
        Aux::String(s) => SerializableAux::String(s.to_owned()),
        Aux::Char(c) => SerializableAux::Char(c),
        Aux::I8(v) => SerializableAux::I8(v),
        Aux::U8(v) => SerializableAux::U8(v),
        Aux::I16(v) => SerializableAux::I16(v),
        Aux::U16(v) => SerializableAux::U16(v),
        Aux::I32(v) => SerializableAux::I32(v),
        Aux::U32(v) => SerializableAux::U32(v),
        Aux::Float(v) => SerializableAux::Float(v),
        Aux::Double(v) => SerializableAux::Double(v),
        Aux::ArrayU8(arr) => SerializableAux::ArrayU8(arr.iter().collect()),
        Aux::ArrayU16(arr) => SerializableAux::ArrayU16(arr.iter().collect()),
        Aux::ArrayU32(arr) => SerializableAux::ArrayU32(arr.iter().collect()),
        Aux::ArrayI8(arr) => SerializableAux::ArrayI8(arr.iter().collect()),
        Aux::ArrayI16(arr) => SerializableAux::ArrayI16(arr.iter().collect()),
        Aux::ArrayI32(arr) => SerializableAux::ArrayI32(arr.iter().collect()),
        Aux::ArrayFloat(arr) => SerializableAux::ArrayFloat(arr.iter().collect()),
        Aux::HexByteArray(s) => {
            let bytes: Vec<u8> = (0..s.len())
                .step_by(2)
                .filter_map(|i| u8::from_str_radix(&s[i..i + 2], 16).ok())
                .collect();
            SerializableAux::HexByteArray(bytes)
        }
    }
}
