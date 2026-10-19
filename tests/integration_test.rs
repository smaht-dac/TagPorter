//! Integration tests using pre-built BAM fixtures in tests/data/.
//!
//! Fixture contents
//! ─────────────────
//! source.bam (unaligned, 4 reads):
//!   read1  MM="C+m,0;"      ML=[200]
//!   read2  MM="C+m,1,2;"    ML=[100,150]
//!   read3  NM=5              (no MM/ML — tests non-requested tag)
//!   read4  MM="C+m,0,0;"    ML=[50,75]  (absent from aligned BAMs)
//!
//! aligned.bam (coordinate-sorted, 4 reads):
//!   read1  chr1:100  (no tags)
//!   read2  chr1:200  (no tags)
//!   read3  chr1:300  (no tags)
//!   read5  chr1:400  (not in source — tests unmatched aligned read)
//!
//! aligned_with_existing_mm.bam (coordinate-sorted, 2 reads):
//!   read1  chr1:100  MM="C+m,0,0,0;"   (conflict test target)
//!   read2  chr1:200  (no tags)

use tagporter::{
    collect_tags, ensure_output_distinct, write_new_bam, ConflictMode, IndexFormat,
    LengthMismatchMode, TransferStats,
};
use rust_htslib::bam::{self, record::Aux, Read};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

// ── helpers ───────────────────────────────────────────────────────────────────

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

/// Write a minimal single-contig reference FASTA (chr1, matching aligned.bam's
/// LN:100000) into `dir` and return its path. htslib auto-builds the .fai.
fn make_reference(dir: &Path) -> PathBuf {
    let path = dir.join("chr1.fa");
    let mut content = String::with_capacity(100_010);
    content.push_str(">chr1\n");
    content.push_str(&"A".repeat(100_000));
    content.push('\n');
    std::fs::write(&path, content).unwrap();
    path
}

/// Transcode a BAM fixture into a CRAM at `out_cram`, using `reference`.
fn make_cram(src_bam: &Path, out_cram: &Path, reference: &Path) {
    let mut reader = bam::Reader::from_path(src_bam).unwrap();
    reader.set_reference(reference).unwrap();
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(out_cram, &header, bam::Format::Cram).unwrap();
    writer.set_reference(reference).unwrap();
    for r in reader.records() {
        writer.write(&r.unwrap()).unwrap();
    }
}

fn run(
    tmp: &TempDir,
    source: &str,
    aligned: &str,
    tags: &[&str],
    on_conflict: ConflictMode,
) -> (bam::Reader, TransferStats) {
    let source_path = fixture(source);
    let aligned_path = fixture(aligned);
    let output_path = tmp.path().join("output.bam");
    let db_path = tmp.path().join("_sled_db");
    let tags_owned: Vec<String> = tags.iter().map(|s| s.to_string()).collect();

    collect_tags(&[source_path], &db_path, 1, &tags_owned, 128, None).unwrap();
    let stats = write_new_bam(
        &db_path,
        1,
        aligned_path,
        &output_path,
        true,
        IndexFormat::Bai,
        on_conflict,
        LengthMismatchMode::Skip,
        128,
        None,
    )
    .unwrap();

    (bam::Reader::from_path(&output_path).unwrap(), stats)
}

fn read_all(reader: &mut bam::Reader) -> Vec<bam::record::Record> {
    reader.records().map(|r| r.unwrap()).collect()
}

// ── tests ─────────────────────────────────────────────────────────────────────

/// MM (String) and ML (ArrayU8) tag values are transferred verbatim.
/// Output record order matches aligned.bam: read1, read2, read3, read5.
#[test]
fn test_mm_ml_values_correct() {
    let tmp = TempDir::new().unwrap();
    let (mut reader, _) = run(
        &tmp,
        "source.bam",
        "aligned.bam",
        &["MM", "ML"],
        ConflictMode::Skip,
    );
    let records = read_all(&mut reader);

    // read1: MM="C+m,0;", ML=[200]
    match records[0].aux(b"MM").unwrap() {
        Aux::String(s) => assert_eq!(s, "C+m,0;"),
        other => panic!("read1 MM: expected String, got {other:?}"),
    }
    match records[0].aux(b"ML").unwrap() {
        Aux::ArrayU8(arr) => assert_eq!(arr.iter().collect::<Vec<_>>(), vec![200u8]),
        other => panic!("read1 ML: expected ArrayU8, got {other:?}"),
    }

    // read2: MM="C+m,1,2;", ML=[100,150]
    match records[1].aux(b"MM").unwrap() {
        Aux::String(s) => assert_eq!(s, "C+m,1,2;"),
        other => panic!("read2 MM: expected String, got {other:?}"),
    }
    match records[1].aux(b"ML").unwrap() {
        Aux::ArrayU8(arr) => assert_eq!(arr.iter().collect::<Vec<_>>(), vec![100u8, 150u8]),
        other => panic!("read2 ML: expected ArrayU8, got {other:?}"),
    }
}

/// Transfer stats correctly count matched, unmatched, and total reads.
/// read1+read2 → tagged; read3 (NM only in source) + read5 (absent) → unmatched.
#[test]
fn test_transfer_stats() {
    let tmp = TempDir::new().unwrap();
    let (_, stats) = run(
        &tmp,
        "source.bam",
        "aligned.bam",
        &["MM", "ML"],
        ConflictMode::Skip,
    );

    assert_eq!(stats.total_aligned, 4);
    assert_eq!(stats.reads_tagged, 2);
    assert_eq!(stats.reads_unmatched, 2); // read3 (not in DB) + read5 (not in source)
    assert_eq!(stats.tags_skipped_conflict, 0);
}

/// read4 is in source but not in aligned — it must not appear in the output.
/// Output record count equals the number of reads in the aligned BAM (4).
#[test]
fn test_source_only_read_absent_from_output() {
    let tmp = TempDir::new().unwrap();
    let (mut reader, _) = run(
        &tmp,
        "source.bam",
        "aligned.bam",
        &["MM", "ML"],
        ConflictMode::Skip,
    );
    assert_eq!(read_all(&mut reader).len(), 4);
}

/// read3 has NM=5 in source.bam but MM/ML are not present.
/// When collecting MM,ML only, read3 is not stored in the DB,
/// so NM must never appear on read3 in the output.
#[test]
fn test_unrequested_tag_not_transferred() {
    let tmp = TempDir::new().unwrap();
    let (mut reader, _) = run(
        &tmp,
        "source.bam",
        "aligned.bam",
        &["MM", "ML"],
        ConflictMode::Skip,
    );
    let records = read_all(&mut reader);
    // records[2] is read3 (chr1:300)
    assert!(
        records[2].aux(b"NM").is_err(),
        "NM from source should not appear when only MM,ML are requested"
    );
}

/// read5 (chr1:400) is in aligned.bam but absent from source.bam.
/// It must be written to the output without any added tags.
#[test]
fn test_unmatched_aligned_read_written_without_tags() {
    let tmp = TempDir::new().unwrap();
    let (mut reader, _) = run(
        &tmp,
        "source.bam",
        "aligned.bam",
        &["MM", "ML"],
        ConflictMode::Skip,
    );
    let records = read_all(&mut reader);
    // records[3] is read5 (chr1:400)
    assert!(records[3].aux(b"MM").is_err(), "read5 should have no MM tag");
    assert!(records[3].aux(b"ML").is_err(), "read5 should have no ML tag");
}

/// With --on-conflict=skip (default), the existing MM tag on the aligned record
/// is preserved and the source value is discarded.
#[test]
fn test_conflict_skip_preserves_existing() {
    let tmp = TempDir::new().unwrap();
    let (mut reader, stats) = run(
        &tmp,
        "source.bam",
        "aligned_with_existing_mm.bam",
        &["MM"],
        ConflictMode::Skip,
    );
    assert_eq!(stats.tags_skipped_conflict, 1);
    // read1 is matched in the DB but its only tag (MM) is skipped → it must land
    // in reads_matched_no_new_tags, not vanish from the accounting. read2 has no
    // existing MM and gets tagged.
    assert_eq!(stats.reads_tagged, 1, "read2 receives its MM");
    assert_eq!(
        stats.reads_matched_no_new_tags, 1,
        "read1 matched but its only tag was conflict-skipped"
    );
    assert_eq!(stats.reads_unmatched, 0);
    // The four outcome buckets must partition total_aligned exactly (#5).
    assert_eq!(
        stats.reads_tagged
            + stats.reads_unmatched
            + stats.reads_matched_no_new_tags
            + stats.reads_secondary_supplementary,
        stats.total_aligned,
        "outcome buckets must sum to total_aligned"
    );

    let records = read_all(&mut reader);
    match records[0].aux(b"MM").unwrap() {
        Aux::String(s) => assert_eq!(s, "C+m,0,0,0;", "existing aligned value should be preserved"),
        other => panic!("Unexpected: {other:?}"),
    }
}

/// With --on-conflict=overwrite, the source value replaces the existing tag.
#[test]
fn test_conflict_overwrite_replaces_existing() {
    let tmp = TempDir::new().unwrap();
    let (mut reader, _) = run(
        &tmp,
        "source.bam",
        "aligned_with_existing_mm.bam",
        &["MM"],
        ConflictMode::Overwrite,
    );
    let records = read_all(&mut reader);
    match records[0].aux(b"MM").unwrap() {
        Aux::String(s) => assert_eq!(s, "C+m,0;", "source value should replace existing"),
        other => panic!("Unexpected: {other:?}"),
    }
}

/// With --on-conflict=error, encountering a pre-existing tag must return Err.
/// The error message should name the conflicting tag.
#[test]
fn test_conflict_error_returns_err() {
    let tmp = TempDir::new().unwrap();
    let source_path = fixture("source.bam");
    let aligned_path = fixture("aligned_with_existing_mm.bam");
    let output_path = tmp.path().join("output.bam");
    let db_path = tmp.path().join("_sled_db");

    collect_tags(&[source_path], &db_path, 1, &["MM".to_string()], 128, None).unwrap();
    let result = write_new_bam(
        &db_path,
        1,
        aligned_path,
        &output_path,
        true,
        IndexFormat::Bai,
        ConflictMode::Error,
        LengthMismatchMode::Skip,
        128,
        None,
    );

    assert!(result.is_err(), "expected Err on conflict");
    assert!(
        format!("{:?}", result.unwrap_err()).contains("MM"),
        "error should mention the conflicting tag name"
    );
}

/// Regression test for CRAM I/O: tags must transfer through a CRAM input and a
/// CRAM output when the reference is supplied only via the `reference` argument
/// (i.e. not discoverable through htslib's md5 cache or an embedded UR path).
/// Before the fix, `--reference` merely set the REF_PATH env var and CRAM I/O
/// failed because the reader/writer never called `set_reference`.
#[test]
fn test_cram_roundtrip_with_reference() {
    let tmp = TempDir::new().unwrap();
    let reference = make_reference(tmp.path());

    // Build a CRAM copy of the aligned BAM fixture.
    let aligned_cram = tmp.path().join("aligned.cram");
    make_cram(&fixture("aligned.bam"), &aligned_cram, &reference);

    let db_path = tmp.path().join("_sled_db");
    let output_cram = tmp.path().join("output.cram");
    let tags = vec!["MM".to_string(), "ML".to_string()];

    collect_tags(&[fixture("source.bam")], &db_path, 1, &tags, 128, Some(&reference)).unwrap();
    let stats = write_new_bam(
        &db_path,
        1,
        aligned_cram,
        &output_cram,
        true, // no_index
        IndexFormat::Bai,
        ConflictMode::Skip,
        LengthMismatchMode::Skip,
        128,
        Some(&reference),
    )
    .unwrap();
    assert_eq!(stats.total_aligned, 4);
    assert_eq!(stats.reads_tagged, 2);

    // Read the CRAM back (reference required) and confirm MM/ML landed on read1.
    let mut reader = bam::Reader::from_path(&output_cram).unwrap();
    reader.set_reference(&reference).unwrap();
    let records = read_all(&mut reader);
    match records[0].aux(b"MM").unwrap() {
        Aux::String(s) => assert_eq!(s, "C+m,0;"),
        other => panic!("read1 MM: expected String, got {other:?}"),
    }
    match records[0].aux(b"ML").unwrap() {
        Aux::ArrayU8(arr) => assert_eq!(arr.iter().collect::<Vec<_>>(), vec![200u8]),
        other => panic!("read1 ML: expected ArrayU8, got {other:?}"),
    }
}

/// Regression test for SAM output: with indexing enabled (the default), a `.sam`
/// output must still succeed — SAM cannot be indexed, so the index step is
/// skipped rather than aborting the whole run. Before the fix this returned Err
/// from `index::build` ("SAM file not BGZF compressed").
#[test]
fn test_sam_output_with_indexing_succeeds() {
    let tmp = TempDir::new().unwrap();
    let db_path = tmp.path().join("_sled_db");
    let output_sam = tmp.path().join("output.sam");
    let tags = vec!["MM".to_string(), "ML".to_string()];

    collect_tags(&[fixture("source.bam")], &db_path, 1, &tags, 128, None).unwrap();
    let stats = write_new_bam(
        &db_path,
        1,
        fixture("aligned.bam"),
        &output_sam,
        false, // indexing ON — must not error for SAM
        IndexFormat::Bai,
        ConflictMode::Skip,
        LengthMismatchMode::Skip,
        128,
        None,
    )
    .unwrap();

    assert_eq!(stats.total_aligned, 4);
    assert!(output_sam.exists(), "SAM output should be written");
    // No .bai index is produced for SAM.
    assert!(!tmp.path().join("output.sam.bai").exists());
}

/// Header with a single chr1 contig, used to build paired/secondary fixtures.
fn single_contig_header(len: i64) -> bam::Header {
    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", &"chr1")
            .push_tag(b"LN", &len),
    );
    header
}

fn mm_string(rec: &bam::record::Record) -> String {
    match rec.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        other => panic!("expected String MM, got {other:?}"),
    }
}

/// Regression test for read-name collision across mates. R1 and R2 of a pair
/// share a QNAME but carry different tags; each aligned mate must receive its own
/// tag, not the other's. Before the fix the DB keyed on QNAME alone, so the
/// second mate collected overwrote the first and both mates got the same value.
#[test]
fn test_paired_end_mates_get_distinct_tags() {
    use rust_htslib::bam::record::{Cigar, CigarString};
    let tmp = TempDir::new().unwrap();
    let header = single_contig_header(1000);
    let qual = [30u8; 4];

    // Source: unaligned pair "pairX", each mate with a distinct MM value.
    let source_path = tmp.path().join("src.bam");
    {
        let mut w = bam::Writer::from_path(&source_path, &header, bam::Format::Bam).unwrap();
        let mut r1 = bam::record::Record::new();
        r1.set(b"pairX", None, b"ACGT", &qual);
        r1.set_flags(0x1 | 0x4 | 0x40); // paired, unmapped, first-in-pair
        r1.push_aux(b"MM", Aux::String("C+m,0;")).unwrap();
        let mut r2 = bam::record::Record::new();
        r2.set(b"pairX", None, b"TTTT", &qual);
        r2.set_flags(0x1 | 0x4 | 0x80); // paired, unmapped, second-in-pair
        r2.push_aux(b"MM", Aux::String("C+m,5;")).unwrap();
        w.write(&r1).unwrap();
        w.write(&r2).unwrap();
    }

    // Aligned: same pair, both mates mapped.
    let aligned_path = tmp.path().join("aln.bam");
    {
        let cigar = CigarString(vec![Cigar::Match(4)]);
        let mut w = bam::Writer::from_path(&aligned_path, &header, bam::Format::Bam).unwrap();
        let mut r1 = bam::record::Record::new();
        r1.set(b"pairX", Some(&cigar), b"ACGT", &qual);
        r1.set_flags(0x1 | 0x2 | 0x40); // first-in-pair
        r1.set_tid(0);
        r1.set_pos(100);
        let mut r2 = bam::record::Record::new();
        r2.set(b"pairX", Some(&cigar), b"TTTT", &qual);
        r2.set_flags(0x1 | 0x2 | 0x80); // second-in-pair
        r2.set_tid(0);
        r2.set_pos(200);
        w.write(&r1).unwrap();
        w.write(&r2).unwrap();
    }

    let db_path = tmp.path().join("_sled_db");
    let output_path = tmp.path().join("out.bam");
    collect_tags(&[source_path], &db_path, 1, &["MM".to_string()], 128, None).unwrap();
    write_new_bam(
        &db_path,
        1,
        aligned_path,
        &output_path,
        true,
        IndexFormat::Bai,
        ConflictMode::Skip,
        LengthMismatchMode::Skip,
        128,
        None,
    )
    .unwrap();

    let mut reader = bam::Reader::from_path(&output_path).unwrap();
    let records = read_all(&mut reader);
    assert_eq!(mm_string(&records[0]), "C+m,0;", "first mate must keep its own MM");
    assert_eq!(mm_string(&records[1]), "C+m,5;", "second mate must keep its own MM");
}

/// Secondary/supplementary alignments share the primary's QNAME. Per the chosen
/// policy they are written through unchanged (no tags) and reported separately;
/// only the primary alignment receives the transferred tags.
#[test]
fn test_secondary_alignment_written_without_tags() {
    use rust_htslib::bam::record::{Cigar, CigarString};
    let tmp = TempDir::new().unwrap();
    let header = single_contig_header(100_000);
    let qual = [30u8; 8];

    // read1 exists in source.bam (flag 4, unpaired) with MM/ML.
    let aligned_path = tmp.path().join("aln.bam");
    {
        let cigar = CigarString(vec![Cigar::Match(8)]);
        let mut w = bam::Writer::from_path(&aligned_path, &header, bam::Format::Bam).unwrap();
        let mut primary = bam::record::Record::new();
        primary.set(b"read1", Some(&cigar), b"ACGTACGT", &qual);
        primary.set_flags(0);
        primary.set_tid(0);
        primary.set_pos(100);
        let mut secondary = bam::record::Record::new();
        secondary.set(b"read1", Some(&cigar), b"ACGTACGT", &qual);
        secondary.set_flags(0x100); // secondary
        secondary.set_tid(0);
        secondary.set_pos(500);
        w.write(&primary).unwrap();
        w.write(&secondary).unwrap();
    }

    let db_path = tmp.path().join("_sled_db");
    let output_path = tmp.path().join("out.bam");
    collect_tags(
        &[fixture("source.bam")],
        &db_path,
        1,
        &["MM".to_string(), "ML".to_string()],
        128,
        None,
    )
    .unwrap();
    let stats = write_new_bam(
        &db_path,
        1,
        aligned_path,
        &output_path,
        true,
        IndexFormat::Bai,
        ConflictMode::Skip,
        LengthMismatchMode::Skip,
        128,
        None,
    )
    .unwrap();

    assert_eq!(stats.total_aligned, 2);
    assert_eq!(stats.reads_tagged, 1, "only the primary is tag-eligible");
    assert_eq!(stats.reads_secondary_supplementary, 1);

    let mut reader = bam::Reader::from_path(&output_path).unwrap();
    let records = read_all(&mut reader);
    assert_eq!(records.len(), 2, "secondary alignment must still be written");
    assert!(records[0].aux(b"MM").is_ok(), "primary read1 should be tagged");
    assert!(
        records[1].aux(b"MM").is_err(),
        "secondary read1 must not be tagged"
    );
}

/// Regression test: writing the output over the aligned input must be refused
/// before the writer truncates it. Previously this "succeeded" (exit 0) while
/// leaving only the few records htslib had already buffered. Aliases of the same
/// file (`..` components, symlinks) must be caught too, and the input left intact.
#[test]
fn test_output_same_as_aligned_input_is_rejected() {
    let tmp = TempDir::new().unwrap();
    let aligned = tmp.path().join("aln.bam");
    std::fs::copy(fixture("aligned.bam"), &aligned).unwrap();
    let original = std::fs::read(&aligned).unwrap();

    std::fs::create_dir(tmp.path().join("sub")).unwrap();
    let mut aliases = vec![aligned.clone(), tmp.path().join("sub/../aln.bam")];
    #[cfg(unix)]
    {
        let link = tmp.path().join("link.bam");
        std::os::unix::fs::symlink(&aligned, &link).unwrap();
        aliases.push(link);
    }

    let db_path = tmp.path().join("_sled_db");
    collect_tags(&[fixture("source.bam")], &db_path, 1, &["MM".to_string()], 128, None).unwrap();
    for output in &aliases {
        let err = write_new_bam(
            &db_path,
            1,
            aligned.clone(),
            output,
            true,
            IndexFormat::Bai,
            ConflictMode::Skip,
            LengthMismatchMode::Skip,
            128,
            None,
        )
        .expect_err(&format!("output '{}' aliases the input", output.display()));
        assert!(err.to_string().contains("same file"), "unexpected error: {err}");
        assert_eq!(std::fs::read(&aligned).unwrap(), original, "input must be untouched");
    }
}

/// The source BAMs are protected as well, and a distinct (not yet existing)
/// output path passes the check.
#[test]
fn test_ensure_output_distinct_checks_every_input() {
    let tmp = TempDir::new().unwrap();
    let aligned = fixture("aligned.bam");
    let source = fixture("source.bam");
    let inputs = [aligned.as_path(), source.as_path()];

    assert!(ensure_output_distinct(&source, &inputs).is_err());
    assert!(ensure_output_distinct(&tmp.path().join("new.bam"), &inputs).is_ok());
}

// ── sequence-length validation of base-modification tags ─────────────────────

/// Write an unaligned source BAM holding one read, `read1`, with `seq` and `tags`.
fn write_source_read(path: &Path, seq: &[u8], tags: Vec<(&[u8; 2], Aux)>) {
    let header = single_contig_header(100_000);
    let mut w = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
    let mut rec = bam::record::Record::new();
    rec.set(b"read1", None, seq, &vec![30u8; seq.len()]);
    rec.set_flags(0x4);
    for (tag, value) in tags {
        rec.push_aux(tag, value).unwrap();
    }
    w.write(&rec).unwrap();
}

/// Write an aligned BAM holding one primary alignment of `read1`.
fn write_aligned_read(path: &Path, flags: u16, cigar: Vec<bam::record::Cigar>, seq: &[u8]) {
    let header = single_contig_header(100_000);
    let mut w = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
    let mut rec = bam::record::Record::new();
    rec.set(
        b"read1",
        Some(&bam::record::CigarString(cigar)),
        seq,
        &vec![30u8; seq.len()],
    );
    rec.set_flags(flags);
    rec.set_tid(0);
    rec.set_pos(100);
    w.write(&rec).unwrap();
}

/// Transfer `tags` from tmp/src.bam to tmp/aln.bam and return the output record.
fn transfer_one(
    tmp: &TempDir,
    tags: &[&str],
    on_length_mismatch: LengthMismatchMode,
) -> anyhow::Result<(bam::record::Record, TransferStats)> {
    let db_path = tmp.path().join("_sled_db");
    let output_path = tmp.path().join("out.bam");
    let tags_owned: Vec<String> = tags.iter().map(|s| s.to_string()).collect();
    collect_tags(&[tmp.path().join("src.bam")], &db_path, 1, &tags_owned, 128, None)?;
    let stats = write_new_bam(
        &db_path,
        1,
        tmp.path().join("aln.bam"),
        &output_path,
        true,
        IndexFormat::Bai,
        ConflictMode::Skip,
        on_length_mismatch,
        128,
        None,
    )?;
    let mut reader = bam::Reader::from_path(&output_path)?;
    let mut records = read_all(&mut reader);
    assert_eq!(records.len(), 1, "the read itself must always be written");
    Ok((records.remove(0), stats))
}

/// Source read1: ACGTACGT with one call on its first C, plus MN and optionally RG.
fn methylated_source(tmp: &TempDir, mn: i32, rg: Option<&str>) {
    let ml = [200u8];
    let mut tags: Vec<(&[u8; 2], Aux)> = vec![
        (b"MM", Aux::String("C+m,0;")),
        (b"ML", Aux::ArrayU8((&ml[..]).into())),
        (b"MN", Aux::I32(mn)),
    ];
    if let Some(rg) = rg {
        tags.push((b"RG", Aux::String(rg)));
    }
    write_source_read(&tmp.path().join("src.bam"), b"ACGTACGT", tags);
}

fn assert_no_mod_tags(rec: &bam::record::Record) {
    for tag in [b"MM", b"ML", b"MN"] {
        assert!(
            rec.aux(tag).is_err(),
            "{} must not be transferred onto a length-mismatched read",
            String::from_utf8_lossy(tag)
        );
    }
}

fn assert_buckets_partition(stats: &TransferStats) {
    assert_eq!(
        stats.reads_tagged
            + stats.reads_unmatched
            + stats.reads_matched_no_new_tags
            + stats.reads_secondary_supplementary,
        stats.total_aligned,
        "outcome buckets must sum to total_aligned"
    );
}

/// Regression test for review finding #3: a hard-clipped primary (2H6M) must not
/// receive the source's MM/ML, whose positions count from the clipped-off bases.
/// Previously the tags were copied verbatim and the call shifted to another C.
#[test]
fn test_hard_clipped_primary_drops_mod_tags() {
    use rust_htslib::bam::record::Cigar;
    let tmp = TempDir::new().unwrap();
    methylated_source(&tmp, 8, None);
    write_aligned_read(
        &tmp.path().join("aln.bam"),
        0,
        vec![Cigar::HardClip(2), Cigar::Match(6)],
        b"GTACGT",
    );

    let (rec, stats) = transfer_one(&tmp, &["MM", "ML", "MN"], LengthMismatchMode::Skip).unwrap();
    assert_no_mod_tags(&rec);
    assert_eq!(stats.reads_length_mismatch, 1);
    assert_eq!(stats.reads_length_mismatch_hard_clipped, 1);
    assert_eq!(stats.reads_tagged, 0);
    assert_eq!(stats.reads_matched_no_new_tags, 1);
    assert_buckets_partition(&stats);
}

/// A read trimmed before alignment (shorter SEQ, no hard clip) is a mismatch
/// that hard clipping does not explain.
#[test]
fn test_trimmed_read_drops_mod_tags() {
    use rust_htslib::bam::record::Cigar;
    let tmp = TempDir::new().unwrap();
    methylated_source(&tmp, 8, None);
    write_aligned_read(&tmp.path().join("aln.bam"), 0, vec![Cigar::Match(6)], b"ACGTAC");

    let (rec, stats) = transfer_one(&tmp, &["MM", "ML", "MN"], LengthMismatchMode::Skip).unwrap();
    assert_no_mod_tags(&rec);
    assert_eq!(stats.reads_length_mismatch, 1);
    assert_eq!(stats.reads_length_mismatch_hard_clipped, 0);
    assert_buckets_partition(&stats);
}

/// Only base-modification tags depend on SEQ; other requested tags (here RG)
/// are still transferred to a length-mismatched read.
#[test]
fn test_length_mismatch_keeps_other_tags() {
    use rust_htslib::bam::record::Cigar;
    let tmp = TempDir::new().unwrap();
    methylated_source(&tmp, 8, Some("grp1"));
    write_aligned_read(
        &tmp.path().join("aln.bam"),
        0,
        vec![Cigar::HardClip(2), Cigar::Match(6)],
        b"GTACGT",
    );

    let (rec, stats) =
        transfer_one(&tmp, &["MM", "ML", "MN", "RG"], LengthMismatchMode::Skip).unwrap();
    assert_no_mod_tags(&rec);
    match rec.aux(b"RG").unwrap() {
        Aux::String(s) => assert_eq!(s, "grp1"),
        other => panic!("RG: expected String, got {other:?}"),
    }
    assert_eq!(stats.reads_tagged, 1, "RG still counts as a transferred tag");
    assert_eq!(stats.reads_length_mismatch, 1);
    assert_buckets_partition(&stats);
}

/// With --on-length-mismatch=error the run fails and names the read.
#[test]
fn test_length_mismatch_error_mode_returns_err() {
    use rust_htslib::bam::record::Cigar;
    let tmp = TempDir::new().unwrap();
    methylated_source(&tmp, 8, None);
    write_aligned_read(
        &tmp.path().join("aln.bam"),
        0,
        vec![Cigar::HardClip(2), Cigar::Match(6)],
        b"GTACGT",
    );

    let err = transfer_one(&tmp, &["MM", "ML", "MN"], LengthMismatchMode::Error)
        .expect_err("length mismatch must abort in error mode");
    let msg = err.to_string();
    assert!(msg.contains("read1"), "error should name the read: {msg}");
    assert!(msg.contains("--on-length-mismatch"), "error should hint at the flag: {msg}");
}

/// Equal lengths transfer normally, including MN, and a reverse-strand alignment
/// (SEQ stored reverse-complemented, same length) is not a mismatch.
#[test]
fn test_equal_length_reverse_strand_transfers_mod_tags_and_mn() {
    use rust_htslib::bam::record::Cigar;
    let tmp = TempDir::new().unwrap();
    methylated_source(&tmp, 8, None);
    // ACGTACGT is its own reverse complement.
    write_aligned_read(&tmp.path().join("aln.bam"), 0x10, vec![Cigar::Match(8)], b"ACGTACGT");

    let (rec, stats) = transfer_one(&tmp, &["MM", "ML", "MN"], LengthMismatchMode::Skip).unwrap();
    assert_eq!(mm_string(&rec), "C+m,0;");
    let mn = match rec.aux(b"MN").unwrap() {
        Aux::I8(v) => v as i64,
        Aux::U8(v) => v as i64,
        Aux::I16(v) => v as i64,
        Aux::U16(v) => v as i64,
        Aux::I32(v) => v as i64,
        Aux::U32(v) => v as i64,
        other => panic!("MN: expected integer, got {other:?}"),
    };
    assert_eq!(mn, 8);
    assert_eq!(stats.reads_length_mismatch, 0);
    assert_eq!(stats.reads_tagged, 1);
}

/// A source whose own MN disagrees with its SEQ already carries invalid MM/ML;
/// they are dropped even though the aligned and source lengths agree.
#[test]
fn test_inconsistent_source_mn_drops_mod_tags() {
    use rust_htslib::bam::record::Cigar;
    let tmp = TempDir::new().unwrap();
    methylated_source(&tmp, 9, None);
    write_aligned_read(&tmp.path().join("aln.bam"), 0, vec![Cigar::Match(8)], b"ACGTACGT");

    let (rec, stats) = transfer_one(&tmp, &["MM", "ML", "MN"], LengthMismatchMode::Skip).unwrap();
    assert_no_mod_tags(&rec);
    assert_eq!(stats.reads_length_mismatch, 1);
    assert_eq!(stats.reads_length_mismatch_hard_clipped, 0);
    assert_buckets_partition(&stats);
}
