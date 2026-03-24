# BV-BRC Module

## Overview

This repository holds a BV-BRC module.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

## Viral assembly entrypoint (`scripts/run_viral_assembly.py`)

`scripts/run_viral_assembly.py` is the service entrypoint that takes a **single job JSON** and runs one of:

- **Reference-guided assembly** (strategy `reference_guided`) via `scripts/reference_guided_assembly.py`
- **IRMA-based assembly** (other strategies; existing behavior)

### How it runs reference-guided assembly

For reference-guided jobs, the script:

- **Obtains input reads** from exactly one source (staging uses `p3-cp` for workspace paths in production; local paths are supported for dev):
  - **Paired-end reads**: `paired_end_lib.read1` + `paired_end_lib.read2`
  - **Single-end reads**: `single_end_lib.read`
  - **SRA**: `srr_id` — `prefetch` / `fasterq-dump` run **in this script** (`scripts/sra_staging.py`) by default; assembly receives local FASTQs only.
- **Stages reference FASTA** (`reference_type: fasta`) into `reference_inputs/` via the same fetch helper as reads, then **resolves GenBank** accessions or local FASTAs and runs the science pipeline in `scripts/reference_guided_assembly.py`:
  - `fastp` trimming
  - `bwa mem` alignment
  - `samtools` sort/index + depth
  - `bcftools mpileup/call` variant calling
  - `bcftools consensus` consensus FASTA (with low-coverage masking)
  - `quast.py` report
- **Generates** `AssemblyReport.html` using the template under `lib/viral_assembly_report_template.html`.

### Environment / paths

- **`KB_TOP`**: In the BV-BRC runtime, `KB_TOP` points at the module root. For local/dev runs, the script now falls back to the repo root automatically if `KB_TOP` is unset.
- **Workspace paths**: In BV-BRC production, read paths and `reference_fasta_file` are workspace object paths and are copied with `p3-cp` (`fetch_file_from_ws` in `run_viral_assembly.py`).
  - **Dev convenience**: if the path is an existing absolute or repo-relative file, it is copied without `p3-cp`.
  - For parity testing without the full service, use `run_reference_guided(...)` with the same job keys as production, or `test_files/compare_reference_guided.py`.

### Local setup and run (quickstart)

Required tools in your env:

- `python3`, `fastp`, `bwa`, `samtools`, `bcftools`, `quast.py`
- for SRA jobs: `prefetch`, `fasterq-dump`
- Python package: `biopython`

Recommended local invocation:

```bash
/home/<user>/miniconda3/bin/conda run -n bvbrc python3 scripts/run_viral_assembly.py -j <job.json> -o <out_dir>
```

Example:

```bash
/home/<user>/miniconda3/bin/conda run -n bvbrc python3 scripts/run_viral_assembly.py \
  -j test_files/job_reference_guided_genbank_nc001474.json \
  -o runs/local_one_job
```

## UI → service JSON contract (reference-guided)

At minimum, the UI must provide:

- `strategy`: `"reference_guided"`
- `reference_type`: `"fasta"` or `"genbank"`
- Reference value (depends on `reference_type` only):
  - **GenBank:** `reference_input` only — one accession, `"acc1;acc2"` for multiple, or a JSON list of accessions.
  - **FASTA:** `reference_fasta_file` only — one path, `"path1;path2"` for multiple files, or a JSON list of paths.
- Exactly one read source:
  - `paired_end_lib` **OR** `single_end_lib` **OR** `srr_id`
- A sample/output name:
  - `output_file` (preferred) or `srr_id`

### Supported optional keys (reference-guided)

- `email`: required for GenBank fetching if `ENTREZ_EMAIL` isn’t set.
- `align_threads`, `fastp_threads`: tool thread knobs.
- `region`: restricts `bcftools mpileup` to a contig/region.
- `depth_cutoff`: integer depth cutoff for masking low-coverage bases with `N` (default `10`, `0` disables).
- `download_sra_from_prefetch` / `download_fastqs_from_sra` (bool, optional): override SRA staging. **Defaults are `true`** when `srr_id` is set (UI does not send these). Set to `false` for offline/CI if you must skip `prefetch` / `fasterq-dump`.

### Segmented references / multi-contig consensus

The reference-guided implementation supports segmented references in two ways:

- **Single multi-record FASTA**: `reference_fasta_file` points at one FASTA with multiple records.
  - If `segment_names` is provided, the FASTA headers are rewritten to those names.
  - By default this emits a combined multi-consensus FASTA; per-segment emission can be controlled via `per_segment_consensus`.
- **Multiple references**: use multiple tokens in the same field — `reference_input` (GenBank) or `reference_fasta_file` (paths) as a list or semicolon-separated string; tokens are concatenated into one multi-FASTA reference.

Optional segmented keys:

- `segment_names`: list or semicolon-separated string. Must match the number of FASTA records after reference resolution.
- `per_segment_consensus`: boolean. If the input is a single multi-record FASTA, default behavior matches the original pipeline (combined-only). Set `true` to also emit per-segment FASTAs.

### JSON examples

**FASTA (paired-end, segmented multi-record FASTA):**

```json
{
  "strategy": "reference_guided",
  "reference_type": "fasta",
  "reference_fasta_file": "test_files/Influenza_Genome.fasta",
  "output_file": "SRR28752452",
  "paired_end_lib": { "read1": "ws:/path/read1.fastq", "read2": "ws:/path/read2.fastq" },
  "segment_names": ["A_HA","A_NA","A_PB1","A_MP","A_NP","A_NS","A_PA","A_PB2"],
  "per_segment_consensus": true
}
```

**GenBank (SRA reads):**

```json
{
  "strategy": "reference_guided",
  "reference_type": "genbank",
  "reference_input": "NC_001474.2",
  "output_file": "SRR27422853",
  "srr_id": "SRR27422853",
  "email": "you@org.org"
}
```

**Multiple references (segmented tokens):**

```json
{
  "strategy": "reference_guided",
  "reference_type": "genbank",
  "reference_input": "NC_045512.2;MN908947.3",
  "segment_names": "seg1;seg2;seg3",
  "output_file": "sampleX",
  "single_end_lib": { "read": "ws:/path/reads.fastq" }
}
```

### Expected outputs

Top-level output folder contains final human-consumable artifacts; intermediate/tool files are under `output_files/`.

Top-level (typical):

- `AssemblyReport.html`
- `report.html` (QUAST moved up from quast subdir)
- consensus FASTA(s):
  - `<sample>.consensus.<accession_or_segment>.fasta`
  - `<sample>.consensus.multi.fasta`
  - FASTA single-reference alias (when applicable): `<sample>.consensus.<ref_label>.fasta`

`output_files/` (typical):

- staged reads (`reads_1.fastq`, `reads_2.fastq` or `reads.fastq`)
- trimmed reads
- reference index files (`.amb/.ann/.bwt/.pac/.sa/.fai`)
- alignment/variant files (`.sam`, `.sorted.bam`, `.bai`, `.vcf.gz`, `.csi`)
- low-coverage BED (`*.lowcov_dp<cutoff>.bed`)
- `ref_cache/` with fetched GenBank FASTAs for accession-based jobs

SRA note:

- For `srr_id` jobs, `.sra` and FASTQs are staged under `<out_dir>/` (and `<out_dir>/<SRR>/` layouts as produced by `prefetch` / `fasterq-dump`).

## Local parity testing and comparison (old vs new)

The repo includes a local comparison runner that:

- runs the **original CSV-driven pipeline** (`test_files/original_reference_guided_assembly.py`)
- runs the **new JSON-driven pipeline** locally (direct `run_reference_guided(...)`, no workspace tooling)
- then compares consensus FASTA plus a “brutal” SHA-256 comparison across multiple artifact patterns.

### Output folders

All local runs should go under `runs/` (ignored by git). Example:

- `runs/old_results/` (old pipeline)
- `runs/new_results/` (new pipeline)

### Run the comparison

From the repo root:

```bash
python test_files/compare_reference_guided.py \
  --csv test_files/sample.csv \
  --email you@org.org \
  --job-json test_files/job_reference_guided_fasta_segmented.json \
  --old-out runs/old_results \
  --new-out runs/new_results
```

If you already ran both and only want to compare:

```bash
python test_files/compare_reference_guided.py \
  --csv test_files/sample.csv \
  --email you@org.org \
  --job-json test_files/job_reference_guided_fasta_segmented.json \
  --old-out runs/old_results \
  --new-out runs/new_results \
  --skip-run
```

For more notes and examples, see `test_files/reference_guided_parity_test.md`.
