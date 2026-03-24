#!/usr/bin/env python3
"""
SRA input staging (prefetch + fasterq-dump) for reference-guided and legacy callers.

Kept separate from reference_guided_assembly so the assembly module can focus on
alignment/variant/consensus work; the run script imports this for production jobs.
"""

import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def run_command(cmd: str, desc: str, capture_stdout: bool = False):
    print(f"[Running] {desc}\n$ {cmd}")
    try:
        if capture_stdout:
            out = subprocess.check_output(cmd, shell=True)
            print(f"[Done] {desc}\n")
            return out
        subprocess.check_call(cmd, shell=True)
        print(f"[Done] {desc}\n")
        return None
    except subprocess.CalledProcessError as e:
        print(f"[Failed] {desc}\nExit code: {e.returncode}")
        raise


def ensure_sra_fastqs(
    output_dir: str,
    srr_id: str,
    job_data: Optional[Dict[str, Any]] = None,
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Ensure local FASTQs exist for SRR *srr_id* under *output_dir*.

    Default: prefetch .sra if missing, then fasterq-dump if FASTQs missing.
    Set job_data.download_sra_from_prefetch / download_fastqs_from_sra to False to disable.

    Returns (read1, read2, read_single) — paired uses read1+read2, single uses read_single.
    """
    job_data = job_data or {}
    outdir = ensure_dir(Path(output_dir))
    srr_out_dir = outdir
    sra_folder = ensure_dir(outdir / str(srr_id))
    sra_path = sra_folder / f"{srr_id}.sra"

    dl_prefetch = job_data.get("download_sra_from_prefetch")
    if dl_prefetch is None:
        dl_prefetch = True
    else:
        dl_prefetch = bool(dl_prefetch)

    dl_fastq = job_data.get("download_fastqs_from_sra")
    if dl_fastq is None:
        dl_fastq = True
    else:
        dl_fastq = bool(dl_fastq)

    if not sra_path.exists() or sra_path.stat().st_size == 0:
        repo_root = Path(__file__).resolve().parents[1]
        local_candidates = [
            Path.cwd() / str(srr_id) / f"{srr_id}.sra",
            repo_root / str(srr_id) / f"{srr_id}.sra",
        ]
        for cand in local_candidates:
            if cand.exists() and cand.stat().st_size > 0:
                ensure_dir(sra_folder)
                shutil.copy2(str(cand), str(sra_path))
                break
        else:
            if not dl_prefetch:
                raise RuntimeError(
                    f"SRR SRA file not found at {sra_path} and local copy candidates are missing. "
                    "Set download_sra_from_prefetch=true (default) to enable prefetch, "
                    "or provide read1/read2/read_single."
                )
            try:
                run_command(
                    f"prefetch -O {srr_out_dir} {srr_id}",
                    f"prefetch {srr_id}",
                )
            except Exception as e:
                raise RuntimeError(f"prefetch failed for {srr_id}: {e}") from e

    if not sra_path.exists() or sra_path.stat().st_size == 0:
        candidates = sorted(
            outdir.glob(f"**/{srr_id}.sra"),
            key=lambda p: len(p.parts),
        )
        if candidates:
            sra_path = candidates[0]
        else:
            raise RuntimeError(f"Expected SRA file not found after prefetch: {sra_path}")

    fq1 = outdir / f"{srr_id}_1.fastq"
    fq2 = outdir / f"{srr_id}_2.fastq"
    fqs = outdir / f"{srr_id}.fastq"
    alt_fq1 = outdir / str(srr_id) / f"{srr_id}_1.fastq"
    alt_fq2 = outdir / str(srr_id) / f"{srr_id}_2.fastq"
    alt_fqs = outdir / str(srr_id) / f"{srr_id}.fastq"

    have_fastqs = (fq1.exists() and fq2.exists()) or fqs.exists() or (
        alt_fq1.exists() and alt_fq2.exists()
    ) or alt_fqs.exists()

    if not have_fastqs:
        if not dl_fastq:
            raise RuntimeError(
                f"Reads were not provided for SRR job {srr_id}. "
                "Set download_fastqs_from_sra=true (default) to run fasterq-dump, "
                "or provide read paths."
            )
        fasterq_cmd = f"fasterq-dump -O {srr_out_dir} {srr_id}"
        disk_limit = job_data.get("fasterq_dump_disk_limit")
        disk_limit_tmp = job_data.get("fasterq_dump_disk_limit_tmp")
        fasterq_temp_dir = job_data.get("fasterq_dump_temp_dir")
        if disk_limit is not None:
            fasterq_cmd += f" --disk-limit {disk_limit}"
        if disk_limit_tmp is not None:
            fasterq_cmd += f" --disk-limit-tmp {disk_limit_tmp}"
        if fasterq_temp_dir:
            fasterq_cmd += f" --temp {fasterq_temp_dir}"
        try:
            run_command(
                fasterq_cmd,
                f"fasterq-dump {srr_id}",
            )
        except Exception as e:
            raise RuntimeError(
                "SRR-based job requires fasterq-dump when reads are not provided. "
                f"fasterq-dump failed for {srr_id}: {e}"
            ) from e

    if fq1.exists() and fq2.exists():
        return str(fq1), str(fq2), None
    if fqs.exists():
        return None, None, str(fqs)
    if alt_fq1.exists() and alt_fq2.exists():
        return str(alt_fq1), str(alt_fq2), None
    if alt_fqs.exists():
        return None, None, str(alt_fqs)
    raise RuntimeError(f"FASTQ files not found after fasterq-dump for {srr_id}.")
