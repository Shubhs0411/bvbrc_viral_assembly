#!/usr/bin/env python3

import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional, List, Dict

from Bio import Entrez, SeqIO


def run_command(cmd: str, desc: str, capture_stdout: bool = False) -> Optional[bytes]:
    print(f"[Running] {desc}\n$ {cmd}")
    try:
        if capture_stdout:
            out = subprocess.check_output(cmd, shell=True)
            print(f"[Done] {desc}\n")
            return out
        else:
            subprocess.check_call(cmd, shell=True)
            print(f"[Done] {desc}\n")
            return None
    except subprocess.CalledProcessError as e:
        print(f"[Failed] {desc}\nExit code: {e.returncode}")
        raise


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def list_contigs_from_fasta(fa_path: Path) -> List[str]:
    contigs: List[str] = []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                contigs.append(line[1:].strip().split()[0])
    return contigs


def bwa_index_present(ref_fa: Path) -> bool:
    suffixes = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    return all((ref_fa.with_suffix(ref_fa.suffix + s)).exists() for s in suffixes)


def is_probable_accession(token: str) -> bool:
    """Heuristic: if token doesn't look like a local path and contains no path separators, treat as accession."""
    p = Path(os.path.expanduser(token))
    if p.exists():
        return False
    return ("/" not in token) and ("\\" not in token)


def fetch_genbank_to_fasta(accession: str, email: str, out_fa: Path) -> Path:
    if not email:
        raise ValueError("Entrez email is required when using GenBank references for reference-guided assembly.")
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    if not records:
        raise ValueError(f"No FASTA sequence returned for accession {accession}")
    with open(out_fa, "w") as fh:
        SeqIO.write(records, fh, "fasta")
    return out_fa


def resolve_reference_fasta(token: str, cache_dir: Path, email: Optional[str]) -> Path:
    """
    Resolve a single reference token that is either a local path or a GenBank accession into a local FASTA path.
    """
    ensure_dir(cache_dir)
    tok = token.strip()
    p = Path(os.path.expanduser(tok)).resolve()
    if p.exists():
        return p

    if not is_probable_accession(tok):
        raise FileNotFoundError(f"Reference path not found: {tok}")

    out_fa = cache_dir / f"{tok}.fasta"
    if not out_fa.exists() or out_fa.stat().st_size == 0:
        fetch_genbank_to_fasta(tok, email, out_fa)  # type: ignore[arg-type]
    return out_fa


def run_reference_guided(
    job_data: Dict,
    output_dir: str,
    read1: Optional[str] = None,
    read2: Optional[str] = None,
    read_single: Optional[str] = None,
) -> Dict:
    """
    Run a single-sample reference-guided assembly using bwa/bcftools.

    Inputs:
      - job_data: job JSON dict, expected to include:
          strategy: \"reference_guided\"
          reference_type: \"genbank\" or \"fasta\"
          reference_assembly / genbank_accession / fasta_file
          email (optional; required for genbank)
      - read1/read2/read_single: local FASTQ paths prepared by the caller
    Outputs (in output_dir):
      - <sample>.consensus.fasta
      - <sample>.vcf.gz (+ index)
      - QUAST output directory with report.txt and report.html
    Returns:
      - summary dict with keys: consensus_fasta, vcf, quast_txt, quast_html
    """
    outdir = ensure_dir(Path(output_dir))

    sample_name = job_data.get("output_file") or job_data.get("srr_id") or "sample"
    sample_name = str(sample_name)

    reference_type = (job_data.get("reference_type") or "").lower()
    if reference_type not in {"genbank", "fasta"}:
        raise ValueError("reference_type must be 'genbank' or 'fasta' for reference-guided strategy.")

    email = job_data.get("email") or os.environ.get("ENTREZ_EMAIL")

    # Resolve reference FASTA
    ref_cache = ensure_dir(outdir / "ref_cache")
    ref_token = None
    if reference_type == "genbank":
        ref_token = job_data.get("genbank_accession") or job_data.get("reference_assembly")
        if not ref_token:
            raise ValueError("genbank_accession or reference_assembly must be provided for GenBank references.")
    elif reference_type == "fasta":
        ref_token = job_data.get("fasta_file") or job_data.get("reference_assembly")
        if not ref_token:
            raise ValueError("fasta_file or reference_assembly must be provided for FASTA references.")

    ref_fasta = resolve_reference_fasta(str(ref_token), ref_cache, email)

    contigs = list_contigs_from_fasta(ref_fasta)
    if not contigs:
        raise ValueError(f"No sequences detected in reference FASTA {ref_fasta}")

    is_paired = bool(read1 and read2)

    # Trimming with fastp
    fastp_threads = job_data.get("fastp_threads", 0)
    fastp_thr_opt = f"--thread {fastp_threads}" if fastp_threads and fastp_threads > 0 else ""

    trim1 = outdir / f"{sample_name}_1.trim.fastq"
    trim2 = outdir / f"{sample_name}_2.trim.fastq"
    trim_single = outdir / f"{sample_name}.trim.fastq"

    if is_paired:
        if not (read1 and read2):
            raise ValueError("Both read1 and read2 must be provided for paired-end reference-guided assembly.")
        run_command(
            f"fastp -i {read1} -I {read2} -o {trim1} -O {trim2} {fastp_thr_opt} "
            f"-h {outdir / (sample_name + '_fastp.html')} -j {outdir / (sample_name + '_fastp.json')}",
            f"Trim paired-end {sample_name}",
        )
    else:
        if not read_single and not read1:
            raise ValueError("A single-end read file must be provided for reference-guided assembly.")
        src = read_single or read1  # type: ignore[assignment]
        run_command(
            f"fastp -i {src} -o {trim_single} {fastp_thr_opt} "
            f"-h {outdir / (sample_name + '_fastp.html')} -j {outdir / (sample_name + '_fastp.json')}",
            f"Trim single-end {sample_name}",
        )

    # Align with bwa
    if not bwa_index_present(ref_fasta):
        run_command(f"bwa index {ref_fasta}", f"BWA index {ref_fasta.name}")

    sam = outdir / f"{sample_name}.sam"
    align_threads = job_data.get("align_threads", 0)
    bwa_thr_opt = f"-t {align_threads}" if align_threads and align_threads > 0 else ""

    if is_paired:
        run_command(
            f"bwa mem {bwa_thr_opt} {ref_fasta} {trim1} {trim2} > {sam}",
            f"Align {sample_name} (PE)",
        )
    else:
        run_command(
            f"bwa mem {bwa_thr_opt} {ref_fasta} {trim_single} > {sam}",
            f"Align {sample_name} (SE)",
        )

    bam_sorted = outdir / f"{sample_name}.sorted.bam"
    run_command(
        f"samtools view -Sb {sam} | samtools sort -o {bam_sorted}",
        f"Sort BAM {sample_name}",
    )
    run_command(f"samtools index {bam_sorted}", f"Index BAM {sample_name}")

    # Variant calling with bcftools
    vcf_gz = outdir / f"{sample_name}.vcf.gz"
    region = job_data.get("region")
    region_opt = f"-r {region} " if region else ""
    run_command(
        f"bcftools mpileup -Ou -f {ref_fasta} {region_opt}{bam_sorted} | "
        f"bcftools call -mv -Oz -o {vcf_gz}",
        f"Variant calling {sample_name}",
    )
    run_command(f"bcftools index {vcf_gz}", f"Index VCF {sample_name}")

    # Consensus
    consensus = outdir / f"{sample_name}.consensus.fasta"
    run_command(
        f"bcftools consensus -f {ref_fasta} {vcf_gz} > {consensus}",
        f"Consensus {sample_name}",
    )

    # QUAST on consensus
    quast_dir = outdir / "quast_reference_guided"
    ensure_dir(quast_dir)
    quast_cmd = (
        f"quast.py -o {quast_dir} -t 4 --min-contig 200 {consensus}"
    )
    try:
        run_command(quast_cmd, f"QUAST reference-guided {sample_name}")
        quast_txt = quast_dir / "report.txt"
        quast_html_rel = "quast_reference_guided/report.html"
    except Exception:
        quast_txt = None
        quast_html_rel = ""

    return {
        "sample": sample_name,
        "consensus_fasta": str(consensus),
        "vcf": str(vcf_gz),
        "quast_txt": str(quast_txt) if quast_txt and quast_txt.exists() else "",
        "quast_html": quast_html_rel,
    }

