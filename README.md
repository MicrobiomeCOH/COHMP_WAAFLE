# COHMP WAAFLE — Multi-Sample LGT Profiling Pipeline

A [Snakemake](https://snakemake.readthedocs.io/) pipeline that runs **WAAFLE** (Workflow to Annotate Assemblies and Find LGT Events) across a cohort of metagenome samples on the COH HPC/SLURM cluster. It wraps the standard five-step WAAFLE workflow, auto-discovers all samples in your input directory, and fans them out as independent cluster jobs to produce per-sample lateral gene transfer (LGT) calls.

---

## Introduction

**Lateral gene transfer (LGT)** — also called horizontal gene transfer — is the movement of genetic material between organisms outside of vertical (parent-to-offspring) inheritance. In microbial communities it is a major driver of genomic diversification, spreading traits such as mobile elements, restriction-modification systems, and antibiotic-resistance genes across species. Detecting LGT directly from assembled metagenomes is difficult: a transferred gene looks like a normal gene, and the signal must be inferred from the taxonomic *disagreement* between neighboring genes on the same contig.

**WAAFLE** is a computational algorithm developed by the Huttenhower lab to profile LGT from assembled metagenomes. It performs a nucleotide homology search of contigs against a curated reference database, calls genes from the hits, then taxonomically scores each gene to find contigs whose genes split cleanly into two different clades — the signature of a putative LGT event, complete with a donor→recipient direction. An optional read-mapping step validates the gene-gene junctions to reject assembly artifacts. WAAFLE prioritizes specificity while retaining high sensitivity for inter-genus LGT (Hsu et al., *Nature Microbiology* 2025).

**This pipeline** adapts WAAFLE for the analytical core's standard practice. The upstream tool processes one contig file at a time; here, Snakemake discovers every sample in a project automatically, builds the dependency graph linking the five WAAFLE steps, and submits each step as a right-sized SLURM job — giving reproducible, parallel, multi-sample runs from a single command.

> For the full algorithm, benchmarks, and biological findings, check the paper (Hsu et al., *Nature Microbiology* 2025) and the upstream repo: <https://github.com/biobakery/waafle>.

---

## Pipeline overview

For each sample, five rules run in sequence (the junction track runs in parallel with gene scoring):

```
{sample}.assembly.fasta
        │
        ▼
 ┌──────────────────┐
 │  waafle_search   │  BLAST contigs against the ChocoPhlAn2 database
 └──────────────────┘
        │ {sample}.blastout
        ├───────────────────────────────────────────────┐
        ▼                                                ▼
 ┌──────────────────┐                          (also uses contigs + reads)
 │ waafle_genecaller│  call genes from BLAST hits
 └──────────────────┘
        │ {sample}.gff
        ├──────────────────────────────┐
        ▼                              ▼
 ┌──────────────────┐        ┌──────────────────┐
 │  waafle_orgscorer│        │ waafle_junctions │  map reads, score
 │  taxonomic score │        │ gene-gene support│  gene-gene junctions
 │  → find LGT      │        └──────────────────┘
 └──────────────────┘                │ {sample}.junctions.tsv
        │ {sample}.lgt.tsv           │
        └──────────────┬─────────────┘
                       ▼
              ┌──────────────────┐
              │    waafle_qc     │  keep LGT calls with read support
              └──────────────────┘
                       │
                       ▼
        {sample}.lgt.tsv.qc_pass   ← final high-confidence LGT calls
```

| Step | What it does |
|------|--------------|
| **waafle_search** | A light wrapper around `blastn` that searches each contig against the WAAFLE-formatted ChocoPhlAn2 nucleotide database. |
| **waafle_genecaller** | Calls gene coordinates directly from the BLAST hits and writes them as a GFF. |
| **waafle_orgscorer** | Taxonomically scores each gene, then evaluates whether a contig's genes resolve to one clade (no LGT) or two clades (putative LGT), assigning donor/recipient direction. |
| **waafle_junctions** | Maps the sample's paired-end reads back to the contigs and quantifies read support across each gene-gene junction. |
| **waafle_qc** | Filters the LGT calls, keeping only those whose junctions have sufficient read support — removing likely misassemblies. |

---

## Requirements

This pipeline assumes the standard COH microbiome program snakemake environment:

- **Snakemake** and **WAAFLE** available on `PATH` (via the cluster module system or a shared conda environment). The five `waafle_*` commands must be callable.
- **BLAST+** (pulled in as a WAAFLE dependency) for the search step.
- **The ChocoPhlAn2 reference database and taxonomy file**, already provisioned on the shared filesystem:
  - DB: `/coh_labs/microbiome/apollo/resources/dbs/chocophlan2/chocophlan2`
  - Taxonomy: `/coh_labs/microbiome/apollo/resources/dbs/chocophlan2/chocophlan2_taxonomy.tsv`
- **SLURM** access on the `compute` partition.

> Installing WAAFLE or building the database from scratch is out of scope here — the shared paths above are the supported defaults. If you need to do that yourself, follow the upstream instructions: <https://github.com/biobakery/waafle#installation>.

---

## Input data layout

The pipeline discovers samples by globbing the contig directory. Each sample needs:

| File | Location (config key) | Naming convention | Used by |
|------|-----------------------|-------------------|---------|
| Assembled contigs (FASTA) | `contig_dir` | `{sample}.assembly.fasta` | search, orgscorer, junctions |
| Paired-end reads (gzipped FASTQ) | `reads_dir` | `{sample}_R1.fastq.gz`, `{sample}_R2.fastq.gz` | junctions only |

The **sample name is derived from the contig filename**: `18193-HTB19-002.assembly.fasta` → sample `18193-HTB19-002`. The reads must use the *same* sample name with the `_R1` / `_R2` suffixes, so the junction step can pair them automatically.

Example:

```
contigs/
├── 18193-HTB19-002.assembly.fasta
├── 18193-HTB19-006.assembly.fasta
└── ...
hostdepleted/
├── 18193-HTB19-002_R1.fastq.gz
├── 18193-HTB19-002_R2.fastq.gz
├── 18193-HTB19-006_R1.fastq.gz
├── 18193-HTB19-006_R2.fastq.gz
└── ...
```

> Reads are only consumed by `waafle_junctions` (and therefore `waafle_qc`). The search → genecaller → orgscorer track runs on contigs alone.

---

## Configuration

All paths live in `config.yaml` — four keys, no code changes needed to point the pipeline at a new project:

```yaml
contig_dir: /scratch/dichen/Ryo_IRB18197-CBM588/contigs            # where {sample}.assembly.fasta files live
reads_dir:  /scratch/dichen/Ryo_IRB18197-CBM588/hostdepleted       # where {sample}_R1/_R2.fastq.gz live
waafle_db:  /coh_labs/microbiome/apollo/resources/dbs/chocophlan2/chocophlan2
taxonomy:   /coh_labs/microbiome/apollo/resources/dbs/chocophlan2/chocophlan2_taxonomy.tsv
```

| Key | Meaning |
|-----|---------|
| `contig_dir` | Directory of assembled contigs; also the source of the sample list. |
| `reads_dir` | Directory of host-depleted paired-end reads for junction validation. |
| `waafle_db` | Path prefix of the WAAFLE-formatted ChocoPhlAn2 BLAST database. |
| `taxonomy` | ChocoPhlAn2 taxonomy table used by `waafle_orgscorer`. |

To run a new project, point `contig_dir` and `reads_dir` at that project's data; the database keys normally stay on the shared defaults.

---

## Running the pipeline

### 1. Preview the plan (recommended)

From the pipeline directory, do a dry run to confirm the sample list and the DAG before submitting anything:

```bash
snakemake -n
```

This prints every job Snakemake intends to run without executing them. If the sample count is wrong, check your `contig_dir` path and the `{sample}.assembly.fasta` naming.

### 2. Submit to SLURM

Edit `config.yaml` for your project, then submit the controller job:

```bash
sbatch waafle.sh
```

`waafle.sh` launches a lightweight Snakemake controller that submits each rule as its own SLURM job:

```bash
snakemake -j 50 \
    --cluster "mkdir -p logs/{rule} && sbatch \
        --partition=compute \
        --mem={resources.mem_mb}M \
        --cpus-per-task={threads} \
        --time={resources.runtime} \
        --output=logs/{rule}/{rule}.%j.out \
        --error=logs/{rule}/{rule}.%j.err" \
    --latency-wait 60
```

- **`-j 50`** — up to 50 jobs in flight at once across all samples/steps.
- **`--cluster ...`** — each rule's `threads` / `resources` block is translated into the matching `sbatch` request, so every step gets right-sized CPU, memory, and walltime.
- **Per-rule logs** — written to `logs/{rule}/{rule}.<jobid>.out` and `.err`.
- **`--latency-wait 60`** — tolerates shared-filesystem lag before a rule's outputs appear.
- **Email** — `waafle.sh` notifies on `BEGIN,END,FAIL`; update `--mail-user` to your own address.

> The workflow file is named `Snakefile`, so plain `snakemake` finds it automatically — no `-s` flag needed.

Snakemake resumes where it left off: re-submitting after a partial or failed run only re-runs the missing steps.

---

## Outputs

Each sample produces the following files in the working directory:

| File | Description |
|------|-------------|
| `{sample}.blastout` | Raw `blastn` hits of contigs against ChocoPhlAn2. |
| `{sample}.gff` | Called gene coordinates per contig. |
| `{sample}.lgt.tsv` | Contigs called as **putative LGT**, with donor/recipient clades and direction. |
| `{sample}.no_lgt.tsv` | Single-clade contigs (no LGT detected). *Written by WAAFLE alongside `.lgt.tsv`.* |
| `{sample}.unclassified.tsv` | Contigs that could not be confidently classified. *Also written alongside `.lgt.tsv`.* |
| `{sample}.junctions.tsv` | Read-support metrics for each gene-gene junction. |
| **`{sample}.lgt.tsv.qc_pass`** | **Final deliverable** — the subset of `.lgt.tsv` LGT calls with adequate junction read support. |

> **Note on side outputs:** `waafle_orgscorer` writes three files (`.lgt.tsv`, `.no_lgt.tsv`, `.unclassified.tsv`), but the pipeline only tracks `.lgt.tsv` as the rule's formal output (it is the sole input to QC). The `.no_lgt.tsv` and `.unclassified.tsv` files are still produced on disk — look for them if you want the full per-contig classification, not just the LGT calls.

### Reading `.lgt.tsv` / `.lgt.tsv.qc_pass`

The most useful columns for interpretation:

| Column | Meaning |
|--------|---------|
| `CONTIG_NAME` | Contig identifier from the input FASTA. |
| `CALL` | Classification of the contig (e.g. `lgt`). |
| `DIRECTION` | Inferred transfer direction between the two clades (e.g. `B>A`). |
| `CLADE_A` / `CLADE_B` | The two organisms involved (donor / recipient). |
| `TAXONOMY_A` / `TAXONOMY_B` | Full taxonomic lineages of each clade. |
| `LOCI` | Per-gene coordinates (`START:STOP:STRAND`). |
| `ANNOTATIONS:UNIPROT` | UniProt annotation(s) per gene. |

For the complete column dictionary, the synteny notation, and the score thresholds, see the upstream documentation: <https://github.com/biobakery/waafle#output-files>.

### Tuning QC stringency

`waafle_qc` exposes two thresholds (defaults shown) if you need to make junction filtering stricter or looser:

- `--min-junction-hits` (default `2`) — minimum spanning mate-pairs to support a junction.
- `--min-junction-ratio` (default `0.5`) — minimum junction-to-flanking-gene coverage ratio.

To change them, add the flag to the `waafle_qc` rule's `shell` command in the `Snakefile`.

---

## Resource profile

Each rule requests resources independently (defined in the `Snakefile`). Tune these for larger cohorts or larger assemblies:

| Rule | Threads | Memory | Walltime |
|------|--------:|-------:|---------:|
| `waafle_search` | 6 | 16 GB | 720 min |
| `waafle_genecaller` | 1 | 8 GB | 180 min |
| `waafle_orgscorer` | 1 | 10 GB | 480 min |
| `waafle_junctions` | 12 | 40 GB | 720 min |
| `waafle_qc` | 1 | 4 GB | 30 min |

`waafle_junctions` is the heaviest step (read mapping, 40 GB) and is usually the limiting factor for throughput.

---

## Troubleshooting

| Symptom | Likely cause / fix |
|---------|--------------------|
| `snakemake -n` lists **0 samples** / "Nothing to be done" | `contig_dir` is wrong, or contig files don't match `{sample}.assembly.fasta` exactly. |
| `MissingInputException` for `_R1`/`_R2` | A sample has contigs but no matching reads in `reads_dir`, or the read suffix isn't `_R1.fastq.gz` / `_R2.fastq.gz` with the *same* sample name. |
| BLAST / database errors in `waafle_search` | Check `waafle_db` points at the database **prefix** (not a single file); confirm the shared path is mounted on the compute nodes. |
| Jobs queue but never start | Confirm partition `compute` is correct for your allocation and that the requested walltime/memory are within limits. |
| A step failed midway | Inspect `logs/{rule}/{rule}.<jobid>.err`, fix the cause, then re-run `sbatch waafle.sh` — completed steps are not redone. |

---

## Citation & references

If you use results from this pipeline, cite the WAAFLE paper:

> Hsu TY, Nzabarushimana E, Wong D, Luo C, Beiko RG, Langille M, Huttenhower C, Nguyen LH, Franzosa EA. **Profiling lateral gene transfer events in the human microbiome using WAAFLE.** *Nature Microbiology* 10:94–111 (2025). doi:[10.1038/s41564-024-01881-w](https://doi.org/10.1038/s41564-024-01881-w)

- Upstream WAAFLE repository: <https://github.com/biobakery/waafle>
- WAAFLE on PyPI: <https://pypi.org/project/waafle/>
- Paper PDF: `waafle_paper.pdf` (in this repository)

---

*Pipeline maintained by the COH Microbiome analytical core. For pipeline issues (sample discovery, SLURM, config), contact the maintainer; for WAAFLE algorithm questions, see the upstream repository.*
