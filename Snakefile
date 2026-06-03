# ============================================================================
# PART 1: CONFIGURATION - Where are your files?
# ============================================================================

import os

# Read configuration
configfile: "config.yaml"

# Get values from config
CONTIG_DIR = config["contig_dir"]
READS_DIR = config["reads_dir"]
WAAFLE_DB = config["waafle_db"]
TAXONOMY = config["taxonomy"]

# ============================================================================
# PART 2: SAMPLE DETECTION - What samples do you have?
# ============================================================================

# Look in CONTIG_DIR and find all files like: *.assembly.fasta
# Extract the sample name (the * part)
# Example: 18193-HTB19-002.assembly.fasta → 18193-HTB19-002

SAMPLES = glob_wildcards(
    os.path.join(CONTIG_DIR, "{sample}.assembly.fasta")
).sample

# This creates a list like: ["18193-HTB19-002", "18193-HTB19-006", ...]

# ============================================================================
# PART 3: THE GOAL - What do you ultimately want?
# ============================================================================

rule all:
    input:
        # For EVERY sample, I want a .lgt.tsv.qc_pass file
        # expand() creates: ["18193-HTB19-002.lgt.tsv.qc_pass", "18193-HTB19-006.lgt.tsv.qc_pass", ...]
        expand("{sample}.lgt.tsv.qc_pass", sample=SAMPLES)

# ============================================================================
# PART 4: THE RULES - How to create each file type
# ============================================================================

# RULE 1: Create .blastout from .fasta
# -----------------------------------------------------------------------------
rule waafle_search:
    input:
        # Input: The contig file
        os.path.join(CONTIG_DIR, "{sample}.assembly.fasta")
    output:
        # Output: Create a .blastout file (in current directory)
        "{sample}.blastout"
    threads: 
        6  # Use 6 CPU cores
    resources:
        mem_mb=16000,   # Request 16GB memory
        runtime=720     # Request 720 minutes (12 hours)
    shell:
        # The command to run
        # {input} gets replaced with the input file path
        # {WAAFLE_DB} gets replaced with the database path
        # {threads} gets replaced with 6
        "waafle_search {input} {WAAFLE_DB} --threads {threads}"

# RULE 2: Create .gff from .blastout
# -----------------------------------------------------------------------------
rule waafle_genecaller:
    input:
        "{sample}.blastout"  # Need the blastout we just created
    output:
        "{sample}.gff"       # Create a gff file
    resources:
        mem_mb=8000,
        runtime=180
    shell:
        "waafle_genecaller {input}"

# RULE 3: Create .lgt.tsv from .fasta + .blastout + .gff
# -----------------------------------------------------------------------------
rule waafle_orgscorer:
    input:
        fasta=os.path.join(CONTIG_DIR, "{sample}.assembly.fasta"),
        blastout="{sample}.blastout",
        gff="{sample}.gff"
    output:
        "{sample}.lgt.tsv"
    resources:
        mem_mb=10000,
        runtime=480
    shell:
        # Use named inputs: {input.fasta}, {input.blastout}, {input.gff}
        "waafle_orgscorer {input.fasta} {input.blastout} {input.gff} {TAXONOMY}"

# RULE 4: Create .junctions.tsv from .fasta + .gff + reads
# -----------------------------------------------------------------------------
rule waafle_junctions:
    input:
        fasta=os.path.join(CONTIG_DIR, "{sample}.assembly.fasta"),
        gff="{sample}.gff",
        r1=os.path.join(READS_DIR, "{sample}_R1.fastq.gz"),
        r2=os.path.join(READS_DIR, "{sample}_R2.fastq.gz")
    output:
        "{sample}.junctions.tsv"
    threads: 12 
    resources:
        mem_mb=40000,
        runtime=720
    shell:
        """
        waafle_junctions {input.fasta} {input.gff} \
            --reads1 {input.r1} \
            --reads2 {input.r2} \
            --threads {threads}
        """

# RULE 5: Create .lgt.tsv.qc_pass from .lgt.tsv + .junctions.tsv
# -----------------------------------------------------------------------------
rule waafle_qc:
    input:
        lgt="{sample}.lgt.tsv",
        junctions="{sample}.junctions.tsv"
    output:
        "{sample}.lgt.tsv.qc_pass"  # THE FINAL OUTPUT!
    resources:
        mem_mb=4000,
        runtime=30
    shell:
        "waafle_qc {input.lgt} {input.junctions}"
