# PacBio HIV/SIV Hypermutation Pipeline

A reproducible HPC workflow for PacBio HiFi amplicon sequencing data from HIV, HIVΔvif, and SIV experiments. The pipeline demultiplexes and primer-trims reads, assigns samples to viral references using barcode metadata, maps reads, quantifies APOBEC-associated hypermutation, calls mutations in multiple sequence contexts, classifies reads by hypermutation status, and performs CDS-aware stop-gain and frameshift analyses.

> This repository provides analysis code, not experimental data. Supply the raw PacBio BAM, barcode and primer FASTAs, metadata CSV, matched HIV/SIV reference FASTAs, and CDS BED annotations before running the workflow.

## Pipeline

```text
PacBio HiFi barcoded BAM
        |
        v
Lima barcode demultiplexing
        |
        v
Lima primer trimming
        |
        v
Metadata validation and HIV/SIV reference assignment
        |
        v
BAM-to-FASTQ conversion and minimap2 mapping
        |
        v
Primary-read BAM filtering
        |
        +-----------------------------+
        |                             |
        v                             v
Hypermut3: RD, GD, AD          Mutation calling: 1-nt, 2-nt, 3-nt
        |                             |
        +--------------+--------------+
                       |
                       v
Read-level hypermutation classification
                       |
                       v
Hypermutated and nonhypermutated BAM subsets
                       |
                       v
APOBEC-compatible stop-gain and frameshift analysis
                       |
                       v
R-based quality control, statistics, and figures
```

The provided shell scripts target SLURM-managed Linux clusters. Python scripts may be run locally after adapting project paths, environment activation, and scheduler-specific commands. Most wrappers are resumable and avoid rerunning stages with valid nonempty outputs.

## Contents

| File | Purpose |
|---|---|
| `pacbio-lima-HIV-SIV.sh` | Demultiplexes PacBio HiFi BAM input by barcode and trims primers using Lima. |
| `HIV-SIV-mapping-2.sh` | Validates barcode metadata, assigns HIV/SIV references, converts BAM to FASTQ, maps reads with minimap2, and produces mapping QC. |
| `BAM_HYPMUT3-3.py` | Hypermut3 implementation for BAM alignments, mutation-context counting, and Fisher exact testing. |
| `hypermut3-all-4.sh` | Runs the Hypermut3 program in RD, GD, and AD contexts and builds merged metadata-aware output. |
| `mutation_counting2-6.py` | Per-read single-base substitution counting. |
| `dinucleotide_mutation_analysis-5.py` | Per-read dinucleotide-context mutation analysis, including directional and opportunity-normalized measures. |
| `trinucleotide_mutation_analysis-8.py` | Per-read trinucleotide-context mutation analysis with raw and normalized measures. |
| `mutational-calling-all-contexts-7.sh` | SLURM-array submitter for 1-nt, 2-nt, and 3-nt mutation calling and table combination. |
| `apobec_stopgain_pipeline-9.py` | CDS-aware analysis of G-to-A events, APOBEC-compatible stop gains, premature stops, indels, and frameshifts. |
| `collapse_stopgain_master-10.py` | Reduces overlapping-CDS output to one row per `sample_name` and `read_id`. |
| `stop-codon-hypermut-nonhypermut-11.sh` | Partitions reads by Hypermut3 status, filters BAMs, runs stop-gain analysis, and produces class-specific master tables. |
| `PacBio-HIV-SIV-2026-08-31-12.R` | Exploratory statistics, threshold sensitivity analyses, mutation-ratio QC, and publication-ready plots. |

## Inputs

### Metadata

The metadata CSV must contain these exact headers:

```csv
barcode,Virus
```

Example barcode:

```text
bc1004_F2--bc1057_R2
```

| `Virus` value | Analysis group | Selected reference |
|---|---|---|
| `HIV` | HIV | HIV near-full-length FASTA |
| `HIVDvif` | HIV | HIV near-full-length FASTA |
| `SIV` | SIV | SIV near-full-length FASTA |

Additional metadata columns are retained and merged into downstream result tables. Barcode assignment is treated as authoritative; scripts flag missing barcodes, unsupported virus labels, ambiguous entries, malformed metadata headers, UTF-8 byte-order marks, and Windows line-ending artifacts.

### Required files

| Input | Used by | Notes |
|---|---|---|
| PacBio HiFi BAM | Lima preprocessing | BAM input is expected rather than raw FASTQ. |
| Barcode FASTA | Lima preprocessing | Must match the barcode design used for the library. |
| Primer FASTA | Lima preprocessing | Used in post-demultiplex primer trimming. |
| Metadata CSV | Mapping and downstream analyses | Requires exact `barcode` and `Virus` columns. |
| HIV FASTA | HIV/HIVΔvif mapping and stop-gain analysis | Must match the HIV CDS BED annotation. |
| SIV FASTA | SIV mapping and stop-gain analysis | Must match the SIV CDS BED annotation. |
| HIV CDS BED | Stop-gain analysis | Must contain validated CDS positions, strand, and frame. |
| SIV CDS BED | Stop-gain analysis | Must contain validated CDS positions, strand, and frame. |

The CDS BED annotations should be validated against their exact plasmid-derived reference sequences. In the source workflow, HIV features were transferred from `KJ925006.1` and SIV features from `M33262.1`, then checked using protein similarity searches and six-frame translation inspection in Geneious Prime.

### File naming

Downstream wrappers expect filenames that preserve barcode identity:

```text
Demultiplexed BAM:  DNA_NFL.demux.<barcode>.bam
Trimmed BAM:        DNA_NFL.demux.<barcode>.trimmed.bam
Mapped BAM:         DNA_NFL.demux.<barcode>.aligned.bam
```

Example:

```text
DNA_NFL.demux.bc1004_F2--bc1057_R2.trimmed.bam
DNA_NFL.demux.bc1004_F2--bc1057_R2.aligned.bam
```

Barcode components may include numeric suffixes such as `F2` and `R2`. These distinguish reused barcode pairs that correspond to different biological samples across sequencing runs or modules.

## Installation

### Requirements

- Linux or a Linux-compatible HPC environment
- Bash 4+
- SLURM for the supplied submission wrappers
- Conda, Miniconda, or Mambaforge
- Sufficient storage for BAMs, FASTQs, BAM indices, intermediate tables, and logs

Core executables:

```text
lima
bam2fastq
pbindex
minimap2
samtools
python3
Rscript
sbatch
```

The shell workflow also uses standard Unix programs including `awk`, `sed`, `tr`, `grep`, `find`, `sort`, `wc`, `head`, `basename`, `dirname`, and `gzip`. `pbindex` is optional in the mapping workflow.

### Conda environments

Separate environments are recommended to avoid conflicts between PacBio utilities and scientific Python dependencies.

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

Create the preprocessing/mapping environment:

```bash
conda create -y -n pacbio-isoseq \
  -c conda-forge -c bioconda \
  python=3.11 lima pbtk minimap2 samtools pysam pandas numpy
```

Create the Hypermut3 and mutation-analysis environment:

```bash
conda create -y -n hypermut3 \
  -c conda-forge -c bioconda \
  python=3.11 samtools pysam pandas scipy numpy
```

Create the R environment:

```bash
conda create -y -n pacbio-r \
  -c conda-forge -c bioconda \
  r-base=4.3 r-readr r-dplyr r-ggplot2 r-scales r-svglite \
  r-tidyr r-viridis r-patchwork
```

Required Python packages:

```text
python >=3.10,<3.13
numpy
pandas
pysam
scipy
```

Required R packages:

```text
readr
dplyr
ggplot2
scales
svglite
tidyr
viridis
patchwork
```

### SLURM Conda initialization

Batch shells do not always load Conda automatically. Source `conda.sh` before activation:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hypermut3
```

On module-based clusters, load the site-specific Conda module first:

```bash
module load miniconda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hypermut3
```

### Validate software

```bash
conda activate pacbio-isoseq
which lima bam2fastq minimap2 samtools
lima --version
bam2fastq --help | head
minimap2 --version
samtools --version
command -v pbindex && pbindex --help | head || true
```

```bash
conda activate hypermut3
python - <<'PY'
import numpy
import pandas
import pysam
import scipy
from scipy.stats import fisher_exact
print("Scientific Python environment: OK")
PY
```

## Configuration

Update the configuration blocks at the beginning of each shell wrapper before submission. At minimum, set:

```bash
PROJECT_DIR="/path/to/HIV-SIV-Hypermutation-Assay-PacBio"
ANALYSIS_DIR="${PROJECT_DIR}/analysis-<run-date>"
```

Also configure paths to raw data, barcode and primer FASTAs, metadata, both reference FASTAs, CDS BED annotations, Python scripts, SLURM logs, and your Conda installation.

Use a separate dated `ANALYSIS_DIR` for each independent run.

Some wrappers refer to stable, unversioned Python script names. Either update the wrapper variables to use versioned filenames directly or create symlinks:

```bash
mkdir -p "${PROJECT_DIR}/scripts"
ln -sfn "$(pwd)/BAM_HYPMUT3-3.py" "${PROJECT_DIR}/scripts/BAM_HYPMUT3.py"
ln -sfn "$(pwd)/mutation_counting2-6.py" "${PROJECT_DIR}/scripts/mutation_counting2.py"
ln -sfn "$(pwd)/dinucleotide_mutation_analysis-5.py" "${PROJECT_DIR}/scripts/dinucleotide_mutation_analysis.py"
ln -sfn "$(pwd)/trinucleotide_mutation_analysis-8.py" "${PROJECT_DIR}/scripts/trinucleotide_mutation_analysis.py"
ln -sfn "$(pwd)/apobec_stopgain_pipeline-9.py" "${PROJECT_DIR}/scripts/apobec_stopgain_pipeline.py"
ln -sfn "$(pwd)/collapse_stopgain_master-10.py" "${PROJECT_DIR}/scripts/collapse_stopgain_master.py"
```

Create FASTA indices required by `pysam.FastaFile` before stop-gain analysis:

```bash
conda activate hypermut3
samtools faidx reference-genomes/pJH048_i3_HIV-NFL.fasta
samtools faidx reference-genomes/SIV239-SpX_i22-NFL.fasta
```

## Usage

### 1. Demultiplex and trim

```bash
conda activate pacbio-isoseq
sbatch scripts/pacbio-lima-HIV-SIV.sh
```

This stage uses Lima to split reads by barcode and trim primers. Trimmed BAMs are written under:

```text
analysis-<run-date>/lima/trimmed/
```

### 2. Map reads

```bash
conda activate pacbio-isoseq
sbatch scripts/HIV-SIV-mapping-2.sh
```

For each primer-trimmed BAM, this workflow extracts the barcode from the filename, validates it against metadata, selects the appropriate viral reference, converts BAM to compressed FASTQ, maps HiFi reads with minimap2, and retains primary mapped reads using:

```bash
samtools view -F 2308
```

Mapped BAMs are routed to `HIV/mapped/` or `SIV/mapped/`.

> Compute mapping efficiency from pre-filter FASTQ read counts and retained mapped-read counts. Do not infer mapping efficiency from a filtered mapped BAM alone, because unmapped reads have already been removed.

### 3. Run Hypermut3

```bash
conda activate hypermut3
sbatch scripts/hypermut3-all-4.sh
```

Hypermut3 analyzes each mapped BAM in three contexts:

| Context | Interpretation |
|---|---|
| `RD` | Broad G-to-A context used for the primary Fisher-test classification. |
| `GD` | GG-associated G-to-A mutation context. |
| `AD` | GA-associated G-to-A mutation context. |

Key output:

```text
analysis-<run-date>/hypermut3/combined/master_hypermut3_RD_GD_AD_metadata.csv
```

Read classification for downstream stop-gain analysis is based on `fisher_p_RD`:

- **Hypermutated:** `fisher_p_RD < 0.05`
- **Nonhypermutated:** `fisher_p_RD >= 0.05`
- **Unclassified:** missing, invalid, or nonnumeric p-value; reported in QC and excluded from both classes

### 4. Call mutation contexts

```bash
conda activate hypermut3
sbatch scripts/mutational-calling-all-contexts-7.sh
```

The submitter creates a metadata-validated BAM manifest, launches SLURM arrays for each context size, and combines results after successful completion.

- **1-nt:** Raw counts of all 12 nonidentity single-base substitutions per read
- **2-nt:** Full dinucleotide context substitutions, directional mutation classes, collapsed classes, and opportunity-normalized values
- **3-nt:** Trinucleotide context counts and opportunity-normalized measures

Duplicate and unmapped reads are excluded. Incomplete or gapped motif windows are not counted. Normalization denominators correspond to valid, ungapped reference windows available per read.

### 5. Analyze stop gains

```bash
conda activate hypermut3
sbatch scripts/stop-codon-hypermut-nonhypermut-11.sh
```

This stage creates class-specific read lists, filters BAMs with `samtools view -N`, applies the relevant HIV or SIV FASTA and CDS BED annotation, and reports coding-sequence mutation consequences.

The analysis identifies:

- Observed CDS-oriented G-to-A events
- APOBEC-compatible G-to-A stop gains
- All premature stop codons
- CDS-overlapping indels
- Frameshift indels
- Frameshift-associated premature stops

A potential APOBEC-compatible stop gain is defined as a G at codon position 1 or 2 in a `GG` or `GA` motif where G-to-A produces `TAA`, `TAG`, or `TGA`.

### 6. Generate statistics and figures

```bash
conda activate pacbio-r
Rscript scripts/PacBio-HIV-SIV-2026-08-31-12.R
```

Update the R script's working directory and input locations as needed. The analysis includes RD p-value threshold sensitivity, donor-level hypermutation summaries, RD/GD/AD ratio QC, GD-versus-AD scatterplots, and pooled SBS distributions for HIV versus HIVΔvif reads.

## Outputs

| Location or file | Description |
|---|---|
| `HIV/*.fastq.gz`, `SIV/*.fastq.gz` | BAM-derived reads routed by metadata-defined virus group. |
| `HIV/mapped/*.aligned.bam`, `SIV/mapped/*.aligned.bam` | Primary mapped reads after alignment and filtering. |
| `mapping_summary.<jobid>.tsv` | Per-sample read counts and mapping percentages. |
| `unmatched_barcodes.<jobid>.tsv` | Inputs skipped due to barcode/metadata matching failures. |
| `ambiguous_barcodes.<jobid>.tsv` | Metadata barcodes with multiple viral assignments. |
| `hypermut3/combined/master_hypermut3_all_contexts_long.csv` | Long-format Hypermut3 results. |
| `hypermut3/combined/master_hypermut3_RD_GD_AD_metadata.csv` | Wide per-read Hypermut3 table with metadata. |
| `1nt/*_mutation_counts.csv` | Per-read single-base substitution counts. |
| `2nt/*_mutation_counts.csv` | Raw dinucleotide-context counts. |
| `2nt/*_normalized_mutation_counts.csv` | Opportunity-normalized dinucleotide values. |
| `3nt/*_mutation_counts.csv` | Raw trinucleotide-context counts. |
| `3nt/*_normalized_mutation_counts.csv` | Opportunity-normalized trinucleotide values. |
| `*.read_summary.tsv` | Per-read/per-CDS stop-gain and indel summary. |
| `*.g2a_events.tsv` | Observed CDS-oriented G-to-A events. |
| `*.stop_gain_events.tsv` | APOBEC-compatible stop-gain events. |
| `*.indel_events.tsv` | CDS-overlapping indels and frame status. |
| `*.frameshift_stop_events.tsv` | Premature stops observed in frameshifted reads. |
| `stopgain_sequence_level_master.tsv` | Class-specific, gene-level master result table. |
| `stopgain_sequence_level_master.collapsed_by_read.tsv` | One row per `sample_name`/`read_id` for read-level summaries. |

Because HIV and SIV have overlapping coding regions, gene-level stop-gain outputs may have more than one row per read. Use the collapsed table for analyses requiring independent read-level observations.

## Quality control

Check that expected outputs exist and BAMs are valid:

```bash
find analysis-<run-date>/lima/trimmed -name '*.trimmed.bam' | wc -l
find analysis-<run-date>/HIV/mapped -name '*.aligned.bam' | wc -l
find analysis-<run-date>/SIV/mapped -name '*.aligned.bam' | wc -l

for bam in analysis-<run-date>/HIV/mapped/*.aligned.bam analysis-<run-date>/SIV/mapped/*.aligned.bam; do
  [ -e "$bam" ] || continue
  samtools quickcheck -v "$bam"
done
```

Inspect the wide Hypermut3 master table and RD p-value availability:

```bash
python - <<'PY'
import pandas as pd
p = 'analysis-<run-date>/hypermut3/combined/master_hypermut3_RD_GD_AD_metadata.csv'
df = pd.read_csv(p)
print(df['fisher_p_RD'].describe())
print('Missing fisher_p_RD:', df['fisher_p_RD'].isna().sum())
print(df.groupby('Virus', dropna=False).size())
PY
```

Check read classification and BAM filtering summaries:

```bash
column -t -s $'\t' analysis-<run-date>/stop-codon/logs/read_classification_summary.tsv | less -S
column -t -s $'\t' analysis-<run-date>/stop-codon/logs/bam_filtering_summary.tsv | less -S
```

Verify that the collapsed stop-gain table has unique sample/read keys:

```bash
python - <<'PY'
import pandas as pd
p = 'analysis-<run-date>/stop-codon/hypermutated/stopgain/stopgain_sequence_level_master.collapsed_by_read.tsv'
df = pd.read_csv(p, sep='\t', dtype={'sample_name': str, 'read_id': str})
print('Rows:', len(df))
print('Duplicate sample/read keys:', df.duplicated(['sample_name', 'read_id']).sum())
PY
```

## Troubleshooting

| Problem | Likely cause | Resolution |
|---|---|---|
| `conda activate` fails in SLURM | Conda is not initialized in the batch shell. | Source `$(conda info --base)/etc/profile.d/conda.sh` before activation; load the cluster's Conda module if needed. |
| `bam2fastq` or `pbindex` is missing | PacBio Toolkit utilities are unavailable in the active environment. | Install or update `pbtk`; verify executable availability with `command -v`. |
| Barcode is absent from metadata | Filename-derived barcode differs from metadata because of suffixes, punctuation, whitespace, or BOM artifacts. | Compare exact strings, clean metadata, and ensure one unambiguous viral assignment per barcode. |
| No mapped BAMs are found downstream | Mapping failed, paths are inconsistent, or BAM naming does not match expected `*.aligned.bam`. | Check mapping logs and verify all wrapper configuration blocks. |
| `pysam.FastaFile` cannot open the reference | FASTA index is missing or unreadable. | Run `samtools faidx` and verify compatible FASTA/BED contig names. |
| Hypermut3 or mutation arrays run out of memory | Per-read alignment reconstruction or wide count tables exceed task memory. | Reduce concurrency, profile representative samples, increase memory per task, or split sample sets. |
| Multiple rows appear per read in stop-gain output | The read overlaps multiple CDS features or reading frames. | Use `stopgain_sequence_level_master.collapsed_by_read.tsv` for one-row-per-read summaries. |

## Reproducibility

- Retain exact versions of scripts, FASTAs, CDS BED files, barcode/primer FASTAs, and metadata used for each run.
- Store SLURM job IDs, stdout/stderr logs, manifests, mapping summaries, and diagnostic outputs with the run directory.
- Use a unique dated analysis directory for every independent analysis.
- Preserve `sample_name` together with `read_id` in all joins because read IDs may recur across independent BAM files.
- Treat metadata corrections as analysis-affecting changes: they may change viral reference assignment, mutation calls, and downstream stop-gain interpretation.
- Export Conda environments after a successful run:

```bash
conda env export --from-history -n pacbio-isoseq > environment-pacbio-isoseq.history.yml
conda env export -n pacbio-isoseq > environment-pacbio-isoseq.lock.yml
conda env export --from-history -n hypermut3 > environment-hypermut3.history.yml
conda env export -n hypermut3 > environment-hypermut3.lock.yml
conda env export --from-history -n pacbio-r > environment-r.history.yml
conda env export -n pacbio-r > environment-r.lock.yml
```

## Citation and license

When using this workflow in a manuscript, cite the applicable versions and methods for PacBio Lima, PacBio Toolkit utilities, minimap2, samtools, pysam, pandas, SciPy, R, the R visualization packages, and the Hypermut method/software appropriate to the implementation used.

No license file is included. Add an explicit license before distribution or reuse, and ensure all raw sequence data, donor information, metadata, and data-sharing practices comply with relevant approvals and institutional policies.

## Contact

Author listed in the supplied scripts: Armando Mendez.
