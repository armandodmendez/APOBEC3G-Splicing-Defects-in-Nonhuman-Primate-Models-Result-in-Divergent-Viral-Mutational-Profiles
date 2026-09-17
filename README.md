# PacBio HIV/SIV Hypermutation and APOBEC Stop-Gain Analysis

## Overview

This repository contains a high-performance-computing (HPC) workflow for analyzing PacBio HiFi amplicon reads from HIV, HIVΔvif, and SIV experiments. The workflow begins with PacBio barcode demultiplexing and primer trimming, assigns each barcode-defined sample to the appropriate viral reference using a metadata table, maps reads, quantifies APOBEC-relevant hypermutation in multiple sequence contexts, calls mutations in 1-, 2-, and 3-nucleotide windows, partitions reads into hypermutated and nonhypermutated groups, and evaluates APOBEC-compatible stop-gain mutations in annotated coding sequences.

The provided scripts are designed for a SLURM-managed Linux cluster, but the Python programs can also be run locally after paths and scheduler-specific commands are adapted. The shell workflows are resumable: completed nonempty output files are generally reused rather than recomputed.

> **Important:** This repository contains analysis scripts, not the raw PacBio data, reference FASTA files, barcode/primer FASTAs, metadata CSV, or CDS BED annotations. Those inputs must be supplied separately and paths in the scripts must be updated for a new project location.

## Workflow at a glance

```text
PacBio HiFi barcode BAM
        |
        v
1. Lima barcode demultiplexing
        |
        v
2. Lima primer trimming
        |
        v
3. Metadata-aware HIV/SIV assignment + BAM-to-FASTQ conversion
        |
        v
4. minimap2 mapping + primary-read BAM filtering
        |
        +-------------------------------+
        |                               |
        v                               v
5. Hypermut3 RD/GD/AD analysis      6. 1-nt, 2-nt, and 3-nt mutation calling
        |                               |
        v                               v
7. Read-level hypermutation classification and merged mutation tables
        |
        v
8. Hypermutated/nonhypermutated BAM subsets
        |
        v
9. APOBEC-compatible stop-gain and frameshift-aware CDS analysis
        |
        v
10. R-based exploratory statistics and publication-ready figures
```

## Repository contents

| File | Role in the workflow |
|---|---|
| `pacbio-lima-HIV-SIV.sh` | SLURM script that demultiplexes a PacBio HiFi BAM by barcode and trims primers with PacBio Lima. |
| `HIV-SIV-mapping-2.sh` | SLURM script that parses barcode pairs from Lima output, validates the barcode against metadata, routes samples to the HIV or SIV reference, converts BAM to FASTQ, maps with minimap2, and writes mapping QC. |
| `BAM_HYPMUT3-3.py` | Python implementation of Hypermut 3.0 for identifying context-defined mutations and computing Fisher exact-test statistics from BAM alignments. |
| `hypermut3-all-4.sh` | SLURM wrapper that runs `BAM_HYPMUT3-3.py` for RD, GD, and AD contexts, merges results, and attaches metadata. |
| `mutation_counting2-6.py` | Per-read single-nucleotide substitution (SBS) counter. It also records selected `G>A` dinucleotide-derived classes. |
| `dinucleotide_mutation_analysis-5.py` | Per-read 2-nt mutation-context analysis, including directional upstream/downstream and collapsed dinucleotide classes, with normalization by available sites. |
| `trinucleotide_mutation_analysis-8.py` | Per-read 3-nt mutation-context analysis, with raw and opportunity-normalized counts. |
| `mutational-calling-all-contexts-7.sh` | Submitter workflow that generates a metadata-aware manifest, launches SLURM arrays for 1-nt, 2-nt, and 3-nt analyses, and combines outputs with metadata and Hypermut3 calls. |
| `apobec_stopgain_pipeline-9.py` | Unified CDS-aware analysis of G-to-A events, APOBEC-compatible stop gains, all premature stops, indels, and frameshift-associated premature stops. |
| `collapse_stopgain_master-10.py` | Collapses overlapping-CDS rows to one row per `sample_name`/`read_id`, summing stop-gain counts and denominators while preserving sample metadata. |
| `stop-codon-hypermut-nonhypermut-11.sh` | SLURM workflow that classifies reads by Hypermut3 RD Fisher p-value, filters BAMs by class, runs stop-gain analysis, and produces collapsed read-level outputs. |
| `PacBio-HIV-SIV-2026-08-31-12.R` | R analysis script for threshold sensitivity, mutation-ratio QC, tables, scatterplots, and pooled SBS summaries. |

## Analysis logic

### Viral-reference assignment

Sample identity is inferred from the barcode pair embedded in each Lima-generated filename. The workflows require a metadata CSV containing at least these exact column names:

```text
barcode,Virus
```

The barcode value must match the barcode pair parsed from the filename, for example:

```text
bc1004_F2--bc1057_R2
```

The `Virus` value determines reference selection:

| Metadata `Virus` value | Analysis group | Reference used |
|---|---|---|
| `HIV` | HIV | HIV near-full-length reference FASTA |
| `HIVDvif` | HIV | HIV near-full-length reference FASTA |
| `SIV` | SIV | SIV near-full-length reference FASTA |

Metadata is treated as authoritative. Scripts explicitly check for missing barcodes, unsupported virus labels, ambiguous barcode assignments, UTF-8 byte-order marks (BOM), Windows line endings, and missing required metadata columns. Diagnostic tables are written for unmatched and ambiguous barcodes.

### Hypermut3 contexts

The Hypermut3 wrapper evaluates each mapped sample in three contexts:

| Context | Intended use |
|---|---|
| `RD` | Broad context used for the primary Fisher-test hypermutation classification. |
| `GD` | GG-associated G-to-A context analysis. |
| `AD` | GA-associated G-to-A context analysis. |

The exact command-line context definitions are configured in `hypermut3-all-4.sh`. Review and preserve those settings when reproducing the analysis because a changed context definition changes the biological interpretation.

For downstream stop-codon partitioning, a read is classified as:

- **Hypermutated:** `fisher_p_RD < 0.05`
- **Nonhypermutated:** `fisher_p_RD >= 0.05`
- **Unclassified:** missing, invalid, or nonnumeric `fisher_p_RD`; these reads are reported in QC output and are not included in either subset.

### Mutation calling

The mutational-calling workflow analyzes primary mapped reads independently in three modes:

- **1-nt analysis:** raw per-read counts for the 12 possible nonidentity single-base substitutions, with a stable output schema.
- **2-nt analysis:** evaluates both `minus1_plus1` and `zero_plus2` windows; reports full context substitutions, directional mutation classes, collapsed dinucleotide classes, and opportunity-normalized values.
- **3-nt analysis:** evaluates trinucleotide windows surrounding single-base differences; writes raw counts and counts normalized to available valid sites in each read.

Reads marked duplicate or unmapped are excluded by the Python mutation scripts. Gapped/incomplete motif windows are not counted. For normalized outputs, the denominator is the number of valid, ungapped reference windows available in that read for the relevant analysis window.

### Stop-gain analysis

The stop-gain workflow performs coding-sequence-aware analysis after Hypermut3 read classification. It:

1. Builds read-name lists separately for every sample using `sample_name` plus read ID as the logical key.
2. Filters each aligned BAM using `samtools view -N` into hypermutated and nonhypermutated subsets.
3. Uses the appropriate HIV or SIV reference FASTA and CDS BED annotation.
4. Handles both plus- and minus-strand CDS features and respects reading frame.
5. Identifies observed G-to-A events in CDS coordinates.
6. Identifies APOBEC-compatible stop-gain opportunities: G positions at codon position 1 or 2 in a `GG` or `GA` motif for which a G-to-A substitution produces `TAA`, `TAG`, or `TGA`.
7. Records all premature stops, indels, frameshift indels, and frameshift-associated premature stops.
8. Writes per-gene/per-read sequence tables and collapses them to one row per sample/read for analyses requiring non-overlapping read-level observations.

Because HIV and SIV genomes contain overlapping reading frames, the initial sequence-level table can contain more than one row for a read. Use the collapsed table for one-row-per-read summaries, while retaining the gene-level table for feature-specific analyses.

## Requirements

### Platform requirements

- Linux or a Linux-compatible HPC environment.
- Bash 4+.
- SLURM for the provided `.sh` submission scripts (`sbatch`, job arrays, and job dependencies are used).
- Conda, Miniconda, or Mambaforge.
- A working C/C++ runtime compatible with the installed bioinformatics packages.
- Sufficient storage for BAM, FASTQ, sorted BAM, index, intermediate CSV/TSV, and logs.

The scripts request substantial resources for full-scale runs. The shipped settings range from 16 GB to 800 GB memory and from 1 to 192 CPU cores depending on the step. These are cluster-specific starting values, not universal requirements. Benchmark a representative sample and adjust `#SBATCH` settings, array concurrency, and `--threads` to match your data volume and scheduler policy.

### Command-line tools

The core workflow expects these programs in `PATH`:

```text
lima
bam2fastq
pbindex
minimap2
samtools
python3
Rscript
awk
sed
tr
grep
find
sort
wc
head
basename
dirname
gzip
sbatch
```

`pbindex` is treated as optional by the mapping script; it warns if unavailable. `awk`, `sed`, `tr`, `grep`, `find`, `sort`, `wc`, `head`, `basename`, `dirname`, and `gzip` are normally provided by the operating system or a standard HPC environment. `sbatch` is only required for the SLURM orchestration scripts.

### Python packages

The supplied Python scripts require:

```text
python >=3.10,<3.13
pandas
pysam
scipy
numpy
```

`BAM_HYPMUT3-3.py` specifically imports `pysam` and `scipy.stats.fisher_exact`. The mutation-calling scripts and stop-gain scripts use `pandas` and `pysam`. `numpy` is included explicitly for a robust scientific Python installation and future extensions.

### R packages

The R script imports or uses the following packages:

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

The script uses `tidyr::crossing()` and `pivot_wider()`, so `tidyr` must be installed even where it is called with an explicit namespace rather than loaded with `library(tidyr)`.

## Conda environment setup

### Required environment separation

Users must create and activate Conda environments containing the exact tool families required by the workflow **before** launching any pipeline script. The original shell scripts activate environments named `pacbio-isoseq` and `hypermut3`; those names are used below to minimize modifications.

A two-environment design is recommended:

- `pacbio-isoseq`: PacBio Lima, PacBio BAM conversion/indexing utilities, minimap2, samtools, and a basic Python runtime for preprocessing/mapping.
- `hypermut3`: Python analysis libraries plus samtools for Hypermut3, mutation calling, and stop-gain analysis.

This separation reduces dependency conflicts between PacBio tooling and Python/scientific packages. If preferred, a single unified environment can be created, but then every shell script must activate that same environment name.

### Channel configuration

Create environments with strict channel priority to reduce cross-channel binary incompatibilities:

```bash
conda config --set channel_priority strict
```

If you have not already configured channels, either include `-c conda-forge -c bioconda` in every command below or configure them once:

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

### Option A: create environments from commands

Create the PacBio/mapping environment:

```bash
conda create -y -n pacbio-isoseq \
  -c conda-forge -c bioconda \
  python=3.11 \
  lima \
  pbtk \
  minimap2 \
  samtools \
  pysam \
  pandas \
  numpy
```

Create the Hypermut3/analysis environment:

```bash
conda create -y -n hypermut3 \
  -c conda-forge -c bioconda \
  python=3.11 \
  samtools \
  pysam \
  pandas \
  scipy \
  numpy
```

Create the R plotting environment:

```bash
conda create -y -n pacbio-r \
  -c conda-forge -c bioconda \
  r-base=4.3 \
  r-readr \
  r-dplyr \
  r-ggplot2 \
  r-scales \
  r-svglite \
  r-tidyr \
  r-viridis \
  r-patchwork
```

Activate an environment as needed:

```bash
conda activate pacbio-isoseq
# or
conda activate hypermut3
# or
conda activate pacbio-r
```

`pbtk` provides PacBio Toolkit utilities including `bam2fastq` and typically `pbindex`. Package availability can vary across channels and platforms; verify executable availability after installation rather than assuming a package name provides a particular binary.

### Option B: reproducible YAML files

Save the following files in the project root and create the environments from them.

#### `environment-pacbio-isoseq.yml`

```yaml
name: pacbio-isoseq
channels:
  - conda-forge
  - bioconda
channel_priority: strict
dependencies:
  - python=3.11
  - lima
  - pbtk
  - minimap2
  - samtools
  - pysam
  - pandas
  - numpy
```

Create it with:

```bash
conda env create -f environment-pacbio-isoseq.yml
```

#### `environment-hypermut3.yml`

```yaml
name: hypermut3
channels:
  - conda-forge
  - bioconda
channel_priority: strict
dependencies:
  - python=3.11
  - samtools
  - pysam
  - pandas
  - scipy
  - numpy
```

Create it with:

```bash
conda env create -f environment-hypermut3.yml
```

#### `environment-r.yml`

```yaml
name: pacbio-r
channels:
  - conda-forge
  - bioconda
channel_priority: strict
dependencies:
  - r-base=4.3
  - r-readr
  - r-dplyr
  - r-ggplot2
  - r-scales
  - r-svglite
  - r-tidyr
  - r-viridis
  - r-patchwork
```

Create it with:

```bash
conda env create -f environment-r.yml
```

### Validate installations

Run the following after environment creation. Each command should resolve to the intended Conda environment and report a version.

```bash
conda activate pacbio-isoseq
which lima bam2fastq minimap2 samtools
lima --version
bam2fastq --help | head
minimap2 --version
samtools --version
command -v pbindex && pbindex --help | head || true

conda activate hypermut3
python --version
python - <<'PY'
import numpy
import pandas
import pysam
import scipy
from scipy.stats import fisher_exact

print("numpy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("pysam:", pysam.__version__)
print("scipy:", scipy.__version__)
print("Fisher exact test import: OK")
PY
samtools --version

conda activate pacbio-r
Rscript -e 'pkgs <- c("readr","dplyr","ggplot2","scales","svglite","tidyr","viridis","patchwork"); stopifnot(all(sapply(pkgs, requireNamespace, quietly=TRUE))); sessionInfo()'
```

### Conda initialization in batch jobs

Noninteractive SLURM jobs frequently do not initialize Conda automatically. At the start of each script, replace the original user-specific initialization path with the correct installation path for your system, for example:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate pacbio-isoseq
```

For scripts that run Hypermut3, mutation calling, or stop-gain analysis, use:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hypermut3
```

On clusters that provide Conda as an environment module, the required sequence may be:

```bash
module load miniconda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hypermut3
```

Do not assume that `conda activate` works in an `sbatch` job unless `conda.sh` has been sourced successfully.

## Input files and layout

The original scripts use absolute project paths and expect a layout similar to the following. This tree is a recommended organization, not a requirement; update all script configuration blocks consistently if you use another layout.

```text
HIV-SIV-Hypermutation-Assay-PacBio/
├── RAW-DATA/
│   └── 2026-08-31/
│       └── <PacBio_HiFi_barcode_BAM>.bam
├── reference-genomes/
│   ├── pJH048_i3_HIV-NFL.fasta
│   ├── SIV239-SpX_i22-NFL.fasta
│   ├── pJH048_i3_-CR_extraction-HIV-NFL-CDS.bed
│   └── RSR_Exp103_SIV239-SpX_i22_extraction-NFL-CDS.bed
├── analysis-2026-08-31/
│   ├── primers-barcodes/
│   │   ├── Aug2026_barcodes.fasta
│   │   └── primers_unique.fasta
│   ├── metadata/
│   │   └── Aug2026_barcode-pcr_type-layout.csv
│   ├── lima/
│   ├── HIV/
│   │   └── mapped/
│   ├── SIV/
│   │   └── mapped/
│   ├── hypermut3/
│   ├── mutational-calling/
│   └── stop-codon/
├── scripts/
│   ├── pacbio-lima-HIV-SIV.sh
│   ├── HIV-SIV-mapping-2.sh
│   ├── hypermut3-all-4.sh
│   ├── mutational-calling-all-contexts-7.sh
│   ├── stop-codon-hypermut-nonhypermut-11.sh
│   ├── BAM_HYPMUT3-3.py
│   ├── mutation_counting2-6.py
│   ├── dinucleotide_mutation_analysis-5.py
│   ├── trinucleotide_mutation_analysis-8.py
│   ├── apobec_stopgain_pipeline-9.py
│   └── collapse_stopgain_master-10.py
└── README.md
```

### Required biological inputs

| Input | Required by | Notes |
|---|---|---|
| PacBio HiFi BAM | Lima preprocessing | The scripts expect a BAM rather than raw FASTQ. |
| Barcode FASTA | Lima preprocessing | Supply barcode sequences compatible with the library design. |
| Primer FASTA | Lima preprocessing | Supply unique primer sequences for post-demultiplex primer trimming. |
| Metadata CSV | Mapping and all downstream workflows | Must include exact `barcode` and `Virus` headers; additional columns are retained and merged into final outputs. |
| HIV FASTA | HIV and HIVΔvif reads | Must match the reference used by minimap2 and associated CDS BED. |
| SIV FASTA | SIV reads | Must match the reference used by minimap2 and associated CDS BED. |
| HIV CDS BED | Stop-gain workflow | Requires CDS coordinates, strand, and frame in a supported BED-like format. |
| SIV CDS BED | Stop-gain workflow | Requires CDS coordinates, strand, and frame in a supported BED-like format. |

### Naming conventions

The shell workflows depend on consistent filenames. The expected lineage is:

```text
Lima demultiplexed BAM:
DNA_NFL.demux.<barcode>.bam

Lima primer-trimmed BAM:
DNA_NFL.demux.<barcode>.trimmed.bam

Mapped BAM:
DNA_NFL.demux.<barcode>.aligned.bam
```

Example:

```text
DNA_NFL.demux.bc1004_F2--bc1057_R2.trimmed.bam
DNA_NFL.demux.bc1004_F2--bc1057_R2.aligned.bam
```

The current barcode pattern includes an optional numeric suffix after `F` and `R`, such as `F2` and `R2`. The stop-gain Python script also supports the older style without those suffixes.

## Before running

### 1. Make paths portable

Every shell script contains a configuration block near the beginning. Update, at minimum:

```bash
PROJECT_DIR="/path/to/HIV-SIV-Hypermutation-Assay-PacBio"
ANALYSIS_DIR="${PROJECT_DIR}/analysis-<run-date>"
```

Also update paths to:

- Raw PacBio BAM input.
- Barcode and primer FASTA files.
- Metadata CSV.
- HIV and SIV reference FASTAs.
- HIV and SIV CDS BED files.
- The locations of the Python scripts.
- SLURM working directory and log paths.
- Conda initialization path.

Avoid mixing files from different runs by giving each analysis a distinct `ANALYSIS_DIR`.

### 2. Place scripts where wrappers expect them

The submitted shell wrappers refer to generic script names such as:

```text
BAM_HYPMUT3.py
mutation_counting2.py
dinucleotide_mutation_analysis.py
trinucleotide_mutation_analysis.py
apobec_stopgain_pipeline.py
collapse_stopgain_master.py
```

The files distributed here have versioned filenames. Either update the wrapper variables to point directly at the versioned files, or copy/symlink them to the names expected by the shell scripts. For example:

```bash
mkdir -p "${PROJECT_DIR}/scripts"
cp BAM_HYPMUT3-3.py "${PROJECT_DIR}/scripts/BAM_HYPMUT3.py"
cp mutation_counting2-6.py "${PROJECT_DIR}/scripts/mutation_counting2.py"
cp dinucleotide_mutation_analysis-5.py "${PROJECT_DIR}/scripts/dinucleotide_mutation_analysis.py"
cp trinucleotide_mutation_analysis-8.py "${PROJECT_DIR}/scripts/trinucleotide_mutation_analysis.py"
cp apobec_stopgain_pipeline-9.py "${PROJECT_DIR}/scripts/apobec_stopgain_pipeline.py"
cp collapse_stopgain_master-10.py "${PROJECT_DIR}/scripts/collapse_stopgain_master.py"
```

Or use symlinks to preserve provenance:

```bash
ln -sfn "$(pwd)/BAM_HYPMUT3-3.py" "${PROJECT_DIR}/scripts/BAM_HYPMUT3.py"
ln -sfn "$(pwd)/mutation_counting2-6.py" "${PROJECT_DIR}/scripts/mutation_counting2.py"
ln -sfn "$(pwd)/dinucleotide_mutation_analysis-5.py" "${PROJECT_DIR}/scripts/dinucleotide_mutation_analysis.py"
ln -sfn "$(pwd)/trinucleotide_mutation_analysis-8.py" "${PROJECT_DIR}/scripts/trinucleotide_mutation_analysis.py"
ln -sfn "$(pwd)/apobec_stopgain_pipeline-9.py" "${PROJECT_DIR}/scripts/apobec_stopgain_pipeline.py"
ln -sfn "$(pwd)/collapse_stopgain_master-10.py" "${PROJECT_DIR}/scripts/collapse_stopgain_master.py"
```

Then set wrapper variables such as:

```bash
PYTHON_SCRIPT_DIR="${PROJECT_DIR}/scripts"
HYPMUT3_SCRIPT="${PYTHON_SCRIPT_DIR}/BAM_HYPMUT3.py"
```

### 3. Index references where required

`pysam.FastaFile` requires a FASTA index. Create `.fai` files before running the stop-gain workflow:

```bash
conda activate hypermut3
samtools faidx reference-genomes/pJH048_i3_HIV-NFL.fasta
samtools faidx reference-genomes/SIV239-SpX_i22-NFL.fasta
```

The mapping process will create sorted and indexed BAMs; do not manually edit BAM filenames after downstream manifests have been created.

### 4. Validate metadata

Before a full run, verify exact headers and barcode uniqueness:

```bash
python - <<'PY'
import pandas as pd

metadata = pd.read_csv(
    "analysis-2026-08-31/metadata/Aug2026_barcode-pcr_type-layout.csv",
    dtype=str,
    encoding="utf-8-sig",
)
metadata.columns = metadata.columns.str.strip()
required = {"barcode", "Virus"}
missing = required - set(metadata.columns)
if missing:
    raise SystemExit(f"Missing required metadata columns: {sorted(missing)}")

metadata["barcode"] = metadata["barcode"].astype(str).str.strip()
metadata["Virus"] = metadata["Virus"].astype(str).str.strip()
print(metadata[["barcode", "Virus"]].head())
print("Rows:", len(metadata))
print("Unique barcodes:", metadata["barcode"].nunique())
print("Duplicate barcode rows:", metadata["barcode"].duplicated(keep=False).sum())
print("Virus labels:", sorted(metadata["Virus"].dropna().unique()))
PY
```

## Running the pipeline

### Step 1: Demultiplex and trim PacBio reads

Activate the PacBio environment and submit the Lima workflow:

```bash
conda activate pacbio-isoseq
sbatch scripts/pacbio-lima-HIV-SIV.sh
```

This script:

- Runs Lima with `--hifi-preset ASYMMETRIC` and `--split-named` to demultiplex barcoded HiFi reads.
- Discovers barcode-split BAM files from the Lima prefix.
- Runs a second Lima pass to trim primers from each barcode-defined BAM.
- Writes trimmed BAM files beneath `analysis-<run-date>/lima/trimmed/`.

Inspect the SLURM `.out` and `.err` logs before proceeding. The script has `set -euo pipefail`, so an unhandled failure should stop the job.

### Step 2: Route, convert, map, and QC

Submit the metadata-aware mapping workflow:

```bash
conda activate pacbio-isoseq
sbatch scripts/HIV-SIV-mapping-2.sh
```

For each primer-trimmed BAM, this workflow:

1. Extracts the barcode pair from the filename.
2. Finds the matching metadata row using the `barcode` field.
3. Uses the metadata `Virus` field to select HIV or SIV routing and reference sequence.
4. Converts the PacBio BAM to compressed FASTQ using `bam2fastq`.
5. Maps HiFi reads with minimap2.
6. Retains primary mapped reads using `samtools view -F 2308`.
7. Produces aligned BAMs in `HIV/mapped/` or `SIV/mapped/`.
8. Writes a mapping summary using pre-filter FASTQ read counts and final primary-mapped BAM read counts.

Do **not** calculate mapping efficiency by running `samtools flagstat` only on the final filtered BAM. Unmapped reads have already been removed, so the mapped percentage is expected to be artificially close to 100%.

### Step 3: Run Hypermut3 in RD, GD, and AD contexts

Update `HYPMUT3_SCRIPT` in the wrapper to the location of `BAM_HYPMUT3-3.py` or its standardized copy, then submit:

```bash
conda activate hypermut3
sbatch scripts/hypermut3-all-4.sh
```

The wrapper scans `HIV/mapped/*.aligned.bam` and `SIV/mapped/*.aligned.bam`, runs each required context unless a nonempty summary sentinel already exists, and creates:

```text
analysis-<run-date>/hypermut3/
├── RD/
├── GD/
├── AD/
├── combined/
│   ├── master_hypermut3_all_contexts_long.csv
│   ├── master_hypermut3_RD_GD_AD_metadata.csv
│   └── hypermut3_run_summary.tsv
├── logs/
└── manifests/
```

The wide master table is the key input for downstream mutation-calling joins, stop-codon classification, and R visualization.

### Step 4: Call 1-nt, 2-nt, and 3-nt mutations

After Hypermut3 completes successfully, submit:

```bash
conda activate hypermut3
sbatch scripts/mutational-calling-all-contexts-7.sh
```

This is a submitter job: it constructs a manifest of metadata-validated BAMs, dynamically writes worker and concatenation scripts, launches independent SLURM arrays for the three mutation-context analyses, and schedules the combiner after successful arrays.

Expected root output:

```text
analysis-<run-date>/mutational-calling/
├── 1nt/
├── 2nt/
├── 3nt/
├── combined/
├── logs/
└── manifests/
```

The workflow is resumable at a sample/context level. It skips a sample/context when its expected completion sentinel is already nonempty.

### Step 5: Classify reads and analyze stop gains

After the wide Hypermut3 master CSV exists, submit:

```bash
conda activate hypermut3
sbatch scripts/stop-codon-hypermut-nonhypermut-11.sh
```

This workflow uses `fisher_p_RD` to make per-sample read lists, filters source BAMs into class-specific BAMs, runs the unified stop-gain program for each class, and collapses overlapping-CDS results.

Expected output layout:

```text
analysis-<run-date>/stop-codon/
├── hypermutated/
│   ├── filtered-bams/
│   ├── read-lists/
│   └── stopgain/
│       ├── stopgain_sequence_level_master.tsv
│       └── stopgain_sequence_level_master.collapsed_by_read.tsv
├── nonhypermutated/
│   ├── filtered-bams/
│   ├── read-lists/
│   └── stopgain/
│       ├── stopgain_sequence_level_master.tsv
│       └── stopgain_sequence_level_master.collapsed_by_read.tsv
├── logs/
└── manifests/
```

### Step 6: R analysis and figures

Copy or link the relevant combined CSV files to the location expected by the R script, or edit the `setwd()` and `input_file` definitions at the top of `PacBio-HIV-SIV-2026-08-31-12.R`.

Run interactively in RStudio or noninteractively:

```bash
conda activate pacbio-r
Rscript scripts/PacBio-HIV-SIV-2026-08-31-12.R
```

The script includes, among other analyses:

- RD p-value threshold sensitivity from 0.01 through 0.10.
- Counts and percentages of hypermutated reads by donor and by donor/virus.
- RD, GD, and AD mutation-ratio calculation and quality control.
- Filtering to selected donors (`N03`, `N04`, `N06`, `N10`), `DPI <= 7`, RD-significant reads, and valid GD/AD denominators.
- Combined and stratified GD-versus-AD scatterplots in PNG and SVG formats.
- Pooled raw single-base-substitution distributions for HIV versus HIVΔvif DNA reads.

## Key outputs

### Mapping outputs

| Output | Description |
|---|---|
| `HIV/*.fastq.gz`, `SIV/*.fastq.gz` | BAM-derived PacBio read sequences routed by metadata-defined virus group. |
| `HIV/mapped/*.aligned.bam`, `SIV/mapped/*.aligned.bam` | Primary mapped reads retained after minimap2 alignment and SAM flag filtering. |
| `mapping_summary.<jobid>.tsv` | Per-sample pre-filter FASTQ reads, retained primary mapped reads, and mapping percentage. |
| `unmatched_barcodes.<jobid>.tsv` | Inputs skipped because filename/barcode/metadata matching failed. |
| `ambiguous_barcodes.<jobid>.tsv` | Barcodes assigned to more than one `Virus` category in metadata. |

### Hypermut3 outputs

| Output | Description |
|---|---|
| `RD/`, `GD/`, `AD/` | Per-context, per-sample Hypermut3 outputs. |
| `master_hypermut3_all_contexts_long.csv` | Long-format merged result table across contexts. |
| `master_hypermut3_RD_GD_AD_metadata.csv` | Wide per-read master table with context metrics and all metadata columns. |
| `hypermut3_run_summary.tsv` | Execution status and summary for processed/skipped samples and contexts. |

### Mutation-calling outputs

| Output | Description |
|---|---|
| `1nt/*_mutation_counts.csv` | Raw per-read counts of nonidentity single-base substitutions. |
| `2nt/*_mutation_counts.csv` | Raw full-context dinucleotide mutation counts. |
| `2nt/*_normalized_mutation_counts.csv` | Opportunity-normalized dinucleotide context counts. |
| `2nt/*_dinucleotide_mutation_class_counts.csv` | Directional upstream/downstream dinucleotide mutation-class counts. |
| `2nt/*_collapsed_dinucleotide_mutation_class_counts.csv` | Context-collapsed dinucleotide class counts. |
| `3nt/*_mutation_counts.csv` | Raw per-read trinucleotide context counts. |
| `3nt/*_normalized_mutation_counts.csv` | Opportunity-normalized trinucleotide counts. |
| `combined/` | Merged mutation tables with metadata and Hypermut3 RD p-values. |

### Stop-gain outputs

| Output | Description |
|---|---|
| `*.read_summary.tsv` | Per-read, per-CDS analysis summary including G-to-A count, stop gains, indels, and frameshift status. |
| `*.g2a_events.tsv` | All observed CDS-oriented G-to-A events. |
| `*.stop_gain_events.tsv` | APOBEC-compatible G-to-A stop-gain events. |
| `*.indel_events.tsv` | CDS-overlapping insertion/deletion events and frameshift flags. |
| `*.frameshift_stop_events.tsv` | Premature stops observed in reads carrying frameshift indels. |
| `*.sequence_stopgain_master.tsv` | Per-sample/per-read/per-gene master table. |
| `stopgain_sequence_level_master.tsv` | Concatenated sequence-level master table for each Hypermut3 class. |
| `stopgain_sequence_level_master.collapsed_by_read.tsv` | One row per `sample_name` and `read_id`, appropriate for read-level summaries. |

## Quality-control checkpoints

Run the following checks before treating an analysis as complete.

### Preprocessing and mapping

```bash
# Number of primer-trimmed barcode BAMs
find analysis-2026-08-31/lima/trimmed -name '*.trimmed.bam' | wc -l

# Number of aligned HIV and SIV sample BAMs
find analysis-2026-08-31/HIV/mapped -name '*.aligned.bam' | wc -l
find analysis-2026-08-31/SIV/mapped -name '*.aligned.bam' | wc -l

# Confirm BAM integrity/index readability
for bam in analysis-2026-08-31/HIV/mapped/*.aligned.bam analysis-2026-08-31/SIV/mapped/*.aligned.bam; do
  [ -e "$bam" ] || continue
  samtools quickcheck -v "$bam"
done
```

Review unmatched and ambiguous barcode logs. No sample should silently proceed with a mismatched barcode/reference assignment.

### Hypermut3 and mutation calling

```bash
# Inspect the first columns and records of the wide master file
head -n 3 analysis-2026-08-31/hypermut3/combined/master_hypermut3_RD_GD_AD_metadata.csv

# Ensure RD p-values are parseable and assess missingness
python - <<'PY'
import pandas as pd
p = 'analysis-2026-08-31/hypermut3/combined/master_hypermut3_RD_GD_AD_metadata.csv'
df = pd.read_csv(p)
print(df['fisher_p_RD'].describe())
print('Missing fisher_p_RD:', df['fisher_p_RD'].isna().sum())
print(df.groupby('Virus', dropna=False).size())
PY
```

Confirm that 1-nt, 2-nt, and 3-nt array status logs show successful completion for all expected samples. Re-run only incomplete work after resolving errors; the wrappers are built to reuse valid output sentinels.

### Stop-gain analysis

```bash
# Read-classification counts by sample
column -t -s $'\t' analysis-2026-08-31/stop-codon/logs/read_classification_summary.tsv | less -S

# Filtering status
column -t -s $'\t' analysis-2026-08-31/stop-codon/logs/bam_filtering_summary.tsv | less -S

# Verify one row per sample/read in collapsed table
python - <<'PY'
import pandas as pd
p = 'analysis-2026-08-31/stop-codon/hypermutated/stopgain/stopgain_sequence_level_master.collapsed_by_read.tsv'
df = pd.read_csv(p, sep='\t', dtype={'sample_name': str, 'read_id': str})
print('Rows:', len(df))
print('Duplicate sample/read keys:', df.duplicated(['sample_name', 'read_id']).sum())
PY
```

## Reproducibility recommendations

- Preserve the exact reference FASTA, CDS BED, barcode FASTA, primer FASTA, metadata CSV, and script version used for each analysis run.
- Record Conda package versions after a successful run:

```bash
conda env export --from-history -n pacbio-isoseq > environment-pacbio-isoseq.history.yml
conda env export -n pacbio-isoseq > environment-pacbio-isoseq.lock.yml
conda env export --from-history -n hypermut3 > environment-hypermut3.history.yml
conda env export -n hypermut3 > environment-hypermut3.lock.yml
conda env export --from-history -n pacbio-r > environment-r.history.yml
conda env export -n pacbio-r > environment-r.lock.yml
```

- Store SLURM job IDs, stdout/stderr files, mapping summaries, manifests, and diagnostic logs with the analysis directory.
- Use a new dated analysis directory for each independent run rather than overwriting outputs.
- Keep `sample_name` together with `read_id` in all merges and downstream analyses. Read IDs are not guaranteed to be globally unique across independent BAM files.
- Do not reinterpret a filtered primary-mapped BAM as representing all sequenced reads. Compute mapping efficiency from counts before unmapped reads are discarded.
- Treat metadata changes as analysis-affecting changes. A corrected barcode or virus label can alter reference assignment, all downstream mutation calls, and stop-gain interpretation.

## Common troubleshooting

### `conda activate` fails in SLURM

Cause: Conda was not initialized in the noninteractive batch shell.

Fix:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hypermut3
```

If `conda` itself is not found, load your site-specific Conda module first or use the absolute path to its installation.

### `bam2fastq` or `pbindex` is not found

Cause: PacBio Toolkit utilities are absent or a different package build is installed.

Fix:

```bash
conda activate pacbio-isoseq
conda install -y -c conda-forge -c bioconda pbtk
command -v bam2fastq
command -v pbindex
```

If your environment provides `bam2fastq` but not `pbindex`, mapping can still continue because the mapping script treats `pbindex` as optional.

### Metadata says barcode not found

Cause: The trimmed/mapped filename barcode does not exactly match metadata, including `F2`/`R2` suffixes, punctuation, or whitespace.

Fix: Compare the filename-derived barcode against the `barcode` column, remove accidental whitespace or BOM artifacts from metadata, and ensure each barcode has one unambiguous virus assignment.

### No mapped BAMs found

Cause: The mapping step did not complete, outputs are in a different directory, or filename suffixes differ.

Fix: Confirm the paths and `*.aligned.bam` naming convention used by `hypermut3-all-4.sh`, `mutational-calling-all-contexts-7.sh`, and `stop-codon-hypermut-nonhypermut-11.sh`.

### `pysam.FastaFile` cannot open a FASTA

Cause: The `.fai` index is missing or the FASTA is unreadable.

Fix:

```bash
samtools faidx /path/to/reference.fasta
```

Ensure the BED contig names are compatible with the corresponding FASTA contig names. The stop-gain script has limited alias handling, but a reference/annotation mismatch should be corrected at the source.

### Hypermut3 or mutation arrays exhaust memory

Cause: Per-read alignment reconstruction and wide count matrices can be memory intensive for high-depth amplicon data.

Fix: Reduce `MAX_CONCURRENT`, lower SLURM array concurrency, increase memory per task after profiling, or split the sample set. Do not simply request more CPUs; these Python loops may not scale with CPU count unless the wrapper explicitly parallelizes work.

### Stop-gain master table has multiple rows per read

This is expected when a read overlaps multiple CDS features. Use:

```text
stopgain_sequence_level_master.collapsed_by_read.tsv
```

for one row per `sample_name`/`read_id`, and use the noncollapsed master when the gene/reading-frame-specific context is biologically important.

## Citation and attribution

If this workflow is used in a manuscript, cite the relevant underlying software and methods, including PacBio Lima, PacBio Toolkit utilities, minimap2, samtools, pysam, pandas, SciPy, R, and the R packages used for visualization. Cite the Hypermut method/software source appropriate to the version and implementation used in your study.

## License and data use

No license file is included with the supplied scripts. Before redistribution, collaboration, or reuse, add an explicit license and ensure that raw sequence data, sample metadata, donor information, and institutional data-sharing requirements are handled according to applicable approvals and policies.

## Contact

Author listed in the provided scripts: Armando Mendez.
