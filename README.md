# APOBEC3G-Splicing-Defects-in-Nonhuman-Primate-Models-Result-in-Divergent-Viral-Mutational-Profiles

## Overview

This repository contains scripts for a PacBio HiFi sequencing pipeline designed to study HIV and SIV hypermutation. The workflow covers barcode demultiplexing, primer trimming, reference genome alignment, and mutational analysis using APOBEC3-context mutation counting and Hypermut 3.0.

A companion Ribo-seq analysis pipeline is also included for APOBEC3G translational profiling in rhesus macaque and human samples.

---

## Repository Structure

```
├── scripts/
│   ├── pacbio-lima-HIV-SIV.sh                    # Step 1: Barcode demultiplexing
│   ├── pacbio-lima-primer-trim-HIV-SIV.sh        # Step 2: Primer trimming
│   ├── pacbio-HIV-SIV-preprocessing-mapping.sh   # Step 3: Preprocessing & alignment (RNA + DNA)
│   ├── pacbio-DNA-HIV-SIV-preprocessing-mapping.sh  # Step 3 (DNA only variant)
│   ├── mutation-counting.sh                      # Step 4: SLURM array for mutation counting
│   ├── mutation_counting2.py                     # Python: per-read mutation counting
│   ├── BAM_HYPMUT3.py                            # Python: Hypermut 3.0 (BAM-based)
│   ├── combine_hypermut3_final.sh                # Step 5: Combine Hypermut3 outputs
│   └── ribo-seq-analysis.R                       # R: Ribo-seQC + ggRibo analysis
├── reference-genomes/
│   ├── pJH048_i3_HIV-NFL.fasta                   # HIV near-full-length reference
│   └── SIV239-SpX_i22-NFL.fasta                  # SIV near-full-length reference
├── LICENSE
└── README.md
```

---

## Dependencies

### Conda Environment: `pacbio-isoseq`
| Tool | Purpose |
|------|---------|
| `lima` | Barcode demultiplexing and primer trimming |
| `pbindex` | PacBio BAM indexing |
| `bam2fastq` | BAM to FASTQ conversion |
| `minimap2` | HiFi long-read alignment |
| `samtools` | BAM filtering, sorting, and indexing |

### Conda Environment: `hypermut3`
| Tool | Purpose |
|------|---------|
| `python3` + `pysam` | BAM parsing |
| `scipy` | Fisher's exact test (Hypermut 3.0) |
| `pandas` | Mutation count aggregation |

### R Packages
| Package | Purpose |
|---------|---------|
| `RiboseQC` | Ribosome profiling QC and P-site analysis |
| `ggRibo` | Visualization of ribosome occupancy |
| `rtracklayer` | GTF import/export |
| `Rsamtools` | BAM header inspection |
| `GenomicRanges` | Genomic interval operations |

---

## Pipeline

### Step 1 — Barcode Demultiplexing (`pacbio-lima-HIV-SIV.sh`)

Demultiplexes raw PacBio HiFi BAM files into RNA and DNA fractions using `lima` with asymmetric barcode mode.

**Inputs:**
- `RNA.m84193_*.hifi_reads.bc2015.bam` — raw RNA HiFi reads
- `DNA.m84193_*.hifi_reads.bc2014.bam` — raw DNA HiFi reads
- `barcodes.fasta` — sample barcodes

**Outputs:**
- `split-barcode/RNA/RNA.demux.*.bam`
- `split-barcode/DNA/DNA.demux.*.bam`

```bash
lima RNA.*.hifi_reads.bc2015.bam barcodes.fasta split-barcode/RNA/RNA.demux.bam \
  --hifi-preset ASYMMETRIC --split-named

lima DNA.*.hifi_reads.bc2014.bam barcodes.fasta split-barcode/DNA/DNA.demux.bam \
  --hifi-preset ASYMMETRIC --split-named
```

---

### Step 2 — Primer Trimming (`pacbio-lima-primer-trim-HIV-SIV.sh`)

Trims amplicon primers from demultiplexed BAM files for both RNA and DNA fractions.

**Inputs:** Demultiplexed BAMs from Step 1, `primers_unique.fasta`

**Outputs:** `split-barcode/{RNA,DNA}/trimmed/*.trimmed.bam`

---

### Step 3 — Preprocessing & Alignment

**Scripts:**
- `pacbio-HIV-SIV-preprocessing-mapping.sh` — processes both RNA and DNA
- `pacbio-DNA-HIV-SIV-preprocessing-mapping.sh` — DNA only

**Per-sample steps:**
1. Index BAM with `pbindex` (skipped if `.pbi` already exists)
2. Convert to FASTQ with `bam2fastq`
3. Align to reference with `minimap2 -ax map-hifi`
4. Filter with `samtools view -F 2308` (removes unmapped, secondary, and supplementary reads)
5. Sort and index aligned BAM
6. Export aligned reads as FASTA for downstream analysis

> **Flag note:** `-F 2308` = `-F 4` (unmapped) + `-F 256` (secondary) + `-F 2048` (supplementary)

**Directory structure expected:**
```
split-barcode/
├── RNA/trimmed/
│   ├── HIV/   ← *.trimmed.bam files go here
│   └── SIV/
└── DNA/trimmed/
    ├── HIV/
    └── SIV/
```

---

### Step 4 — Mutation Counting

#### 4a. Per-read mutation counting (`mutation-counting.sh` + `mutation_counting2.py`)

A SLURM array job that runs `mutation_counting2.py` on each aligned BAM independently.

**Features:**
- Counts all 12 single-nucleotide substitution types (e.g., `G_A`, `C_T`)
- Counts GA→AA dinucleotide context mutations (APOBEC3-relevant `G→A` in `GА` context)
- Output: one CSV per BAM with per-read mutation counts

**Output columns (fixed schema):**

*Single-base substitutions:* `A_C`, `A_G`, `A_T`, `C_A`, `C_G`, `C_T`, `G_A`, `G_C`, `G_T`, `T_A`, `T_C`, `T_G`

*GA-context G→A substitutions:* `GA_AA`, `GC_AA`, `GG_AA`, `GT_AA` (and AC/AG/AT variants)

**SLURM setup:**
```bash
#SBATCH --array=1-20%10   # adjust upper bound to total BAM count
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
```

**Required environment variables (set by the shell script):**
```bash
export REF_PATH=/your/path/to/reference.fasta
export BAM_PATH=/your/path/to/sample.aligned.bam
export RESULT_DIR=/your/path/to/output/
```

#### 4b. Hypermut 3.0 analysis (`BAM_HYPMUT3.py`)

A BAM-native re-implementation of the HYPERMUT algorithm. For each read, counts mutations in a user-defined nucleotide context (e.g., `G→A` in `GA` context) and performs Fisher's exact test comparing primary vs. control sites.

**Usage:**
```bash
python BAM_HYPMUT3.py \
  sample.aligned.bam \
  reference.fasta \
  G A \                        # mutationfrom mutationto
  --upstreamcontext "" \
  --downstreamcontext "A|D" \  # GA context (APOBEC3G)
  --enforce D \
  --prefix output_prefix_
```

**Outputs:**
- `*_args.csv` — run parameters
- `*_summary.csv` — per-read Fisher's exact test results
- `*_positions.csv` — per-position site annotations

---

### Step 5 — Combine Hypermut3 Results (`combine_hypermut3_final.sh`)

Aggregates all `*_hypermut3_*_summary.csv` files into a single master table.

```bash
# Finds all summary CSVs and passes them to the combiner
find /your/project/path -type f -name "*_hypermut3_*_summary.csv" > hypermut3_csv_list.txt

python3 combine_hypermut3_final.py \
  --list hypermut3_csv_list.txt \
  --out  master_hypermut3_final.csv
```

---

## Ribo-seq Analysis (`ribo-seq-analysis.R`)

Ribosome profiling pipeline for APOBEC3G translational analysis in rhesus macaque (*Macaca mulatta*) and human (*Homo sapiens*).

### Workflow

1. **RiboseQC** — Quality control and P-site calling from BAM files
   - Cleans GTF (consistent strand/chromosome per gene)
   - Prepares genome annotation files
   - Extracts TPM/RPKM values per gene

2. **P-site reformatting** — Converts RiboseQC `.bedgraph` outputs (plus/minus strands) into the combined `.ggRibo` format

3. **ggRibo visualization** — Plots ribosome occupancy alongside RNA-seq for specific isoforms

**Key genes visualized:**
- `APOBEC3G` — full transcript, Exon 1, Exon 1b, Exon 2

**Datasets used:**
- Rhesus macaque riboseq: PRJNA292112
- Rhesus macaque RNAseq: PRJNA214703
- Human riboseq: PRJNA292112
- Human RNAseq: PRJNA122271

**Reference genomes:**
- *Macaca mulatta* Mmul_10 (Ensembl release 114)
- *Homo sapiens* GRCh38 (Ensembl release 114)

---

## Configuration

All scripts use hardcoded paths that must be updated before running. Search for `/your/path/here` placeholders and replace with your actual directory structure. The general layout assumed is:

```
/your/project/root/
├── reference-genomes/
├── split-barcode/
│   ├── RNA/trimmed/{HIV,SIV}/
│   └── DNA/trimmed/{HIV,SIV}/
├── mutational-calling/
│   ├── HIV/res/
│   └── SIV/res/
└── hypermut3-combined/
```

---

## Grid Engine / SLURM Notes

- Demultiplexing and trimming scripts are written for **SGE** (`#$ -pe smp N`)
- Mutation counting uses **SLURM** (`#SBATCH --array`)
- Update `--array=1-N` in `mutation-counting.sh` to match your actual BAM count

---

## License

MIT License — see [LICENSE](LICENSE) for details.

## Author

Armando Mendez — Texas Biomedical Research Institute
