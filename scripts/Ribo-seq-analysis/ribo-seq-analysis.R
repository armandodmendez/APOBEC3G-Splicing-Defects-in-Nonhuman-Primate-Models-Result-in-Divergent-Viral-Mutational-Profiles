###############################################################################
# Ribo-seQC on PRJNA292112 (riboseq data) for rhesus macaque

library(devtools)
library(RiboseQC)
library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)

# Set directory and file paths
dir <- "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ribo-seqc/rhesus"

mmul_gtf <- file.path(dir, "Macaca_mulatta.Mmul_10.114.filtered2.gtf")
mmul_fasta <- file.path(dir, "Macaca_mulatta.Mmul_10.dna.toplevel.fa")
bam_filepath <- file.path(dir, "PRJNA292112-merged-mmul_10.bam")

# Clean GTF file
gtf_gr <- import(mmul_gtf)
seqlevels(gtf_gr) <- seqlevels(gtf_gr)[!is.na(seqlevels(gtf_gr))]
seqlevels(gtf_gr) <- unique(seqlevels(gtf_gr))
export(gtf_gr, file.path(dir, "Macaca_mulatta.Mmul_10.114_cleaned.gtf"))


# Check BAM seqlevels
bam_header <- scanBamHeader(bam_filepath)[[1]]
bam_seqlevels <- names(bam_header$targets)


# Filter genes with consistent strand AND chromosome
gene_ids <- mcols(gtf_gr)$gene_id
chrom_uniform <- tapply(as.character(seqnames(gtf_gr)), gene_ids, 
                        function(x) length(unique(x)) == 1)
strand_uniform <- tapply(as.character(strand(gtf_gr)), gene_ids, 
                         function(x) length(unique(x)) == 1)
valid_genes <- names(which(chrom_uniform & strand_uniform))

# Apply filtering
gtf_gr <- gtf_gr[gene_ids %in% valid_genes]

# Export with compression (recommended)
export(gtf_gr, file.path(dir, "Macaca_mulatta.Mmul_10.114_cleaned.gtf"), format = "GTF")


# Rebuild annotation with consistent seqlevels
annot_file <- prepare_annotation_files(
  annotation_directory = dir,
  genome_seq = mmul_fasta,
  gtf_file = file.path(dir, "Macaca_mulatta.Mmul_10.114_cleaned.gtf"),
  scientific_name = "Macaca.mulatta",
  annotation_name = "Mmul_10",
  export_bed_tables_TxDb = TRUE,
  forge_BSgenome = FALSE,
  create_TxDb = TRUE
)

# Run analysis with corrected files
load_annotation(annot_file)
resfile <- RiboseQC_analysis(
  annotation_file = annot_file,
  bam_files = bam_filepath,
  fast_mode = FALSE,
  report_file = file.path(dir, "mmul_PRJNA292112.html"),
  write_tmp_files = TRUE,
  create_report = TRUE,
  extended_report = TRUE,
  readlength_choice_method = "all"
)

##extract RPKM and TPM values

load("/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ribo-seqc/rhesus/PRJNA292112-merged-mmul_10.bam_results_RiboseQC_all")

df <- data.frame(value = res_all$read_stats$counts_all_genes$TPM)
write.csv(df, "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ribo-seqc/rhesus/rhesus-PRJNA292112-TPM.csv", row.names = TRUE)


## Reformatting Ribo-seQC P-site output
## The two P-site files can be converted into the format compatible with ggRibo using the following code in R:
# Set directory and file paths
dir <- "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ribo-seqc/rhesus"

### Define file paths and load the plus and minus strand BedGraph files
in_plus <- file.path(dir, "PRJNA292112-merged-mmul_10.bam_P_sites_plus.bedgraph")
in_minus <- file.path(dir, "PRJNA292112-merged-mmul_10.bam_P_sites_minus.bedgraph")
out_file <- file.path(dir,"PRJNA292112-merged-mmul_10.bam_riboseq_Psites_combined.ggRibo")

plus <- read.table(in_plus, sep="\t")
minus <- read.table(in_minus, sep="\t")

### Reformat dataframe and output table
plus[,ncol(plus)+1] <- "+"  # assign strands
minus[,ncol(minus)+1] <- "-"

comb <- rbind(plus, minus)  # combine dataframes
comb <- comb[,c(4, 1, 3, 5)]  # reorder columns to match ggRibo format

write.table(comb, file=out_file, col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep="\t")

###################################################################################

#Unload the ggRibo first if you have it loaded
detach("package:ggRibo", unload=TRUE)
#Remove the old version of ggRibo if you installed it before
remove.packages("ggRibo")
#Follow the stage below to install ggRibo.
#Install ggRibo
library(devtools)
install_github("hsinyenwu/ggRibo@v2025.9.24")





# ggRibo
#load packages
library(devtools)
library(ggRibo)

ggRibo_dir <- "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ggRibo/rhesus"

#path for annotated gtf
mm_gtf <- file.path(ggRibo_dir, "Macaca_mulatta.Mmul_10.114_cleaned.gtf")

#load transcriptome annotation
gtf_import(annotation = mm_gtf, format ="gtf",organism = "Macaca mulatta")

#load genome FASTA
mm_fasta <- FaFile("/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ggRibo/rhesus/Macaca_mulatta.Mmul_10.dna.toplevel.fa")

# set up variables for the ggRibo function
# input single-end data
inputs_full <- create_seq_input(
  rna_files = "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ggRibo/rhesus/PRJNA214703-merged-mmul_10.bam",
  ribo_files = "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ggRibo/rhesus/PRJNA292112-merged-mmul_10.bam_riboseq_Psites_combined.ggRibo",
  sample_names = "Rhesus LCL",
  rna_paired = "single"
)


### Plot subset of isoforms
tx_subset <- c("ENSMMUT00000101598", "PB.1708.4", "A3G.M1")

ggRibo(tx_id = "PB.1708.4", NAME = "APOBEC3G", FASTA = mm_fasta, show_seq = F, Extend = 50, selected_isoforms = tx_subset, gene_model_coord_font_size =16, 
       axis_label_font_size = 18, axis_title_font_size = 20, transcript_label_font_size = 16, sample_label_font_size = 0, title_font_size = 0)


# Plot with ggRibo to zoom-in into Exon 1b
range <- c(12499165,12498975)
ggRibo(tx_id = "A3G.M1", NAME = "APOBEC3G Exon 1b", FASTA = mm_fasta, show_seq = T, Extend = 25,
       dna_aa_height_ratio = 0.20, plot_range = range,  selected_isoforms = tx_subset, gene_model_coord_font_size =16, 
       axis_label_font_size = 18, axis_title_font_size = 20, transcript_label_font_size = 16, sample_label_font_size = 0, title_font_size = 0)

#plot with ggRibo to zoom-in into Exon1
range1 <- c(12499530,12499360)
ggRibo(tx_id = "A3G.M1", NAME = "APOBEC3G Exon 1", FASTA = mm_fasta, show_seq = T, Extend = 25,
       dna_aa_height_ratio = 0.20, plot_range = range1,  selected_isoforms = tx_subset, gene_model_coord_font_size =16, 
       axis_label_font_size = 18, axis_title_font_size = 20, transcript_label_font_size = 16, sample_label_font_size = 0, title_font_size = 0)

#plot with ggRibo to zoom-in into Exon2
range2 <- c(12497525,12497360)
ggRibo(tx_id = "A3G.M1", NAME = "APOBEC3G Exon 2", FASTA = mm_fasta, show_seq = T, Extend = 25,
       dna_aa_height_ratio = 0.20, plot_range = range2,  selected_isoforms = tx_subset, gene_model_coord_font_size =16, 
       axis_label_font_size = 18, axis_title_font_size = 20, transcript_label_font_size = 16, sample_label_font_size = 0, title_font_size = 0)









#################################################################################
#Homo sapiens Ribo-seq analysis
# Ribo-seQC on PRJNA292112 (riboseq data) for human hg38

library(devtools)
library(RiboseQC)
library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)

# Set directory and file paths
dir <- "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ribo-seqc/human"

hg38_gtf <- file.path(dir, "Homo_sapiens.GRCh38.114.gtf")
hg38_fasta <- file.path(dir, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
bam_human_path <- file.path(dir, "PRJNA122271-merged-hg38.bam")

# Clean GTF file
gtf_gr <- import(hg38_gtf)
seqlevels(gtf_gr) <- seqlevels(gtf_gr)[!is.na(seqlevels(gtf_gr))]
seqlevels(gtf_gr) <- unique(seqlevels(gtf_gr))
export(gtf_gr, file.path(dir, "hg38_cleaned.gtf"))


# Check BAM seqlevels
bam_header <- scanBamHeader(bam_human_path)[[1]]
bam_seqlevels <- names(bam_header$targets)


# Filter genes with consistent strand AND chromosome
gene_ids <- mcols(gtf_gr)$gene_id
chrom_uniform <- tapply(as.character(seqnames(gtf_gr)), gene_ids, 
                        function(x) length(unique(x)) == 1)
strand_uniform <- tapply(as.character(strand(gtf_gr)), gene_ids, 
                         function(x) length(unique(x)) == 1)
valid_genes <- names(which(chrom_uniform & strand_uniform))

# Apply filtering
gtf_gr <- gtf_gr[gene_ids %in% valid_genes]

# Export with compression (recommended)
export(gtf_gr, file.path(dir, "hg38_cleaned.gtf"), format = "GTF")


# Rebuild annotation with consistent seqlevels
annot_file_human <- prepare_annotation_files(
  annotation_directory = dir,
  genome_seq = hg38_fasta,
  gtf_file = file.path(dir, "hg38_cleaned.gtf"),
  scientific_name = "Homo.sapiens",
  annotation_name = "hg38",
  export_bed_tables_TxDb = TRUE,
  forge_BSgenome = FALSE,
  create_TxDb = TRUE
)

# Run analysis with corrected files
load_annotation(annot_file_human)
resfile <- RiboseQC_analysis(
  annotation_file = annot_file_human,
  bam_files = bam_human_path,
  fast_mode = FALSE,
  report_file = file.path(dir, "human_PRJNA292112.html"),
  write_tmp_files = TRUE,
  create_report = TRUE,
  extended_report = TRUE,
  readlength_choice_method = "all"
)


###extract RPKM and TPM

load("/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ribo-seqc/human/PRJNA292112-merged-hg38.bam_results_RiboseQC_all")

df <- data.frame(value = res_all$read_stats$counts_all_genes$TPM)
write.csv(df, "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ribo-seqc/human/human-PRJNA292112-TPM.csv", row.names = TRUE)
df <- data.frame(value = res_all$read_stats$counts_all_genes$RPKM)
write.csv(df, "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ribo-seqc/human/human-PRJNA292112-RPKM.csv", row.names = TRUE)




## Reformatting Ribo-seQC P-site output
## The two P-site files can be converted into the format compatible with ggRibo using the following code in R:
# Set directory and file paths
dir <- "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ribo-seqc/human"

### Define file paths and load the plus and minus strand BedGraph files
in_plus <- file.path(dir, "PRJNA292112-merged-hg38.bam_P_sites_plus.bedgraph")
in_minus <- file.path(dir, "PRJNA292112-merged-hg38.bam_P_sites_minus.bedgraph")
out_file <- file.path(dir, "PRJNA292112-merged-hg38.bam_riboseq_Psites_combined.ggRibo")

plus <- read.table(in_plus, sep="\t")
minus <- read.table(in_minus, sep="\t")

### Reformat dataframe and output table
plus[,ncol(plus)+1] <- "+"  # assign strands
minus[,ncol(minus)+1] <- "-"

comb <- rbind(plus, minus)  # combine dataframes
comb <- comb[,c(4, 1, 3, 5)]  # reorder columns to match ggRibo format

write.table(comb, file=out_file, col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep="\t")

###################################################################################
# ggRibo
#load packages
library(devtools)
library(ggRibo)

ggRibo_dir <- "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ggRibo/human"

#path for annotated gtf
hu_gtf <- file.path(ggRibo_dir, "hg38_cleaned.gtf")

#load transcriptome annotation
gtf_import(annotation = hu_gtf, format ="gtf",organism = "Homo sapiens")

#load genome FASTA
hg38_fasta <- FaFile("/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ggRibo/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa")


# set up variables for the ggRibo function
# input single-end data
inputs_full <- create_seq_input(
  rna_files = "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ggRibo/human/PRJNA122271-merged-hg38.bam",
  ribo_files = "/Users/armandomendez/TX Biomed Dropbox/Armando Mendez/Armando-Mendez/R/ggRibo/human/PRJNA292112-merged-hg38.bam_riboseq_Psites_combined.ggRibo",
  sample_names = "Human LCL",
  rna_paired = "single"
)

tx_subset_hg <- c("ENST00000407997")

# Plot with ggRibo
ggRibo(gene_id = "ENSG00000239713", tx_id = "ENST00000407997", NAME = "APOBEC3G", FASTA = hg38_fasta, show_seq = F, Extend = 25, 
       selected_isoforms = tx_subset_hg, gene_model_coord_font_size =16, 
       axis_label_font_size = 18, axis_title_font_size = 20, transcript_label_font_size = 16, sample_label_font_size = 0, title_font_size = 0)

# Plot with ggRibo to zoom-in into Exon 2
range1 <- c(39078920,39079090)
ggRibo(tx_id = "ENST00000407997", NAME = "APOBEC3G Exon 2", FASTA = hg38_fasta, show_seq = F, Extend = 25,
       dna_aa_height_ratio = 0.15, plot_range = range1, selected_isoforms = tx_subset_hg, gene_model_coord_font_size =16, 
       axis_label_font_size = 18, axis_title_font_size = 20, transcript_label_font_size = 16, sample_label_font_size = 0, title_font_size = 0)

# Plot with ggRibo between Exon1 and 2
range2 <- c(39077000,39079090)
ggRibo(tx_id = "ENST00000407997", NAME = "APOBEC3G Exon 1-2", FASTA = hg38_fasta, show_seq = F, Extend = 25,
       dna_aa_height_ratio = 0.15, plot_range = range2, selected_isoforms = tx_subset_hg, gene_model_coord_font_size =16, 
       axis_label_font_size = 18, axis_title_font_size = 20, transcript_label_font_size = 16, sample_label_font_size = 0, title_font_size = 0)

#Plot with ggRibo Exon 1
range3 <- c(39077080,39077385)
ggRibo(tx_id = "ENST00000407997", NAME = "APOBEC3G Exon 1", FASTA = hg38_fasta, show_seq = F, Extend = 25,
       dna_aa_height_ratio = 0.15, plot_range = range3, selected_isoforms = tx_subset_hg, gene_model_coord_font_size =16, 
       axis_label_font_size = 18, axis_title_font_size = 20, transcript_label_font_size = 16, sample_label_font_size = 0, title_font_size = 0)







