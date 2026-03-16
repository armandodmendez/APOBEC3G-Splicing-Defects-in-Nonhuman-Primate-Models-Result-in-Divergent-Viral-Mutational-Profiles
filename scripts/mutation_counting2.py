from collections import defaultdict
import re
import os
import pandas as pd
import pysam

def reconstruct_pirwise_alignment(query, ref, aligned_pairs):
    reconstructed_query = []
    reconstructed_ref = []
    for Q_pos, R_pos in aligned_pairs:
        if Q_pos is None:
            reconstructed_ref.append(ref[R_pos])
            reconstructed_query.append("-")
        elif R_pos is None:
            reconstructed_ref.append("-")
            reconstructed_query.append(query[Q_pos])
        else:
            reconstructed_ref.append(ref[R_pos])
            reconstructed_query.append(query[Q_pos])
    return "".join(reconstructed_query), "".join(reconstructed_ref)

def read_fasta(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: The file '{file_path}' does not exist.")
    headers = []
    sequences = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        if not lines:
            raise ValueError(f"Error: The file '{file_path}' is empty.")
        sequence_data = ''
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                if sequence_data:
                    sequence_data = re.sub(r"[^TCGAU]", "-", sequence_data.upper()).replace("U", "T")
                    sequences.append(sequence_data)
                sequence_data = ''
                headers.append(line[1:])
            elif line:
                sequence_data += line
        if sequence_data:
            sequence_data = re.sub(r"[^TCGAU]", "-", sequence_data.upper()).replace("U", "T")
            sequences.append(sequence_data)
    if not headers:
        raise ValueError(f"Error: The file '{file_path}' does not contain any headers.")
    if len(headers) != len(sequences):
        raise ValueError(f"Error: Mismatch between headers ({len(headers)}) and sequences ({len(sequences)}).")
    return headers, sequences

def find_mutations(query, reference):
    return [
        i for i in range(len(query))
        if query[i] != reference[i]
        and 1 <= i < len(query) - 1
        and query[i] != '-'
        and reference[i] != '-'
    ]

# fixed column schema (exact order you want in every CSV)
MUTATION_COLUMNS = [
    # 1-bp motifs
    "A_C", "A_G", "A_T",
    "C_A", "C_G", "C_T",
    "G_A", "G_C", "G_T",
    "T_A", "T_C", "T_G",

    # GA contexts
    "GA_AA", "GA_AC", "GA_AG", "GA_AT",
    # GC contexts
    "GC_AA", "GC_AC", "GC_AG", "GC_AT",
    # GG contexts
    "GG_AA", "GG_AC", "GG_AG", "GG_AT",
    # GT contexts
    "GT_AA", "GT_AC", "GT_AG", "GT_AT",
]

# environment variables from SLURM
ref_path   = os.environ["REF_PATH"]
bam_path   = os.environ["BAM_PATH"]
result_dir = os.environ["RESULT_DIR"]

print(f"Loading ref: {ref_path}")
print(f"Loading BAM: {bam_path}")
print(f"Saving to: {result_dir}")

# handle single/multi-seq FASTA
_, ref_list = read_fasta(ref_path)
if len(ref_list) > 1:
    print(f"Warning: Multiple refs found ({len(ref_list)}), using first.")
ref = ref_list[0]

motif_start = 0
motif_end = 1

try:
    name = os.path.splitext(os.path.basename(bam_path))[0]
    print(f"Processing {name}...")
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    mutation_counts = defaultdict(lambda: defaultdict(int))

    j = 0
    for read in bamfile.fetch(until_eof=True):
        if read.is_duplicate:
            continue

        header = read.qname
        query_seq, ref_seq = reconstruct_pirwise_alignment(
            read.query_sequence, ref, read.aligned_pairs
        )
        mutation_list = find_mutations(query_seq, ref_seq)

        for mutation_point in mutation_list:
            motif_ref = ref_seq[(mutation_point + motif_start):(mutation_point + motif_end)]
            motif_query = query_seq[(mutation_point + motif_start):(mutation_point + motif_end)]

            if "-" in motif_ref or "-" in motif_query:
                continue

            from_seq = motif_ref
            to_seq = motif_query
            key = from_seq + "_" + to_seq
            mutation_counts[header][key] += 1

            # GA context (2bp) special case
            if from_seq == "G" and to_seq == "A":
                context_ref = ref_seq[mutation_point:(mutation_point + 2)]
                context_query = query_seq[mutation_point:(mutation_point + 2)]
                context_key = context_ref + "_" + context_query
                mutation_counts[header][context_key] += 1

        j += 1
        if j % 10000 == 0:
            print(f"Processed {j} reads")

    bamfile.close()
    print(f"Total reads processed: {j}")

    # build DataFrame; ensure all rows have all columns, in fixed order
    df = pd.DataFrame.from_dict(mutation_counts, orient="index").fillna(0)

    # add any missing columns from the schema as 0
    for col in MUTATION_COLUMNS:
        if col not in df.columns:
            df[col] = 0

    # keep only schema columns, in this exact order
    df = df[MUTATION_COLUMNS]

    # optional: sort rows by read name
    df = df.sort_index()

    # write output
    output_file = os.path.join(result_dir, name + "_mutation_counts.csv")
    df.to_csv(output_file)
    print(f"Result saved: {output_file}")

except Exception as e:
    print(f"Error: {str(e)}")
    raise
