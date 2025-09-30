from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
from tqdm import tqdm

# Input FASTA file paths
bacts_file = "/cs/usr/stavperez/sp/InSilicoSeq/dataset1/output_reads_mix_256_orgs/assembly_output/128_euks/contigs.fasta"
classified_file = "/cs/usr/stavperez/sp/InSilicoSeq/dataset1/output_reads_mix_256_orgs/assembly_output/classified_euks_256_orgs/contigs.fasta"

# Load contigs from FASTA
print("Loading contigs...")
bacts = list(SeqIO.parse(bacts_file, "fasta"))[:50]
classified = list(SeqIO.parse(classified_file, "fasta"))[:50]

aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = -1
aligner.extend_gap_score = -0.5

results = []

for b in tqdm(bacts, desc="Aligning"):
    for c in classified:
        score = aligner.score(b.seq, c.seq)
        bact_len = len(b.seq)
        class_len = len(c.seq)
        length_diff = abs(bact_len - class_len)
        percent_identity = (score / max(bact_len, class_len)) * 100

        results.append({
            "bact_id": b.id,
            "classified_id": c.id,
            "bact_len": bact_len,
            "classified_len": class_len,
            "alignment_score": round(score, 2),
            "percent_identity": round(percent_identity, 2),
            "length_diff": length_diff
        })

# Convert to DataFrame and export
df = pd.DataFrame(results)
df.to_csv("/cs/usr/stavperez/sp/InSilicoSeq/dataset1/output_reads_mix_256_orgs/assembly_output/contig_euks_pairwise_alignment_results.csv", index=False)
print("Done! Results saved to /cs/usr/stavperez/sp/InSilicoSeq/dataset1/output_reads_mix_256_orgs/assembly_output/contig_euks_pairwise_alignment")
