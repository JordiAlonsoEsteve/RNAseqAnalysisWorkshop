from Bio import SeqIO

# Parse the FASTA file and convert it into a dictionary
# Keys will be the sequence IDs, and values will be SeqRecord objects
fasta_dict = SeqIO.to_dict(SeqIO.parse("genomic_data/col-0_keygene.fasta", "fasta"))

# Example: Print IDs and their sequences
for seq_id, record in fasta_dict.items():
  print(f"ID: {seq_id}")
  print(f"Length: {len(record.seq)}")
  #print(f"Sequence: {record.seq}\n")