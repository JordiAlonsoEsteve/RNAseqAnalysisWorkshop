import re

input_file = "data/genomic_data/TAIR12_raw/col-0_liftoff_TAIR12_polished_tbtools_cleaned_agat.gff"
output_file = "data/genomic_data/TAIR_12/annotation.gff"

# Define the regex pattern for the prefix you want to remove
pattern = re.compile(r"gene-TAIR12_TAIR12_")

with open(input_file, "r") as fin, open(output_file, "w") as fout:
  for line in fin:
    # Keep header/comment lines intact
    if line.startswith("#"):
      fout.write(line)
      continue

    parts = line.strip().split("\t")
    if len(parts) < 9:
      fout.write(line)
      continue

    # Column 9 contains the semicolon-separated attributes (ID=, Parent=, etc.)
    # Clean up the prefix inside the attributes column
    parts[8] = pattern.sub("", parts[8])

    # Write the modified line back out
    fout.write("\t".join(parts) + "\n")

print(f"Cleaned GFF written to {output_file}")