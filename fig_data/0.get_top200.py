import pandas as pd

# Read the CSV file
df = pd.read_csv("HP-select.csv")  # Replace with your file path

# Extract the sequences of the first 200 rows
sequences = df['Seq'][:200]  # Get the first 200 sequences

# Create and write the FASTA file
with open("HP-select.fasta", "w") as f:
    for idx, seq in enumerate(sequences, start=1):
        f.write(f">sequence_{idx}\n")  # Sequence title, modify as needed
        f.write(f"{seq}\n")  # Write the sequence

print("FASTA file has been saved as 'sequences.fasta'.")
