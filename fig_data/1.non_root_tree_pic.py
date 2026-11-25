import os
import subprocess
import numpy as np
import pandas as pd  # Used to save the distance matrix to a CSV file
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio.SeqIO import parse
from tqdm import tqdm  # Used to display a progress bar

def calculate_similarity(alignment):
    """Calculate the similarity between two sequences"""
    seq1 = alignment[0].seq
    seq2 = alignment[1].seq
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return matches / len(seq1)


# Set the query folder path
query_dir = "fig_data/fasta_results/"  # Directory containing the query sequences

# Path to the EMBOSS needle executable
needle_path = "needle"  # Make sure `needle` is installed and callable

# Iterate over every FASTA file in the query directory
for fasta_file in os.listdir(query_dir):
    if fasta_file.endswith(".fasta"):
        fasta_path = os.path.join(query_dir, fasta_file)

        # Read the FASTA file and extract the sequences
        with open(fasta_path, "r") as handle:
            sequences = list(parse(handle, "fasta"))

        # Extract the IDs of all sequences
        sequence_ids = [seq.id for seq in sequences]

        # Initialize the distance matrix filled with zeros
        num_sequences = len(sequences)
        matrix = np.zeros((num_sequences, num_sequences))  # Use numpy to create an all-zero matrix

        # Align all pairs and compute the distances
        total_comparisons = (num_sequences * (num_sequences - 1)) // 2  # Compute the total number of pairs to compare

        count = 0

        # Use tqdm to display the progress bar
        with tqdm(total=total_comparisons, desc=f"Comparing sequences in {fasta_file}") as pbar:
            for i in range(1, num_sequences):
                for j in range(i):
                    seq1 = sequences[i]
                    seq2 = sequences[j]

                    # Save temporary files for the pairwise alignment
                    try:
                        with open("temp1.fasta", "w") as f1, open("temp2.fasta", "w") as f2:
                            f1.write(f">{seq1.id}\n{seq1.seq}\n")
                            f2.write(f">{seq2.id}\n{seq2.seq}\n")
                            f1.flush()  # Force the file to be written
                            f2.flush()  # Force the file to be written
                            f1.close()
                            f2.close()
                    except IOError as e:
                        print(f"Error writing to file: {e}")
                        continue  # Skip the current alignment

                    # Call needle for the alignment, adding gapopen and gapextend parameters
                    needle_cmd = [
                        needle_path,
                        "-asequence=temp1.fasta",
                        "-bsequence=temp2.fasta",
                        "-gapopen", "10",  # Set the gap opening penalty
                        "-gapextend", "0.5",  # Set the gap extension penalty
                        "-outfile=temp_aligned.needle"
                    ]

                    try:
                        result = subprocess.run(needle_cmd, capture_output=True, text=True, check=True, timeout=60)  # Add a timeout
                        if result.returncode != 0:
                            print(f"Error in needle command: {result.stderr}")
                            continue  # Skip this alignment and continue with the next pair
                    except subprocess.TimeoutExpired:
                        print("needle command timed out, skipping this pair")
                        continue
                    except subprocess.CalledProcessError as e:
                        print(f"Error running needle: {e.stderr}")
                        continue  # Skip this alignment and continue with the next pair
                    finally:
                        # Delete the temporary files
                        if os.path.exists('temp1.fasta'):
                            os.remove('temp1.fasta')
                        if os.path.exists('temp2.fasta'):
                            os.remove('temp2.fasta')

                    # Parse the alignment result and compute the distance (sequence similarity is used here)
                    alignment = AlignIO.read("temp_aligned.needle", "emboss")
                    score = calculate_similarity(alignment)

                    # Compute the distance and fill the lower triangle of the matrix
                    distance = 1 - score
                    matrix[i, j] = distance  # Fill the lower triangle of the matrix

                    # Update the progress bar after each alignment
                    pbar.update(1)

        # Create the distance matrix object, ensuring it is a lower triangular matrix
        matrix_list = matrix.tolist()
        dm = DistanceMatrix(sequence_ids, [row[:i+1] for i, row in enumerate(matrix_list)])

        # Build the phylogenetic tree
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

        # Save the phylogenetic tree in Newick format
        output_tree = f"fig_data/cal_result/HP-select/{fasta_file}_evolutionary_tree.nwk"
        Phylo.write(tree, output_tree, format="newick")

        print(f"Evolutionary tree saved to: {output_tree}")

        # Save the distance matrix as a CSV file
        distance_matrix_df = pd.DataFrame(matrix, index=sequence_ids, columns=sequence_ids)
        distance_matrix_output_file = f"fig_data/cal_result/HP-select/{fasta_file}_distance_matrix.csv"
        distance_matrix_df.to_csv(distance_matrix_output_file)

        print(f"Distance matrix saved to: {distance_matrix_output_file}")

print("All processes completed.")


