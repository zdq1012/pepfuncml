import os
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import math

# Set the file paths
folder_path = "analysis/predicted/seqs"  # Replace with your actual path
output_folder = "analysis/predicted/seqs_properties"  # Replace with the output path
os.makedirs(output_folder, exist_ok=True)  # Make sure the output folder exists

# Eisenberg hydrophobicity reference scale
eisenberg_scale = {
    'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90, 'C': 0.29,
    'Q': -0.85, 'E': -0.74, 'G': 0.48, 'H': -0.40, 'I': 1.38,
    'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19, 'P': 0.12,
    'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08
}

def calculate_eisenberg_hydrophobicity(seq):
    """Calculate the Eisenberg hydrophobicity"""
    return sum(eisenberg_scale.get(aa, 0) for aa in seq) / len(seq)

def calculate_hydrophobic_moment(seq):
    """Calculate the hydrophobic moment"""
    moment = 0.0
    for i, aa in enumerate(seq):
        hydrophobicity = eisenberg_scale.get(aa, 0)
        moment += hydrophobicity * math.sin((i + 1) * math.pi / len(seq))
    return moment / len(seq)

def calculate_aliphatic_index(seq):
    """Calculate the aliphatic index"""
    count_a = seq.count('A')
    count_v = seq.count('V')
    count_i = seq.count('I')
    count_l = seq.count('L')
    return (count_a + count_v * 2.9 + (count_i + count_l) * 3.9) / len(seq)

def calculate_properties(seq):
    """Calculate all physicochemical properties of a peptide sequence"""
    analysis = ProteinAnalysis(seq)

    # Compute the properties in the order of Fig. 4
    gravy = analysis.gravy()
    hydrophobic_moment = calculate_hydrophobic_moment(seq)
    aliphatic_index = calculate_aliphatic_index(seq)
    aromaticity = analysis.aromaticity()
    instability_index = analysis.instability_index()
    isoelectric_point = analysis.isoelectric_point()
    charge_at_ph7 = analysis.charge_at_pH(7.0)
    charge_density = charge_at_ph7 / len(seq)
    hydrophobic_ratio = sum(1 for aa in seq if aa in "AILMFWV") / len(seq)
    eisenberg_hydrophobicity = calculate_eisenberg_hydrophobicity(seq)

    return {
        "Sequence": seq,
        "GRAVY": gravy,
        "Hydrophobic_Moment": hydrophobic_moment,
        "Aliphatic_Index": aliphatic_index,
        "Aromaticity": aromaticity,
        "Instability_Index": instability_index,
        "Eisenberg_Hydrophobicity": eisenberg_hydrophobicity,
        "Isoelectric_Point": isoelectric_point,
        "Hydrophobic_Ratio": hydrophobic_ratio,
        "Charge_Density": charge_density,
        "Charge": charge_at_ph7,
        
        
    }

# Iterate over the CSV files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".csv"):
        file_path = os.path.join(folder_path, filename)
        df = pd.read_csv(file_path)

        # Store the computed results
        results = []

        # Compute the properties of each sequence
        for seq in df['seq'].dropna():
            properties = calculate_properties(seq)
            results.append(properties)

        # Save the results as a CSV file
        result_df = pd.DataFrame(results)
        output_filename = f"{os.path.splitext(filename)[0]}.csv"
        output_path = os.path.join(output_folder, output_filename)
        result_df.to_csv(output_path, index=False)

        print(f"Results saved to: {output_path}")
