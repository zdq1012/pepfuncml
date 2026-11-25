import pandas as pd
from Bio.Align import PairwiseAligner
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

# Initialize the PairwiseAligner
aligner = PairwiseAligner()
aligner.mode = 'global'  # Use global alignment mode (similar to the Needleman-Wunsch algorithm)

# Define the name list
names = [ 'Anti-inflammatory']

# Define the colors (based on the uploaded image)
colors = [ '#ADD8E6']  # Use a HEX color code

# Used to store the match percentage data of all names
all_match_percentages = []

# Iterate over each name, align the sequences and save the results
for name in names:
    print(f"Processing: {name}")

    # Read the two CSV files
    csv1 = pd.read_csv(f'analysis/predicted/seqs/{name}.csv')
    csv2 = pd.read_csv(f'analysis/trained/seqs/{name}.csv')

    # Assume the 'seq' column stores the peptide sequences
    sequences1 = csv1['seq'].tolist()
    sequences2 = csv2['seq'].tolist()

    # Store the alignment results
    for seq1 in sequences1:
        for seq2 in sequences2:
            # Compute the alignment score
            score = aligner.score(seq1, seq2)
            match_percentage = score / max(len(seq1), len(seq2))  # Compute the match ratio
            all_match_percentages.append({'Name': name, 'Match_Percentage': match_percentage})

# Aggregate all the match percentage data into a DataFrame
df_all = pd.DataFrame(all_match_percentages)

# Create the chart
fig, ax = plt.subplots(figsize=(10, 6))

# Generate the x-axis for plotting (from 0 to 1, step of 0.01)
x = np.linspace(0, 1, 10)  # Use more points to make the curve smoother

# Draw a curve for each name using the corresponding color
for name, color in zip(names, colors):
    # Get the match percentage data of the current name
    subset = df_all[df_all['Name'] == name]['Match_Percentage']

    # Use KDE to compute a smooth curve
    kde = gaussian_kde(subset, bw_method=0.0000001)  # Set the bandwidth to control the smoothness

    # Compute the y values at the x points and normalize the density to a proportion
    y = kde(x)
    y_normalized = y / y.sum()  # Normalize the density curve to a proportion

    # Draw a smooth density plot and fill the area under the curve
    plt.fill_between(x, y_normalized, alpha=0.5, label=name, color=color)

# Set the font sizes and labels
plt.xticks(fontsize=20, ha='center')  # Center-align the x-axis labels
plt.xlabel('Match Percentage', fontsize=24)
plt.yticks(fontsize=20)  # Set the y-axis ticks to 0, 0.1, 0.2, ..., 1.0
plt.ylabel('Proportion', fontsize=24)

# Set the spine styles
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')

# Set the tick mark width to 2
ax.tick_params(axis='both', width=2, length=12)

ax.legend(prop={'size': 24})  # Set the legend font size to 24

# Adjust the layout to avoid overlapping labels
plt.tight_layout()

# Save and display the chart
plt.savefig('analysis/pic/match_percentage_density_colored.png')
plt.savefig("analysis/pic/match_percentage_density_colored.svg", format='svg')
plt.show()
