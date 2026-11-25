import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Define a function to read the data and label the source
def load_data(name, source):
    ori_source = source
    # Map 'trained' to 'Training' and 'predicted' to 'Predicted'
    source_map = {'trained': 'Training', 'predicted': 'Predicted'}
    source = source_map.get(source, source)  # Keep the original value if the source is neither 'trained' nor 'predicted'

    # Build the file path
    file_path = f"analysis/{ori_source}/seqs/{name}.csv" 
    
    # Read the data
    df = pd.read_csv(file_path)
    
    # Add the source information
    df['Source'] = f"{source} {name}"  # Changed to 'Training' or 'Predicted' + ' Anti-oxidant' or ' Non-cytotoxic'
    
    # Initialize the amino acid count dictionary
    aa_counts = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}
    
    # Iterate over each sequence and count the amino acid frequencies
    for seq in df['seq']:
        for aa in seq:  # Iterate over each amino acid in the sequence
            if aa in aa_counts:
                aa_counts[aa] += 1
    
    # Compute the total number of amino acids
    total_aa = sum(aa_counts.values())
    
    # Compute the proportion of each amino acid
    aa_proportion = {aa: count / total_aa for aa, count in aa_counts.items()}
    
    # Convert the proportions into a DataFrame
    aa_proportion_df = pd.DataFrame(list(aa_proportion.items()), columns=['Amino_Acid', 'Proportion'])
    aa_proportion_df['Source'] = f"{source} {name}"  # Add the source information
    
    return aa_proportion_df

# Merge all the datasets
names = ['anti-inflammatory']
all_data = []

# Load the data and merge it
for name in names:
    all_data.append(load_data(name, 'trained'))
    all_data.append(load_data(name, 'predicted'))

# Merge all the data into a single DataFrame
merged_df = pd.concat(all_data)

# Define the amino acid ordering list
amino_acid_order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 
                    'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 
                    'V', 'W', 'Y']

# Set the Amino_Acid column as a categorical type sorted by the defined order
merged_df['Amino_Acid'] = pd.Categorical(merged_df['Amino_Acid'], categories=amino_acid_order, ordered=True)

# Set the plot size
plt.figure(figsize=(16, 6))

# Define a custom palette: 'Training' and 'Predicted' are distinguished by different color schemes
custom_palette = {
    'Training anti-inflammatory': '#B0B0B0',  # Darker gray
    'Predicted anti-inflammatory': '#ADD8E6'   # More saturated blue
}

# Use Seaborn to draw the bar plot with the custom palette
ax = sns.barplot(
    data=merged_df, 
    x='Amino_Acid', 
    y='Proportion', 
    hue='Source', 
    dodge=True, 
    ci=None,
    palette=custom_palette
)

# Set the font sizes and labels
plt.xticks(fontsize=26, ha='right')  # Rotate the x-axis labels to avoid overlapping
plt.xlabel('')
plt.yticks(fontsize=26)
plt.ylabel('Proportion', fontsize=28)

# Set the spine styles
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')

# Set the tick mark width to 2
ax.tick_params(axis='both', width=2, length=12)

# Set the legend font size
ax.legend(prop={'size': 24})  # Set the font size to 24

# Adjust the layout to avoid overlapping labels
plt.tight_layout()

# Save the chart in PNG and SVG formats
plt.savefig("analysis/pic/propor_pic.png")
plt.savefig("analysis/pic/propor_pic.svg", format='svg')

# Display the figure
plt.show()
