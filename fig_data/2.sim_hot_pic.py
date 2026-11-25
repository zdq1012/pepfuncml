import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Read the saved CSV file
csv_file = "fig_data/cal_result/HP-select/HP-select.fasta_distance_matrix.csv"  # Replace with your actual CSV file path
distance_matrix_df = pd.read_csv(csv_file, index_col=0)  # index_col=0 means the first column holds the row labels (sequence IDs)

# Get the size of the matrix
num_sequences = distance_matrix_df.shape[0]

# Symmetric filling: mirror the lower-left values to the upper-right
for i in range(1, num_sequences):
    for j in range(i):
        distance_matrix_df.iloc[j, i] = distance_matrix_df.iloc[i, j]

# Set all the values on the diagonal to 1
for i in range(num_sequences):
    distance_matrix_df.iloc[i, i] = 1

# Create a custom color map (from sky blue to pink)
colors = ["skyblue", "pink"]
cmap = LinearSegmentedColormap.from_list("skyblue_to_pink", colors)

# Use seaborn's clustermap to cluster and draw the heatmap, removing the frame and applying the custom color map
g = sns.clustermap(distance_matrix_df, annot=False, cmap=cmap, fmt='.2f', linewidths=0, figsize=(6, 6), cbar_pos=(0, .2, .03, .4))

# Remove the ticks and labels of the x- and y-axes
g.ax_heatmap.set_xticklabels([])  # Remove the x-axis labels
g.ax_heatmap.set_yticklabels([])  # Remove the y-axis labels
g.ax_heatmap.set_xticks([])  # Remove the x-axis ticks
g.ax_heatmap.set_yticks([])  # Remove the y-axis ticks

# Remove the dendrogram labels
g.ax_row_dendrogram.set_xticks([])  # Remove the x-axis ticks of the row dendrogram
g.ax_col_dendrogram.set_yticks([])  # Remove the y-axis ticks of the column dendrogram

# Remove the color bar labels
g.cax.set_visible(False)  # Hide the color bar labels

# Display the heatmap
plt.tight_layout()  # Adjust the layout to avoid overlapping labels
plt.show()

# Save the heatmap as an SVG file
output_image_svg = "fig_data/cal_result/HP-select/distance_matrix_heatmap_clustermap_no_text.svg"  # Save in SVG format
plt.savefig(output_image_svg, format='svg', bbox_inches='tight', pad_inches=0)

# Save the heatmap as a high-resolution PNG file
output_image_png = "fig_data/cal_result/HP-select/distance_matrix_heatmap_clustermap_no_text.png"  # Save in PNG format
plt.savefig(output_image_png, format='png', dpi=300, bbox_inches='tight', pad_inches=0)  # Set dpi=300 and save as a high-resolution PNG
