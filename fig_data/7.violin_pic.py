import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Set folder paths
predicted_dir = 'analysis/predicted/seqs_properties'
trained_dir = 'analysis/trained/seqs_properties'

# Function to read a file
def read_file(file_path):
    try:
        df = pd.read_csv(file_path)
        return df
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

# Function to extract the specified columns
def extract_columns(df, columns):
    return df[columns]

# Merge the data and add labels
def prepare_data(predicted_df, trained_df, label):
    predicted_df['Source'] = 'Predicted'
    trained_df['Source'] = 'Trained'
    predicted_df['Category'] = label
    trained_df['Category'] = label
    return pd.concat([predicted_df, trained_df], ignore_index=True)

# Columns to extract
columns = [
    'GRAVY', 'Hydrophobic_Moment', 'Aliphatic_Index', 'Aromaticity',
    'Instability_Index', 'Eisenberg_Hydrophobicity', 'Isoelectric_Point',
    'Hydrophobic_Ratio', 'Charge_Density', 'Charge'
]

# Files in the folder
files = ['Anti-inflammatory.csv']

# Process and merge the data
def combine_data():
    combined_data = pd.DataFrame()

    for file_name in files:
        predicted_file_path = os.path.join(predicted_dir, file_name)
        trained_file_path = os.path.join(trained_dir, file_name)

        predicted_df = read_file(predicted_file_path)
        trained_df = read_file(trained_file_path)

        if predicted_df is not None and trained_df is not None:
            # Extract the required columns
            predicted_df = extract_columns(predicted_df, columns)
            trained_df = extract_columns(trained_df, columns)
            
            # Merge the data
            merged_data = prepare_data(predicted_df, trained_df, file_name.split('.')[0])
            combined_data = pd.concat([combined_data, merged_data], ignore_index=True)

    return combined_data

# Function to draw the violin plot
def plot_violin(data, category_col, value_col, save_path):
    fig, ax = plt.subplots(figsize=(6, 6))
    
    # Replace underscores in the Y-axis label with spaces
    value_col_formatted = value_col.replace('_', ' ')
    
    # Define a custom color palette: 'Trained' is light gray, 'Predicted' is light blue
    custom_palette = {'Trained': 'lightgray', 'Predicted': 'lightblue'}
    
    # Draw the violin plot and apply the custom color palette
    sns.violinplot(x=category_col, y=value_col, hue="Source", data=data, split=True, inner=None, ax=ax,
                   linewidth=0, hue_order=['Trained', 'Predicted'], scale='width', width=0.7, palette=custom_palette)

    # Set font sizes and styles
    plt.xticks(fontsize=26, ha='center')  # Center-align the x-axis labels
    plt.xlabel('')
    plt.yticks(fontsize=26)  # Set the font size of the y-axis ticks
    plt.ylabel(value_col_formatted, fontsize=28)  # Replace the Y-axis label

    # Set the spine styles
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')

    # Set the width and length of the tick marks
    ax.tick_params(axis='both', width=2, length=12)

    # Set the legend font size
    ax.legend(prop={'size': 24})
    plt.tight_layout()

    # Save in PNG and SVG formats
    plt.savefig(save_path.replace(".png", ".svg"))
    plt.savefig(save_path)
    print(f"Saved plot for {value_col} to {save_path}")
    plt.close()

# Output save directory
output_dir = "analysis/pic/violin_plots"
os.makedirs(output_dir, exist_ok=True)

# Merge all the data and draw the violin plots
combined_data = combine_data()

# Draw a violin plot for each column and save it as a separate image
for column in columns:
    save_path = os.path.join(output_dir, f"{column}_violin_plot.png")
    plot_violin(combined_data, 'Category', column, save_path)

# Indicate completion
print("All plots saved.")
