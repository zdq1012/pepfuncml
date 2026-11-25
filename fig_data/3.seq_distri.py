import pandas as pd

# 1. Read the original CSV file
df = pd.read_csv('HP-select.csv')

# 2. Extract the rows where Anti-inflammatory == 1
filtered_df = df[df['Anti-inflammatory'] == 1]

# 3. Generate seq_names: name the sequences after sorting them alphabetically
filtered_df = filtered_df.copy()
filtered_df['seq_names'] = [f"seq_{i+1}" for i in range(len(filtered_df))]

# 4. Compute the length of each sequence
filtered_df['seq_length'] = filtered_df['Seq'].apply(len)

# 5. Count the length distribution
length_counts = filtered_df['seq_length'].value_counts().sort_index()

# 6. Save the length distribution to a CSV file
length_counts_df = length_counts.reset_index()
length_counts_df.columns = ['Sequence Length', 'Frequency']
length_counts_df.to_csv('sequence_length_distribution.csv', index=False)

# 7. Keep only the required columns
result = filtered_df[['Seq', 'seq_names']]

# 8. Export to a new CSV file
result.to_csv('filtered_sequences.csv', index=False)

print("✅ Generated filtered_sequences.csv and saved the sequence length distribution to sequence_length_distribution.csv")
