import pandas as pd
import matplotlib.pyplot as plt

pd.options.display.max_rows = 9999

# Read the CSV or as a dataframe and sort the dataframe
# Change to path to the CSV file, or put a copy of the CSV file in the same folder as this python script
# Comment: You can add or remove which columns to include, but 'start', 'end', 'protein.id' and 'is.neighbour' are necessary for the calculations
df = pd.read_csv('df_300_15.csv', usecols=['start', 'end', 'PIGI', 'ID', 'is.neighbour'])
sort_df = df.sort_values(by=['PIGI', 'start'])

# Calculate the intergenic space. The "last" gene will not have a gap. 
sort_df['gap_to_next'] = (
    sort_df.groupby('PIGI')['start']
      .shift(-1) - sort_df['end']
)

# Optional: Save the gaps for all proteins of interest and their neighbours in a CSV file (the columns chosen in the CSV reading + a column with gap sizes will be included)
# Includes genes with "no gap" 
# sort_df.to_csv('all_gene_gaps.csv', index=False)

# Save only the neighbours in a new dataframe 
neighbour_df = sort_df[sort_df['is.neighbour'] == True].copy()

# Include only neighbours that has a gap
neighbour_df = neighbour_df[neighbour_df['gap_to_next'].notna()]

# Optional: Remove negative gaps (overlapping genes) or zero gaps
#neighbour_df = neighbour_df[neighbour_df['gap_to_next'] > 0]

# Optional: Convert gap sizes to integers (will look prettier in the CSV file) 
neighbour_df['gap_to_next'] = neighbour_df['gap_to_next'].astype(int)

# Sort by gap size and compute cumulative count
neighbour_df = neighbour_df.sort_values('gap_to_next')
neighbour_df['cumulative_neighbours'] = range(1, len(neighbour_df) + 1)

# Plot the gap sizes
plt.figure(figsize=(10, 6))
plt.plot(neighbour_df['gap_to_next'], neighbour_df['cumulative_neighbours'], linestyle='-', color='blue')
plt.xlabel('Gap Size', fontsize=15)
plt.ylabel('Cumulative Number of Neighbours', fontsize=15)
#plt.title('Cumulative Distribution of Neighbouring Proteins by Gap Size', fontsize=15)
plt.grid(True)
plt.show()

# Optional: Save the gaps for the neighbours in a CSV file (the columns chosen in the CSV reading + a column with gap sizes + a column with the cumulative count will be included)
# Does not include proteins of interest or genes with "no gap"
#neighbour_df.to_csv('neighbour_gaps.csv', index=False)