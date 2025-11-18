import pandas as pd
import matplotlib.pyplot as plt

pd.options.display.max_rows = 9999

# Read the CSV as a dataframe and sort the dataframe
# Change to path to the CSV file, or put a copy of the CSV file in the same folder as this python script
# Comment: You can add or remove which columns to include, but 'start', 'end', 'protein.id', 'is.neighbour', 'strand' and 'COG_category' are necessary for the calculations

file_name = "50_15" # Change accoring to your filename
df = pd.read_csv('df_' + file_name + '.csv', usecols=['seqid', 'start', 'end', 'protein.id', 'ID', 'is.neighbour', 'clade', 'COG_category', 'PFAMs', 'product', 'Description', 'strand'])
grouped = df.sort_values(by=['protein.id']).groupby('protein.id')

# Function that ranks the neighbour protein genes based on their distance to the protein of interest (the hydrogenase genes) and calculates the distance on the genome
def neighbour_rank(group_df):

    protein_of_interest = None
    neighbour_poi = []

    for index, row in group_df.iterrows():
        if row['is.neighbour'] is False:
            protein_of_interest = row
        else: 
            neighbour_poi.append(row)

    poi_df = pd.DataFrame([protein_of_interest])
    ne_poi_df = pd.DataFrame(neighbour_poi)

    ne_upstream = []
    ne_downstream = []
    dist_down = []
    dist_up = []

    poi_start = protein_of_interest['start']
    poi_end = protein_of_interest['end']
    poi_strand = protein_of_interest['strand']

    for index, row in ne_poi_df.iterrows():
        ne_start = row['start']
        ne_end = row['end']

        if poi_strand == '+':
            if ne_start < poi_start: 
                dist_poi = poi_start - ne_end
                dist_up.append(dist_poi)
                ne_upstream.append(row)
            else: 
                dist_poi = ne_start - poi_end
                dist_down.append(dist_poi)
                ne_downstream.append(row)
        else:
            if ne_start > poi_start:
                dist_poi = ne_start - poi_end
                dist_up.append(dist_poi)
                ne_upstream.append(row)
            else:
                dist_poi = poi_start - ne_end
                dist_down.append(dist_poi)
                ne_downstream.append(row)

    df_upstream = pd.DataFrame(ne_upstream)
    df_downstream = pd.DataFrame(ne_downstream)

    if not df_upstream.empty:
        df_upstream['distance'] = dist_up
        df_upstream = df_upstream.sort_values(by='start', ascending=False)
        df_upstream['rank'] = [-i for i in range(1, len(df_upstream) + 1)]
    
    if not df_downstream.empty:
        df_downstream['distance'] = dist_down
        df_downstream = df_downstream.sort_values(by='start', ascending=True)
        df_downstream['rank'] = [i for i in range(1, len(df_downstream) + 1)]

    if poi_strand == '-':
        if not df_upstream.empty:
            df_upstream = df_upstream.sort_values(by='distance', ascending=True)
            df_upstream['rank'] = [-i for i in range(1, len(df_upstream) + 1)]

        if not df_downstream.empty:
            df_downstream = df_downstream.sort_values(by='distance', ascending=True)
            df_downstream['rank'] = [i for i in range(1, len(df_downstream) + 1)]

    poi_df['distance'] = 0
    poi_df['rank'] = 0
    
    df_ranked = pd.concat([df_upstream, df_downstream], ignore_index=True)

    return df_ranked

list_of_ranked_neighbours = []

for protein_id, group_df in grouped:
    ranked_neighbours = neighbour_rank(group_df.reset_index(drop=True))
    list_of_ranked_neighbours.append(ranked_neighbours)
    
df_combined = pd.concat(list_of_ranked_neighbours, ignore_index=True)
sorted_by_rank = df_combined.sort_values(by='rank')

# Creates a CSV file with the ranks, distances and previously read columns (can be skipped)
sorted_by_rank.to_csv('ranks_and_distance_df_' + file_name + '.csv', index=False)

# Creates a dataframe to plot
plot_data = sorted_by_rank.groupby(['rank', 'COG_category']).size().unstack(fill_value=0)

cog_colors = {
    'J': "#9fb90a", 'A': '#aec7e8', 'K': '#98df8a', 'L': '#ffbb78',
    'B': '#2ca02c', 'D': '#ff7f0e', 'Y': '#d62728', 'V': '#ff9896',
    'T': '#c5b0d5', 'M': '#9467bd', 'N': '#8c564b', 'Z': '#c49c94',
    'W': "#c958a7", 'U': '#f7b6d2', 'O': '#17becf', 'X': '#9edae5',
    'C': '#1f77b4', 'G': '#dbdb8d', 'E': "#e7208a", 'F': '#66c2a5',
    'H': "#e793c7", 'I': '#a6d854', 'P': "#ffd622", 'Q': "#d15e5e",
    'R': '#1a55FF', 'S': '#8da0cb', 'unknown': "#616060", 
    'No neighbour': "#e4e4e4"
}

# Plots the dataframe (not clade-specific)
plot_data.plot(kind='bar', stacked=True, figsize=(10, 6), color=[cog_colors[c] for c in plot_data.columns])

custom_labels = {
    'J': 'J: Translation, ribosomal structure and biogenesis',
    'A': 'A: RNA processing and modification',
    'K': 'K: Transcription',
    'L': 'L: Replication, recombination and repair',
    'B': 'B: Chromatin structure and dynamics',
    'D': 'D: Cell cycle control, cell division, chromosome partitioning',
    'Y': 'Y: Nuclear structure',
    'V': 'V: Defense mechanisms',
    'T': 'T: Signal transduction mechanisms',
    'M': 'M: Cell wall/membrane/envelope biogenesis',
    'N': 'N: Cell motility',
    'Z': 'Z: Cytoskeleton',
    'W': 'W: Extracellular structures',
    'U': 'U: Intracellular trafficking, secretion, and vesicular transport',
    'O': 'O: Posttranslational modification, protein turnover, chaperones',
    'X': 'X: Mobilome: prophages, transposons',
    'C': 'C: Energy production and conversion',
    'G': 'G: Carbohydrate transport and metabolism',
    'E': 'E: Amino acid transport and metabolism',
    'F': 'F: Nucleotide transport and metabolism',
    'H': 'H: Coenzyme transport and metabolism',
    'I': 'I: Lipid transport and metabolism',
    'P': 'P: Inorganic ion transport and metabolism',
    'Q': 'Q: Secondary metabolites biosynthesis, transport and metabolism',
    'R': 'R: General function prediction only',
    'S': 'S: Function unknown'
}

handles, labels = plt.gca().get_legend_handles_labels()
new_labels = [custom_labels.get(lbl, lbl) for lbl in labels]
plt.legend(handles, new_labels, title='COG Category', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
plt.xlabel('Ranking', fontsize=15)
plt.ylabel('Amount of Neighbours', fontsize=15)
plt.title('All Groups', fontsize=15)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
#plt.legend(title='COG Category', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save figure in pre-existing folder called "plots" (can be skipped)
#plt.savefig("plots/" + file_name + ".png")
plt.show()

group_by_clade = sorted_by_rank.groupby(['clade', 'rank', 'COG_category']).size().unstack(fill_value=0)

# Makes clade-specific plots 
for clade, clade_group in group_by_clade.groupby(level=0):
    subset = clade_group.droplevel(0)
    subset.plot(kind='bar', stacked=True, figsize=(10, 6), color=[cog_colors[c] for c in plot_data.columns])
    plt.xlabel('Ranking', fontsize=23)
    plt.ylabel('Amount of Neighbours', fontsize=23)
    #plt.legend(title='COG Category', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    #plt.legend(handles, new_labels, title='COG Category', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.legend('', frameon=False)
    plt.title(f'Group {clade}', fontsize=25)
    plt.tight_layout()
    
    # Save figure in pre-existing folder called "plots" (can be skipped)
    #plt.savefig("plots/" + clade + ".png")
    plt.show()

