import pandas as pd
import matplotlib.pyplot as plt

pd.options.display.max_rows = 9999

# Read the XLSX file as a dataframe and sort the dataframe
# Change to path to the XLSX file, or put a copy of the XLSX file in the same folder as this python script
file_name = "D_domains"
df = pd.read_excel(file_name + '.xlsx')

domain_colors = {
    'PHP': "#f03d76",
    'CBS,HATPase_c,HATPase_c_2': "#55d847",
    'DRTGG': "#5652b0",
    'SPOR': "#e4df3a",
    'Other': "#616060", 
    'No neighbour': "#e4e4e4"
}

# Normalises the numbers of neighbours against the number of hydrogenase genes
df_plot = df[df['rank'] != 0]
ranks = sorted(df_plot['rank'].unique())
plot_data = {}
ref_genes = df['protein.id'].unique()

for r in ranks:
    counts = df_plot[df_plot['rank'] == r].groupby('Domains')['protein.id'].nunique().to_dict()
    n_absent = max(len(ref_genes) - sum(counts.values()), 0)
    counts['No neighbour'] = n_absent
    plot_data[r] = counts

plot_df = pd.DataFrame.from_dict(plot_data, orient='index').fillna(0)
plot_df_norm = plot_df.div(plot_df.sum(axis=1), axis=0)
cols = [c for c in plot_df_norm.columns if c != 'No neighbour'] + ['No neighbour']
plot_df_norm = plot_df_norm[cols]
plot_df_norm = plot_df_norm.sort_index()

# Plots the dataframe (not clade-specific)
plot_df_norm.plot(kind='bar', stacked=True, figsize=(10,6), color=[domain_colors[c] for c in plot_df_norm.columns])
plt.legend(title='Domain', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
plt.xlabel('Ranking', fontsize=15)
plt.ylabel('Amount of Neighbours', fontsize=15)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.title('Group D', fontsize=15)
plt.tight_layout()

# Save figure in pre-existing folder called "D_domainplots" (can be skipped)
#plt.savefig("D_domainplots/" + file_name + ".png")
plt.show()

# Makes clade-specific plots 
for clade, clade_group in df_plot.groupby('clade'):

    ranks = sorted(clade_group['rank'].unique())
    plot_data = {}
    ref_genes = df[df['clade'] == clade]['protein.id'].unique()

    for r in ranks:
        counts = clade_group[clade_group['rank'] == r].groupby('Domains')['protein.id'].nunique().to_dict()
        n_absent = max(len(ref_genes) - sum(counts.values()), 0)
        counts['No neighbour'] = n_absent
        plot_data[r] = counts

    plot_df = pd.DataFrame.from_dict(plot_data, orient='index').fillna(0)
    plot_df_norm = plot_df.div(plot_df.sum(axis=1), axis=0)
    cols = [c for c in plot_df_norm.columns if c != 'No neighbour'] + ['No neighbour']
    plot_df_norm = plot_df_norm[cols]
    plot_df_norm = plot_df_norm.sort_index()

    plot_df_norm.plot(kind='bar', stacked=True, figsize=(10, 6), color=[domain_colors[c] for c in plot_df_norm.columns])
    plt.xlabel('Ranking', fontsize=23)
    plt.ylabel('Amount of Neighbours', fontsize=23)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title(f'Group {clade}', fontsize=25)
    plt.legend(title='Domain', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    plt.tight_layout()
    
    # Save figures in pre-existing folder called "D_domainplots" (can be skipped)
    #plt.savefig("D_domainplots/" + clade + file_name + ".png")
    plt.show()