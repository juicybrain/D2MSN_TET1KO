

import pandas as pd

import sys

if len(sys.argv) < 2:
    print("Usage: python example.py inputfile")
    sys.exit(1)

        
file_path = sys.argv[1]
data = pd.read_csv(file_path, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'diff'])

# Sort the data by chromosome and start position to ensure correct distance calculation
data.sort_values(by=['Chromosome', 'Start'], inplace=True)

# Calculate the distance to the previous mutation (upstream)
data['Distance to Previous'] = data.groupby('Chromosome')['Start'].diff().fillna(0).astype(int)

# Calculate the distance to the next mutation (downstream)
data['Distance to Next'] = data.groupby('Chromosome')['Start'].shift(-1) - data['Start']
data['Distance to Next'].fillna(0, inplace=True)  # Handle the last mutation which has no next mutation

# Calculate the average of the distances to the previous and next mutations
data['Intermutational Distance'] = ((data['Distance to Previous'] + data['Distance to Next']) / 2).astype(int)

# Drop the individual distance columns if they are not needed
data.drop(['Distance to Previous', 'Distance to Next'], axis=1, inplace=True)

# Save the updated dataframe to a new file
output_file_path = file_path + 'output_file_with_distances.txt'    # Update with your desired output file path
data.to_csv(output_file_path, sep='\t', index=False, header=False)

















