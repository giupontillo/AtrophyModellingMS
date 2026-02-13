#!/data/anw/knw-work/g.pontillo/predictive_modelling/python-env/atrophymodelling/bin/python

#authors: m.barrantescepas@amsterdamumc.nl
#updates: 22 August 2024
#description: read epicenter individuals and select significant values

import os 
import re
import numpy as np 
import pandas as pd
from enigmatoolbox.plotting import plot_cortical, plot_subcortical
from enigmatoolbox.utils.parcellation import parcel_to_surface

#%%

def process_file(file_path, top_percent=5):
    #function to select top 5 and binarise    
    filename = file_path.split("/")[-1] 
    match = re.match(r".*_(epi)_.*\.csv", filename) #rvalues
    if not match:
        #print(f"Skipping {filename} (not epi_p)")
        return None
    
    df = pd.read_csv(file_path, header=None)  #no header
    #threshold = np.percentile(df[0], 5)
    #1 if value <= threshold AND p <= 0.05
    #binary_df = ((df[0] <= threshold) & (df[0] <= 0.05)).astype(int).to_frame()
    
    return df

def parse_filename(filename):
    # read filename information 
    pattern = r"(sub-[\w\d]+)_(fc|sc)_(ctx|sctx)_(epi(?:_p)?)_(d100|d400)\.csv" # sub-037595_fc_sctx_epi_d100.csv
    match = re.match(pattern, filename)
    if match:
        subject, modality, region, typ, atl = match.groups()
        return subject, modality, region, typ, atl
    else:
        return None, None, None, None, None
    
#%%
#plotting functions
def plot_cortical_top_epi(data, atlas, cmap, color_range, output_folder, filename):
    surf = parcel_to_surface(data, atlas) 
    plot_cortical(array_name=surf, surface_name="fsa5", size=(2000,1500), 
                  cmap=cmap, color_bar=True, color_range=color_range, 
                  screenshot=True, transparent_bg=False, 
                  filename=f"{output_folder}/{filename}.png")
    
def plot_subcortical_top_epi(data, cmap, color_range, output_folder, filename):
    plot_subcortical(array_name=data, ventricles=False, 
                     size=(2000,1500), cmap=cmap, color_bar=True, color_range=color_range, 
                     screenshot=True, transparent_bg=False,
                     filename=f"{output_folder}/{filename}.png")

#%%

# import epicenter data
dir="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling" 

folder_path=f"{dir}/output/june/cross_indv_epi_results"
output_folder = f"{dir}/output/oct/cross_indv_epi_plots_v2"
os.makedirs(output_folder, exist_ok=True)

csv_files = [f for f in os.listdir(folder_path) if f.endswith(".csv")]

#%%
#import atrophy data and match to check whether subjects missing
atrophydir = f'{dir}/data/atrophy_cross'
atrophy_csv= pd.read_csv(f'{atrophydir}/atrophy_baseline_SchaeferSubcortex414Parcels_individuals_wscore.csv')
missing_subjects=atrophy_csv[['id']]

#extract subjects from list of strings 
df_subj = pd.DataFrame({"file_name": csv_files})
df_subj[["id", "modality"]] = df_subj["file_name"].str.split("_", n=1, expand=True)

#check d100 
df100 = df_subj[df_subj["modality"] == "sc_ctx_epi_d400.csv"]
merged = df100.merge(missing_subjects, on=["id", "id"], how="outer", indicator=True)
different = merged[merged["_merge"] != "both"]

#%%
# Group files by modality and region type
groups = {}
for csv_file in csv_files:
    subject, modality, region, typ, atl = parse_filename(csv_file)
    if subject is None:
        continue  # skip files that don't match the pattern
    key = f"{modality}_{region}_{typ}_{atl}"
    if key not in groups:
        groups[key] = []
    groups[key].append(csv_file)

#%%
# Process each group
group_dfs = {}
for key, files in groups.items():
    all_subjects_data = []
    for f in files:
        full_path = os.path.join(folder_path, f)
        binary_df = process_file(full_path)
                
        if binary_df is None: 
            continue 
        
        parsed = parse_filename(f)
        if parsed[0] is None:
            print(f"Skipping {f}: filename pattern did not match")
            continue
        subject = parsed[0]
 
        binary_df = binary_df.T
        binary_df.columns = [f"V{i+1}" for i in range(binary_df.shape[1])]
        binary_df['subject'] = subject
        all_subjects_data.append(binary_df)
        
    # Skip the group if no valid files
    if not all_subjects_data:
        print(f"No valid files for group {key}, skipping.")
        continue

    #concatenate
    group_df = pd.concat(all_subjects_data)
    group_df.set_index('subject', inplace=True)

    output_file = os.path.join(output_folder, f"{key}_summary.csv")
    group_df.to_csv(output_file)
    print(f"Saved {output_file}")
    group_dfs[key] = group_df
    
#%%
group_sums = {}

for key, group_df in group_dfs.items():
    
    column_sum = group_df.sum(axis=0)  
    sum_df = pd.DataFrame(column_sum).T  
    group_sums[key] = sum_df

    output_file = os.path.join(output_folder, f"{key}_total.csv")
    sum_df.to_csv(output_file)
    print(f"Saved sum CSV for group {key}: {output_file}")

#%%

for key, sum_df in group_sums.items():
   
    data = sum_df.iloc[0].values  
    
    if "100" in key:
        atlas = 'schaefer_100_fsa5'
    elif "400" in key:
        atlas = 'schaefer_400_fsa5'
    else:
        print(f"Unknown resolution in key {key}, skipping...")
        continue
     
    cmap = "hot_r"  
    color_range = (0, np.max(data))  
    filename = f"{key}_top_epi"
    
    if "sctx" in key:
        plot_subcortical_top_epi(data, cmap, color_range, output_folder, filename)
    elif "ctx" in key:
        plot_cortical_top_epi(data, atlas, cmap, color_range, output_folder, filename)
    else:
        print(f"Unknown type in key {key}, skipping...")
        continue
    
    print(f"Plotted column sums for group {key}")