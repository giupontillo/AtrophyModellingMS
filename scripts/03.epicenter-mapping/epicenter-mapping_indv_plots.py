#!/data/anw/knw-work/g.pontillo/predictive_modelling/python-env/atrophymodelling/bin/python

#authors: m.barrantescepas@amsterdamumc.nl
#updates: 22 August 2024
#description: read epicenter individuals and plot top 5% epicenters 

import os 
import re
import numpy as np 
import pandas as pd
from enigmatoolbox.plotting import plot_cortical, plot_subcortical
from enigmatoolbox.utils.parcellation import parcel_to_surface
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import fdrcorrection

#%%

def select_top(file_path, top_percent=5):
    #function to select top 5 and binarise    
    filename = file_path.split("/")[-1] 
    match = re.match(r".*_(epi_p)_.*\.csv", filename)
    if not match:
        #print(f"Skipping {filename} (not epi_p)")
        return None
    
    df = pd.read_csv(file_path, header=None)  #no header
    threshold = np.percentile(df[0], 5)
    #1 if value <= threshold AND p <= 0.05
    binary_df = ((df[0] <= threshold) & (df[0] <= 0.05)).astype(int).to_frame()
    
    return binary_df

def collect_r_values_for_group(key, files, folder_path, output_folder):
    """
    For a given group (key), reads all _epi_ csv, concatenates r-values per region,
    and saves a combined dataframe with subjects as rows and regions as columns.
    """
    all_subjects_data = []

    for f in files:
        # Only process files containing "_epi_"
        if "_epi_p_" in f:
            continue
        
        full_path = os.path.join(folder_path, f)
        filename = os.path.basename(f)

        match = re.match(r"(sub-[A-Za-z0-9]+)", filename)
        if not match:
            print(f"Skipping {filename} (subject ID not found)")
            continue
        subject = match.group(1)

        try:
            df = pd.read_csv(full_path, header=None)
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            continue

        # Transpose so regions become columns
        df = df.T
        df.columns = [f"region_{i+1}" for i in range(df.shape[1])]
        df["subject"] = subject

        all_subjects_data.append(df)

    if not all_subjects_data:
        print(f"No valid _epi_ files for {key}, skipping.")
        return None

    # Combine all subjects
    group_df = pd.concat(all_subjects_data, ignore_index=True)
    group_df.set_index("subject", inplace=True)

    # Save output
    os.makedirs(output_folder, exist_ok=True)
    output_file = os.path.join(output_folder, f"{key}_r_values_summary.csv")
    group_df.to_csv(output_file)
    print(f"Saved combined r-values for {key} → {output_file}")

    return group_df

def fisher_r_to_z(r):
    r = np.clip(r, -0.999999, 0.999999)  # avoid infinities
    return np.arctanh(r)

def parse_filename(filename):
    # read filename information 
    pattern = r"(sub-[\w\d]+)_(fc|sc)_(ctx|sctx)_(epi(?:_p)?)_(d100|d400)\.csv" # sub-037595_fc_sctx_epi_d100.csv
    match = re.match(pattern, filename)
    if match:
        subject, modality, region, typ, atl = match.groups()
        return subject, modality, region, typ, atl
    else:
        return None, None, None, None, None
    
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
output_folder = f"{dir}/output/oct/cross_indv_epi_plots_v3"
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
# group_dfs = {}
# for key, files in groups.items():
#     all_subjects_data = []
#     for f in files:
        
#         full_path = os.path.join(folder_path, f)
#         binary_df = select_top(full_path)
                
#         if binary_df is None: 
#             continue 
        
#         parsed = parse_filename(f)
#         if parsed[0] is None:
#             print(f"Skipping {f}: filename pattern did not match")
#             continue
#         subject = parsed[0]
 
#         binary_df = binary_df.T
#         binary_df.columns = [f"V{i+1}" for i in range(binary_df.shape[1])]
#         binary_df['subject'] = subject
#         all_subjects_data.append(binary_df)
        
#     # Skip the group if no valid files
#     if not all_subjects_data:
#         print(f"No valid files for group {key}, skipping.")
#         continue

#     #concatenate
#     group_df = pd.concat(all_subjects_data)
#     group_df.set_index('subject', inplace=True)

#     output_file = os.path.join(output_folder, f"{key}_top5_summary.csv")
#     group_df.to_csv(output_file)
#     print(f"Saved {output_file}")
#     group_dfs[key] = group_df
    

        
        
#%%
# group_sums = {}

# for key, group_df in group_dfs.items():
    
#     column_sum = group_df.sum(axis=0)  
#     sum_df = pd.DataFrame(column_sum).T  
#     group_sums[key] = sum_df

#     output_file = os.path.join(output_folder, f"{key}_total.csv")
#     sum_df.to_csv(output_file)
#     print(f"Saved sum CSV for group {key}: {output_file}")

#%%

# for key, sum_df in group_sums.items():
   
#     data = sum_df.iloc[0].values  
    
#     if "100" in key:
#         atlas = 'schaefer_100_fsa5'
#     elif "400" in key:
#         atlas = 'schaefer_400_fsa5'
#     else:
#         print(f"Unknown resolution in key {key}, skipping...")
#         continue
     
#     cmap = "hot_r"  
#     color_range = (0, np.max(data))  
#     filename = f"{key}_top_epi"
    
#     if "sctx" in key:
#         plot_subcortical_top_epi(data, cmap, color_range, output_folder, filename)
#     elif "ctx" in key:
#         plot_cortical_top_epi(data, atlas, cmap, color_range, output_folder, filename)
#     else:
#         print(f"Unknown type in key {key}, skipping...")
#         continue
    
#     print(f"Plotted column sums for group {key}")
    
    
#%%

group_dfs = {}

# concatenate r_values for all the subjects
for key, files in groups.items(): 
    group_df=collect_r_values_for_group(key, files, folder_path, output_folder)
    if group_df is not None:
        group_dfs[key] = group_df
        
#%%

group_stats = {}

for key, df_rvals in group_dfs.items():
    print(f"Processing group: {key}")

    # Fisher transform r → z
    z_vals = df_rvals.apply(fisher_r_to_z)

    # One-sample t-test vs 0 (per region)
    t_vals, p_vals = ttest_1samp(z_vals, popmean=0, axis=0, nan_policy='omit')
    _, p_fdr = fdrcorrection(p_vals) #fdr correction
    
    # Alternative: Bonferroni correction
    n_tests = len(df_rvals.columns) #num of regions 
    p_bonf = np.minimum(p_vals * n_tests, 1.0)
    
    # One-tailed test for r > 0
    p_one_tail = p_vals / 2
    p_one_tail[t_vals < 0] = 1 
    _, p_fdr_one_tail = fdrcorrection(p_one_tail)

    # Combine results
    results_df = pd.DataFrame({
        'region': df_rvals.columns,
        't_value': t_vals,
        'p_value': p_vals,
        'p_fdr': p_fdr, 
        'p_bonf': p_bonf,
        'p_value_one_tail': p_one_tail,
        'p_fdr_one_tail': p_fdr_one_tail,
        'n_subjects': z_vals.notna().sum().values
    })

    # Save
    output_file = os.path.join(output_folder, f"{key}_r_to_z_ttest_results.csv")
    results_df.to_csv(output_file, index=False)
    
    group_stats[key] = results_df

    print(f"Saved results for {key} → {output_file}")

print("All groups processed successfully.")

#%%

#plot epicenters
for key, df in group_stats.items():
    print(f"Plotting for group: {key}")
    
    # Cortical plotting example
    # Separate t-values by significance (p < 0.05)
    correction = 'p_fdr_one_tail'
    sig_mask = df[correction] < 0.05
    t_vals = df['t_value'].values
    t_vals_sig = t_vals.copy()
    t_vals_sig[~sig_mask] = 0  # 0 non-significant regions

    t_vals_nonsig = t_vals.copy()
    t_vals_nonsig[sig_mask] = 0  # 0 significant regions
    
    if "100" in key:
        atlas = 'schaefer_100_fsa5'
    elif "400" in key: 
        atlas = 'schaefer_400_fsa5'
    else: 
        print(f"Unknown resolution in key {key}, skipping...")
        continue

    if "sctx" in key.lower():
        #print("sctx")
        # Plot significant regions
        plot_subcortical_top_epi(
            data=t_vals_sig,
            cmap="coolwarm",
            #color_range=(-np.max(np.abs(t_vals)), np.max(np.abs(t_vals))),
            color_range=(-20, 20),
            output_folder=output_folder,
            filename=f"{key}_subcortical_{correction}_tvals_sig"
        )

        # Plot non-significant regions
        plot_subcortical_top_epi(
            data=t_vals_nonsig,
            cmap="coolwarm",
            #color_range=(-np.max(np.abs(t_vals_nonsig)), np.max(np.abs(t_vals_nonsig))),
            color_range=(-20, 20),
            output_folder=output_folder,
            filename=f"{key}_subcortical_{correction}_tvals_nonsig"
        )
        
    elif "ctx" in key.lower():
        #print("ctx")
        # Plot significant regions
        plot_cortical_top_epi(
            data=t_vals_sig,
            atlas=atlas,   
            cmap="coolwarm",
            #color_range=(-np.max(np.abs(t_vals)), np.max(np.abs(t_vals))),
            color_range=(-20, 20),
            output_folder=output_folder,
            filename=f"{key}_cortical_{correction}_tvals_sig"
        )

        # Plot non-significant regions
        plot_cortical_top_epi(
            data=t_vals_nonsig,
            atlas=atlas,
            cmap="coolwarm",
            #color_range=(-np.max(np.abs(t_vals_nonsig)), np.max(np.abs(t_vals_nonsig))),
            color_range=(-20, 20),
            output_folder=output_folder,
            filename=f"{key}_cortical_{correction}_tvals_nonsig"
        )
    else:
         print(f"Unknown type in key {key}, skipping...")
    continue