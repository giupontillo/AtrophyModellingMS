#!/data/anw/knw-work/g.pontillo/predictive_modelling/python-env/atrophymodelling/bin/python

#authors: m.barrantescepas@amsterdamumc.nl
#updates: 22 Jan 2024

import pandas as pd

def load_csv_to_numpy(file_mapping, base_path):
    "reads multiple csv files into numpy arrays"
    arrays = {}
    for key, filename in file_mapping.items():
        arrays[key] = pd.read_csv(f"{base_path}/{filename}", header=None).to_numpy()
    return arrays