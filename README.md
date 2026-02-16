# The Role of Network Connectivity and Transcriptomic Vulnerability in Shaping Grey Matter Atrophy in Multiple Sclerosis

[![Paper Status](https://img.shields.io/badge/Paper-Under%20Review-yellow)]() 
[![Preprint DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15745212.svg)](https://doi.org/10.64898/2026.02.13.26346243) 
[![Github repo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15745212.svg)](https://zenodo.org/records/18632256)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)]()

This repository contains the code associated with our research paper to reproduce the results presented in the manuscript: 

**The Role of Network Connectivity and Transcriptomic Vulnerability in Shaping Grey Matter Atrophy in Multiple Sclerosis.** 

**Authors:** Mar Barrantes-Cepas, Mario Tranfa, David R. van Nederpelt, Ismail Koubiyr, Luigi Lorenzini, Birgit Helmlinger, Stefan Ropele, Daniela Pinter, Christian Enzinger, Tomas Uher, Manuela Vaneckova, Joep Killestein, Eva M.M. Strijbis, Martijn D. Steenwijk, Hugo Vrenken, Frederik Barkhof, Menno M. Schoonheim†, Giuseppe Pontillo†

**Correspondence:** M. Barrantes-Cepas, m.barrantescepas@amsterdamumc.nl 

† These authors contributed equally to this work. 

## Abstract

Clinical progression is strongly linked to grey matter atrophy in multiple sclerosis (MS), detectable early on MRI and progressing non-randomly across the brain. However, the mechanisms driving its spatio-temporal progression and individual variability remain unclear. Using MRIs from 2,187 participants, alongside normative data, we systematically investigated network-based mechanisms underlying MS-related atrophy. Regional atrophy colocalised with functional cortical hubs, supporting the *nodal stress* hypothesis, and propagated along anatomical and functional connections, consistent with *transneuronal degeneration*. *Lesional disconnection* and *transcriptomic vulnerability* played marginal roles. Patient- and subgroup-level analyses revealed that network-based mechanisms are specifically linked to MS-related neurodegeneration and may operate differently in distinct subtypes or disease phases. Atrophy patterns were anchored to the connectivity profiles of disease epicentres involving the visual, sensorimotor, and temporal cortices, and the hippocampi and thalami. Network-based measures enhanced the prediction of future atrophy progression in individual with MS, providing a mechanistic framework to understand neurodegeneration in MS.

![Fig. 1. ](figures/Fig1.png)
<small>**Fig. 1. Study workflow.** (a) The study integrated three datasets: MS Cohort and healthy controls, including T1-weighted (T1w) and FLAIR scans; the Human Connectome Project (HCP), providing structural and functional connectomes (SC/FC); and the Allen Human Brain Atlas (AHBA), providing gene expression data. (b) In the MS cohort, lesions were segmented using LST-AI and subsequently T1w scans were filled with NiftySeg. Disconnectome matrices were computed using the Lesion Quantification Toolkit (LQT), while cortical thickness and subcortical volumes were extracted with FreeSurfer. Atrophy maps were generated using Cohen’s d, and connectome and gene expression data were preprocessed using the ENIGMA toolbox. Neighbouring atrophy maps were obtained as described in Methods (A. Vo et *al.*). (c) Statistical analyses included spin permutation tests (n=10,000) to assess spatial correspondence. Created using https://BioRender.com and adapted from (S. Lariviere et *al.*, S. Gabarino et *al.*, and A. Vo et *al.*).</small>

## Repository Structure

```
.
├── figures/                # Figures 
├── LICENSE                 # License 
├── README.md               # This file
└── scripts/                # Analysis scripts and models
    ├── 01.group            	# group-level analyses
    ├── 02.indv		    	    # individual-level analyses
    ├── 03.epicenter-mapping	# epicenter mapping analyses
    ├── 04.GAM			        # generalised additive modelling
    ├── 05.hydra		        # HYDRA clustering
    └── 06.stats		        # data org and statistical analyses 


```

### Data Availability

- MS imaging data is subject to ethical restrictions and available upon reasonable request to m.barrantescepas@amsterdamumc.nl
- Normative datasets were used from the ENIGMA toolbox but are available: 
	- Human Connectome Project (HCP): https://www.humanconnectome.org/ 
	- Allen Human Brain Atlas (AHBA): https://human.brain-map.org/

## Requirements

- R >= 4.1.1
- Python >= 3.8
- ENIGMA Toolbox (S. Lariviere et *al.*): https://enigma-toolbox.readthedocs.io/en/latest/
- HYDRA (E. Varol, A. Sotiras, C. Davatzikos): https://github.com/evarol/HYDRA

## Bibliography

- S. Larivière et *al.*, The ENIGMA Toolbox: multiscale neural contextualization of multisite neuroimaging datasets. Nature Methods 18, 698–700 (2021).
- S. Garbarino et *al.*, Differences in topological progression profile among neurodegenerative diseases from imaging data. Elife 8,  (2019).
- A. Vo et *al.*, Network connectivity and local transcriptomic vulnerability underpin cortical atrophy progression in Parkinson's disease. Neuroimage Clin 40, 103523 (2023).
- E. Varol, A. Sotiras, C. Davatzikos, HYDRA: Revealing heterogeneity of imaging and genetic patterns through a multiple max-margin discriminative analysis framework. Neuroimage 145, 346–364 (2017).

## Citation

If you use this code or data in your research, please cite our paper:

```bibtex
@article{YourName2026,
  title={TITLE},
  author={[Author 1] and [Author 2] and [Author 3]},
  journal={[Journal Name]},
  year={2026},
  note={Under Review}
}
```

**Note:** This citation format will be updated once the paper is published.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

