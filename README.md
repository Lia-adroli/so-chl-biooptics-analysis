# Southern Ocean Phytoplankton PSC & PFT Reconstruction (1998–2024)

This repository contains the full data processing, model training, validation, and reconstruction workflow used to estimate phytoplankton size classes (PSCs) and phytoplankton functional types (PFTs) across the Southern Ocean (40–90°S) for the spring and summer period 1998–2024.

## Approach
The approach integrates:
1. HPLC pigment observations (SeaBASS + Hayward 2024)
2. Satellite chlorophyll from the Ocean Colour Climate Change Initiative (OC-CCI) Level-3 daily
3. NNLS pigment linear mixing
4. Gradient Boosting Machine (GBM) regression
5. Bin-based empirical reconstruction
6. Independent OC-CCI-based validation

---

### Core Variables

| Variable | Description |
|----------|-------------|
| `Tchla`  | In situ total chlorophyll (HPLC), mg m⁻³ |
| `OCChla` | Satellite chlorophyll from Ocean Colour Climate Change Initiative (OC-CCI), mg m⁻³ |

---

### Pigment Variables (8)

| Abbreviation | Full Name |
|--------------|-----------|
| Fuco     | Fucoxanthin |
| Per      | Peridinin |
| X19hex   | 19′-Hexanoyloxyfucoxanthin |
| Allo     | Alloxanthin |
| X19but   | 19′-Butanoyloxyfucoxanthin |
| Chl_b    | Chlorophyll b |
| Zea      | Zeaxanthin |
| DVChla   | Divinyl Chlorophyll a |

---

### Metadata

- `lat`, `lon` — geographic coordinates  
- `time` — sampling date (dd/mm/yyyy)
- **Source information**
  - For **SeaBASS**, full metadata are provided (cruise and station information).
  - For **AlexiHayw (2024)**, only the citation is used as the data source. Detailed cruise or time metadata should be obtained directly from the original Hayward (2024) reference.

---

### Spatial–Temporal Coverage

- **Region:** Southern Ocean (40–90°S, 180°W–180°E)  
- **Depth range:** 0–20 m (surface layer only)  
- **Period:** 1998–2024  

---

## Data Processing Workflow (Analysis Code Only)

1. **Initial Data Filtering**
   - ≤ 50% relative error retained from OCChla and TChla

2. **Year-Stratified Split**
   - 70% training per year  
   - 30% independent validation per year  

3. **NNLS Pigment Decomposition**
   - Ecologically constrained non-negative unmixing  
   - Linear diagnostic pigment sum (DP_lin)

4. **GBM Total Chlorophyll Reconstruction**
   - Input features:
     - Raw pigments  
     - Log-transformed pigments  
     - DP_lin and log(DP_lin)

5. **GBM-Anchored PSC & PFT Fractions**

6. **Bin-Based Empirical Model**
   - 100 log-spaced chlorophyll bins
   - 4th-degree polynomial fits
   - Weighted by bin population

7. **Independent OC-CCI Validation**
   - Apply bin model to satellite chlorophyll
   - Compare against GBM-derived PSC/PFT
---

## Software Requirements

- MATLAB R2024a or newer  
- Statistics and Machine Learning Toolbox  
- Optimization Toolbox  
---

## Data Source and DOI

The full dataset associated with this repository is archived at:

**Zenodo DOI:**  
https://doi.org/10.5281/zenodo.17875100

---

## Contact

**Nurmalia Adroli**  
PhD Candidate, Research School of Earth Science,
Australian National University  
nurmalia.adroli@anu.edu.au
