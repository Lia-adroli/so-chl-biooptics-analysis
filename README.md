# Southern Ocean Phytoplankton Functional Types and Size-Class Diversity, Phenology, and Ecosystem Regime Shifts from Ocean Colour and Abundance-Based Diagnostics

This repository contains the full data processing, model training, validation, and reconstruction workflow used to estimate phytoplankton size classes (PSCs) and phytoplankton functional types (PFTs) across the Southern Ocean (40–90°S) for the spring and summer period 1998–2024.

## Approach
The approach integrates:
1. HPLC pigment observations (SeaBASS + Hayward 2024)  
2. Satellite chlorophyll from the Ocean Colour Climate Change Initiative (OC-CCI) Level-3 daily  
3. NNLS pigment linear mixing to derive ecologically constrained pigment contributions  
4. Gradient Boosting Machine (GBM) regression with K-fold cross-validation for total chlorophyll (Tchla) reconstruction  
5. Bin-based empirical reconstruction of PSCs and PFTs as functions of log₁₀(Tchla)  
6. Independent OC-CCI-based validation using the bin model applied to satellite chlorophyll

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
  - `DP_lin` and `log(DP_lin)`
- Gradient Boosting Machine (GBM) regression
- **K-fold cross-validation (default K = 5)** applied to the training set to obtain out-of-fold predictions of Tchla and robust skill metrics (r, R², RMSE)

5. **GBM-Anchored PSC & PFT Fractions**
- GBM-derived total chlorophyll (`Tchla_GBM`) used as the denominator
- NNLS pigment contributions converted to:
  - Phytoplankton size classes (PSC: micro, nano, pico)
  - Phytoplankton functional types (PFT: diatoms, dinoflagellates, prymnesiophytes, cryptophytes, Prok+Prochlorococcus, PicoEuk+Green)
- Fractions constrained to 0–1 and then expressed as concentrations (mg m⁻³)

6.  **Bin-Based Empirical Model**
- 100 log-spaced Tchla bins
- For each PSC/PFT:
  - Mean and standard deviation of fractional contribution (% of total) computed per bin
  - 4th-degree polynomial fits of PSC/PFT (%) as a function of log₁₀(Tchla)
  - Fits weighted by bin population (weighted R², r, and approximate p-values)

7. **Independent OC-CCI Validation**
- Reconstructed PSC and PFT concentrations from OCChla compared against GBM-derived PSC/PFT for the independent validation subset
- Bin-model polynomial relationships applied to satellite OC-CCI chlorophyll (OCChla)

---
## Software Requirements

- MATLAB R2024a or newer  
- Statistics and Machine Learning Toolbox  
- Optimization Toolbox  

---
## Citation

If you use this code or dataset, please cite:

Adroli, N., Hayward, A., Strutton, P., & Ellwood, M. (2025).
*Southern Ocean phytoplankton functional types and size-class diversity, phenology, and ecosystem regime shifts from ocean colour and abundance-based diagnostics (1998–2024).*  
Zenodo. https://doi.org/10.5281/zenodo.17875100

---
## Contact

Nurmalia (Lia) Adroli  
PhD Candidate, Earth Science  
Remote Sensing of Ocean Biogeochemistry  
Research School of Earth Sciences  
The Australian National University   
https://orcid.org/0009-0000-2982-1900  
