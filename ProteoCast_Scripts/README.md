# ProteocastScripts

A collection of Python and R scripts for protein sequence analysis, variant effect prediction, and structural analysis. These tools are designed for analyzing protein mutations using GEMME predictions, AlphaFold structures, and multiple sequence alignments (MSA).

## Scripts Overview

### Data Retrieval

| Script | Description |
|--------|-------------|
| `Retreive_AFstructure.py` | Downloads AlphaFold 3D structure (PDB format) for a given UniProt ID from the EBI AlphaFold database. |
| `Retreive_AFMSA.py` | Downloads AlphaFold multiple sequence alignment (A3M format) for a given UniProt ID. |

### MSA Processing

| Script | Description |
|--------|-------------|
| `AF_plot.py` | Generates a visual representation of a multiple sequence alignment showing sequence coverage and identity to the query sequence. |
| `RemoveXfromMSA.py` | Removes undefined amino acids ('X') from the query sequence in an MSA by either replacing them with the most common residue at that position or removing the column if it's mostly gaps. |

### Sequence Analysis

| Script | Description |
|--------|-------------|
| `msa_pdb_seqAlignment.py` | Performs pairwise alignment between a FASTA sequence and a PDB structure sequence, generating a position mapping file between the two. |
| `sequence.py` | Provides a `Sequence` class for handling protein sequences with support for FASTA I/O, gap analysis, and mutation compatibility checking. |
| `amino_acid.py` | Provides an `AminoAcid` class with mappings between amino acid representations (one-letter, three-letter codes, IDs) and handles non-standard amino acids. |

### Variant Effect Prediction

| Script | Description |
|--------|-------------|
| `local_confidence.py` | Calculates local confidence scores for mutational landscapesf predictions using weighted smoothing with a Gaussian kernel. Identifies low-confidence regions based on coverage, conservation, and score variability. |
| `GMM.py` | Applies a 3-component Gaussian Mixture Model to classify variant scores into neutral, mild, and impactful categories. Generates classification thresholds and visualization plots. |

### Structural Analysis

| Script | Description |
|--------|-------------|
| `run_rsa.py` | Calculates Relative Solvent Accessibility (RSA) for each residue in a PDB structure using the Shrake-Rupley algorithm. Handles multi-chain structures and non-standard amino acids. |
| `bfactors_change.py` | Maps variant scores and RSA values onto PDB B-factors for 3D visualization. Generates multiple PDB files with different metrics (residue class, sensitivity, RSA-weighted scores). |
| `prep4seg.py` | Prepares data for segmentation analysis by extracting mean mutational sensitivity scores and pLDDT values (if available from AlphaFold structures). |

### Segmentation Analysis

| Script | Description |
|--------|-------------|
| `main.r` | R script that runs protein segmentation analysis using the FPOP algorithm to identify regions of varying mutational sensitivity. |
| `segmentation_plot.py` | Generates visualization plots showing segmentation results with mutational sensitivity profiles and pLDDT heatmaps. |
| `utils.r` | R utility functions containing the core segmentation algorithm (`seg`) and wrapper function (`segmentation`) using the FPOP changepoint detection method. |

### Utilities

| Script | Description |
|--------|-------------|
| `utils.py` | Python utility functions including mutational landscape matrix reading, wild-type value removal, residue sensitivity analysis, PDB B-factor extraction, and sequence manipulation. |



## Usage Examples

### Download AlphaFold Structure
```bash
python Retreive_AFstructure.py P12345 -o ./structures/
```

### Generate MSA Plot
```bash
python AF_plot.py --msa path/to/alignment.fasta
```

### Calculate RSA
```bash
python run_rsa.py structure.pdb A
```

### Run GMM Classification
```bash
python GMM.py path/to/ProteoCast.csv true
```

### Run Segmentation (R)
```bash
Rscript main.r path/to/ProteoCast.csv
```

### Generate Segmentation Plot
```bash
python segmentation_plot.py --seg-file segmentation.csv --input-dir ./results/
```

## Dependencies

### Python
- numpy
- pandas
- matplotlib
- seaborn
- scikit-learn
- biopython
- gemmi
- prody

### R
- fpopw

## File Formats

- **FASTA**: Multiple sequence alignment files
- **A3M**: AlphaFold MSA format
- **PDB**: Protein structure files
- **CSV**: Tabular data for variant scores, mappings, and results
