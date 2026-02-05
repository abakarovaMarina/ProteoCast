# ProteoCast

**Proteome-wide Prediction of the Functional Impact of Missense Variants**

ProteoCast is a scalable and interpretable computational framework for proteome-wide classification of genetic variants and functional protein site identification. Leveraging evolutionary information from protein sequences across organisms, ProteoCast enables researchers to:

- **Predict variant effects** across entire proteomes with high accuracy
- **Classify mutations** as functionally neutral, mild, or impactful
- **Identify sensitive residues** critical for protein function
- **Discover regulatory sites** in unstructured protein regions
- **Prioritize targets** for CRISPR genome editing experiments

Built on the GEMME evolutionary model (Laine et *al.*), ProteoCast provides comprehensive mutational landscapes for all protein isoforms, complete with confidence metrics and 3D structural mapping. Originally developed and validated on *Drosophila melanogaster*, ProteoCast is applicable to any organism and has been benchmarked against human clinical variants (ClinVar), demonstrating strong performance in distinguishing pathogenic from benign mutations.

**Publication:** Abakarova M, Freiberger MI, Liehrmann A, Rera M*, Laine E*. "Proteome-wide Prediction of the Functional Impact of Missense Variants with ProteoCast." *Nature Communications* (2026). [Link to paper coming soon]

**Key Features:**
- 🔬 Experimentally validated through CRISPR genome editing
- 📊 293+ million variant predictions for *D. melanogaster*
- 🎯 85% accuracy in detecting developmentally lethal mutations
- 🧬 Identifies functional sites in disordered protein regions
- ⚡ Fast and scalable
- 🆓 Fully open-source with Docker deployment

---

**Web Server:** [https://proteocast.ijm.fr/](https://proteocast.ijm.fr/)  
**Drosophila Database:** [https://proteocast.ijm.fr/drosophiladb/](https://proteocast.ijm.fr/drosophiladb/)  
**Docker Image:** [https://hub.docker.com/r/marinaabakarova/proteocast](https://hub.docker.com/r/marinaabakarova/proteocast)


## Usage

### Docker (Recommended)

The easiest way to use ProteoCast is via [https://hub.docker.com/r/marinaabakarova/proteocast](Docker):

```bash
# Pull the latest image
docker pull marinaabakarova/proteocast
```


### Web Server

The [https://proteocast.ijm.fr/](web server) provides an interactive interface for running predictions by only providing a Multiple Sequence Alignment or simply a UniProt code if the reference exists on [https://alphafold.ebi.ac.uk](AlphaFoldDB)


## Key Components of the repository

- **web_server/**: Full-featured Django web application for interactive predictions and result visualization
- **ProteoCast_scripts/**: Standalone Python scripts for batch processing, benchmarking, and research analysis
- **model/**: Complete model implementation including architecture, training procedures, and inference pipelines (extracted from Docker image)
- **browser/**: Integration files for visualizing predictions in genome browsers

## Data and Predictions

### Pre-computed Predictions

Pre-computed predictions for all *Drosophila melanogaster* proteins are available through the web interface or can be downloaded directly:

- **Web Access:** [https://proteocast.ijm.fr/drosophiladb/](https://proteocast.ijm.fr/drosophiladb/)
- **Bulk Download:** Contact us for access to complete datasets

### Input Data Format

ProteoCast accepts protein sequences in FASTA format:

```
>protein_id|organism|description
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQV
GKIPLGTVNNKVYTDVGPKGPPPPKAAQAAAKAKKAAKASVKAAKPKKVKPKPVKKAE
```

### Output Format

Predictions are provided in CSV format with the following columns:

- `protein_id`: Protein identifier
- `position`: Amino acid position (1-indexed)
- `amino_acid`: Amino acid at position
- `ptm_type`: Type of PTM (e.g., phosphorylation, acetylation)
- `score`: Prediction confidence score (0-1)
- `prediction`: Binary prediction (modified/unmodified)

## Model Information

While the complete trained model is distributed via Docker for reproducibility, this repository provides:

- Model architecture implementation
- Training loss functions and procedures
- Data preprocessing pipelines
- Inference scripts

**Note:** For production use, we strongly recommend using the Docker image to ensure consistent predictions and proper dependency management.

## Development

### Running Tests

```bash
# Run all tests
python manage.py test

# Run specific test module
python manage.py test ProteoCast_scripts.tests.test_predictions
```

### Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Code Style

We follow PEP 8 style guidelines. Please run code formatting before submitting:

```bash
# Install formatting tools
pip install black flake8

# Format code
black .

# Check style
flake8 .
```

## Citation

If you use ProteoCast in your research, please cite:

```bibtex
@article{ProteoCast2024,
  author    = {[Authors]},
  title     = {ProteoCast: Deep Learning-based Prediction of Protein Post-translational Modifications},
  journal   = {[Journal Name]},
  year      = {2024},
  volume    = {[Volume]},
  pages     = {[Pages]},
  doi       = {[DOI]}
}
```

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Support

For questions, issues, or feature requests:

- **Issues:** [GitHub Issues](https://github.com/abakarovaMarina/ProteoCast/issues)
- **Email:** abakamarina@gmail.com
- **Web:** [https://proteocast.ijm.fr/](https://proteocast.ijm.fr/)

---

**Last Updated:** February 2026  
**Version:** 1.0.0
