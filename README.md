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

**Web Server:** [proteocast.ijm.fr/](https://proteocast.ijm.fr/)  
**Drosophila Database:** [proteocast.ijm.fr/drosophiladb/](https://proteocast.ijm.fr/drosophiladb/)  
**Docker Image:** [marinaabakarova/proteocast](https://hub.docker.com/r/marinaabakarova/proteocast)
**Bulk Download:** [Zenodo repository](https://zenodo.org/records/14871341))


## Usage

### Web Server

The [web server](https://proteocast.ijm.fr/) provides an interactive interface for proteome-wide variant effect predictions. Simply provide:
- A **UniProt accession code** (e.g., P30542) for proteins available in [AlphaFold DB](https://alphafold.ebi.ac.uk/), or
- A **Multiple Sequence Alignment** in FASTA or A3M format for any protein sequence


### Docker

For large-scale predictions, use the [Docker image](https://hub.docker.com/r/marinaabakarova/proteocast):

```bash
# Pull the latest image
docker pull marinaabakarova/proteocast
```


## Key Components of the repository

- **web_server/**: Full-featured Django web application for interactive predictions and result visualization
- **ProteoCast_scripts/**: Standalone Python scripts for batch processing, benchmarking, and research analysis


## Citation

If you use ProteoCast in your research, please cite:

```bibtex
@article{ProteoCast2026,
  author    = {[Abakarova M, Freiberger MI, Liehrmann A, Rera M, Laine E]},
  title     = {Proteome-wide Prediction of the Functional Impact of Missense Variants},
  journal   = {[Nature Communication]},
  year      = {2026},
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
