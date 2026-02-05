# ProteoCast


### You can find our web server at https://proteocast.ijm.fr/

### And the database for *Drosophila melanogaster* at  https://proteocast.ijm.fr/drosophiladb/ 

### Docker image at https://hub.docker.com/r/marinaabakarova/proteocast

#### License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.

# ProteoCast

This repository provides the implementation and analysis tools for ProteoCast, a deep learning-based method for predicting protein post-translational modifications (PTMs) and their functional impacts.

**Web Server:** [https://proteocast.ijm.fr/](https://proteocast.ijm.fr/)  
**Drosophila Database:** [https://proteocast.ijm.fr/drosophiladb/](https://proteocast.ijm.fr/drosophiladb/)  
**Docker Image:** [https://hub.docker.com/r/marinaabakarova/proteocast](https://hub.docker.com/r/marinaabakarova/proteocast)

## Overview

ProteoCast is a computational framework for predicting and analyzing protein post-translational modifications across multiple species. The repository includes:

* A web server implementation for interactive predictions
* Analysis scripts for PTM prediction evaluation and benchmarking
* Data processing pipelines for proteomics data
* Pre-computed predictions for *Drosophila melanogaster*

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  - [Docker (Recommended)](#docker-recommended)
  - [Web Server](#web-server)
  - [Analysis Scripts](#analysis-scripts)
- [Repository Structure](#repository-structure)
- [Data and Predictions](#data-and-predictions)
- [Citation](#citation)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Installation

### Docker (Recommended)

The easiest way to use ProteoCast is via Docker:

```bash
# Pull the latest image
docker pull marinaabakarova/proteocast:latest

# Run the container
docker run -p 8000:8000 marinaabakarova/proteocast:latest
```

The web interface will be available at `http://localhost:8000`

### Local Installation

For development or custom deployments:

1. **Clone the repository:**

```bash
git clone https://github.com/abakarovaMarina/ProteoCast.git
cd ProteoCast
```

2. **Create a virtual environment:**

```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install dependencies:**

```bash
pip install -r requirements.txt
```

4. **Run migrations (for web server):**

```bash
python manage.py migrate
```

5. **Start the development server:**

```bash
python manage.py runserver
```

## Usage

### Docker (Recommended)

For production deployment or reproducible environments:

```bash
# Build the image locally (optional)
docker build -t proteocast:local .

# Run with custom port
docker run -p 8080:8000 marinaabakarova/proteocast:latest

# Run with volume mounting for data persistence
docker run -v $(pwd)/data:/app/data -p 8000:8000 marinaabakarova/proteocast:latest
```

### Web Server

The web server provides an interactive interface for PTM predictions:

1. Navigate to `http://localhost:8000` (or your deployment URL)
2. Upload protein sequences in FASTA format or enter sequences directly
3. Select the organism and PTM types of interest
4. Submit for prediction and view results

#### API Access

The web server also provides a REST API for programmatic access:

```python
import requests

# Example: Submit a prediction job
url = "http://localhost:8000/api/predict"
data = {
    "sequence": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQV",
    "organism": "drosophila",
    "ptm_types": ["phosphorylation", "acetylation"]
}

response = requests.post(url, json=data)
result = response.json()
```

### Analysis Scripts

The `ProteoCast_scripts/` directory contains tools for data analysis and benchmarking:

#### PTM Prediction Analysis

```bash
python ProteoCast_scripts/analyze_ptm_predictions.py \
    --input predictions.csv \
    --output results/ \
    --threshold 0.5
```

#### Benchmarking

```bash
python ProteoCast_scripts/benchmark.py \
    --predictions predictions.csv \
    --ground_truth validation_set.csv \
    --output benchmark_results/
```

#### Data Collection and Preprocessing

```bash
python ProteoCast_scripts/collect_ptm_data.py \
    --source uniprot \
    --organism drosophila \
    --output processed_data/
```

## Repository Structure

```
ProteoCast/
├── web_server/              # Django-based web server implementation
│   ├── templates/           # HTML templates for web interface
│   ├── static/              # CSS, JavaScript, and assets
│   ├── views.py             # View controllers
│   ├── models.py            # Database models
│   └── urls.py              # URL routing
│
├── ProteoCast_scripts/      # Analysis and data processing scripts
│   ├── ptm_prediction/      # PTM prediction analysis tools
│   ├── benchmarking/        # Benchmarking and validation scripts
│   ├── data_collection/     # Data collection and preprocessing
│   └── visualization/       # Result visualization tools
│
├── model/                   # Model architecture and weights (from Docker)
│   ├── architecture/        # Neural network architecture definitions
│   ├── training/            # Training scripts and configurations
│   └── inference/           # Inference pipeline
│
├── browser/                 # Genome browser integration files
│   └── tracks/              # Browser track configurations
│
├── csv/                     # Sample data and prediction outputs
│   ├── examples/            # Example input files
│   └── predictions/         # Pre-computed predictions
│
├── manage.py                # Django management script
├── requirements.txt         # Python dependencies
├── requirements-min.txt     # Minimal dependencies for core functionality
├── Dockerfile               # Docker container configuration
├── docker-compose.yml       # Docker Compose configuration
└── README.md                # This file
```

### Key Components

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
