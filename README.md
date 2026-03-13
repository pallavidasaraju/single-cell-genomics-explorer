# Single-Cell Genomics Explorer

A Streamlit-based proof-of-concept demonstrating single-cell RNA-seq data exploration: UMAP clustering, gene expression analysis, and quality control metrics.

Built for the **Bioinformatician** role at **Wellcome Sanger Institute** (Cellular Genomics Informatics).

## Features

- **UMAP Clustering** — Interactive cell type visualization with synthetic scRNA-seq data
- **Gene Expression** — Per-gene UMAP overlay, violin plots by cell type, and marker gene heatmaps
- **Quality Control** — Mitochondrial read filtering, genes-per-cell distribution, doublet score analysis
- **Interactive Filters** — Cell type selection, gene picker, QC threshold sliders

## Tech Stack

- **Python** — Core language
- **Streamlit** — Interactive dashboard framework
- **Pandas / NumPy** — Data processing and analysis
- **Plotly** — Interactive visualizations
- **Synthetic Data** — Realistic scRNA-seq dataset with 8 cell types and 40 marker genes

## Quick Start

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Running Tests

```bash
pytest tests/ -v
```

## Project Structure

```
single-cell-genomics-explorer/
├── app.py              # Streamlit dashboard application
├── sample_data.py      # Synthetic scRNA-seq data generator
├── requirements.txt    # Python dependencies
├── README.md           # This file
└── tests/
    └── test_genomics.py  # Unit tests for data pipeline
```

## Author

**Pallavi Dasaraju** — [GitHub](https://github.com/pallavidasaraju) | [LinkedIn](https://www.linkedin.com/in/pallavi-dasaraju-4442a73b6/)
