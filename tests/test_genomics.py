"""Unit tests for the single-cell genomics data pipeline."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from sample_data import (
    generate_expression_matrix, generate_qc_metrics,
    generate_umap_coordinates, ALL_GENES, CELL_TYPES, MARKER_GENES,
)


def test_expression_matrix_shape():
    df = generate_expression_matrix(n_cells=200)
    assert df.shape[0] == 200
    assert "cell_id" in df.columns
    assert "cell_type" in df.columns
    assert all(g in df.columns for g in ALL_GENES)


def test_expression_values_non_negative():
    df = generate_expression_matrix(n_cells=200)
    numeric = df.select_dtypes(include=[np.number])
    assert (numeric >= 0).all().all()


def test_cell_types_present():
    df = generate_expression_matrix(n_cells=400)
    found_types = set(df["cell_type"].unique())
    assert found_types == set(CELL_TYPES)


def test_marker_gene_upregulation():
    df = generate_expression_matrix(n_cells=1000)
    for ct in CELL_TYPES:
        ct_mask = df["cell_type"] == ct
        other_mask = ~ct_mask
        for gene in MARKER_GENES[ct][:2]:
            if gene in ALL_GENES:
                ct_mean = df.loc[ct_mask, gene].mean()
                other_mean = df.loc[other_mask, gene].mean()
                assert ct_mean > other_mean, f"{gene} should be upregulated in {ct}"


def test_qc_metrics_columns():
    qc = generate_qc_metrics(n_cells=100)
    expected = ["cell_id", "n_genes_detected", "total_counts",
                "pct_mitochondrial", "pct_ribosomal", "doublet_score"]
    assert all(c in qc.columns for c in expected)
    assert len(qc) == 100


def test_qc_metrics_ranges():
    qc = generate_qc_metrics(n_cells=500)
    assert (qc["pct_mitochondrial"] >= 0).all()
    assert (qc["pct_mitochondrial"] <= 30).all()
    assert (qc["doublet_score"] >= 0).all()
    assert (qc["doublet_score"] <= 1).all()
    assert (qc["n_genes_detected"] > 0).all()
    assert (qc["total_counts"] > 0).all()


def test_umap_coordinates():
    expr = generate_expression_matrix(n_cells=200)
    umap = generate_umap_coordinates(expr)
    assert len(umap) == 200
    assert "UMAP_1" in umap.columns
    assert "UMAP_2" in umap.columns
    assert "cell_type" in umap.columns


def test_umap_clusters_separated():
    expr = generate_expression_matrix(n_cells=800)
    umap = generate_umap_coordinates(expr)
    centroids = umap.groupby("cell_type")[["UMAP_1", "UMAP_2"]].mean()
    # Check that at least some cluster centroids are well-separated
    for i, ct1 in enumerate(centroids.index):
        for ct2 in centroids.index[i+1:]:
            dist = np.sqrt(
                (centroids.loc[ct1, "UMAP_1"] - centroids.loc[ct2, "UMAP_1"])**2 +
                (centroids.loc[ct1, "UMAP_2"] - centroids.loc[ct2, "UMAP_2"])**2
            )
            assert dist > 1.0, f"{ct1} and {ct2} clusters too close"
