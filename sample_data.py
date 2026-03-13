"""Synthetic single-cell RNA-seq dataset generator for the POC dashboard."""
import numpy as np
import pandas as pd

CELL_TYPES = [
    "T cell", "B cell", "NK cell", "Monocyte", "Dendritic cell",
    "Macrophage", "Epithelial", "Fibroblast"
]

MARKER_GENES = {
    "T cell": ["CD3D", "CD3E", "CD4", "CD8A", "TRAC"],
    "B cell": ["CD19", "MS4A1", "CD79A", "CD79B", "BANK1"],
    "NK cell": ["NKG7", "GNLY", "KLRD1", "NCAM1", "PRF1"],
    "Monocyte": ["CD14", "LYZ", "FCGR3A", "S100A8", "S100A9"],
    "Dendritic cell": ["FCER1A", "CLEC10A", "CD1C", "HLA-DRA", "ITGAX"],
    "Macrophage": ["CD68", "CSF1R", "MSR1", "MRC1", "MARCO"],
    "Epithelial": ["EPCAM", "KRT18", "KRT19", "CDH1", "MUC1"],
    "Fibroblast": ["COL1A1", "COL3A1", "VIM", "DCN", "THY1"],
}

ALL_GENES = sorted(set(g for genes in MARKER_GENES.values() for g in genes))


def generate_expression_matrix(n_cells=2000, seed=42):
    """Generate a synthetic single-cell expression matrix with known cell type structure."""
    rng = np.random.default_rng(seed)
    n_genes = len(ALL_GENES)
    cells_per_type = n_cells // len(CELL_TYPES)

    cell_labels = []
    expression = np.zeros((n_cells, n_genes))

    for i, ct in enumerate(CELL_TYPES):
        start = i * cells_per_type
        end = start + cells_per_type if i < len(CELL_TYPES) - 1 else n_cells
        n = end - start
        cell_labels.extend([ct] * n)

        # Background expression
        expression[start:end, :] = rng.exponential(0.3, size=(n, n_genes))

        # Marker gene upregulation
        for gene in MARKER_GENES[ct]:
            if gene in ALL_GENES:
                gene_idx = ALL_GENES.index(gene)
                expression[start:end, gene_idx] += rng.normal(4.0, 0.8, size=n)

    # Clip negatives
    expression = np.clip(expression, 0, None)

    df = pd.DataFrame(expression, columns=ALL_GENES)
    df.insert(0, "cell_id", [f"cell_{i:04d}" for i in range(n_cells)])
    df.insert(1, "cell_type", cell_labels)
    return df


def generate_qc_metrics(n_cells=2000, seed=42):
    """Generate quality control metrics per cell."""
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "cell_id": [f"cell_{i:04d}" for i in range(n_cells)],
        "n_genes_detected": rng.integers(500, 5000, size=n_cells),
        "total_counts": rng.integers(1000, 50000, size=n_cells),
        "pct_mitochondrial": np.clip(rng.exponential(3.0, size=n_cells), 0, 30).round(2),
        "pct_ribosomal": np.clip(rng.normal(15, 5, size=n_cells), 0, 40).round(2),
        "doublet_score": np.clip(rng.beta(2, 10, size=n_cells), 0, 1).round(3),
    })


def generate_umap_coordinates(expression_df, seed=42):
    """Generate synthetic UMAP coordinates that cluster by cell type."""
    rng = np.random.default_rng(seed)
    n = len(expression_df)
    centers = {
        "T cell": (2, 5), "B cell": (-4, 3), "NK cell": (5, 2),
        "Monocyte": (-2, -3), "Dendritic cell": (0, -5),
        "Macrophage": (-5, -1), "Epithelial": (4, -4), "Fibroblast": (-1, 4),
    }
    umap1, umap2 = [], []
    for _, row in expression_df.iterrows():
        cx, cy = centers[row["cell_type"]]
        umap1.append(cx + rng.normal(0, 0.8))
        umap2.append(cy + rng.normal(0, 0.8))
    return pd.DataFrame({
        "cell_id": expression_df["cell_id"],
        "cell_type": expression_df["cell_type"],
        "UMAP_1": umap1,
        "UMAP_2": umap2,
    })
