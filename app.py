"""Single-Cell Genomics Explorer — Interactive Streamlit dashboard for scRNA-seq data."""
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from sample_data import (
    generate_expression_matrix, generate_qc_metrics,
    generate_umap_coordinates, ALL_GENES, CELL_TYPES, MARKER_GENES,
)

st.set_page_config(page_title="Single-Cell Genomics Explorer", layout="wide")
st.title("Single-Cell Genomics Explorer")
st.caption("Interactive scRNA-seq data exploration — gene expression, clustering, and QC metrics")

# Load data
@st.cache_data
def load_data():
    expr = generate_expression_matrix()
    qc = generate_qc_metrics()
    umap = generate_umap_coordinates(expr)
    return expr, qc, umap

expr_df, qc_df, umap_df = load_data()

# Sidebar filters
st.sidebar.header("Filters")
selected_types = st.sidebar.multiselect(
    "Cell Types", CELL_TYPES, default=CELL_TYPES
)
selected_gene = st.sidebar.selectbox("Gene of Interest", ALL_GENES, index=ALL_GENES.index("CD3D"))
mito_threshold = st.sidebar.slider("Max % Mitochondrial", 0.0, 30.0, 15.0, 0.5)

# Filter data
mask = expr_df["cell_type"].isin(selected_types)
filtered_expr = expr_df[mask]
filtered_umap = umap_df[mask]
filtered_qc = qc_df[qc_df["cell_id"].isin(filtered_expr["cell_id"])]
filtered_qc_clean = filtered_qc[filtered_qc["pct_mitochondrial"] <= mito_threshold]

# KPIs
col1, col2, col3, col4 = st.columns(4)
col1.metric("Total Cells", f"{len(filtered_expr):,}")
col2.metric("After QC Filter", f"{len(filtered_qc_clean):,}")
col3.metric("Cell Types", len(selected_types))
col4.metric("Genes Profiled", len(ALL_GENES))

# Tabs
tab1, tab2, tab3 = st.tabs(["UMAP Clustering", "Gene Expression", "Quality Control"])

with tab1:
    st.subheader("UMAP Cell Type Clustering")
    fig_umap = px.scatter(
        filtered_umap, x="UMAP_1", y="UMAP_2", color="cell_type",
        hover_data=["cell_id"],
        color_discrete_sequence=px.colors.qualitative.Set2,
        title="Cell Type Clusters (UMAP Projection)"
    )
    fig_umap.update_traces(marker=dict(size=4, opacity=0.7))
    fig_umap.update_layout(height=550, legend_title="Cell Type")
    st.plotly_chart(fig_umap, use_container_width=True)

    st.subheader("Cell Type Distribution")
    counts = filtered_expr["cell_type"].value_counts().reset_index()
    counts.columns = ["Cell Type", "Count"]
    fig_bar = px.bar(counts, x="Cell Type", y="Count", color="Cell Type",
                     color_discrete_sequence=px.colors.qualitative.Set2)
    fig_bar.update_layout(height=350, showlegend=False)
    st.plotly_chart(fig_bar, use_container_width=True)

with tab2:
    st.subheader(f"Gene Expression: {selected_gene}")

    col_a, col_b = st.columns(2)
    with col_a:
        # UMAP colored by gene expression
        gene_vals = filtered_expr[selected_gene].values
        plot_df = filtered_umap.copy()
        plot_df["expression"] = gene_vals
        fig_gene = px.scatter(
            plot_df, x="UMAP_1", y="UMAP_2", color="expression",
            color_continuous_scale="Viridis",
            title=f"{selected_gene} Expression on UMAP"
        )
        fig_gene.update_traces(marker=dict(size=4, opacity=0.7))
        fig_gene.update_layout(height=450)
        st.plotly_chart(fig_gene, use_container_width=True)

    with col_b:
        # Violin plot of expression by cell type
        violin_df = filtered_expr[["cell_type", selected_gene]].copy()
        fig_violin = px.violin(
            violin_df, x="cell_type", y=selected_gene, color="cell_type",
            box=True, points=False,
            color_discrete_sequence=px.colors.qualitative.Set2,
            title=f"{selected_gene} Distribution by Cell Type"
        )
        fig_violin.update_layout(height=450, showlegend=False, xaxis_title="Cell Type")
        st.plotly_chart(fig_violin, use_container_width=True)

    # Marker gene heatmap
    st.subheader("Marker Gene Heatmap")
    mean_expr = filtered_expr.groupby("cell_type")[ALL_GENES].mean()
    fig_heat = px.imshow(
        mean_expr, aspect="auto", color_continuous_scale="RdYlBu_r",
        title="Mean Expression per Cell Type"
    )
    fig_heat.update_layout(height=400)
    st.plotly_chart(fig_heat, use_container_width=True)

with tab3:
    st.subheader("Quality Control Metrics")

    col_q1, col_q2 = st.columns(2)
    with col_q1:
        fig_mito = px.histogram(
            filtered_qc, x="pct_mitochondrial", nbins=50,
            title="% Mitochondrial Reads Distribution"
        )
        fig_mito.add_vline(x=mito_threshold, line_dash="dash", line_color="red",
                           annotation_text=f"Threshold: {mito_threshold}%")
        fig_mito.update_layout(height=350)
        st.plotly_chart(fig_mito, use_container_width=True)

    with col_q2:
        fig_genes = px.histogram(
            filtered_qc, x="n_genes_detected", nbins=50,
            title="Genes Detected per Cell"
        )
        fig_genes.update_layout(height=350)
        st.plotly_chart(fig_genes, use_container_width=True)

    fig_scatter_qc = px.scatter(
        filtered_qc, x="total_counts", y="n_genes_detected",
        color="pct_mitochondrial", color_continuous_scale="Reds",
        title="Total Counts vs Genes Detected (colored by % Mito)",
        opacity=0.5
    )
    fig_scatter_qc.update_traces(marker=dict(size=3))
    fig_scatter_qc.update_layout(height=400)
    st.plotly_chart(fig_scatter_qc, use_container_width=True)

    # Doublet scores
    fig_doublet = px.histogram(
        filtered_qc, x="doublet_score", nbins=50,
        title="Doublet Score Distribution"
    )
    fig_doublet.update_layout(height=300)
    st.plotly_chart(fig_doublet, use_container_width=True)
