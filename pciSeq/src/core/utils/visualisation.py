import plotly.express as px


def heatmap_counts_per_class(obj):
    data = obj.cells.mean_gene_reads_per_class()
    gene_names = obj.genes.gene_panel
    class_names = obj.cells.class_names

    fig = px.imshow(
        data,
        labels=dict(x="Class", y="Gene", color="Mean Reads"),
        x=class_names,
        y=gene_names,
        color_continuous_scale='inferno',
        aspect="auto"
    )

    fig.update_traces(
        hovertemplate="Gene: %{y}<br>Class: %{x}<br>Mean Reads: %{z:.2f}<extra></extra>",
        colorbar=dict(
            len=0.5,  # Colorbar is 50% of the plot height.
            y=0.5,  # Centered vertically.
            yanchor="middle"  # Anchor the center at y=0.5.
        )
    )

    # Build a list of shapes for grid lines.
    shapes = []

    # Add vertical lines between each class.
    for i in range(1, len(class_names)):
        shapes.append({
            "type": "line",
            "x0": i - 0.5,
            "y0": -0.5,
            "x1": i - 0.5,
            "y1": data.shape[0] - 0.5,
            "line": {"color": "white", "width": 1},
            "opacity": 0.3,  # Set the opacity here, not inside the "line" dict.
            "xref": "x",
            "yref": "y"
        })

    # Add horizontal lines between each gene.
    for j in range(1, len(gene_names)):
        shapes.append({
            "type": "line",
            "x0": -0.5,
            "y0": j - 0.5,
            "x1": data.shape[1] - 0.5,
            "y1": j - 0.5,
            "line": {"color": "white", "width": 1},
            "opacity": 0.3,  # Set the opacity as a top-level attribute.
            "xref": "x",
            "yref": "y"
        })

    fig.update_layout(
        title="pciSeq: Mean Gene Reads per Class",
        height=2000,
        xaxis_nticks=len(class_names),
        yaxis_nticks=150,
        xaxis_title="Cell Classes",
        yaxis_title="Genes",
        title_x=0.5,
        shapes=shapes
    )

    fig.update_xaxes(
        tickmode='array',
        tickvals=list(range(len(class_names))),
        ticktext=class_names,
        tickangle=90
    )

    # Update the color axis settings for the colorbar
    fig.update_layout(
        coloraxis_colorbar=dict(
            len=0.5,  # Set to 25% of the plot height (i.e. half of current)
            y=0.5,  # Center vertically
            yanchor="middle"  # Anchor at the middle
        )
    )

    fig.show()
