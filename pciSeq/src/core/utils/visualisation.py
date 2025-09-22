import plotly.express as px
import pandas as pd
import numpy as np
from scipy.special import softmax
import plotly.graph_objects as go
from ..utils.io_utils import get_out_dir
import os


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


def check_spot(self, spot_id):
    """
    Analyze a spot by creating visualization charts and returning score/probability arrays.

    Parameters:
    spot_id (int): The ID of the spot to analyze

    Returns:
    tuple: (scores_array, probabilities_array)s
    """
    # Get data for the specified spot
    # First find the row position of the spot_id
    row_pos = self.spots.data.index.get_loc(spot_id)

    gene_name = self.spots.data.iloc[row_pos].gene_name # I could have used loc[spot_id] here too
    x = self.spots.data.iloc[row_pos].x.astype(np.int32).tolist()
    y = self.spots.data.iloc[row_pos].y.astype(np.int32).tolist()
    z = self.spots.data.iloc[row_pos].z.astype(np.int32).tolist()
    n_cells = len(self.spots.parent_cell_id[row_pos]) - 1  # Exclude background
    cell_ids = self.spots.parent_cell_id[row_pos][:-1]
    mvn_loglik = self.spots.mvn_loglik_arr[row_pos][:-1]
    attention = self.spots.attention[row_pos][:-1]
    expr_fluct = self.spots.expr_fluctuations[row_pos][:-1]
    misread = np.log(self.genes.misread_density[gene_name])

    # Calculate scores and probabilities
    scores = mvn_loglik + attention + expr_fluct
    scores = np.append(scores, misread)
    probabilities = softmax(scores)

    # Create labels. If the segmentation has been relabelled, map the labels back to the original ones.
    if self.config['label_map']:
        reverse_map = {v:k for k, v in self.config['label_map'].items()}
        cell_ids = [reverse_map[d] for d in cell_ids]

    labels = [f'Cell {cid}' for cid in cell_ids] + ['Misread']

    datadict = {
        'spot_id': spot_id,
        'gene_name': gene_name,
        'x': x,  # Already converted to list of int32
        'y': y,  # (same as above)
        'z': z,  # (same as above)
        'n_cells': n_cells,
        'cell_ids': cell_ids,
        'mvn_loglik': mvn_loglik,
        'attention': attention,
        'expr_fluct': expr_fluct,
        'misread': float(misread),  # Convert numpy float to native Python float
        'score': scores,
        'prob': probabilities,
        'labels': labels
    }

    df = pd.DataFrame({
        'Name': labels[:-1],
        # 'internal_tag':self.spots.parent_cell_id[row_pos][:-1],
        'mvn_loglik': mvn_loglik,
        'attention': attention,
        'expr_fluct': expr_fluct}).set_index(['Name'])
    df['sum'] = df[['mvn_loglik', 'attention', 'expr_fluct']].sum(axis=1)
    df.loc['misread'] = [np.nan, np.nan, np.nan, misread]

    spot_to_cell_score_plot(datadict)
    spot_to_cell_prob_plot(datadict)
    return df


def spot_to_cell_prob_plot(data):
    gene_name = data['gene_name']
    spot_id = data['spot_id']
    x = data['x']
    y = data['y']
    z = data['z']
    n_cells = data['n_cells']
    cell_ids = data['cell_ids']
    mvn_loglik = data['mvn_loglik']
    attention = data['attention']
    expr_fluct = data['expr_fluct']
    misread = data['misread']

    # Calculate scores and probabilities
    scores = mvn_loglik + attention + expr_fluct
    scores = np.append(scores, misread)
    prob = softmax(scores)

    # Labels (cells + misread)
    labels = [f'Cell {cid}' for cid in cell_ids] + ['Misread']

    # Create bar chart with consistent styling
    fig = go.Figure()

    # Bar trace (blue for cells, red for misread)
    fig.add_trace(go.Bar(
        x=labels,
        y=prob,
        marker_color=['#1f77b4'] * len(labels[:-1]) + ['#d62728'],  # Blue for cells, red for misread
        width=0.7,
        hovertemplate="<b>%{x}</b><br>Probability: %{y:.4f}<extra></extra>",
        marker_line=dict(width=0)  # No border on bars
    ))

    # Layout (matches previous plot exactly)
    fig.update_layout(
        title={
            'text': f'Spot {spot_id} ({gene_name}, x={x}, y={y}, z={z}) - Assignment Probabilities',
            'y': 0.95,
            'x': 0.02,
            'xanchor': 'left',
            'yanchor': 'top',
            'pad': {'b': 30}  # Add padding below the title (adjust as needed)
        },
        xaxis=dict(
            title='Candidate Assignment',
            tickangle=45,
            tickfont=dict(size=12),
            showline=True,
            linecolor='black',
            gridcolor='rgba(0,0,0,0.05)'
        ),
        yaxis=dict(
            title='Probability',
            gridcolor='rgba(0,0,0,0.1)',
            showline=True,
            linecolor='black',
            range=[0, min(1.1, max(prob) * 1.2)]  # Add headroom
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=50, r=50, b=80, t=70),
        height=500,
        width=max(800, len(labels) * 100),
        hoverlabel=dict(
            bgcolor="white",
            font_size=12,
            font_family="Arial"
        ),
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02
        )
    )

    # Highlight max probability
    max_prob_idx = np.argmax(prob)
    fig.add_annotation(
        x=labels[max_prob_idx],
        y=prob[max_prob_idx] + 0.02,
        text="Most probable",
        showarrow=True,
        arrowhead=1,
        font=dict(size=12)
    )

    fig.show()


def spot_to_cell_score_plot(my_dict):
    gene_name = my_dict['gene_name']
    my_spot = my_dict['spot_id']
    x = my_dict['x']
    y = my_dict['y']
    z = my_dict['z']
    n_cells = my_dict['n_cells']
    cell_ids = my_dict['cell_ids']
    mvn_loglik = my_dict['mvn_loglik']
    attention = my_dict['attention']
    expr_fluct = my_dict['expr_fluct']
    misread = my_dict['misread']

    # Labels
    labels = [f'Cell {cid}' for cid in cell_ids] + ['Misread']

    # Create figure with Matplotlib-like aesthetics
    fig = go.Figure()

    # Add stacked bars (with individual hover for each component)
    fig.add_trace(go.Bar(
        x=labels[:-1],
        y=mvn_loglik,
        name='MVN Log-Likelihood',
        marker_color='#1f77b4',
        hovertemplate="<b>%{x}</b><br>MVN: %{y:.2f}<extra></extra>",
        width=0.7  # Matplotlib-like bar width
    ))

    fig.add_trace(go.Bar(
        x=labels[:-1],
        y=attention,
        name='Attention',
        marker_color='#ff7f0e',
        hovertemplate="<b>%{x}</b><br>Attention: %{y:.2f}<extra></extra>",
        width=0.7
    ))

    fig.add_trace(go.Bar(
        x=labels[:-1],
        y=expr_fluct,
        name='Expr Fluctuations',
        marker_color='#2ca02c',
        hovertemplate="<b>%{x}</b><br>Expr Fluct: %{y:.2f}<extra></extra>",
        width=0.7
    ))

    # Misread bar (standalone)
    fig.add_trace(go.Bar(
        x=[labels[-1]],
        y=[misread],
        name='Misread Density (log)',
        marker_color='#d62728',
        hovertemplate="<b>Misread</b><br>Value: %{y:.2f}<extra></extra>",
        width=0.7
    ))

    # Update layout to mimic Matplotlib
    fig.update_layout(
        title={
            'text': f'Spot {my_spot} {gene_name} - Score Decomposition<br><span style="font-size:12px; color:gray">The higher the better</span>',
            'y': 0.95,
            'x': 0.02,
            'xanchor': 'left',
            'yanchor': 'top',
            'pad': {'b': 30}  # Add padding below the title (adjust as needed)
        },
        yaxis_title='Log-Likelihood Score',
        barmode='stack',
        hovermode='closest',  # Tooltip shows only the hovered segment
        plot_bgcolor='white',
        font=dict(size=12),
        margin=dict(l=50, r=50, b=100, t=60),
        height=500,
        width=max(800, n_cells * 100),  # Dynamic width
        xaxis=dict(
            tickangle=45,
            tickfont=dict(size=12),
            title_standoff=25
        ),
        yaxis=dict(
            gridcolor='rgba(0,0,0,0.1)',
            showline=True,
            linecolor='black'
        ),
        legend=dict(
            orientation="v",  # Vertical layout
            yanchor="top",  # Anchor to top of legend
            y=1,  # Position at top of plot area
            xanchor="left",  # Anchor to left of legend
            x=1.02,  # Push right (into whitespace)
            bgcolor="rgba(255,255,255,0.8)",  # Optional: semi-transparent white
            bordercolor="rgba(0,0,0,0.2)",  # Optional: subtle border
            borderwidth=1
        )
    )

    fig.show()


def make_trellis_enh3(df, highlight_label=None):
    """
    Create a trellis plot with:
    - Different colors for each gene
    - Optional circles highlighting specific neighbors

    Parameters:
    - df: DataFrame containing spatial gene data
    - highlight_label: (optional) Draw circles around markers containing this neighbor ID,
                      with circle size proportional to the neighbor's probability
    """

    # Sort data by plane_id and gene_name for consistent coloring
    df = df.sort_values(['plane_id', 'gene_name'])

    # Calculate grid dimensions
    n_planes = df['plane_id'].nunique()
    n_cols = 6  # Number of columns in grid
    subplot_size = 250  # Size for each subplot

    # Create base scatter plot with color by gene_name
    fig = px.scatter(
        df,
        x="x",
        y="y",
        color="gene_name",
        facet_col="plane_id",
        facet_col_wrap=n_cols,
        hover_data={
            "gene_name": True,
            "x": ":.2f",
            "y": ":.2f",
            "z": ":.2f",
            "neighbour_array": True,
            "neighbour_prob": True,
            "plane_id": False
        },
        height=max(600, (n_planes // n_cols + 1) * subplot_size),
        width=n_cols * subplot_size,
        category_orders={"gene_name": sorted(df['gene_name'].unique())}  # Consistent color mapping
    )

    # Add highlighting circles if a label is specified
    if highlight_label is not None:
        shapes = []
        for plane_id in df['plane_id'].unique():
            plane_df = df[df['plane_id'] == plane_id]
            for _, row in plane_df.iterrows():
                if highlight_label in row['neighbour_array']:
                    # Get the probability for this label
                    idx = row['neighbour_array'].index(highlight_label)
                    prob = row['neighbour_prob'][idx]

                    # Calculate circle radius proportional to probability (range 10-50 units)
                    radius = 10 + 40 * prob

                    shapes.append({
                        'type': 'circle',
                        'xref': f'x{plane_id + 1}',
                        'yref': f'y{plane_id + 1}',
                        'x0': row['x'] - radius,
                        'y0': row['y'] - radius,
                        'x1': row['x'] + radius,
                        'y1': row['y'] + radius,
                        'line': {
                            'color': 'black' if prob > 0.5 else 'darkgray',  # High prob = black border
                            'width': 1 + 2 * prob,
                            'dash': 'dot' if prob < 0.3 else 'solid'  # Low prob = dotted line
                        },
                        'opacity': 0.7
                    })

        fig.update_layout(shapes=shapes)

    # Custom hover template
    hover_template = (
        "<b>%{customdata[0]}</b><br>"
        "Coord: (%{x:.2f}, %{y:.2f}, %{customdata[1]:.2f})<br>"
        "Neighbors: %{customdata[2]}<br>"
        "Probs: %{customdata[3]}"
    )

    # Add highlight info to hover if specified
    if highlight_label is not None:
        hover_template = (
                "<b>%{customdata[0]}</b><br>"
                "Coord_xyz: (%{x:.2f}, %{y:.2f}, %{customdata[1]:.2f})<br>"
                f"Neighbor {highlight_label} prob: " +
                ("%.3f<br>" % df.loc[df.index, 'neighbour_prob'].apply(
                    lambda probs, arr=df.loc[df.index, 'neighbour_array']:
                    probs[arr.index(highlight_label)] if highlight_label in arr else 'N/A'
                )) +
                "All neighbors: %{customdata[2]}<br>"
                "All probs: %{customdata[3]}"
        )

    # Visual enhancements
    fig.update_traces(
        marker=dict(
            size=8,  # Slightly larger for better color visibility
            opacity=0.9,
            line=dict(width=1, color='black')  # Dark outline for contrast
        ),
        hovertemplate=hover_template + "<extra></extra>"
    )

    # Layout adjustments
    fig.update_layout(
        margin=dict(l=5, r=5, t=25, b=5),
        grid=dict(rows=None, columns=n_cols, xgap=0.01, ygap=0.01),
        # plot_bgcolor='white',
        # paper_bgcolor='white',

        legend=dict(
            title_text='Gene',
            orientation='v',
            yanchor='middle',
            y=0.5,  # Centered vertically
            xanchor='left',
            x=1.02,  # Closer to plot
            entrywidthmode='pixels',
            entrywidth=40  # Fixed width for alignment
        )
    )

    # Equal axes and clean annotations
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    for annotation in fig.layout.annotations:
        annotation.text = f"Plane {annotation.text.split('=')[1]}"
        annotation.font.size = 9

    fig.show()




def cell_class_stacked_bar(obj, class_col='top_class'):
    """
    Create a stacked bar chart where:
    - X-axis: Integer Z values
    - Y-axis: Count of cells
    - Stacking: Classes ordered by count (highest at bottom, lowest at top)
    - Zero class is colored black
    - Tooltips show individual counts (not cumulative)
    """

    from .geometry import anisotropy_calc

    centroids = obj.cells.centroid.values
    voxel_size = obj.config['voxel_size']
    data = anisotropy_calc(centroids, voxel_size, inverse=True)
    plane_id = data[:, -1].astype(int)

    idx = np.argmax(obj.cells.classProb, axis=1)
    cell_class = obj.cells.class_names[idx]

    df = pd.DataFrame({'plane_id': plane_id,
                       'top_class': cell_class})

    # Get counts for each class at each Z level
    z_class_counts = {}
    for z in sorted(df['plane_id'].unique()):
        z_data = df[df['plane_id'] == z]
        class_counts = z_data[class_col].value_counts().to_dict()
        z_class_counts[z] = class_counts

    # Get all unique classes
    all_classes = set()
    for counts in z_class_counts.values():
        all_classes.update(counts.keys())
    all_classes = list(all_classes)

    # Sort classes globally by total frequency for consistent colors
    global_totals = {}
    for cls in all_classes:
        global_totals[cls] = sum(z_counts.get(cls, 0) for z_counts in z_class_counts.values())

    classes_by_total = sorted(all_classes, key=lambda x: global_totals[x], reverse=True)

    # Create the plot
    fig = go.Figure()
    z_values = sorted(z_class_counts.keys())

    # Calculate cumulative heights for proper stacking
    cumulative_data = {z: {} for z in z_values}

    for z in z_values:
        # Sort classes by count for this Z value (highest first)
        z_counts = z_class_counts[z]
        sorted_classes = sorted(z_counts.items(), key=lambda x: x[1], reverse=True)

        cumulative = 0
        for cls, count in sorted_classes:
            cumulative_data[z][cls] = {
                'bottom': cumulative,
                'height': count
            }
            cumulative += count

    # Add traces in the order that ensures proper stacking
    for cls in classes_by_total:
        y_values = []
        base_values = []
        hover_texts = []

        for z in z_values:
            if cls in cumulative_data[z]:
                height = cumulative_data[z][cls]['height']
                bottom = cumulative_data[z][cls]['bottom']
            else:
                height = 0
                bottom = 0

            y_values.append(height)
            base_values.append(bottom)

            # Custom hover text showing individual count
            hover_texts.append(f"Z-Plane: {z}<br>Class: {cls}<br>Count: {height}")

        if any(y > 0 for y in y_values):
            # Set color for Zero class to black
            color = 'black' if cls == 'Zero' else None

            fig.add_trace(go.Bar(
                x=z_values,
                y=y_values,
                base=base_values,
                name=cls,
                offsetgroup=1,
                marker_color=color,
                hovertemplate='%{hovertext}<extra></extra>',
                hovertext=hover_texts
            ))

    fig.update_layout(
        title='Cell Class Distribution by Z-Plane',
        xaxis_title='Z-Plane (Integer)',
        yaxis_title='Cell Count',
        barmode='group',
        showlegend=True,
        height = 700
    )

    return fig

# def trellis_plot(self, label, flatfile_folder):
#
#     cellBoundaries_tsv = os.path.join(flatfile_folder, 'cellBoundaries.tsv')
#     cell_boundaries = self.read_tsv(cellBoundaries_tsv)
#     target_cell = cell_boundaries[cell_boundaries.cell_id == label]
#
#     coords = target_cell.coords.squeeze()
#     min_x = min(x for x, y in coords)
#     min_y = min(y for x, y in coords)
#     max_x = max(x for x, y in coords)
#     max_y = max(y for x, y in coords)
#
#     bbox = [min_x, min_y, max_x, max_y]
#
#     # geneData = self.read_tsv('/tmp/pciSeq/data/geneData.tsv')
#     geneData_tsv = os.path.join(flatfile_folder, 'geneData.tsv')
#     geneData = self.read_tsv(geneData_tsv)
#     mask = (
#             (geneData['x'] >= bbox[0]) &  # x >= x_min
#             (geneData['x'] <= bbox[2]) &  # x <= x_max
#             (geneData['y'] >= bbox[1]) &  # y >= y_min
#             (geneData['y'] <= bbox[3])  # y <= y_max
#     )
#     df = geneData[mask]
#     make_trellis_enh3(df)
