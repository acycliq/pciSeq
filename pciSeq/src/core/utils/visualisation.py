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
    Show the spot-to-cell score breakdown for one spot.

    Draws the score-decomposition and assignment-probability charts, then returns the
    breakdown as a table.

    Parameters:
    spot_id (int): The ID of the spot to analyze

    Returns:
    pd.DataFrame: One row per candidate cell (plus a background row), with the score
        terms, the misread value, and their sum.
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
    cell_inefficiency = self.spots.cell_inefficiency[row_pos][:-1]
    gene_inefficiency = self.spots.gene_inefficiency[row_pos][:-1]
    gene_idx = np.where(self.genes.gene_panel == gene_name)[0][0]
    misread = self.genes.log_rho_bar[gene_idx]
    # the inside-cell bonus the model adds before the softmax in spots_to_cell. it is
    # nonzero only for the cell whose boundary the spot sits in, and zero for background.
    bonus = self.spots.bonus_mask[row_pos][:-1] * self.config['InsideCellBonus']

    # Calculate scores and probabilities
    scores = mvn_loglik + attention + expr_fluct + cell_inefficiency + gene_inefficiency + bonus
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
        'cell_inefficiency': cell_inefficiency,
        'gene_inefficiency': gene_inefficiency,
        'bonus': bonus,
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
        'expr_fluct': expr_fluct,
        'cell_inefficiency': cell_inefficiency,
        'gene_inefficiency': gene_inefficiency,
        'bonus': bonus}).set_index(['Name'])
    df['misread'] = np.nan
    df['sum'] = df[['mvn_loglik', 'attention', 'expr_fluct', 'cell_inefficiency', 'gene_inefficiency', 'bonus']].sum(axis=1)
    df.loc['background'] = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, misread, misread]

    spot_to_cell_score_plot(datadict)
    spot_to_cell_prob_plot(datadict)
    return df


def spot_to_cell_prob_plot(data):
    gene_name = data['gene_name']
    spot_id = data['spot_id']
    x = data['x']
    y = data['y']
    z = data['z']
    cell_ids = data['cell_ids']
    # check_spot already computed this (with the inside-cell bonus); reuse it so the two
    # never drift apart.
    prob = data['prob']

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
    cell_inefficiency = my_dict['cell_inefficiency']
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

    fig.add_trace(go.Bar(
        x=labels[:-1],
        y=cell_inefficiency,
        name='Cell Inefficiency',
        marker_color='#9467bd',
        hovertemplate="<b>%{x}</b><br>Cell Inefficiency"
                      ": %{y:.2f}<extra></extra>",
        width=0.7
    ))

    fig.add_trace(go.Bar(
        x=labels[:-1],
        y=my_dict['gene_inefficiency'],
        name='Gene Inefficiency',
        marker_color='#8c564b',
        hovertemplate="<b>%{x}</b><br>Gene Inefficiency"
                      ": %{y:.2f}<extra></extra>",
        width=0.7
    ))

    fig.add_trace(go.Bar(
        x=labels[:-1],
        y=my_dict['bonus'],
        name='Inside-cell Bonus',
        marker_color='#e377c2',
        hovertemplate="<b>%{x}</b><br>Inside-cell Bonus"
                      ": %{y:.2f}<extra></extra>",
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

