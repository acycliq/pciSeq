import plotly.express as px
import numpy as np
from scipy.special import softmax
import plotly.graph_objects as go


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
    gene_name = self.spots.data.iloc[spot_id].gene_name
    x = self.spots.data.iloc[spot_id].x.astype(np.int32).tolist()
    y = self.spots.data.iloc[spot_id].y.astype(np.int32).tolist()
    z = self.spots.data.iloc[spot_id].z.astype(np.int32).tolist()
    n_cells = len(self.spots.parent_cell_id[spot_id]) - 1  # Exclude background
    cell_ids = self.spots.parent_cell_id[spot_id][:-1]
    mvn_loglik = self.spots.mvn_loglik_arr[spot_id][:-1]
    attention = self.spots.attention[spot_id][:-1]
    expr_fluct = self.spots.expr_fluctuations[spot_id][:-1]
    misread = np.log(self.genes.misread_density[gene_name])

    # Calculate scores and probabilities
    scores = mvn_loglik + attention + expr_fluct
    scores = np.append(scores, misread)
    probabilities = softmax(scores)

    # Create labels
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

    spot_to_cell_score_plot(datadict)
    spot_to_cell_prob_plot(datadict)


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
