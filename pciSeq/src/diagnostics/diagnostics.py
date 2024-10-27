"""
This module represents the View component in the Model-View-Controller (MVC) pattern.
It is responsible for presenting data to the user through a Streamlit dashboard.

In the MVC pattern:
- Model: Handled by DiagnosticsModel (imported from model.py)
- View: This file (diagnostics.py)
- Controller: Handled by launch_diagnostics.py

The View (this file) focuses on data presentation and user interface. It receives data
from the Model and displays it, but does not handle data processing or application logic.
"""

import pandas as pd
import streamlit as st
import altair as alt
import logging
from pciSeq.src.diagnostics.model import DiagnosticsModel

diagnostics_logger = logging.getLogger(__name__)


def barchart(df: pd.DataFrame, nominal_col: str, val_col: str) -> alt.Chart:
    """
    Create a bar chart using Altair.

    This function is part of the View, responsible for data visualization.

    Args:
        df (pd.DataFrame): The input DataFrame.
        nominal_col (str): The column name for nominal (categorical) data.
        val_col (str): The column name for quantitative data.

    Returns:
        alt.Chart: An Altair chart object representing the bar chart.
    """
    chart = alt.Chart(df).mark_bar().encode(
        y=alt.Y(f'{nominal_col}:N', title=nominal_col),
        x=f'{val_col}:Q',
        color=alt.Color(f'{nominal_col}:N', legend=None),
        tooltip=[
            alt.Tooltip(f'{nominal_col}:N'),
            alt.Tooltip(f'{val_col}:Q')
        ]
    ).properties(height=1200)
    return chart


def main():
    """
    Main function to run the Streamlit dashboard.

    This function sets up the page, initializes the model, and continuously
    updates the dashboard based on incoming Redis messages.

    In the MVC pattern, this function acts as the main entry point for the View,
    coordinating the display of data received from the Model.
    """
    # Set up the Streamlit page configuration
    st.set_page_config(
        page_title="Diagnostics: pciSeq",
        page_icon="✅",
        layout="wide",
    )

    # Initialize dashboard title and placeholder
    title = st.title("Convergence screen. (Waiting for data....)")
    placeholder = st.empty()

    # Initialize the model and subscribe to Redis channels
    # Note: While this initializes the model, the View doesn't process data directly
    model = DiagnosticsModel()
    pubsub = model.subscribe_to_channels()
    diagnostics_logger.info("Subscribed to channels: 'gene_efficiency', 'cell_type_posterior'")

    # Main loop to listen for messages and update the dashboard
    for message in pubsub.listen():
        if message['type'] == 'pmessage':
            channel = message['channel'].decode('UTF-8')
            
            with placeholder.container():
                # Create two columns for charts
                fig_col1, fig_col2 = st.columns(2)
                
                # Update gene efficiency chart
                with fig_col1:
                    if channel == 'gene_efficiency':
                        update_gene_efficiency_chart(model, title)

                # Update cell type posterior chart
                with fig_col2:
                    if channel == 'cell_type_posterior':
                        update_cell_type_posterior_chart(model, title)


def update_gene_efficiency_chart(model: DiagnosticsModel, title: st.delta_generator.DeltaGenerator):
    """
    Update the gene efficiency chart in the Streamlit dashboard.

    This function is part of the View, responsible for updating a specific chart.
    It retrieves data from the Model and updates the UI accordingly.

    Args:
        model (DiagnosticsModel): The data model instance.
        title (st.delta_generator.DeltaGenerator): The Streamlit title object to update.
    """
    gene_efficiency = model.get_gene_efficiency()
    if isinstance(gene_efficiency, pd.DataFrame):
        title.title("Convergence screen")
        st.markdown(f"#### Gene efficiency after iteration {gene_efficiency.iteration.max()}")
        bar_chart_1 = barchart(gene_efficiency, nominal_col='gene', val_col='gene_efficiency')
        st.altair_chart(bar_chart_1, use_container_width=True)


def update_cell_type_posterior_chart(model: DiagnosticsModel, title: st.delta_generator.DeltaGenerator):
    """
    Update the cell type posterior chart in the Streamlit dashboard.

    This function is part of the View, responsible for updating a specific chart.
    It retrieves data from the Model and updates the UI accordingly.

    Args:
        model (DiagnosticsModel): The data model instance.
        title (st.delta_generator.DeltaGenerator): The Streamlit title object to update.
    """
    cell_type_posterior = model.get_cell_type_posterior()
    if isinstance(cell_type_posterior, pd.DataFrame):
        title.title("Convergence screen")
        st.markdown(f"#### Cell counts per cell class after iteration {cell_type_posterior.iteration.max()}")
        bar_chart_2 = barchart(cell_type_posterior, nominal_col='class_name', val_col='counts')
        bar_chart_2 = bar_chart_2.properties(
            title=alt.TitleParams(
                [f'#cells: {cell_type_posterior.counts.sum()}'],
                baseline='bottom',
                orient='bottom',
                anchor='end',
                fontWeight='normal',
                fontSize=12,
            ),
            height=1200
        )
        st.altair_chart(bar_chart_2, use_container_width=True)


if __name__ == "__main__":
    main()
