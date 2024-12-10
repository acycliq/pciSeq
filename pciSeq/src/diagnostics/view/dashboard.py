"""
View component of the MVC pattern for pciSeq diagnostics.

This module is responsible for presenting data to the user through a Streamlit dashboard.
It maintains the original visualization style while working with the new MVC structure.
"""

import streamlit as st
import pandas as pd
import json
from io import StringIO
import altair as alt
import logging
from pciSeq.src.diagnostics.constants import DiagnosticKeys
from pciSeq.src.diagnostics.model.diagnostic_model import DiagnosticModel

dashboard_logger = logging.getLogger(__name__)


class DiagnosticDashboard:
    def __init__(self, model: DiagnosticModel):
        """
        Initialize dashboard with model connection.

        The dashboard subscribes to Redis channels that the Model publishes to:
        - gene_efficiency: Model publishes gene efficiency metrics
        - cell_type_posterior: Model publishes cell type distribution data

        This pub/sub mechanism enables real-time updates as the Model
        publishes new diagnostic data to Redis.
        """
        self.model = model
        self._setup_page()

        # Subscribe to the Redis channels that the Model publishes to
        self.pubsub = self.model.redis_client.pubsub()

        # These channels match the Model's publishing channels defined in DiagnosticKeys
        # Model -> Redis -> View data flow
        self.pubsub.subscribe(
            DiagnosticKeys.GENE_EFFICIENCY.value,  # Model publishes gene metrics here
            DiagnosticKeys.CELL_TYPE_POSTERIOR.value  # Model publishes cell data here
        )

    def _setup_page(self) -> None:
        """Configure Streamlit page settings."""
        st.set_page_config(
            page_title="Diagnostics: pciSeq",
            page_icon="✅",
            layout="wide",
        )

    @staticmethod
    def _create_barchart(df: pd.DataFrame, nominal_col: str, val_col: str) -> alt.Chart:
        """Create a bar chart using Altair."""
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

    def render(self) -> None:
        """
        Main rendering loop for the dashboard.
        Listens for Model's publications to Redis and updates visualizations.
        """
        title = st.title("Convergence screen. (Waiting for data....)")
        placeholder = st.empty()

        while True:
            with placeholder.container():
                fig_col1, fig_col2 = st.columns(2)

                with fig_col1:
                    self._render_gene_efficiency(title)  # Visualize Model's gene efficiency data

                with fig_col2:
                    self._render_cell_distribution(title)  # Visualize Model's cell type data

            # Check for new messages from Model's publications. This is the heart of the LIVE UPDATE MECHANISM
            message = self.pubsub.get_message()  # Check if Model published new data
            if message and message['type'] == 'message':  # If it's a real message (not control message)
                st.rerun()  # Tell Streamlit to refresh the dashboard

    def _render_gene_efficiency(self, title) -> None:
        """Render gene efficiency visualization."""
        data = self.model.get_diagnostic_data(DiagnosticKeys.GENE_EFFICIENCY)
        if data:
            try:
                df = pd.read_json(StringIO(data['data']))
                metadata = data['metadata']

                title.title("Convergence screen")
                st.markdown(f"#### Gene efficiency after iteration {metadata['iteration']}")

                bar_chart = self._create_barchart(df, nominal_col='gene', val_col='gene_efficiency')
                st.altair_chart(bar_chart, use_container_width=True)
            except json.JSONDecodeError as e:
                dashboard_logger.info('guru meditation....')

    def _render_cell_distribution(self, title) -> None:
        """Render cell type distribution visualization."""
        data = self.model.get_diagnostic_data(DiagnosticKeys.CELL_TYPE_POSTERIOR)
        if data:
            try:
                df = pd.read_json(StringIO(data['data']))
                metadata = data['metadata']

                title.title("Convergence screen")
                st.markdown(f"#### Cell counts per cell class after iteration {metadata['iteration']}")

                bar_chart = self._create_barchart(df, nominal_col='class_name', val_col='counts')
                bar_chart = bar_chart.properties(
                    title=alt.TitleParams(
                        [f'#cells: {df.counts.sum()}'],
                        baseline='bottom',
                        orient='bottom',
                        anchor='end',
                        fontWeight='normal',
                        fontSize=12,
                    ),
                    height=1200
                )
                st.altair_chart(bar_chart, use_container_width=True)
            except json.JSONDecodeError as e:
                dashboard_logger.info('guru meditation....')


def main():
    """Entry point for the Streamlit dashboard."""
    try:
        model = DiagnosticModel()
        dashboard = DiagnosticDashboard(model)
        dashboard.render()
    except Exception as e:
        st.error(f"Dashboard initialization failed: {e}")
        dashboard_logger.error(f"Dashboard error: {e}", exc_info=True)


if __name__ == "__main__":
    main()
