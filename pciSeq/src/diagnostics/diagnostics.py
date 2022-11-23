import pandas as pd
import streamlit as st  # 🎈 data web app development
from pciSeq.src.cell_call.log_config import logger
from pciSeq.src.diagnostics.utils import redis_db
import altair as alt
import pickle


def barchart(df, nominal_col, val_col):
    chart = alt.Chart(df).mark_bar().encode(
        y=alt.Y('%s:N' % nominal_col, title=nominal_col),
        x='%s:Q' % val_col,
        # color=alt.Color('%s:N' % nominal_col, legend=None),
        tooltip=[
            alt.Tooltip('%s:N' % nominal_col),
            alt.Tooltip('%s:Q' % val_col)
        ]
    ).properties(height=1200)
    return chart


def parse_msg(message):
    gene_efficiency = None
    cell_type_posterior = None
    if message['type'] == 'pmessage' and isinstance(message, dict):
        data = message['data']
        channel = message['channel'].decode('UTF-8')
        if channel == "gene_efficiency":
            gene_efficiency = pickle.loads(data)
        elif channel == "cell_type_posterior":
            cell_type_posterior = pickle.loads(data)
    return gene_efficiency, cell_type_posterior


def main():
    st.set_page_config(
        page_title="Diagnostics: pciSeq",
        page_icon="✅",
        layout="wide",
    )

    # dashboard title
    title = st.title("Convergence monitor.")

    # creating a single-element container
    placeholder = st.empty()

    logger.info('Getting redis_db')
    redis = redis_db(flush=False)
    logger.info('subscribing')
    p = redis.redis_client.pubsub()
    p.psubscribe('gene_efficiency')
    p.psubscribe('cell_type_posterior')
    logger.info('subscribed to two channels')
    for message in p.listen():
        gene_efficiency, cell_type_posterior = parse_msg(message)
        with placeholder.container():
            # create two columns for charts
            fig_col1, fig_col2 = st.columns(2)
            with fig_col1:
                if isinstance(gene_efficiency, pd.DataFrame):
                    st.markdown("### Gene efficiency after iteration %d" % gene_efficiency.iteration.max())
                    bar_chart_1 = barchart(gene_efficiency, nominal_col='gene', val_col='gene_efficiency')
                    fig1 = st.altair_chart(bar_chart_1, use_container_width=True)

            with fig_col2:
                if isinstance(cell_type_posterior, pd.DataFrame):
                    st.markdown("### Posterior cell class weight after iteration %d" % cell_type_posterior.iteration.max())
                    bar_chart_2 = barchart(cell_type_posterior, nominal_col='class_name', val_col='prob')
                    fig2 = st.altair_chart(bar_chart_2, use_container_width=True)


if __name__ == "__main__":
    logger.info('in diagnostics_listen, main()')
    main()
