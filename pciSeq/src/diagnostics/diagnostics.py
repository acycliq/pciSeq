import streamlit as st  # 🎈 data web app development
import numpy as np
import pandas as pd
import time
from pciSeq.src.cell_call.log_config import logger
# from pciSeq.src.cell_call.utils import db_connect, get_db_tables
import pciSeq.src.diagnostics.config as diagnostics_cfg
from pciSeq.src.diagnostics.utils import redis_db
import altair as alt
import sys
import os


# DB_FILE = pathlib.Path(__file__).resolve().parent.parent.parent.joinpath("pciSeq.db").resolve()
# DB_FILE = r"D:\Home\Dimitris\OneDrive - University College London\dev\Python\pciSeq\pciSeq\pciSeq.db"
# DB_FILE = "file:memdb1?mode=memory&cache=shared"
# DB_FILE = diagnostics_cfg.SETTINGS['DB_URL']

conn = None
checked_tables = False

st.set_page_config(
    page_title="Diagnostics: pciSeq",
    page_icon="✅",
    layout="wide",
)

# def sql_query(table_name):
#     sql_str = "select * from %s where iteration = (" \
#           "select max(iteration) from %s" \
#           ")" % (table_name, table_name)
#     return sql_str


# dashboard title
title = st.title("Convergence monitor.")


# creating a single-element container
placeholder = st.empty()

step = 1
previous_iteration = -1


# def query_tables(tables, con):
#     db_tables = get_db_tables(con)
#     if set(tables) != set(db_tables):
#         sys.exit(1)


redis = redis_db()
redis.attach_client()
while True:
    try:
        # print(redis.from_redis('foo'))
        # print(pd.DataFrame(redis.from_redis("gene_efficiency")))
        # if conn is None:
        #     assert os.path.isfile(DB_FILE)
        #     conn = db_connect(DB_FILE, remove_if_exists=False)
        #     logger.info('Connection to %s established' % DB_FILE)
        # if not checked_tables:
        #     set_tables = {"gene_efficiency", "cell_type_prior", "cell_type_posterior"}
        #     set_db = set(get_db_tables(conn))
        #     assert set_tables.issubset(set_db)
        #     checked_tables = True

        # sql_str = sql_query("gene_efficiency")
        gene_efficiency = redis.from_redis("gene_efficiency")

        # sql_str = sql_query("cell_type_prior")
        # cell_type_prior = pd.read_sql(sql_str, conn)
        cell_type_prior = redis.from_redis("cell_type_prior")

        # sql_str = sql_query("cell_type_posterior")
        # cell_type_posterior = pd.read_sql(sql_str, conn)
        cell_type_posterior = redis.from_redis("cell_type_posterior")

        # print(data)
        iter = gene_efficiency.iteration
        assert len(np.unique(iter)) == 1
        i = gene_efficiency.iteration.max()
        if previous_iteration == i:
            pass
            # logger.info('do nothing')
        elif i % step == 0:
            title.title("Convergence monitor")
            # logger.info('iteration: %d' % i)
            with placeholder.container():
                # create two columns for charts
                fig_col1, fig_col2 = st.columns(2)
                with fig_col1:
                    st.markdown("### Gene efficiency after iteration %d" % i)
                    bar_chart_1 = alt.Chart(gene_efficiency).mark_bar().encode(
                        y=alt.Y('gene:N', title='Gene'),
                        x='gene_efficiency:Q',
                        color=alt.Color('class:N', legend=None),
                        tooltip=[
                            alt.Tooltip('gene:N', title='gene'),
                            alt.Tooltip('gene_efficiency:Q', title='efficiency')
                        ]
                    ).properties(height=1200)
                    fig1 = st.altair_chart(bar_chart_1, use_container_width=True)
                    # fig2 = px.histogram(data_frame=df, x="age_new")
                    # st.write(fig1)

                with fig_col2:
                    st.markdown("### Posterior cell class weight after iteration %d" % i)
                    bar_chart_2 = alt.Chart(cell_type_posterior).mark_bar().encode(
                        y=alt.Y('class_name:N', title='Cell Type'),
                        x='prob:Q',
                        color=alt.Color('class_name:N', legend=None),
                        tooltip=[
                            alt.Tooltip('class_name:N'),
                            alt.Tooltip('prob:Q')
                        ]
                    ).properties(height=1200)
                    fig2 = st.altair_chart(bar_chart_2, use_container_width=True)

                # with fig_col3:
                #     st.markdown("### Prior cell class weight at iteration %d" % i)
                #     bar_chart_3 = alt.Chart(cell_type_prior).mark_bar().encode(
                #         y=alt.Y('class:N', title='Cell Type'),
                #         x='weight:Q',
                #         color=alt.Color('class:N', legend=None),
                #         tooltip=[
                #             alt.Tooltip('class:N'),
                #             alt.Tooltip('weight:Q')
                #         ]
                #     ).properties(height=1200)
                #     fig3 = st.altair_chart(bar_chart_3, use_container_width=True)
            previous_iteration = i

        # wait 1 sec before pingin the db again
        time.sleep(1)

    except AssertionError:
        logger.info('..waiting...')
        time.sleep(1)
        pass
    except Exception as error:
        logger.info(error)

