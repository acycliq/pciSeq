import plotly.express as px  # interactive charts
import streamlit as st  # 🎈 data web app development
import numpy as np
import pandas as pd
import time
import pathlib
from pciSeq.src.cell_call.utils import db_connect
import altair as alt
import asyncio
import datetime
from pciSeq.app import run_me

DB_FILE = pathlib.Path(__file__).resolve().parent.parent.parent.joinpath("pciSeq.db").resolve()


conn = db_connect(DB_FILE, remove_if_exists=False)
cell_type_prior = pd.read_sql("select * from cell_type_prior", conn)
cell_types = cell_type_prior[cell_type_prior.iteration==0]['class'].values

gene_efficiency = pd.read_sql("select * from gene_efficiency", conn)

class_prob = pd.read_sql("select * from classProb", conn)
class_prob = class_prob[class_prob.iteration == 0]
class_prob = class_prob.set_index('cell_label')
class_prob = class_prob[cell_types]
class_prob = class_prob.stack().reset_index()
class_prob.columns = ['cell_label', 'cell_type', 'prob']
class_prob = class_prob[class_prob.cell_label > 0]
temp = class_prob.prob.values * 10000
class_prob.prob = temp.astype(np.int)/100
st.set_page_config(
    page_title="Real-Time Data Science Dashboard",
    page_icon="✅",
    layout="wide",
)

def sql_query(table_name):
    sql_str = "select * from %s where iteration = (" \
          "select max(iteration) from %s" \
          ")" % (table_name, table_name)
    return sql_str


df = px.data.tips()
# dashboard title
title = st.title("Convergence monitor.")


# creating a single-element container
placeholder = st.empty()

source = pd.DataFrame({
    "Price ($)": [10, 15, 20],
    "Month": ["January", "February", "March"]
})
step = 3
previous_iteration = -1
# st.header("refreshing every %d iterations" % step)

while True:
    try:
        sql_str = sql_query("gene_efficiency")
        gene_efficiency = pd.read_sql(sql_str, conn)

        sql_str = sql_query("cell_type_prior")
        cell_type_prior = pd.read_sql(sql_str, conn)

        # print(data)
        iter = gene_efficiency.iteration
        assert len(np.unique(iter)) == 1
        i = gene_efficiency.iteration.max()
        if previous_iteration == i:
            print('do nothing')
        elif i % step == 0:
            title.title("Convergence monitor: iteration %d, step: %d" % (i, step))
            print('iteration: %d' % i)
            with placeholder.container():
                # create two columns for charts
                fig_col1, fig_col2 = st.columns(2)
                with fig_col1:
                    st.markdown("### Gene efficiency at iteration %d" % i)
                    bar_chart_1 = alt.Chart(gene_efficiency).mark_bar().encode(
                        y=alt.Y('gene:N', title='Cell Type'),
                        x='gene_efficiency:Q',
                        color=alt.Color('class:N', legend=None),
                        tooltip=[
                            alt.Tooltip('gene:N', title='Date'),
                            alt.Tooltip('gene_efficiency:Q', title='Max Temp')
                        ]
                    ).properties(height=1200)
                    fig1 = st.altair_chart(bar_chart_1, use_container_width=True)
                    # fig2 = px.histogram(data_frame=df, x="age_new")
                    # st.write(fig1)

                with fig_col2:
                    st.markdown("### Cell class weight at iteration %d" % i)
                    bar_chart_2 = alt.Chart(cell_type_prior).mark_bar().encode(
                        y=alt.Y('class:N', title='Cell Type'),
                        x='weight:Q',
                        color=alt.Color('class:N', legend=None),
                        tooltip=[
                            alt.Tooltip('class:N'),
                            alt.Tooltip('weight:Q')
                        ]
                    ).properties(height=1200)
                    fig2 = st.altair_chart(bar_chart_2, use_container_width=True)
            previous_iteration = i


        # wait 1 sec before pingin the db again
        time.sleep(1)

    except Exception as error:
        print(error)



# for seconds in range(200):
#     df['random_bill'] = 100 * np.random.random_sample(df.shape[0])
#     with placeholder.container():
#         # create two columns for charts
#         fig_col1, fig_col2 = st.columns(2)
#         with fig_col1:
#             st.markdown("### Chart Num %d" % seconds)
#             fig1 = px.histogram(df, x="random_bill")
#             # fig2 = px.histogram(data_frame=df, x="age_new")
#             st.write(fig1)
#
#
#     # st.dataframe(df)
#     time.sleep(1)
