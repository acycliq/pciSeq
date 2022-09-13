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
st.title("pciSeq: Diagnostics.")

# creating a single-element container
placeholder = st.empty()

source = pd.DataFrame({
    "Price ($)": [10, 15, 20],
    "Month": ["January", "February", "March"]
})
step = 1
while True:
    try:
        sql_str = sql_query("gene_efficiency")
        data = pd.read_sql(sql_str, conn)
        data_2 = pd.read_sql(sql_query("classProb"), conn)
        data_2 = data_2.set_index('cell_label')
        cell_types = [col for col in data_2.columns if data_2[col].dtype in ['float64']]
        data_2 = data_2[cell_types]
        mask = data_2.drop(['Zero'], axis=1).max(axis=1) > 0.02
        data_2 = data_2[mask]
        data_2 = data_2.stack().reset_index()
        data_2.columns = ['cell_label', 'cell_type', 'prob']
        # print(data)
        iter = data.iteration
        assert len(np.unique(iter)) == 1
        i = data.iteration.max()
        if i % step == 0:
            print('iteration: %d' % i)
            with placeholder.container():
                # create two columns for charts
                fig_col1, fig_col2 = st.columns(2)
                with fig_col1:
                    st.markdown("### Chart Num %d" % i)
                    # bar_chart = alt.Chart(source).mark_bar().encode(
                    #     x="sum(Price ($)):Q",
                    #     y=alt.Y("Month:N", sort="-x")
                    # )
                    bar_chart = alt.Chart(data).mark_bar().encode(
                        x=alt.X('gene:N', title='Cell Type'),
                        y='gene_efficiency:Q',
                        color=alt.Color('class:N', legend=None),
                        tooltip=[
                            alt.Tooltip('gene:N', title='Date'),
                            alt.Tooltip('gene_efficiency:Q', title='Max Temp')
                        ]
                    ).properties(width=400, height=550)
                    fig1 = st.altair_chart(bar_chart, use_container_width=True)
                    # fig2 = px.histogram(data_frame=df, x="age_new")
                    # st.write(fig1)

                    heatmap = alt.Chart(
                        data_2,
                        title="2010 Daily High Temperature (F) in Seattle, WA"
                    ).mark_rect().encode(
                        x='cell_label:N',
                        y='cell_type:N',
                        color=alt.Color('prob:Q', scale=alt.Scale(scheme="inferno")),
                        tooltip=[
                            alt.Tooltip('cell_label:N', title='Cell Label'),
                            alt.Tooltip('cell_type:N', title='Cell Type'),
                            alt.Tooltip('prob:Q', title='Prob'),
                        ]
                    ).properties(width=550)
                    fig2 = st.altair_chart(heatmap, use_container_width=True)


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
