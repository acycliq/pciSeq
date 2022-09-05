import plotly.express as px  # interactive charts
import streamlit as st  # 🎈 data web app development
import numpy as np
import time
from pciSeq.app import run_me

st.set_page_config(
    page_title="Real-Time Data Science Dashboard",
    page_icon="✅",
    layout="wide",
)

df = px.data.tips()
# dashboard title
st.title("pciSeq: Diagnostics.")

# creating a single-element container
placeholder = st.empty()

for seconds in range(200):
    df['random_bill'] = 100 * np.random.random_sample(df.shape[0])
    with placeholder.container():
        # create two columns for charts
        fig_col1, fig_col2 = st.columns(2)
        with fig_col1:
            st.markdown("### Chart Num %d" % seconds)
            fig1 = px.histogram(df, x="random_bill")
            # fig2 = px.histogram(data_frame=df, x="age_new")
            st.write(fig1)


    # st.dataframe(df)
    time.sleep(1)


