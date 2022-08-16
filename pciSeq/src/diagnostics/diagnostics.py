import plotly.express as px  # interactive charts
import streamlit as st  # 🎈 data web app development

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

with placeholder.container():
    # create two columns for charts
    fig_col1, fig_col2 = st.columns(2)
    with fig_col1:
        st.markdown("### Second Chart")
        fig1 = px.histogram(df, x="total_bill")
        # fig2 = px.histogram(data_frame=df, x="age_new")
        st.write(fig1)


st.dataframe(df)


