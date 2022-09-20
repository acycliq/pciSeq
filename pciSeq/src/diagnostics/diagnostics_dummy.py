import plotly.express as px  # interactive charts
import pandas as pd
import sqlite3
import os
import time
import streamlit as st  # 🎈 data web app development

title = st.title("Convergence monitor.")
while True:
    try:
        con = sqlite3.connect("file:memdb1?mode=memory&cache=shared")
        df = pd.read_sql_query("SELECT * FROM spots ", con)
        # print(df)

        st.dataframe(df)
        time.sleep(1)
    except RuntimeError as e:
        if str(e) == "Event loop is closed":
            pass
