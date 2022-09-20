import plotly.express as px  # interactive charts
import pandas as pd
import sqlite3
import os
import time
import streamlit as st  # 🎈 data web app development

title = st.title("Convergence monitor.")
while True:
    if os.path.exists('my_db.db'):
        con = sqlite3.connect(r'my_db.db')
        df = pd.read_sql_query("SELECT * FROM spots ", con)
        print(df)

        st.dataframe(df)
        time.sleep(1)
