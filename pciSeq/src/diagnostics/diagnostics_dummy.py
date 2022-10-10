import plotly.express as px  # interactive charts
import pandas as pd
import sqlite3
import os
from pciSeq.src.cell_call.utils import get_db_tables
import time
import streamlit as st  # 🎈 data web app development

title = st.title("Convergence monitor (dummy).")
while True:
    try:
        # con = sqlite3.connect("file:memdb1?mode=memory&cache=shared")
        con = sqlite3.connect(r"D:\Home\Dimitris\OneDrive - University College London\dev\Python\pciSeq\pciSeq\pciSeq.db")
        # con = sqlite3.connect('my_db.db')
        # tables = {"spots"}
        # db_tables = set(get_db_tables(con))
        # if tables.issubset(db_tables):
        df = pd.read_sql_query("SELECT * FROM geneCount", con)
        st.dataframe(df)
        time.sleep(1)
    except RuntimeError as e:
        if str(e) == "Event loop is closed":
            pass
