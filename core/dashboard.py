# core/dashboard.py
import streamlit as st
import numpy as np
import time
import os
import re
import pandas as pd

st.set_page_config(page_title="Euclid Sentinel | Control Room", layout="wide")

# --- CSS STYLING FOR TERMINAL LOOK ---
st.markdown("""
    <style>
    .reportview-container { background: #0e1117; }
    .stCodeBlock { border: 1px solid #4a4a4a; }
    </style>
    """, unsafe_allow_html=True)

st.title("🛰️ Euclid Sentinel: Mimetic-Conformal Dashboard")

# --- SIDEBAR: CONTROL & CONFIG ---
st.sidebar.header("🕹️ Mission Control")
target = st.sidebar.selectbox("Active Target", 
                              ["Bullet Cluster", "Abell 370", "El Gordo", "HUDF"])

st.sidebar.divider()
st.sidebar.markdown(f"**Engine Version:** V3.1.3")
st.sidebar.markdown(f"**Sentinel Constant ($\kappa$):** 0.80")

# Logic to check if engine is actually writing to the log
log_file = "sentinel_engine.log"
engine_status = "ACTIVE" if os.path.exists(log_file) and (time.time() - os.path.getmtime(log_file) < 30) else "IDLE"
st.sidebar.markdown(f"**Engine Status:** `{engine_status}`")

# --- MAIN LAYOUT: VISUALIZATION ---
col1, col2 = st.columns(2)

with col1:
    st.subheader("📡 Input Source")
    # In a future update, we can extract a thumbnail from the FITS here
    st.image("https://via.placeholder.com/600x400/1a1a1a/ffffff?text=FITS+Stream+Ready", use_column_width=True)

with col2:
    st.subheader("🔮 Emergent Potential")
    if os.path.exists("Sentinel_FINAL_REPORT.png"):
        st.image("Sentinel_FINAL_REPORT.png", caption=f"Last Validation: {target}")
    else:
        st.info("Awaiting Analysis Locked signal...")

st.divider()

# --- LIVE TELEMETRY & STABILITY ---
col_log, col_graph = st.columns([1, 1])

def parse_logs(filename):
    """Extracts stability metrics for real-time graphing."""
    ratios = []
    if os.path.exists(filename):
        with open(filename, "r") as f:
            lines = f.readlines()[-50:] # Look at last 50 entries
            for line in lines:
                match = re.search(r"Ratio=(\d+\.\d+)x", line)
                if match:
                    ratios.append(float(match.group(1)))
    return ratios if ratios else [0.0]

with col_log:
    st.subheader("📜 Live Engine Logs")
    log_display = st.empty()

with col_graph:
    st.subheader("📈 DM Ratio Stability (Universality Tracker)")
    stability_chart = st.empty()

# --- REFRESH LOOP ---
while True:
    # 1. Update Logs
    if os.path.exists(log_file):
        with open(log_file, "r") as f:
            raw_logs = "".join(f.readlines()[-15:])
            log_display.code(raw_logs, language="bash")
    
    # 2. Update Stability Graph
    ratio_data = parse_logs(log_file)
    stability_chart.line_chart(pd.DataFrame(ratio_data, columns=["DM Ratio"]))
    
    time.sleep(1) # Frequency: 1Hz