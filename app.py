import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def calculate_msc_growth(N0, doubling_time, passages):
    time_per_passage = 24  # hours per passage
    times = np.arange(0, passages * time_per_passage, time_per_passage)
    N = N0 * np.exp((np.log(2) / doubling_time) * times)
    return times, N

def calculate_flask_usage(N0, final_cells, flask_capacity):
    N = N0
    passage = 0
    flask_usage = []
    while N < final_cells and passage < 10:  # Limit to 10 passages to avoid infinite loops
        N *= 2  # Assume one doubling per passage
        passage += 1
        required_flasks = np.ceil(N / flask_capacity)
        flask_usage.append((passage, N, required_flasks))
    return pd.DataFrame(flask_usage, columns=["Passage", "Cell Count", "Flasks Needed"])

def plot_growth_curve(times, N):
    plt.figure(figsize=(8, 5))
    plt.plot(times, N, marker='o', linestyle='-', color='b')
    plt.xlabel("Time (hours)")
    plt.ylabel("Cell Count")
    plt.title("MSC Growth Curve")
    plt.grid(True)
    st.pyplot(plt)

def main():
    st.title("MSC Therapy Calculator")
    
    # User Inputs
    N0 = st.number_input("Initial Cell Count (N0)", min_value=1, value=1000000, step=1000000)
    doubling_time = st.slider("Doubling Time (hours)", min_value=10, max_value=72, value=24)
    passages = st.slider("Number of Passages", min_value=1, max_value=10, value=5)
    final_cells = st.number_input("Final Cell Requirement", min_value=1e6, value=1e9, step=1e7)
    flask_capacity = st.number_input("Flask Capacity (cells)", min_value=1e6, value=1e7, step=1e6)
    
    if st.button("Calculate MSC Growth"):
        times, N = calculate_msc_growth(N0, doubling_time, passages)
        plot_growth_curve(times, N)
    
    if st.button("Calculate Flask Usage"):
        df = calculate_flask_usage(N0, final_cells, flask_capacity)
        st.write(df)

if __name__ == "__main__":
    main()
