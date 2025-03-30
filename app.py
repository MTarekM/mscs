import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math

# EXACT VALUES FROM REFERENCE TABLE
FLASK_TYPES = {
    'T-25': {
        'surface': 25,       # cm²
        'seeding': 0.7e6,    # cells
        'confluent': 2.8e6,  # cells
        'media': 5           # mL
    },
    'T-75': {
        'surface': 75,
        'seeding': 2.1e6,
        'confluent': 8.4e6,
        'media': 15
    },
    'T-160': {
        'surface': 160,
        'seeding': 4.6e6,
        'confluent': 18.4e6,
        'media': 30
    }
}

GVHD_RESPONSE = {
    'Grade I': {'min_dose': 0.5, 'max_dose': 1.0, 'response': [60, 80]},
    'Grade II': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]},
    'Grade III-IV': {'min_dose': 1.5, 'max_dose': 2.0, 'response': [40, 60]}
}

SEPARATOR_YIELD = {
    'Haemonetics': 1.0e6,
    'Spectra Optia': 2.0e6
}

# Constants
SAFETY_FACTOR = 1.2
MAX_PASSAGES = 3
PASSAGE1_DAYS = 14
PASSAGE_DAYS = 7
MAX_PBSC = 50  # mL
PLASMA_PERCENT = 0.15

def calculate_therapy(weight, dose, separator, flask_type, plasma_priming):
    target = dose * weight * 1e6
    pbsc = min(target / SEPARATOR_YIELD[separator], MAX_PBSC)
    
    # Flask calculations
    flasks = []
    n_flasks = math.ceil(target / FLASK_TYPES[flask_type]['confluent'])
    
    for passage in range(MAX_PASSAGES + 1):
        if passage == 0:
            flasks.append({
                'passage': 0,
                'count': n_flasks,
                'yield': n_flasks * FLASK_TYPES[flask_type]['confluent'],
                'days': 0
            })
        else:
            prev_yield = flasks[-1]['yield']
            if prev_yield >= target:
                break
                
            n_flasks = math.ceil((prev_yield / FLASK_TYPES[flask_type]['seeding']) * SAFETY_FACTOR)
            flasks.append({
                'passage': passage,
                'count': n_flasks,
                'yield': n_flasks * FLASK_TYPES[flask_type]['confluent'],
                'days': PASSAGE1_DAYS if passage == 1 else PASSAGE_DAYS
            })
    
    total_days = sum(f['days'] for f in flasks[1:])
    media_changes = 4 + (len(flasks)-2)*2
    total_media = flasks[0]['count'] * FLASK_TYPES[flask_type]['media'] * media_changes
    plasma = flasks[0]['count'] * FLASK_TYPES[flask_type]['media'] * PLASMA_PERCENT if plasma_priming else 0
    
    return {
        'pbsc': pbsc,
        'flasks': flasks,
        'total_days': total_days,
        'total_media': total_media,
        'plasma': plasma,
        'target': target,
        'final_yield': flasks[-1]['yield']
    }

def plot_growth(results, flask_type):
    fig, ax = plt.subplots(figsize=(10, 6))
    days = [0]
    cells = [0]
    
    for i, f in enumerate(results['flasks'][1:]):
        days.append(days[-1] + f['days'])
        cells.append(f['yield'])
    
    ax.plot(days, [x/1e6 for x in cells], 'go-', linewidth=2)
    ax.axhline(results['target']/1e6, color='r', linestyle='--')
    ax.set_ylabel('Cells (×10⁶)', fontsize=12)
    ax.set_xlabel('Culture Days', fontsize=12)
    ax.grid(True, alpha=0.3)
    return fig

def plot_gvhd_probability(grade, dose):
    data = GVHD_RESPONSE[grade]
    min_d = data['min_dose']
    max_d = data['max_dose']
    resp_min, resp_max = data['response']
    
    x = np.linspace(0.5, 2.5, 100)
    opt = (min_d + max_d)/2
    y = resp_min + (resp_max - resp_min) * np.exp(-((x - opt)/0.3)**2)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x, y, 'b-', linewidth=2)
    ax.axvline(dose, color='r', linestyle='--')
    ax.set_ylim(0, 100)
    ax.set_xlabel('Dose (×10⁶ cells/kg)', fontsize=12)
    ax.set_ylabel('Remission Probability (%)', fontsize=12)
    ax.grid(True, alpha=0.3)
    return fig

def main():
    st.set_page_config(page_title="MSC Therapy Calculator", layout="wide")
    st.title("MSC Therapy Calculator (8kg+)")
    
    with st.sidebar:
        st.header("Patient Parameters")
        weight = st.slider("Weight (kg)", 8.0, 120.0, 8.0, 0.1)
        grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        dose = st.slider("Dose (×10⁶/kg)", 0.5, 2.0, 1.0, 0.1)
        
        st.header("Lab Parameters")
        separator = st.selectbox("Cell Separator", list(SEPARATOR_YIELD.keys()))
        flask_type = st.selectbox("Flask Type", list(FLASK_TYPES.keys()))
        plasma = st.checkbox("Plasma Priming")
    
    # Calculations
    results = calculate_therapy(weight, dose, separator, flask_type, plasma)
    grade_data = GVHD_RESPONSE[grade]
    
    # Results display
    cols = st.columns(4)
    cols[0].metric("PBSC Volume", f"{results['pbsc']:.1f} mL")
    cols[1].metric("Initial Flasks", results['flasks'][0]['count'])
    cols[2].metric("Total Media", f"{results['total_media']:.0f} mL")
    cols[3].metric("Culture Days", results['total_days'])
    
    # Plots
    col1, col2 = st.columns(2)
    with col1:
        st.header("Expansion Progress")
        st.pyplot(plot_growth(results, flask_type))
    
    with col2:
        st.header("GVHD Remission Probability")
        st.pyplot(plot_gvhd_probability(grade, dose))
    
    # Protocol details
    st.header("Protocol Summary")
    st.markdown(f"""
    **For {grade} GVHD:**
    - Recommended dose: {grade_data['min_dose']}-{grade_data['max_dose']}×10⁶/kg
    - Expected response: {grade_data['response'][0]}–{grade_data['response'][1]}%
    
    **Culture Parameters:**
    - Flask: {flask_type} ({FLASK_TYPES[flask_type]['surface']} cm²)
    - Seeding: {FLASK_TYPES[flask_type]['seeding']/1e6:.1f}×10⁶ cells/flask
    - Confluent yield: {FLASK_TYPES[flask_type]['confluent']/1e6:.1f}×10⁶ cells/flask
    
    **Final Yield:** {results['final_yield']/1e6:.1f}×10⁶ cells
    """)

if __name__ == "__main__":
    main()
