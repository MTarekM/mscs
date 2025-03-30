import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math

# Constants
SEEDING_DENSITY = {
    'T25': 5000, 'T75': 5000, 'T175': 5000  # cells/cm²
}
GROWTH_RATE = 0.5  # doublings per day
MIN_CULTURE_DAYS = 3
MAX_PASSAGES = 3  # Limit culture to 3 passages
MEDIA_CHANGE_DAYS = [1, 4, 7, 11]

# Flask properties
FLASK_TYPES = {
    'T25': {'surface': 25, 'media': 5, 'max_cells': 2.5e5},
    'T75': {'surface': 75, 'media': 15, 'max_cells': 2.0e6},
    'T175': {'surface': 175, 'media': 30, 'max_cells': 5.0e6}
}

# GVHD Dose Response
GVHD_RESPONSE = {
    'Grade I': {'min_dose': 0.5, 'max_dose': 1.0, 'response': [60, 80]},
    'Grade II': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]},
    'Grade III-IV': {'min_dose': 1.5, 'max_dose': 2.0, 'response': [40, 60]}
}

def calculate_therapy(weight, dose, flask_type, plasma_priming):
    total_cells = dose * weight * 1e6  # Total cells needed
    flask = FLASK_TYPES[flask_type]
    seeding_density = SEEDING_DENSITY[flask_type]
    
    initial_cells_per_flask = flask['surface'] * seeding_density
    expansion_factor = flask['max_cells'] / initial_cells_per_flask
    
    passages = min(MAX_PASSAGES, math.ceil(np.log(total_cells / initial_cells_per_flask) / np.log(expansion_factor)))
    culture_days = passages * 14  # Each passage takes ~14 days
    
    flasks_needed = math.ceil(total_cells / flask['max_cells'])
    total_media = flasks_needed * flask['media'] * len(MEDIA_CHANGE_DAYS)
    plasma_volume = 0.15 * flask['media'] * flasks_needed if plasma_priming else 0
    
    return {
        'flasks': flasks_needed,
        'passages': passages,
        'culture_days': culture_days,
        'total_media': total_media,
        'plasma_volume': plasma_volume
    }

def plot_gvhd_response(dose, gvhd_grade):
    min_dose, max_dose = GVHD_RESPONSE[gvhd_grade]['min_dose'], GVHD_RESPONSE[gvhd_grade]['max_dose']
    response_range = GVHD_RESPONSE[gvhd_grade]['response']
    
    doses = np.linspace(min_dose, max_dose, 100)
    responses = np.interp(doses, [min_dose, max_dose], response_range)
    
    fig, ax = plt.subplots()
    ax.plot(doses, responses, 'g-', label='GVHD Remission Probability')
    ax.axvline(dose, color='r', linestyle='--', label=f'Chosen Dose: {dose}×10⁶/kg')
    ax.set_xlabel('MSC Dose (×10⁶/kg)')
    ax.set_ylabel('Remission Probability (%)')
    ax.legend()
    ax.grid(True)
    return fig

def main():
    st.set_page_config(page_title="MSC Therapy Calculator", layout="wide")
    st.title("MSC Therapy Optimization Tool")
    
    with st.sidebar:
        weight = st.slider("Weight (kg)", 30, 120, 70)
        gvhd_grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        dose = st.slider("MSC Dose (×10⁶/kg)", 0.5, 2.0, 1.0, 0.1)
        flask_type = st.selectbox("Flask Type", list(FLASK_TYPES.keys()))
        plasma_priming = st.checkbox("Use Plasma Priming?")
    
    results = calculate_therapy(weight, dose, flask_type, plasma_priming)
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Flasks Needed", results['flasks'])
        st.metric("Passages", results['passages'])
    with col2:
        st.metric("Culture Days", results['culture_days'])
        st.metric("Total Media (mL)", f"{results['total_media']:.0f}")
    with col3:
        st.metric("Plasma Volume (mL)", f"{results['plasma_volume']:.0f}")
    
    st.pyplot(plot_gvhd_response(dose, gvhd_grade))
    
if __name__ == "__main__":
    main()
