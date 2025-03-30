import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math

# Constants from reference tables
FLASK_TYPES = {
    'T-25': {'seeding': 0.7e6, 'confluent': 2.8e6, 'media': 5},
    'T-75': {'seeding': 2.1e6, 'confluent': 8.4e6, 'media': 15},
    'T-160': {'seeding': 4.6e6, 'confluent': 18.4e6, 'media': 30}
}

GVHD_RESPONSE = {
    'Grade I': {'min_dose': 0.5, 'max_dose': 1.0, 'response': [60, 80]},
    'Grade II': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]},
    'Grade III-IV': {'min_dose': 1.5, 'max_dose': 2.0, 'response': [40, 60]}
}

SEPARATOR_YIELD = {'Haemonetics': 1e6, 'Spectra Optia': 2e6}

# Protocol constants
SAFETY_FACTOR = 1.2
MAX_PASSAGES = 3
PASSAGE1_DAYS = 14
PASSAGE_DAYS = 7
PLASMA_PERCENT = 0.15
MAX_PBSC = 50  # mL

def calculate_therapy(weight, dose, separator, flask_type, plasma_priming, media_freq):
    # Validate inputs
    weight = max(8.0, min(weight, 120.0))
    dose = max(0.5, min(dose, 2.0))
    
    # Calculate target cells and PBSC volume
    target_cells = dose * weight * 1e6
    pbsc_ml = min(target_cells / SEPARATOR_YIELD[separator], MAX_PBSC)
    
    # Initialize passage tracking
    passages = []
    total_days = 0
    total_media = 0
    
    # Initial seeding (Passage 0)
    initial_flasks = math.ceil(target_cells / FLASK_TYPES[flask_type]['confluent'])
    passages.append({
        'passage_num': 0,
        'flasks': initial_flasks,
        'input': initial_flasks * FLASK_TYPES[flask_type]['seeding'],
        'output': initial_flasks * FLASK_TYPES[flask_type]['confluent'],
        'days': 0,
        'media_changes': 0
    })
    
    # Subsequent passages
    for passage_num in range(1, MAX_PASSAGES + 1):
        if passages[-1]['output'] >= target_cells:
            break
            
        current_flasks = math.ceil(
            (passages[-1]['output'] / FLASK_TYPES[flask_type]['seeding']) * SAFETY_FACTOR
        )
        
        days = PASSAGE1_DAYS if passage_num == 1 else PASSAGE_DAYS
        media_changes = (days // media_freq) + 1  # +1 for initial media
        
        passages.append({
            'passage_num': passage_num,
            'flasks': current_flasks,
            'input': current_flasks * FLASK_TYPES[flask_type]['seeding'],
            'output': current_flasks * FLASK_TYPES[flask_type]['confluent'],
            'days': days,
            'media_changes': media_changes
        })
        
        total_days += days
        total_media += current_flasks * FLASK_TYPES[flask_type]['media'] * media_changes
    
    plasma_vol = passages[0]['flasks'] * FLASK_TYPES[flask_type]['media'] * PLASMA_PERCENT if plasma_priming else 0
    
    return {
        'pbsc_ml': pbsc_ml,
        'passages': passages,
        'total_days': total_days,
        'total_media': total_media,
        'plasma_vol': plasma_vol,
        'target_cells': target_cells,
        'final_yield': passages[-1]['output']
    }

def main():
    st.set_page_config(page_title="MSC Therapy Calculator", layout="wide")
    st.title("MSC Therapy Calculator")
    
    with st.sidebar:
        st.header("Patient Parameters")
        weight = st.slider("Weight (kg)", 8.0, 120.0, 70.0, 0.1, format="%.1f")
        grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        dose = st.slider("Dose (×10⁶/kg)", 0.5, 2.0, 1.0, 0.1)
        
        st.header("Lab Parameters")
        separator = st.selectbox("Cell Separator", list(SEPARATOR_YIELD.keys()))
        flask_type = st.selectbox("Flask Type", list(FLASK_TYPES.keys()))
        plasma_priming = st.checkbox("Plasma Priming (15% of initial media)")
        media_freq = st.slider("Media Change Frequency (days)", 1, 7, 2)
    
    results = calculate_therapy(weight, dose, separator, flask_type, plasma_priming, media_freq)
    
    st.header("Therapy Parameters")
    cols = st.columns(4)
    cols[0].metric("PBSC Volume", f"{results['pbsc_ml']:.1f} mL")
    cols[1].metric("Initial Flasks", results['passages'][0]['flasks'])
    cols[2].metric("Total Media", f"{results['total_media']} mL")
    cols[3].metric("Culture Days", results['total_days'])
    
    if plasma_priming:
        st.info(f"**Plasma Volume Needed:** {results['plasma_vol']:.1f} mL")
    
    st.header("GVHD Remission Probability")
    grade_data = GVHD_RESPONSE[grade]
    st.markdown(f"**For {grade} GVHD:** Recommended dose: {grade_data['min_dose']}-{grade_data['max_dose']}×10⁶/kg, Expected response: {grade_data['response'][0]}–{grade_data['response'][1]}%")
    
if __name__ == "__main__":
    main()
