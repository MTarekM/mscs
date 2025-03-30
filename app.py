import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math

# Constants from reference table
FLASK_DATA = {
    'T-25': {
        'surface_cm2': 25,
        'seeding_cells': 0.7e6,
        'max_cells': 2.8e6,
        'media_ml': 5  # Using upper end of 3-5ml range
    },
    'T-75': {
        'surface_cm2': 75,
        'seeding_cells': 2.1e6,
        'max_cells': 8.4e6,
        'media_ml': 15  # Using upper end of 8-15ml range
    },
    'T-160': {
        'surface_cm2': 160,
        'seeding_cells': 4.6e6,
        'max_cells': 18.4e6,
        'media_ml': 30  # Using upper end of 15-30ml range
    }
}

# Clinical parameters
GVHD_RESPONSE = {
    'Grade I': {'min_dose': 0.5, 'max_dose': 1.0, 'response': [60, 80]},
    'Grade II': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]},
    'Grade III-IV': {'min_dose': 1.5, 'max_dose': 2.0, 'response': [40, 60]}
}

SEPARATOR_YIELD = {
    'Haemonetics': 1.0e6,
    'Spectra Optia': 2.0e6
}

# Culture parameters
SAFETY_FACTOR = 1.2
MAX_PASSAGES = 3
MIN_PASSAGE1_DAYS = 14
ADDITIONAL_PASSAGE_DAYS = 7
PLASMA_PERCENTAGE = 0.15
MAX_PBSC_VOLUME = 50  # ml

def calculate_therapy(weight, desired_dose, separator, flask_type, plasma_priming):
    # Validate inputs
    weight = max(8.0, min(weight, 120.0))
    desired_dose = max(0.5, min(desired_dose, 2.0))
    
    # Get equipment parameters
    flask = FLASK_DATA.get(flask_type, FLASK_DATA['T-75'])
    sep_yield = SEPARATOR_YIELD.get(separator, 1.5e6)
    
    # Calculate required cells
    target_cells = desired_dose * weight * 1e6
    
    # PBSC volume calculation
    pbsc_ml = min(target_cells / sep_yield, MAX_PBSC_VOLUME)
    
    # Flask calculations
    initial_flasks = None
    for N0 in range(1, 201):
        p1 = N0 * flask['max_cells']
        
        if p1 >= target_cells:
            initial_flasks = N0
            p2, p3 = 0, 0
            break
        
        n2 = math.ceil((p1 / flask['seeding_cells']) * SAFETY_FACTOR)
        p2 = n2 * flask['max_cells']
        
        if p2 >= target_cells:
            initial_flasks = N0
            p3 = 0
            break
        
        n3 = math.ceil((p2 / flask['seeding_cells']) * SAFETY_FACTOR)
        p3 = n3 * flask['max_cells']
        
        if p3 >= target_cells:
            initial_flasks = N0
            break
    
    if not initial_flasks:
        initial_flasks = 200
        n2 = math.ceil((initial_flasks * flask['max_cells'] / flask['seeding_cells']) * SAFETY_FACTOR)
        n3 = math.ceil((n2 * flask['max_cells'] / flask['seeding_cells']) * SAFETY_FACTOR)
        p3 = n3 * flask['max_cells']
    
    # Calculate culture duration
    passages = 3 if p3 < target_cells else 2 if p2 < target_cells else 1
    total_days = MIN_PASSAGE1_DAYS + (passages-1)*ADDITIONAL_PASSAGE_DAYS
    
    # Media calculations
    media_changes = 4 + (passages-1)*2  # 4 changes for P1
    total_media = initial_flasks * flask['media_ml'] * media_changes
    
    # Plasma priming
    plasma_vol = initial_flasks * flask['media_ml'] * PLASMA_PERCENTAGE if plasma_priming else 0
    
    return {
        'pbsc_ml': pbsc_ml,
        'initial_flasks': initial_flasks,
        'passages': passages,
        'total_days': total_days,
        'total_media': total_media,
        'plasma_vol': plasma_vol,
        'target_cells': target_cells,
        'final_cells': max(p1, p2, p3)
    }

def plot_growth(results, flask_type):
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        days = [0, 14, 21, 28]
        cells = [0, 
                results['initial_flasks'] * FLASK_DATA[flask_type]['seeding_cells'],
                results['initial_flasks'] * FLASK_DATA[flask_type]['max_cells'],
                results['final_cells']]
        
        ax.plot(days[:results['passages']+1], 
               [x/1e6 for x in cells[:results['passages']+1]], 
               'g-', marker='o', linewidth=2)
        ax.axhline(results['target_cells']/1e6, color='r', linestyle='--')
        ax.set_ylabel('Cells (×10⁶)', fontsize=12)
        ax.set_xlabel('Culture Days', fontsize=12)
        ax.grid(True, alpha=0.3)
        return fig
    except:
        return None

def main():
    st.set_page_config(page_title="Pediatric MSC Calculator", layout="wide")
    st.header("Infant MSC Therapy Calculator (8kg+)")
    
    with st.sidebar:
        weight = st.slider("Patient Weight (kg)", 8.0, 120.0, 8.0)
        grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        dose = st.slider("Target Dose (×10⁶/kg)", 0.5, 2.0, 1.0)
        separator = st.selectbox("Cell Separator", list(SEPARATOR_YIELD.keys()))
        flask_type = st.selectbox("Flask Type", list(FLASK_DATA.keys()))
        plasma = st.checkbox("Use Plasma Priming")
    
    results = calculate_therapy(weight, dose, separator, flask_type, plasma)
    grade_info = GVHD_RESPONSE[grade]
    
    # Display results
    cols = st.columns(4)
    cols[0].metric("PBSC Volume", f"{results['pbsc_ml']:.1f} mL")
    cols[1].metric("Initial Flasks", results['initial_flasks'])
    cols[2].metric("Total Media", f"{results['total_media']} mL")
    cols[3].metric("Culture Days", results['total_days'])
    
    # Plots
    growth_fig = plot_growth(results, flask_type)
    if growth_fig:
        st.pyplot(growth_fig)
    else:
        st.warning("Could not generate growth curve")
    
    # Protocol details
    st.markdown(f"""
    **Clinical Protocol for {grade} GVHD:**
    - Recommended dose: {grade_info['min_dose']}-{grade_info['max_dose']}×10⁶/kg
    - Expected response: {grade_info['response'][0]}-{grade_info['response'][1]}%
    
    **Culture Parameters:**
    - Flask: {flask_type} ({FLASK_DATA[flask_type]['surface_cm2']} cm²)
    - Seeding density: {FLASK_DATA[flask_type]['seeding_cells']/1e6:.1f}×10⁶ cells/flask
    - Max yield: {FLASK_DATA[flask_type]['max_cells']/1e6:.1f}×10⁶ cells/flask
    - Media per flask: {FLASK_DATA[flask_type]['media_ml']} mL
    """)

if __name__ == "__main__":
    main()
