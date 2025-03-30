import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

# Constants
SEEDING_DENSITY = 5000  # cells/cm²
MAX_CONFLUENCY = 80  # %
HARVEST_CONFLUENCY = 80  # %
GROWTH_RATE = 0.5  # doublings per day

FLASK_TYPES = {
    'T25': {'surface': 25, 'media_min': 5, 'media_max': 10},
    'T75': {'surface': 75, 'media_min': 15, 'media_max': 20},
    'T175': {'surface': 175, 'media_min': 30, 'media_max': 40}
}

GVHD_RESPONSE = {
    'Grade I': {'min_dose': 0.5, 'max_dose': 1.0, 'response': [60, 80]},
    'Grade II': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]},
    'Grade III-IV': {'min_dose': 1.5, 'max_dose': 2.0, 'response': [40, 60]}
}

def calculate_msc_therapy(weight, desired_dose, flask_type, media_volume, media_freq):
    # Calculate total MSCs needed (×10⁶ cells)
    total_cells = desired_dose * weight  # ×10⁶ cells
    
    # Calculate required PBSC volume (adjusts with dose)
    pbsc_volume = max(0.05, 0.05 * (desired_dose / 1.0))  # Minimum 50mL, scales with dose
    
    # Calculate flasks needed based on seeding density and growth
    flask_area = FLASK_TYPES[flask_type]['surface']
    initial_cells = flask_area * SEEDING_DENSITY / 1e6  # ×10⁶ cells
    
    # Calculate passages needed
    passages = 0
    remaining_cells = total_cells
    while remaining_cells > initial_cells:
        passages += 1
        remaining_cells /= 5  # Assume 5-fold expansion per passage
    
    # Calculate culture days
    days_per_passage = np.log(HARVEST_CONFLUENCY/SEEDING_DENSITY) / np.log(2) / GROWTH_RATE
    culture_days = int(np.ceil(days_per_passage * (passages + 1)))
    
    # Calculate total media needed
    media_per_change = media_volume
    total_media = media_per_change * (culture_days // media_freq) * (passages + 1)
    
    return {
        'pbsc_volume': pbsc_volume,
        'flasks': passages + 1,
        'passages': passages,
        'culture_days': culture_days,
        'total_media': total_media,
        'initial_cells': initial_cells,
        'target_cells': total_cells
    }

def plot_growth_curve(initial_cells, target_cells, days):
    x = np.linspace(0, days, 100)
    plateau = min(3, days-2)  # Time to reach plateau
    growth = initial_cells * np.exp(GROWTH_RATE * x)
    growth[x > plateau] = growth[int(plateau/days*100)]  # Plateau effect
    
    fig, ax = plt.subplots()
    ax.plot(x, growth, 'g-', linewidth=2)
    ax.axhline(y=target_cells, color='r', linestyle='--', label='Target Dose')
    ax.set_xlabel('Culture Days')
    ax.set_ylabel('MSC Count (×10⁶)')
    ax.set_title('MSC Growth Curve to Target Dose')
    ax.legend()
    ax.grid(True, alpha=0.3)
    return fig

def plot_dose_response(grade, desired_dose):
    grade_data = GVHD_RESPONSE[grade]
    doses = np.linspace(0.5, 2.0, 100)
    response = np.interp(doses, 
                       [grade_data['min_dose'], grade_data['max_dose']],
                       grade_data['response'])
    
    fig, ax = plt.subplots()
    ax.plot(doses, response, 'b-', linewidth=2)
    ax.axvline(x=desired_dose, color='r', linestyle='--', label='Desired Dose')
    ax.set_xlabel('Dose (×10⁶ MSCs/kg)')
    ax.set_ylabel('Response Probability (%)')
    ax.set_title(f'GVHD {grade} Dose-Response Curve')
    ax.legend()
    ax.grid(True, alpha=0.3)
    return fig

def main():
    st.set_page_config(page_title="MSC Therapy Calculator", layout="wide")
    st.title("Precision MSC Therapy Calculator for GVHD")
    
    # Inputs
    with st.sidebar:
        st.header("Patient Parameters")
        weight = st.number_input("Weight (kg)", min_value=30, max_value=120, value=70)
        gvhd_grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        desired_dose = st.slider("Desired Dose (×10⁶ MSCs/kg)", 
                               min_value=0.5, max_value=2.0, value=1.0, step=0.1)
        
        st.header("Culture Parameters")
        flask_type = st.selectbox("Flask Type", list(FLASK_TYPES.keys()))
        media_volume = st.slider(f"Media Volume per {flask_type} (mL)",
                               min_value=FLASK_TYPES[flask_type]['media_min'],
                               max_value=FLASK_TYPES[flask_type]['media_max'],
                               value=FLASK_TYPES[flask_type]['media_min'])
        media_freq = st.slider("Media Change Frequency (days)", 
                              min_value=1, max_value=4, value=2)
    
    # Calculations
    results = calculate_msc_therapy(weight, desired_dose, flask_type, media_volume, media_freq)
    recommended_dose = GVHD_RESPONSE[gvhd_grade]
    
    # Results
    st.header("Therapy Parameters")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("PBSC Volume Needed", f"{results['pbsc_volume']*1000:.0f} mL")
        st.metric(f"{flask_type} Flasks Needed", results['flasks'])
        
    with col2:
        st.metric("Passages Needed", results['passages'])
        st.metric("Culture Duration", f"{results['culture_days']} days")
        
    with col3:
        st.metric("Total Media Needed", f"{results['total_media']:.0f} mL")
        st.metric("Recommended Dose Range", 
                 f"{recommended_dose['min_dose']}-{recommended_dose['max_dose']} ×10⁶/kg")
    
    # Plots
    st.header("Biological Projections")
    
    col1, col2 = st.columns(2)
    with col1:
        st.pyplot(plot_growth_curve(results['initial_cells'], results['target_cells'], results['culture_days']))
    
    with col2:
        st.pyplot(plot_dose_response(gvhd_grade, desired_dose))
    
    # Protocol Notes
    st.header("Optimized Protocol")
    st.markdown(f"""
    **For GVHD {gvhd_grade}:**
    - Recommended dose: {recommended_dose['min_dose']}-{recommended_dose['max_dose']} ×10⁶ MSCs/kg
    - Expected response: {recommended_dose['response'][0]}–{recommended_dose['response'][1]}%
    
    **Culture Protocol:**
    - Seeding density: {SEEDING_DENSITY} cells/cm²
    - Target confluency: {HARVEST_CONFLUENCY}%
    - Media changes: Every {media_freq} days with {media_volume} mL per {flask_type}
    - Expected expansion: {GROWTH_RATE:.1f} doublings/day
    
    **Quality Control:**
    - Viability >90% (Trypan Blue)
    - Surface markers: CD73+/CD90+/CD105+ >95%
    - Sterility testing mandatory
    """)

if __name__ == "__main__":
    main()
