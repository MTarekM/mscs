import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

# Constants
SEEDING_DENSITY = 5000  # cells/cm²
GROWTH_RATE = 0.5  # doublings per day
MIN_CULTURE_DAYS = 3  # Minimum days per passage

def plot_growth_curve(initial_cells, target_cells, days):
    # Ensure we have valid time points
    days = max(1, days)  # Minimum 1 day
    x = np.linspace(0, days, 100)
    
    # Calculate exponential growth
    growth = initial_cells * np.exp(GROWTH_RATE * x)
    
    # Find when we reach the target (if we do)
    reached_target = np.where(growth >= target_cells)[0]
    
    # Apply plateau if target is reached
    if len(reached_target) > 0:
        plateau_idx = reached_target[0]
        # Ensure we don't exceed array bounds
        if plateau_idx < len(growth):
            growth[plateau_idx:] = target_cells
    
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
    doses = np.linspace(grade_data['min_dose'], grade_data['max_dose'], 100)
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
    - Seeding density: {SEEDING_DENSITY:,} cells/cm²
    - Target confluency: {HARVEST_CONFLUENCY}%
    - Media changes: Every {media_freq} days with {media_volume} mL per {flask_type}
    - Expected expansion rate: {GROWTH_RATE:.1f} doublings/day
    
    **Quality Control:**
    - Viability >90% (Trypan Blue)
    - Surface markers: CD73+/CD90+/CD105+ >95%
    - Sterility testing mandatory
    """)

if __name__ == "__main__":
    main()
