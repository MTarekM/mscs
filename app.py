import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

# Constants from literature
MSC_YIELD = {
    'Haemonetics': (0.8, 1.2),  # (low, high) ×10⁶ MSCs/L due to RBC contamination
    'Spectra Optia': (1.5, 2.5)  # Better MSC yield
}

FLASK_CAPACITY = 175  # mL for T75 flask with 15-20mL media
PASSAGE_EXPANSION = 5  # 5-fold expansion per passage
MEDIA_CHANGE_FREQ = 2  # Recommended days between media changes
MAX_CULTURE_DAYS = 14  # Maximum recommended culture duration

GVHD_DOSING = {
    'Grade I': (0.5, 1.0),
    'Grade II': (1.0, 1.5),
    'Grade III-IV': (1.5, 2.0)
}

PLASMA_PRIMING_RATIO = 0.1  # 10% of PBSC volume needed for plasma priming

def calculate_msc_parameters(weight, dose, separator, gvhd_grade, plasma_priming):
    # Calculate total required MSCs
    total_msc = dose * weight
    
    # Get yield range based on separator
    yield_low, yield_high = MSC_YIELD[separator]
    
    # Calculate required PBSC volume
    pbsc_volume_low = total_msc / yield_high
    pbsc_volume_high = total_msc / yield_low
    
    # Calculate flasks needed
    flasks = int(np.ceil(pbsc_volume_high / (FLASK_CAPACITY/1000)))  # Convert flask capacity to liters
    
    # Calculate culture duration and passages
    passages = max(1, int(np.ceil(np.log(total_msc/1e6)/np.log(PASSAGE_EXPANSION)))
    culture_days = min(passages * MEDIA_CHANGE_FREQ * 2, MAX_CULTURE_DAYS)
    
    # Plasma priming calculation
    plasma_volume = pbsc_volume_high * PLASMA_PRIMING_RATIO if plasma_priming else 0
    
    return {
        'pbsc_volume': (pbsc_volume_low, pbsc_volume_high),
        'flasks': flasks,
        'passages': passages,
        'culture_days': culture_days,
        'plasma_volume': plasma_volume,
        'recommended_dose': GVHD_DOSING[gvhd_grade]
    }

def main():
    st.set_page_config(page_title="MSC Therapy Calculator", layout="wide")
    st.title("Donor MSC Therapy Calculator for GVHD")
    
    # Input parameters
    with st.sidebar:
        st.header("Patient Parameters")
        weight = st.slider("Patient Weight (kg)", 30, 120, 70)
        gvhd_grade = st.selectbox("GVHD Grade", list(GVHD_DOSING.keys()))
        dose = st.slider("Desired Dose (×10⁶ MSCs/kg)", 0.5, 2.0, 1.0, 0.1)
        separator = st.selectbox("Cell Separator", ['Haemonetics', 'Spectra Optia'])
        plasma_priming = st.checkbox("Include Plasma Priming")
    
    # Calculations
    results = calculate_msc_parameters(weight, dose, separator, gvhd_grade, plasma_priming)
    
    # Display results
    st.header("Therapy Parameters")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Required PBSC Volume", 
                 f"{results['pbsc_volume'][0]:.1f}-{results['pbsc_volume'][1]:.1f} L")
        st.metric("T75 Flasks Required", results['flasks'])
        
    with col2:
        st.metric("Estimated Passages", results['passages'])
        st.metric("Culture Duration", f"{results['culture_days']} days")
        
    with col3:
        if plasma_priming:
            st.metric("Plasma Volume Needed", f"{results['plasma_volume']:.1f} L")
        st.metric("Recommended Dose Range", 
                 f"{results['recommended_dose'][0]}-{results['recommended_dose'][1]} ×10⁶/kg")
    
    # GVHD Efficacy Plot
    st.header("Clinical Efficacy by Dose")
    grades = list(GVHD_DOSING.keys())
    response_rates = [60, 75, 85]  # Hypothetical response rates based on literature
    
    fig, ax = plt.subplots()
    bars = ax.bar(grades, response_rates, color=['#4CAF50', '#FFC107', '#F44336'])
    
    ax.set_ylabel('Response Rate (%)')
    ax.set_ylim(0, 100)
    ax.set_title('Estimated Clinical Response by GVHD Grade')
    
    # Add dose ranges to bars
    for bar, grade in zip(bars, grades):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{GVHD_DOSING[grade][0]}-{GVHD_DOSING[grade][1]} ×10⁶/kg',
                ha='center', va='bottom')
    
    st.pyplot(fig)
    
    # Protocol Notes
    st.header("Protocol Considerations")
    st.markdown("""
    - **Cell Separation:** 
      - Spectra Optia recommended for better MSC yield
      - Process within 6 hours of collection
    - **Culture Conditions:**
      - Maintain 37°C with 5% CO₂
      - Media change every 2 days
      - Monitor confluence daily
    - **Quality Control:**
      - CD73+/CD90+/CD105+ > 95%
      - CD45- < 2%
      - Sterility testing required
    """)
    
    # References
    st.caption("References: Frontiers in Immunology 2021; Stem Cell Research & Therapy 2020; Blood 2005")

if __name__ == "__main__":
    main()
