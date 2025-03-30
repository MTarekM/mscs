import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

# Constants
MSC_YIELD = {
    'Haemonetics': (1.0, 1.5),  # ×10⁶ MSCs/50mL
    'Spectra Optia': (1.5, 2.0)  # Better yield
}

FLASK_TYPES = {
    'T25': {'surface': 25, 'media': 5, 'max_cells': 0.25},  # ×10⁶ cells
    'T75': {'surface': 75, 'media': 15, 'max_cells': 0.75},
    'T175': {'surface': 175, 'media': 30, 'max_cells': 1.75}
}

GVHD_RESPONSE = {
    'Grade I': {'base': 70, 'dose_factor': 15},
    'Grade II': {'base': 50, 'dose_factor': 20},
    'Grade III-IV': {'base': 30, 'dose_factor': 25}
}

PLASMA_PRIMING_RATIO = 0.05  # 5% of PBSC volume

def calculate_msc_therapy(weight, dose, separator, gvhd_grade, flask_type, media_freq, plasma_priming):
    # Calculate total MSCs needed (×10⁶ cells)
    total_msc = dose * weight
    
    # Standard 50mL collection volume
    pbsc_volume = 0.05  # 50mL in liters
    yield_low, yield_high = MSC_YIELD[separator]
    
    # Calculate flasks needed (1 passage only)
    max_cells_per_flask = FLASK_TYPES[flask_type]['max_cells']
    flasks = int(np.ceil(total_msc / max_cells_per_flask))
    
    # Culture duration (max 14 days)
    culture_days = min(14, flasks * media_freq)
    
    # Plasma volume
    plasma_volume = pbsc_volume * PLASMA_PRIMING_RATIO if plasma_priming else 0
    
    # Media calculation
    media_per_flask = FLASK_TYPES[flask_type]['media']
    total_media = media_per_flask * flasks * (culture_days/media_freq)
    
    return {
        'pbsc_volume': pbsc_volume,
        'flasks': flasks,
        'passages': 1,
        'culture_days': culture_days,
        'plasma_volume': plasma_volume,
        'total_media': total_media,
        'yield_range': (yield_low, yield_high)
    }

def plot_growth_curve(days):
    x = np.linspace(0, days, 100)
    y = 1 / (1 + np.exp(-0.5*(x-3)))  # Sigmoid growth using numpy only
    fig, ax = plt.subplots()
    ax.plot(x, y, 'b-', linewidth=2)
    ax.set_xlabel('Culture Days')
    ax.set_ylabel('Relative MSC Growth')
    ax.set_title('Expected MSC Growth Curve')
    ax.grid(True, alpha=0.3)
    return fig

def plot_response_curve(grade, dose):
    base = GVHD_RESPONSE[grade]['base']
    factor = GVHD_RESPONSE[grade]['dose_factor']
    doses = np.linspace(0.5, 2.0, 100)
    response = base + factor*(doses-0.5)
    response = np.clip(response, 0, 95)
    
    fig, ax = plt.subplots()
    ax.plot(doses, response, 'r-', linewidth=2)
    ax.axvline(dose, color='k', linestyle='--')
    ax.set_xlabel('Dose (×10⁶ MSCs/kg)')
    ax.set_ylabel('Response Probability (%)')
    ax.set_title(f'GVHD {grade} Response Probability')
    ax.grid(True, alpha=0.3)
    return fig

def main():
    st.set_page_config(page_title="MSC Therapy Calculator", layout="wide")
    st.title("Advanced MSC Therapy Calculator for GVHD")
    
    # Inputs
    with st.sidebar:
        st.header("Patient Parameters")
        weight = st.slider("Weight (kg)", 30, 120, 70)
        gvhd_grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        dose = st.slider("Desired Dose (×10⁶ MSCs/kg)", 0.5, 2.0, 1.0, 0.1)
        
        st.header("Lab Parameters")
        separator = st.selectbox("Cell Separator", list(MSC_YIELD.keys()))
        flask_type = st.selectbox("Flask Type", list(FLASK_TYPES.keys()))
        media_freq = st.slider("Media Change (days)", 1, 4, 2)
        plasma_priming = st.checkbox("Include Plasma Priming")
    
    # Calculations
    results = calculate_msc_therapy(weight, dose, separator, gvhd_grade, flask_type, media_freq, plasma_priming)
    
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
        if plasma_priming:
            st.metric("Plasma Volume Needed", f"{results['plasma_volume']*1000:.0f} mL")
        st.metric("Total Media Needed", f"{results['total_media']:.0f} mL")
    
    # Plots
    st.header("Biological Projections")
    col1, col2 = st.columns(2)
    
    with col1:
        st.pyplot(plot_growth_curve(results['culture_days']))
    
    with col2:
        st.pyplot(plot_response_curve(gvhd_grade, dose))
    
    # Protocol Notes
    st.header("Clinical Protocol")
    st.markdown(f"""
    **For GVHD {gvhd_grade}:**
    - Expected yield: {results['yield_range'][0]}–{results['yield_range'][1]} ×10⁶ MSCs/50mL
    - Recommended media change: Every {media_freq} days
    
    **Culture Setup:**
    - Flask type: {flask_type} ({FLASK_TYPES[flask_type]['media']} mL media/flask)
    - Target confluence: 70-80% before harvest
    - Maximum culture duration: 14 days
    
    **Quality Control:**
    - Viability >90% by trypan blue
    - Surface markers: CD73+/CD90+/CD105+ >95%
    - Hematopoietic markers: CD45- <2%
    - Sterility testing required pre-infusion
    """)

if __name__ == "__main__":
    main()
