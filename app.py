import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

# Constants
SEEDING_DENSITY = 5000  # cells/cm²
GROWTH_RATE = 0.5  # doublings per day
MIN_CULTURE_DAYS = 3  # Minimum days per passage
TARGET_CONFLUENCY = 80  # %

# Clinical parameters
GVHD_RESPONSE = {
    'Grade I': {'min_dose': 0.5, 'max_dose': 1.0, 'response': [60, 80]},
    'Grade II': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]},
    'Grade III-IV': {'min_dose': 1.5, 'max_dose': 2.0, 'response': [40, 60]}
}

# Cell separator yields (cells/mL)
SEPARATOR_YIELD = {
    'Haemonetics': 1.0e6,  # cells/mL
    'Spectra Optia': 2.0e6
}

# Flask parameters
FLASK_TYPES = {
    'T25': {'surface': 25, 'media': 5, 'max_cells': 2.5e5},
    'T75': {'surface': 75, 'media': 15, 'max_cells': 2.0e6},
    'T175': {'surface': 175, 'media': 30, 'max_cells': 5.0e6}
}

def calculate_msc_therapy(weight, desired_dose, separator, flask_type, media_freq):
    # Validate inputs
    weight = max(30, weight)
    desired_dose = max(0.5, min(desired_dose, 2.0))
    media_freq = max(1, media_freq)
    
    # Calculate total cells needed (×10⁶)
    total_cells = desired_dose * weight  # in million cells
    
    # Calculate PBSC volume needed (mL) with validation
    try:
        pbsc_ml = (total_cells * 1e6) / SEPARATOR_YIELD[separator]
    except ZeroDivisionError:
        pbsc_ml = 0
    pbsc_ml = max(50, pbsc_ml)  # Minimum 50mL
    
    # Calculate initial cells per flask
    flask_data = FLASK_TYPES[flask_type]
    initial_cells = (flask_data['surface'] * SEEDING_DENSITY) / 1e6  # million cells
    
    # Calculate required passages safely
    passages = 0
    current_cells = initial_cells
    while current_cells < total_cells:
        passages += 1
        current_cells *= 5  # 5-fold expansion per passage
        if passages > 10:  # Safety break
            break
    
    # Calculate culture duration with validation
    try:
        days_per_passage = np.log((TARGET_CONFLUENCY/100 * SEEDING_DENSITY)/SEEDING_DENSITY) / np.log(2) / GROWTH_RATE
    except:
        days_per_passage = MIN_CULTURE_DAYS
    total_days = max(MIN_CULTURE_DAYS, int(np.ceil(days_per_passage * passages)))
    
    # Calculate flasks needed with validation
    try:
        flasks = int(np.ceil(total_cells / (flask_data['max_cells'] / 1e6))
    except ZeroDivisionError:
        flasks = 0
    flasks = max(1, flasks)
    
    # Calculate media needed safely
    media_changes = max(1, total_days // media_freq)
    total_media = flask_data['media'] * flasks * media_changes
    
    return {
        'pbsc_ml': pbsc_ml,
        'flasks': flasks,
        'passages': passages,
        'total_days': total_days,
        'total_media': total_media,
        'initial_cells': initial_cells,
        'target_cells': total_cells
    }

def plot_growth_curve(initial, target, days):
    # Validate inputs
    days = max(1, days)
    initial = max(0.1, initial)
    target = max(initial, target)
    
    x = np.linspace(0, days, 100)
    growth = initial * np.exp(GROWTH_RATE * x)
    
    # Safe plateau detection
    reached = np.where(growth >= target)[0]
    if reached.size > 0 and reached[0] < len(growth):
        plateau_idx = reached[0]
        growth[plateau_idx:] = target
    
    fig, ax = plt.subplots()
    ax.plot(x, growth, 'g-', label='MSC Growth')
    ax.axhline(target, color='r', linestyle='--', label='Target')
    ax.set_xlabel('Culture Days')
    ax.set_ylabel('Cells (×10⁶)')
    ax.legend()
    ax.grid(True)
    return fig

def main():
    st.set_page_config(page_title="MSC Therapy Calculator", layout="wide")
    st.title("Advanced MSC Therapy Calculator")
    
    with st.sidebar:
        st.header("Patient Parameters")
        weight = st.slider("Weight (kg)", 40, 120, 70)
        gvhd_grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        desired_dose = st.slider("Dose (×10⁶/kg)", 0.5, 2.0, 1.0)
        
        st.header("Lab Parameters")
        separator = st.selectbox("Cell Separator", list(SEPARATOR_YIELD.keys()))
        flask_type = st.selectbox("Flask Type", list(FLASK_TYPES.keys()))
        media_freq = st.slider("Media Change (days)", 1, 4, 2)
    
    results = calculate_msc_therapy(weight, desired_dose, separator, flask_type, media_freq)
    
    st.header("Therapy Parameters")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("PBSC Volume", f"{results['pbsc_ml']:.0f} mL")
        st.metric("Flasks Needed", results['flasks'])
    with col2:
        st.metric("Passages", results['passages'])
        st.metric("Culture Days", results['total_days'])
    with col3:
        st.metric("Total Media", f"{results['total_media']} mL")
        st.metric("Recommended Dose", f"{GVHD_RESPONSE[gvhd_grade]['min_dose']}-{GVHD_RESPONSE[gvhd_grade]['max_dose']} ×10⁶/kg")
    
    st.header("Growth Projection")
    st.pyplot(plot_growth_curve(results['initial_cells'], results['target_cells'], results['total_days']))
    
    st.markdown("""
    **Protocol Notes:**
    - Maintain 37°C with 5% CO₂
    - Media changes every selected interval
    - Harvest at 80% confluency
    - Validate CD73+/CD90+/CD105+ >95%
    """)

if __name__ == "__main__":
    main()
