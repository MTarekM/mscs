import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

# Constants with validation
SEEDING_DENSITY = 5000  # cells/cm²
GROWTH_RATE = 0.5  # doublings per day
MIN_CULTURE_DAYS = 3  # Minimum days per passage
TARGET_CONFLUENCY = 80  # %

# Clinical parameters with fallbacks
GVHD_RESPONSE = {
    'Grade I': {'min_dose': 0.5, 'max_dose': 1.0, 'response': [60, 80]},
    'Grade II': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]},
    'Grade III-IV': {'min_dose': 1.5, 'max_dose': 2.0, 'response': [40, 60]},
    'Default': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]}  # Fallback
}

# Cell separator yields (cells/mL) with validation
SEPARATOR_YIELD = {
    'Haemonetics': 1.0e6,
    'Spectra Optia': 2.0e6,
    'Default': 1.5e6  # Fallback
}

# Flask parameters with validation
FLASK_TYPES = {
    'T25': {'surface': 25, 'media': 5, 'max_cells': 2.5e5},
    'T75': {'surface': 75, 'media': 15, 'max_cells': 2.0e6},
    'T175': {'surface': 175, 'media': 30, 'max_cells': 5.0e6},
    'Default': {'surface': 75, 'media': 15, 'max_cells': 2.0e6}  # Fallback
}

def safe_get(dictionary, key, default_key='Default'):
    """Safely get dictionary value with fallback"""
    return dictionary.get(key, dictionary.get(default_key))

def calculate_msc_therapy(weight, desired_dose, separator, flask_type, media_freq):
    # Validate all inputs
    weight = max(30, min(weight, 120))
    desired_dose = max(0.5, min(desired_dose, 2.0))
    media_freq = max(1, min(media_freq, 4))
    
    # Get parameters with fallbacks
    sep_yield = safe_get(SEPARATOR_YIELD, separator)
    flask_data = safe_get(FLASK_TYPES, flask_type)
    
    # Calculate total cells needed (×10⁶)
    total_cells = desired_dose * weight
    
    # Calculate PBSC volume (mL)
    try:
        pbsc_ml = (total_cells * 1e6) / sep_yield
    except:
        pbsc_ml = 50  # Fallback to 50mL
    pbsc_ml = max(50, pbsc_ml)
    
    # Calculate initial cells per flask (×10⁶)
    initial_cells = (flask_data['surface'] * SEEDING_DENSITY) / 1e6
    
    # Calculate passages needed safely
    passages = 0
    current_cells = initial_cells
    while current_cells < total_cells and passages <= 10:  # Max 10 passages
        passages += 1
        current_cells *= 5  # 5-fold expansion
    
    # Calculate culture duration safely
    try:
        days_per_passage = np.log((TARGET_CONFLUENCY/100)) / np.log(2) / GROWTH_RATE
    except:
        days_per_passage = MIN_CULTURE_DAYS
    total_days = max(MIN_CULTURE_DAYS, int(np.ceil(days_per_passage * passages)))
    
    # Calculate flasks needed safely
    try:
        flasks = int(np.ceil(total_cells / (flask_data['max_cells'] / 1e6)))
    except:
        flasks = 1
    flasks = max(1, flasks)
    
    # Calculate total media needed
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
    initial = max(0.1, initial)
    target = max(initial, target)
    days = max(1, days)
    
    x = np.linspace(0, days, 100)
    growth = initial * np.exp(GROWTH_RATE * x)
    
    # Safe plateau detection
    with np.errstate(all='ignore'):  # Ignore numpy warnings
        reached = np.where(growth >= target)[0]
        if reached.size > 0 and reached[0] < len(growth):
            growth[reached[0]:] = target
    
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
        weight = st.slider("Weight (kg)", 30, 120, 70)
        gvhd_grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        desired_dose = st.slider("Dose (×10⁶/kg)", 0.5, 2.0, 1.0, 0.1)
        
        st.header("Lab Parameters")
        separator = st.selectbox("Cell Separator", list(SEPARATOR_YIELD.keys()))
        flask_type = st.selectbox("Flask Type", list(FLASK_TYPES.keys()))
        media_freq = st.slider("Media Change (days)", 1, 4, 2)
    
    # Get grade data with fallback
    grade_data = safe_get(GVHD_RESPONSE, gvhd_grade)
    
    # Calculate therapy parameters
    results = calculate_msc_therapy(weight, desired_dose, separator, flask_type, media_freq)
    
    # Display results
    st.header("Therapy Parameters")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("PBSC Volume", f"{results['pbsc_ml']:.0f} mL")
        st.metric("Flasks Needed", results['flasks'])
    with col2:
        st.metric("Passages", results['passages'])
        st.metric("Culture Days", results['total_days'])
    with col3:
        st.metric("Total Media", f"{results['total_media']:.0f} mL")
        st.metric("Recommended Dose", 
                f"{grade_data['min_dose']}-{grade_data['max_dose']} ×10⁶/kg")
    
    # Growth curve plot
    st.header("Growth Projection")
    try:
        st.pyplot(plot_growth_curve(results['initial_cells'], 
                                  results['target_cells'], 
                                  results['total_days']))
    except:
        st.warning("Could not generate growth curve with current parameters")
    
    # Protocol notes
    st.header("Protocol Notes")
    st.markdown(f"""
    **For GVHD {gvhd_grade}:**
    - Recommended dose: {grade_data['min_dose']}-{grade_data['max_dose']} ×10⁶/kg
    - Expected response: {grade_data['response'][0]}–{grade_data['response'][1]}%
    
    **Culture Protocol:**
    - Seeding density: {SEEDING_DENSITY:,} cells/cm²
    - Target confluency: {TARGET_CONFLUENCY}%
    - Media changes: Every {media_freq} days
    - Expected expansion rate: {GROWTH_RATE:.1f} doublings/day
    
    **Quality Control:**
    - Viability >90% (Trypan Blue)
    - CD73+/CD90+/CD105+ >95%
    - CD45- <2%
    - Sterility testing required
    """)

if __name__ == "__main__":
    main()
