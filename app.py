import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math

# Constants with validation
SEEDING_DENSITY = 5000  # cells/cm² (applied per flask)
GROWTH_RATE = 0.5  # doublings per day
MIN_CULTURE_DAYS = 3  # Minimum days per passage
TARGET_CONFLUENCY = 80  # % (no longer used directly for time)

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
# Note: max_cells is expressed in absolute numbers.
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
    # Validate inputs
    weight = max(30, min(weight, 120))              # in kg
    desired_dose = max(0.5, min(desired_dose, 2.0))   # in ×10⁶/kg
    media_freq = max(1, min(media_freq, 4))           # days between media changes
    
    # Get parameters with fallbacks
    sep_yield = safe_get(SEPARATOR_YIELD, separator)
    flask_data = safe_get(FLASK_TYPES, flask_type)
    
    # Total cell dose needed (absolute number of cells)
    total_cells = desired_dose * weight * 1e6  # Convert from ×10⁶ to cells
    
    # Calculate PBSC volume (mL)
    try:
        pbsc_ml = total_cells / sep_yield
    except Exception:
        pbsc_ml = 50  # Fallback to 50 mL if any issue
    pbsc_ml = max(50, pbsc_ml)
    
    # Determine initial cells in one flask (in cells)
    initial_cells = flask_data['surface'] * SEEDING_DENSITY  # cells seeded in a flask
    
    # Calculate the maximum expansion factor per passage given the flask's capacity
    expansion_factor = flask_data['max_cells'] / initial_cells
    # Ensure expansion_factor is at least a minimal expansion (>1)
    expansion_factor = max(1.1, expansion_factor)
    
    # Calculate number of passages needed using logarithms (round up)
    if initial_cells <= 0:
        passages = 0
    else:
        passages = math.ceil(np.log(total_cells / initial_cells) / np.log(expansion_factor))
        passages = max(1, passages)
    
    # Calculate days per passage from the expansion factor
    # Number of doublings required: log2(expansion_factor)
    doublings_needed = np.log(expansion_factor) / np.log(2)
    days_per_passage = max(MIN_CULTURE_DAYS, doublings_needed / GROWTH_RATE)
    
    # Total culture duration is passages times days per passage
    total_days = int(np.ceil(passages * days_per_passage))
    
    # Calculate number of flasks needed in the final expansion round
    flasks = math.ceil(total_cells / flask_data['max_cells'])
    flasks = max(1, flasks)
    
    # Calculate total media required:
    # Determine the number of media changes (at least one) over the culture period.
    media_changes = max(1, total_days // media_freq)
    total_media = flask_data['media'] * flasks * media_changes
    
    return {
        'pbsc_ml': pbsc_ml,
        'flasks': flasks,
        'passages': passages,
        'total_days': total_days,
        'total_media': total_media,
        'initial_cells': initial_cells,
        'target_cells': total_cells,
        'expansion_factor': expansion_factor,
        'days_per_passage': days_per_passage
    }

def plot_growth_curve(initial, target, days):
    # Validate inputs
    initial = max(0.1, initial)
    target = max(initial, target)
    days = max(1, days)
    
    x = np.linspace(0, days, 100)
    # Simulate exponential growth: initial cells * 2^(doublings)
    # Here we assume a constant doubling rate given by GROWTH_RATE.
    growth = initial * np.power(2, GROWTH_RATE * x)
    
    # Limit growth to the target (simulate plateau)
    with np.errstate(all='ignore'):
        reached = np.where(growth >= target)[0]
        if reached.size > 0:
            growth[reached[0]:] = target
    
    fig, ax = plt.subplots()
    ax.plot(x, growth, 'g-', label='MSC Growth')
    ax.axhline(target, color='r', linestyle='--', label='Target')
    ax.set_xlabel('Culture Days')
    ax.set_ylabel('Cells (absolute count)')
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
        media_freq = st.slider("Media Change Frequency (days)", 1, 4, 2)
    
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
        # Use the initial seeding from one flask and target as the final cell count per flask
        # (for plotting we use the expansion in one flask)
        plot_fig = plot_growth_curve(results['initial_cells'], 
                                  results['initial_cells'] * results['expansion_factor'],
                                  results['days_per_passage'])
        st.pyplot(plot_fig)
    except Exception as e:
        st.warning("Could not generate growth curve with current parameters")
    
    # Protocol notes
    st.header("Protocol Notes")
    st.markdown(f"""
    **For GVHD {gvhd_grade}:**
    - Recommended dose: {grade_data['min_dose']}-{grade_data['max_dose']} ×10⁶/kg
    - Expected response: {grade_data['response'][0]}–{grade_data['response'][1]}%
    
    **Culture Protocol:**
    - Seeding density: {SEEDING_DENSITY:,} cells/cm² (applied per flask)
    - Flask type: {flask_type} with surface area {safe_get(FLASK_TYPES, flask_type)['surface']} cm²
    - Max capacity per flask: {safe_get(FLASK_TYPES, flask_type)['max_cells']:,} cells
    - Media changes every {media_freq} days
    - Estimated expansion per passage: {results['expansion_factor']:.2f}-fold in ~{results['days_per_passage']:.1f} days
    - Total passages: {results['passages']} and culture duration: {results['total_days']} days
    
    **Quality Control:**
    - Viability >90% (Trypan Blue)
    - CD73+/CD90+/CD105+ >95%
    - CD45- <2%
    - Sterility testing required
    """)

if __name__ == "__main__":
    main()
