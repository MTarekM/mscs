import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math

# Constants with validation
SEEDING_DENSITY = 5000  # cells/cm² (applied per flask)
GROWTH_RATE = 0.5       # doublings per day (for intra-passage exponential growth)
MIN_PASSAGE1_DAYS = 14  # Fixed duration for the first passage (includes trypsinization)
ADDITIONAL_PASSAGE_DAYS = 7  # Each subsequent passage lasts 7 days
MAX_PASSAGES = 3       # Maximum number of passages allowed
SAFETY_FACTOR = 1.2    # Safety factor for re-seeding losses
PLASMA_VOLUME_PER_FLASK = 50  # mL plasma required per initial flask for priming

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
# Note: max_cells is the maximum number of cells obtainable in that flask at confluency.
FLASK_TYPES = {
    'T25': {'surface': 25, 'media': 5, 'max_cells': 2.5e5},
    'T75': {'surface': 75, 'media': 15, 'max_cells': 2.0e6},
    'T175': {'surface': 175, 'media': 30, 'max_cells': 5.0e6},
    'Default': {'surface': 75, 'media': 15, 'max_cells': 2.0e6}  # Fallback
}

def safe_get(dictionary, key, default_key='Default'):
    """Safely get dictionary value with fallback"""
    return dictionary.get(key, dictionary.get(default_key))

def calculate_msc_therapy(weight, desired_dose, separator, flask_type, plasma_priming=False):
    # Validate inputs
    weight = max(30, min(weight, 120))              # in kg
    desired_dose = max(0.5, min(desired_dose, 2.0))   # in ×10⁶/kg
    
    # Get parameters with fallbacks
    sep_yield = safe_get(SEPARATOR_YIELD, separator)
    flask_data = safe_get(FLASK_TYPES, flask_type)
    
    # Total cell dose needed (absolute number of cells)
    total_cells = desired_dose * weight * 1e6  # Convert ×10⁶/kg to cells
    
    # Calculate PBSC volume (mL)
    try:
        pbsc_ml = total_cells / sep_yield
    except Exception:
        pbsc_ml = 50  # Fallback
    pbsc_ml = max(50, pbsc_ml)
    
    # Determine initial cells seeded in one flask
    initial_cells = flask_data['surface'] * SEEDING_DENSITY  # cells per flask
    
    # --- Determine the minimal number of initial flasks (N0) needed so that after MAX_PASSAGES yield meets/exceeds target.
    # We now search up to 200 flasks and use a safety factor when re-seeding.
    found = False
    for N0 in range(1, 201):
        # Passage 1:
        passage1_yield = N0 * flask_data['max_cells']
        # Passage 2:
        flasks_p2 = math.ceil((passage1_yield / initial_cells) * SAFETY_FACTOR)
        passage2_yield = flasks_p2 * flask_data['max_cells']
        # Passage 3:
        flasks_p3 = math.ceil((passage2_yield / initial_cells) * SAFETY_FACTOR)
        passage3_yield = flasks_p3 * flask_data['max_cells']
        
        final_yield = passage3_yield  # Maximum obtainable with 3 passages
        if final_yield >= total_cells:
            initial_flasks = N0
            found = True
            break
    if not found:
        # If even 200 flasks aren’t enough, take the best obtainable yield
        initial_flasks = 200
        flasks_p2 = math.ceil((initial_flasks * flask_data['max_cells'] / initial_cells) * SAFETY_FACTOR)
        flasks_p3 = math.ceil((flasks_p2 * flask_data['max_cells'] / initial_cells) * SAFETY_FACTOR)
        passage3_yield = flasks_p3 * flask_data['max_cells']
        final_yield = passage3_yield

    # Total culture duration:
    # Passage 1 takes MIN_PASSAGE1_DAYS and each subsequent passage adds ADDITIONAL_PASSAGE_DAYS.
    total_days = MIN_PASSAGE1_DAYS + (MAX_PASSAGES - 1) * ADDITIONAL_PASSAGE_DAYS

    # Media changes:
    # Passage 1: media changes on days 1,4,7,11  → 4 changes
    # Passage 2 and 3 (each 7 days): assume 2 changes each.
    total_media_changes = 4 + 2 * (MAX_PASSAGES - 1)
    # For media volume, we now base it on the initial flasks used.
    total_media = initial_flasks * flask_data['media'] * total_media_changes

    # Plasma priming volume (if selected) based on initial flasks.
    plasma_volume = initial_flasks * PLASMA_VOLUME_PER_FLASK if plasma_priming else 0

    return {
        'pbsc_ml': pbsc_ml,
        'initial_flasks': initial_flasks,
        'passages': MAX_PASSAGES,
        'total_days': total_days,
        'total_media': total_media,
        'initial_cells': initial_cells,
        'target_cells': total_cells,
        'passage1_yield': initial_flasks * flask_data['max_cells'],
        'flasks_p2': flasks_p2,
        'passage2_yield': passage1_yield if MAX_PASSAGES < 2 else flasks_p2 * flask_data['max_cells'],
        'flasks_p3': flasks_p3,
        'passage3_yield': passage3_yield,
        'plasma_volume': plasma_volume
    }

def plot_growth_curve(results, flask_data):
    """
    Plots a piecewise growth curve showing overall cell yield over time.
    Passage 1: from day 0 to MIN_PASSAGE1_DAYS, exponential growth from (N0*initial_cells) to passage1_yield.
    Passage 2: next ADDITIONAL_PASSAGE_DAYS days, from (flasks_p2*initial_cells) to passage2_yield.
    Passage 3: next ADDITIONAL_PASSAGE_DAYS days, from (flasks_p3*initial_cells) to passage3_yield.
    """
    initial_cells = results['initial_cells']
    N0 = results['initial_flasks']
    
    # Calculate starting yields for each passage:
    start1 = N0 * initial_cells
    end1 = results['passage1_yield']
    
    flasks_p2 = results['flasks_p2']
    start2 = flasks_p2 * initial_cells
    end2 = results['passage2_yield']
    
    flasks_p3 = results['flasks_p3']
    start3 = flasks_p3 * initial_cells
    end3 = results['passage3_yield']
    
    # Define time segments
    t1 = np.linspace(0, MIN_PASSAGE1_DAYS, 100)
    t2 = np.linspace(MIN_PASSAGE1_DAYS, MIN_PASSAGE1_DAYS + ADDITIONAL_PASSAGE_DAYS, 100)
    t3 = np.linspace(MIN_PASSAGE1_DAYS + ADDITIONAL_PASSAGE_DAYS, results['total_days'], 100)
    
    # Exponential growth in each segment (using the appropriate growth rate calculated from endpoints)
    factor1 = end1 / start1 if start1 > 0 else 1
    rate1 = np.log(factor1) / (MIN_PASSAGE1_DAYS)
    growth1 = start1 * np.exp(rate1 * t1)
    
    factor2 = end2 / start2 if start2 > 0 else 1
    rate2 = np.log(factor2) / (ADDITIONAL_PASSAGE_DAYS)
    growth2 = start2 * np.exp(rate2 * (t2 - MIN_PASSAGE1_DAYS))
    
    factor3 = end3 / start3 if start3 > 0 else 1
    rate3 = np.log(factor3) / (ADDITIONAL_PASSAGE_DAYS)
    growth3 = start3 * np.exp(rate3 * (t3 - (MIN_PASSAGE1_DAYS + ADDITIONAL_PASSAGE_DAYS)))
    
    # Combine segments
    t_total = np.concatenate([t1, t2, t3])
    growth_total = np.concatenate([growth1, growth2, growth3])
    
    fig, ax = plt.subplots()
    ax.plot(t_total, growth_total, 'g-', label='Overall MSC Yield')
    # Plot horizontal lines showing target and best obtainable yield
    ax.axhline(results['target_cells'], color='r', linestyle='--', label='Desired Dose')
    ax.axhline(end3, color='b', linestyle='--', label='Max obtainable (3 passages)')
    ax.set_xlabel('Culture Days')
    ax.set_ylabel('Total Cells')
    ax.legend()
    ax.grid(True)
    return fig

def plot_gvhd_probability(grade_data, dose):
    """
    Plots a predicted GVHD remission probability curve based on clinical grade data.
    We assume a Gaussian-like curve peaking at the optimum dose = (min_dose + max_dose)/2.
    """
    min_dose = grade_data['min_dose']
    max_dose = grade_data['max_dose']
    opt_dose = (min_dose + max_dose) / 2
    lower_prob, upper_prob = grade_data['response']  # e.g. 60 and 80%
    sigma = (max_dose - min_dose) / 2  # spread such that min and max roughly span 2σ
    x = np.linspace(0.5, 2.0, 200)
    # Gaussian-like function scaled between lower_prob and upper_prob
    y = lower_prob + (upper_prob - lower_prob) * np.exp(-((x - opt_dose)**2) / (2 * sigma**2))
    
    fig, ax = plt.subplots()
    ax.plot(x, y, label="Predicted GVHD Remission Probability (%)")
    ax.axvline(dose, color='r', linestyle='--', label=f"Selected Dose: {dose:.2f}")
    ax.set_xlabel("Dose (×10⁶/kg)")
    ax.set_ylabel("Probability (%)")
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
        plasma_priming = st.checkbox("Plasma priming from GVHD patient", value=False)
    
    # Get grade data with fallback
    grade_data = safe_get(GVHD_RESPONSE, gvhd_grade)
    
    # Calculate therapy parameters (pass plasma_priming option)
    results = calculate_msc_therapy(weight, desired_dose, separator, flask_type, plasma_priming)
    flask_data = safe_get(FLASK_TYPES, flask_type)
    
    # Display results
    st.header("Therapy Parameters")
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("PBSC Volume", f"{results['pbsc_ml']:.0f} mL")
        st.metric("Initial Flasks Needed", results['initial_flasks'])
    with col2:
        st.metric("Total Passages", results['passages'])
        st.metric("Culture Duration", f"{results['total_days']} days")
    with col3:
        st.metric("Total Media", f"{results['total_media']:.0f} mL")
        st.metric("Recommended Dose", f"{grade_data['min_dose']}-{grade_data['max_dose']} ×10⁶/kg")
    with col4:
        if plasma_priming:
            st.metric("Plasma Volume", f"{results['plasma_volume']:.0f} mL")
        else:
            st.write("Plasma priming not selected")
    
    # Growth curve plot
    st.header("Growth Projection")
    try:
        plot_fig = plot_growth_curve(results, flask_data)
        st.pyplot(plot_fig)
    except Exception as e:
        st.warning("Could not generate growth curve with current parameters")
    
    # GVHD Remission Probability plot
    st.header("GVHD Remission Probability")
    try:
        gvhd_fig = plot_gvhd_probability(grade_data, desired_dose)
        st.pyplot(gvhd_fig)
    except Exception as e:
        st.warning("Could not generate GVHD remission probability plot")
    
    # Protocol notes
    st.header("Protocol Notes")
    st.markdown(f"""
    **For GVHD {gvhd_grade}:**
    - Recommended dose: {grade_data['min_dose']}-{grade_data['max_dose']} ×10⁶/kg
    - Expected remission response: {grade_data['response'][0]}–{grade_data['response'][1]}%
    
    **Culture Protocol:**
    - Seeding density: {SEEDING_DENSITY:,} cells/cm² (per flask)
    - Flask type: {flask_type} with surface area {flask_data['surface']} cm² and max capacity {flask_data['max_cells']:,} cells
    - **Passage Schedule:**  
      - Passage 1: {MIN_PASSAGE1_DAYS} days (media changes on days 1,4,7,11)  
      - Passage 2 & 3: {ADDITIONAL_PASSAGE_DAYS} days each (assumed 2 media changes per passage)  
      - Total maximum passages: {MAX_PASSAGES} (beyond which cryopreservation or infusion is recommended)
    - Initial flasks required: {results['initial_flasks']}
    - Final yield (after 3 passages): {results['passage3_yield']:,} cells
    - Culture duration: {results['total_days']} days
    - Total media used: {results['total_media']:.0f} mL
    - { "Plasma priming selected: " + str(results['plasma_volume']) + " mL plasma required." if results['plasma_volume'] else "No plasma priming." }
    
    **Quality Control:**
    - Viability >90% (Trypan Blue)
    - CD73+/CD90+/CD105+ >95%
    - CD45- <2%
    - Sterility testing required
    """)
    
if __name__ == "__main__":
    main()
