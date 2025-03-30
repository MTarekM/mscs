import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math

# Constants with validation
SEEDING_DENSITY = 5000  # cells/cm² (applied per flask)
GROWTH_RATE = 0.5       # doublings per day (no longer used in discrete passage model)
MIN_CULTURE_DAYS = 14   # fixed culture period per passage (media changes on days 1,4,7,11 and passage at day 14)
# TARGET_CONFLUENCY is no longer used in these calculations

# Clinical parameters with fallbacks
GVHD_RESPONSE = {
    'Grade I': {'min_dose': 0.5, 'max_dose': 1.0, 'response': [60, 80]},
    'Grade II': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]},
    'Grade III-IV': {'min_dose': 1.5, 'max_dose': 2.0, 'response': [40, 60]},
    'Default': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]}  # Fallback
}

# Cell separator yields (cells/mL) with validation (for PBSC volume calculation)
SEPARATOR_YIELD = {
    'Haemonetics': 1.0e6,
    'Spectra Optia': 2.0e6,
    'Default': 1.5e6  # Fallback
}

# Flask parameters with validation
# Note: max_cells is the maximum cell count attainable in a given flask at confluency.
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
    """
    Revised calculation that:
      - Limits expansion to 3 passages (each 14 days)
      - Uses fixed media change schedule (days 1,4,7,11 per passage)
      - Computes yield at each passage based on flask seeding density and maximum cell capacity.
    
    Assumptions:
      * Initial culture is started in one flask with cells seeded = (flask surface * SEEDING_DENSITY).
      * At confluency (day 14), a flask yields up to flask_data['max_cells'] cells.
      * For subsequent passages, cells are split into as many flasks as needed so that each receives the ideal seeding count.
      * Total yield after each passage is computed and compared with the target dose.
      * If the target dose is reached before passage 3, that passage is used; otherwise, the best yield after 3 passages is used.
    """
    # Validate inputs
    weight = max(30, min(weight, 120))              # in kg
    desired_dose = max(0.5, min(desired_dose, 2.0))   # in ×10⁶/kg
    media_freq = max(1, min(media_freq, 4))           # not used in fixed schedule below

    # Get parameters with fallbacks
    sep_yield = safe_get(SEPARATOR_YIELD, separator)
    flask_data = safe_get(FLASK_TYPES, flask_type)
    
    # Total cell dose needed (absolute number of cells)
    target_cells = desired_dose * weight * 1e6  # Convert from ×10⁶ to absolute cells
    
    # Calculate PBSC volume (mL) needed from separator yield
    try:
        pbsc_ml = target_cells / sep_yield
    except Exception:
        pbsc_ml = 50  # Fallback to 50 mL if any issue
    pbsc_ml = max(50, pbsc_ml)
    
    # Determine ideal seeding density (cells per flask) from flask surface
    initial_seeding = flask_data['surface'] * SEEDING_DENSITY  # cells seeded in one flask
    
    # --- Calculate yield per passage (discrete model) ---
    # Passage 1 (14 days): starting with 1 flask
    yield_p1 = flask_data['max_cells']  # one flask grows to max capacity
    
    # Passage 2 (days 15-28):
    # The cells from passage 1 are split into as many flasks as needed,
    # with each flask receiving the ideal seeding count.
    flasks_needed_p2 = math.ceil(yield_p1 / initial_seeding)
    yield_p2 = flasks_needed_p2 * flask_data['max_cells']
    
    # Passage 3 (days 29-42):
    flasks_needed_p3 = math.ceil(yield_p2 / initial_seeding)
    yield_p3 = flasks_needed_p3 * flask_data['max_cells']
    
    # --- Determine which passage yields the target dose (limit to 3 passages) ---
    if yield_p1 >= target_cells:
        passage_needed = 1
        final_yield = target_cells  # assume we harvest exactly target dose
    elif yield_p2 >= target_cells:
        passage_needed = 2
        final_yield = target_cells
    elif yield_p3 >= target_cells:
        passage_needed = 3
        final_yield = target_cells
    else:
        passage_needed = 3
        final_yield = yield_p3  # even after 3 passages, target not met; use best yield
    
    total_days = passage_needed * MIN_CULTURE_DAYS  # 14 days per passage
    
    # Final number of flasks in the final passage required to yield the target (or best yield)
    final_flasks = math.ceil(final_yield / flask_data['max_cells'])
    
    # --- Media calculation ---
    # Instead of using media_freq, we use the fixed schedule:
    # For each passage, media is changed on days 1, 4, 7, 11 (4 changes per 14-day block)
    # The number of flasks used in each passage is:
    # Passage 1: assume 1 flask
    # Passage 2: flasks_needed_p2
    # Passage 3: flasks_needed_p3
    media_p1 = 1 * flask_data['media'] * 4
    media_p2 = flasks_needed_p2 * flask_data['media'] * 4 if passage_needed >= 2 else 0
    media_p3 = flasks_needed_p3 * flask_data['media'] * 4 if passage_needed >= 3 else 0
    total_media = media_p1 + media_p2 + media_p3
    
    return {
        'pbsc_ml': pbsc_ml,
        'final_flasks': final_flasks,
        'passages': passage_needed,
        'total_days': total_days,
        'total_media': total_media,
        'initial_seeding': initial_seeding,
        'target_cells': target_cells,
        'yield_p1': yield_p1,
        'yield_p2': yield_p2,
        'yield_p3': yield_p3,
        'flasks_p2': flasks_needed_p2,
        'flasks_p3': flasks_needed_p3
    }

def plot_growth_curve(results):
    """
    Plot a stepwise growth curve using the discrete passage model.
    
    The x-axis marks the end of each passage (14, 28, 42 days).
    The y-axis shows the cell yield.
    If the target dose is reached before 3 passages, the yield is capped at the target.
    """
    # Define time points (in days) and corresponding yields
    time_points = [0, MIN_CULTURE_DAYS]  # start and end of passage 1
    yields = [results['initial_seeding'], results['yield_p1']]
    
    # Add passage 2 if applicable
    if results['passages'] >= 2:
        time_points.append(MIN_CULTURE_DAYS*2)
        yields.append(results['yield_p2'])
    # Add passage 3 if applicable
    if results['passages'] >= 3:
        time_points.append(MIN_CULTURE_DAYS*3)
        yields.append(results['yield_p3'])
    
    # If target is reached before the last passage, cap the yield
    if yields[-1] > results['target_cells']:
        yields[-1] = results['target_cells']
    
    # Create a stepwise plot
    fig, ax = plt.subplots()
    ax.step(time_points, yields, where='post', label='Cell Yield', linewidth=2)
    ax.plot(time_points, yields, 'o', color='red')
    ax.axhline(results['target_cells'], color='gray', linestyle='--', label='Target Cells')
    ax.set_xlabel("Culture Days")
    ax.set_ylabel("Cell Count (absolute)")
    ax.set_title("Growth Projection (Discrete Passage Model)")
    ax.grid(True)
    ax.legend()
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
        # Note: Media change frequency is fixed by the protocol (days 1,4,7,11), so this slider is not used.
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
        st.metric("Final Flask Count (Passage {0})".format(results['passages']), results['final_flasks'])
    with col2:
        st.metric("Passages", results['passages'])
        st.metric("Culture Duration", f"{results['total_days']} days")
    with col3:
        st.metric("Total Media", f"{results['total_media']:.0f} mL")
        st.metric("Recommended Dose", f"{grade_data['min_dose']}-{grade_data['max_dose']} ×10⁶/kg")
    
    # Growth curve plot
    st.header("Growth Projection")
    try:
        plot_fig = plot_growth_curve(results)
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
    - **Passage Limit:** Maximum of 3 passages. If the target dose is not met by Passage 3, proceed with cryopreservation or infusion of the best‐obtained dose.
    - **Media Schedule:** Media changes on days 1, 4, 7, and 11 of each 14‑day passage.
    - **Passage Timing:** First passage (trypsinization) is done at day 14.
    - **Flask Details:** 
         - Seeding density: {SEEDING_DENSITY:,} cells/cm² (per flask)
         - {flask_type} Flask with surface area of {safe_get(FLASK_TYPES, flask_type)['surface']} cm² and maximum capacity of {safe_get(FLASK_TYPES, flask_type)['max_cells']:,} cells.
    - **Expansion Estimates:**
         - Passage 1 yield: ~{results['yield_p1']:.0f} cells.
         - Passage 2 yield: ~{results['yield_p2']:.0f} cells (using ~{results['flasks_p2']} flask(s)).
         - Passage 3 yield: ~{results['yield_p3']:.0f} cells (using ~{results['flasks_p3']} flask(s)).
    - **Overall:** Culture duration is {results['total_days']} days (14 days per passage), and media consumption is estimated at {results['total_media']:.0f} mL.
    
    **Quality Control:**
    - Viability >90% (Trypan Blue)
    - CD73+/CD90+/CD105+ >95%
    - CD45- <2%
    - Sterility testing required
    """)
    
if __name__ == "__main__":
    main()
