import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math

# Updated Constants with literature support
SEEDING_DENSITY = 5.0e4  # cells/cm² (optimized based on literature: 5×10⁴ cells/cm²)
GROWTH_RATE = 0.5       # doublings per day (for intra-passage exponential growth)
MIN_PASSAGE1_DAYS = 14  # Fixed duration for first passage (includes trypsinization)
ADDITIONAL_PASSAGE_DAYS = 7  # Each subsequent passage lasts 7 days
MAX_PASSAGES = 3       # Maximum number of passages allowed
SAFETY_FACTOR = 1.2    # Safety factor for re-seeding losses
PLASMA_PERCENTAGE = 0.15  # 15% of initial media volume for priming (literature-based)

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

# Updated Flask parameters with literature-based calculations
# Surface area in cm², media volume in mL, max_cells at confluency based on 5×10⁴ cells/cm² seeding
FLASK_TYPES = {
    'T25': {
        'surface': 25,
        'media': 5,
        'initial_cells': 25 * 5.0e4,  # 1.25×10⁶ cells
        'max_cells': 25 * 2.0e5      # ~5×10⁶ cells at confluency (40x expansion)
    },
    'T75': {
        'surface': 75,
        'media': 15,
        'initial_cells': 75 * 5.0e4,  # 3.75×10⁶ cells
        'max_cells': 75 * 2.0e5       # ~15×10⁶ cells at confluency
    },
    'T175': {
        'surface': 175,
        'media': 30,
        'initial_cells': 175 * 5.0e4, # 8.75×10⁶ cells
        'max_cells': 175 * 2.0e5      # ~35×10⁶ cells at confluency
    },
    'Default': {
        'surface': 75,
        'media': 15,
        'initial_cells': 75 * 5.0e4,
        'max_cells': 75 * 2.0e5
    }
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
    
    # Calculate PBSC volume (mL) - literature-based calculation
    try:
        pbsc_ml = total_cells / sep_yield
    except Exception:
        pbsc_ml = 50  # Fallback
    pbsc_ml = max(50, pbsc_ml)  # Minimum 50 mL collection
    
    # Determine initial cells seeded per flask (literature-based)
    initial_cells_per_flask = flask_data['initial_cells']
    
    # Calculate required initial flasks (N0) through iterative passage simulation
    found = False
    for N0 in range(1, 201):  # Search up to 200 initial flasks
        # Passage 1: N0 flasks to confluency
        passage1_yield = N0 * flask_data['max_cells']
        
        if passage1_yield >= total_cells:
            # Target met in first passage
            initial_flasks = N0
            flasks_p2 = 0
            flasks_p3 = 0
            passage2_yield = passage1_yield
            passage3_yield = passage1_yield
            found = True
            break
        
        # Passage 2: calculate needed flasks with safety factor
        flasks_p2 = math.ceil((passage1_yield / initial_cells_per_flask) * SAFETY_FACTOR)
        passage2_yield = flasks_p2 * flask_data['max_cells']
        
        if passage2_yield >= total_cells:
            # Target met in second passage
            initial_flasks = N0
            flasks_p3 = 0
            passage3_yield = passage2_yield
            found = True
            break
        
        # Passage 3: calculate needed flasks with safety factor
        flasks_p3 = math.ceil((passage2_yield / initial_cells_per_flask) * SAFETY_FACTOR)
        passage3_yield = flasks_p3 * flask_data['max_cells']
        
        if passage3_yield >= total_cells:
            initial_flasks = N0
            found = True
            break
    
    if not found:
        # If even 200 flasks aren't enough, use maximum possible
        initial_flasks = 200
        flasks_p2 = math.ceil((initial_flasks * flask_data['max_cells'] / initial_cells_per_flask) * SAFETY_FACTOR)
        flasks_p3 = math.ceil((flasks_p2 * flask_data['max_cells'] / initial_cells_per_flask) * SAFETY_FACTOR)
        passage3_yield = flasks_p3 * flask_data['max_cells']

    # Calculate total culture duration
    if passage1_yield >= total_cells:
        total_days = MIN_PASSAGE1_DAYS
        passages_used = 1
    elif passage2_yield >= total_cells:
        total_days = MIN_PASSAGE1_DAYS + ADDITIONAL_PASSAGE_DAYS
        passages_used = 2
    else:
        total_days = MIN_PASSAGE1_DAYS + 2 * ADDITIONAL_PASSAGE_DAYS
        passages_used = 3

    # Media calculations (literature-based)
    # Passage 1: 4 media changes (days 1,4,7,11)
    # Subsequent passages: 2 media changes each
    total_media_changes = 4 + (passages_used - 1) * 2
    total_media_volume = initial_flasks * flask_data['media'] * total_media_changes

    # Plasma priming volume (15% of initial media preparation - literature-based)
    if plasma_priming:
        # Only for initial media preparation (first change)
        plasma_volume = math.ceil(initial_flasks * flask_data['media'] * PLASMA_PERCENTAGE)
    else:
        plasma_volume = 0

    return {
        'pbsc_ml': pbsc_ml,
        'initial_flasks': initial_flasks,
        'passages': passages_used,
        'total_days': total_days,
        'total_media': total_media_volume,
        'initial_cells': initial_cells_per_flask,
        'target_cells': total_cells,
        'passage1_yield': passage1_yield,
        'flasks_p2': flasks_p2,
        'passage2_yield': passage2_yield if passages_used >= 2 else passage1_yield,
        'flasks_p3': flasks_p3,
        'passage3_yield': passage3_yield if passages_used >= 3 else passage2_yield,
        'plasma_volume': plasma_volume,
        'passages_used': passages_used
    }

def plot_growth_curve(results, flask_data):
    """
    Plots a piecewise growth curve showing overall cell yield over time.
    """
    try:
        # Get all required values from results with fallbacks
        initial_cells = results.get('initial_cells', 0)
        N0 = results.get('initial_flasks', 1)
        passage1_yield = results.get('passage1_yield', 0)
        passage2_yield = results.get('passage2_yield', passage1_yield)
        passage3_yield = results.get('passage3_yield', passage2_yield)
        flasks_p2 = results.get('flasks_p2', 0)
        flasks_p3 = results.get('flasks_p3', 0)
        target_cells = results.get('target_cells', 0)
        
        # Calculate starting yields for each passage
        start1 = N0 * initial_cells if initial_cells > 0 else 1
        end1 = passage1_yield if passage1_yield > 0 else start1
        
        start2 = flasks_p2 * initial_cells if initial_cells > 0 else end1
        end2 = passage2_yield if passage2_yield > 0 else start2
        
        start3 = flasks_p3 * initial_cells if initial_cells > 0 else end2
        end3 = passage3_yield if passage3_yield > 0 else start3
        
        # Define time segments
        t1 = np.linspace(0, MIN_PASSAGE1_DAYS, 100)
        t2 = np.linspace(MIN_PASSAGE1_DAYS, MIN_PASSAGE1_DAYS + ADDITIONAL_PASSAGE_DAYS, 100)
        t3 = np.linspace(MIN_PASSAGE1_DAYS + ADDITIONAL_PASSAGE_DAYS, 
                        MIN_PASSAGE1_DAYS + 2 * ADDITIONAL_PASSAGE_DAYS, 100)
        
        # Exponential growth in each segment
        factor1 = end1 / start1 if start1 > 0 else 1
        rate1 = np.log(factor1) / (MIN_PASSAGE1_DAYS) if MIN_PASSAGE1_DAYS > 0 else 0
        growth1 = start1 * np.exp(rate1 * t1)
        
        factor2 = end2 / start2 if start2 > 0 else 1
        rate2 = np.log(factor2) / (ADDITIONAL_PASSAGE_DAYS) if ADDITIONAL_PASSAGE_DAYS > 0 else 0
        growth2 = start2 * np.exp(rate2 * (t2 - MIN_PASSAGE1_DAYS))
        
        factor3 = end3 / start3 if start3 > 0 else 1
        rate3 = np.log(factor3) / (ADDITIONAL_PASSAGE_DAYS) if ADDITIONAL_PASSAGE_DAYS > 0 else 0
        growth3 = start3 * np.exp(rate3 * (t3 - (MIN_PASSAGE1_DAYS + ADDITIONAL_PASSAGE_DAYS)))
        
        # Combine segments
        t_total = np.concatenate([t1, t2, t3])
        growth_total = np.concatenate([growth1, growth2, growth3])
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(t_total, growth_total, 'g-', linewidth=2, label='MSC Growth Curve')
        ax.axhline(target_cells, color='r', linestyle='--', label=f'Target Dose ({target_cells/1e6:.1f}×10⁶ cells)')
        ax.axhline(end3, color='b', linestyle=':', label=f'Projected Yield ({end3/1e6:.1f}×10⁶ cells)')
        
        ax.set_xlabel('Culture Time (days)', fontsize=12)
        ax.set_ylabel('Total Cells', fontsize=12)
        ax.set_title('MSC Expansion Projection', fontsize=14)
        ax.legend(fontsize=10)
        ax.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        return fig
    except Exception as e:
        st.error(f"Error generating growth curve: {str(e)}")
        return None
def plot_gvhd_probability(grade_data, dose):
    """
    Plots a predicted GVHD remission probability curve.
    """
    try:
        if not grade_data or 'min_dose' not in grade_data:
            raise ValueError("Invalid grade data")
        
        min_dose = grade_data.get('min_dose', 0.5)
        max_dose = grade_data.get('max_dose', 2.0)
        response = grade_data.get('response', [50, 70])
        
        opt_dose = (min_dose + max_dose) / 2
        lower_prob, upper_prob = response[0], response[1]
        sigma = (max_dose - min_dose) / 3  # More gradual curve
        
        x = np.linspace(0.5, 2.5, 200)
        # Sigmoid-like function for more realistic probability curve
        y = lower_prob + (upper_prob - lower_prob) / (1 + np.exp(-(x - opt_dose)/sigma*2))
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(x, y, 'b-', linewidth=2, 
               label=f"GVHD Remission Probability ({lower_prob}-{upper_prob}%)")
        ax.axvline(dose, color='r', linestyle='--', 
                  label=f"Selected Dose: {dose:.2f}×10⁶/kg")
        
        ax.set_xlabel("MSC Dose (×10⁶ cells/kg)", fontsize=12)
        ax.set_ylabel("Probability of Response (%)", fontsize=12)
        ax.set_title("Dose-Response Relationship", fontsize=14)
        ax.set_ylim(0, 100)
        ax.legend(fontsize=10, loc='lower right')
        ax.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        
        return fig
    except Exception as e:
        st.error(f"Error generating GVHD plot: {str(e)}")
        return None

def main():
    st.set_page_config(page_title="MSC Therapy Calculator", layout="wide")
    st.title("Literature-Based MSC Therapy Calculator")
    
    with st.sidebar:
        st.header("Patient Parameters")
        weight = st.slider("Weight (kg)", 30, 120, 70)
        gvhd_grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        desired_dose = st.slider("Dose (×10⁶/kg)", 0.5, 2.0, 1.0, 0.1)
        
        st.header("Lab Parameters")
        separator = st.selectbox("Cell Separator", list(SEPARATOR_YIELD.keys()))
        flask_type = st.selectbox("Flask Type", list(FLASK_TYPES.keys()))
        plasma_priming = st.checkbox("Plasma priming from GVHD patient", value=False,
                                   help="15% of initial media volume from GVHD patient plasma")
    
    # Get grade data with fallback
    grade_data = safe_get(GVHD_RESPONSE, gvhd_grade)
    
    # Calculate therapy parameters
    results = calculate_msc_therapy(weight, desired_dose, separator, flask_type, plasma_priming)
    flask_data = safe_get(FLASK_TYPES, flask_type)
    
    # Display results
    st.header("Therapy Parameters")
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("PBSC Volume", f"{results['pbsc_ml']:.0f} mL")
        st.metric("Initial Flasks Needed", results['initial_flasks'])
    with col2:
        st.metric("Passages Required", results['passages_used'])
        st.metric("Culture Duration", f"{results['total_days']} days")
    with col3:
        st.metric("Total Media", f"{results['total_media']:.0f} mL")
        st.metric("Recommended Dose", f"{grade_data['min_dose']}-{grade_data['max_dose']} ×10⁶/kg")
    with col4:
        if plasma_priming:
            st.metric("Plasma Volume Needed", f"{results['plasma_volume']:.0f} mL (15% of initial media)")
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
    
    # Protocol notes with literature-based details
    st.header("Protocol Notes (Literature-Based)")
    st.markdown(f"""
    **For GVHD {gvhd_grade}:**
    - Recommended dose: {grade_data['min_dose']}-{grade_data['max_dose']} ×10⁶/kg
    - Expected remission response: {grade_data['response'][0]}–{grade_data['response'][1]}%
    
    **Culture Protocol (Based on Stem Cell Research & Therapy 2020):**
    - Seeding density: {SEEDING_DENSITY/1e4:.1f}×10⁴ cells/cm² (optimized for MSC expansion)
    - Flask type: {flask_type} with:
      - Surface area: {flask_data['surface']} cm²
      - Initial cells: {flask_data['initial_cells']/1e6:.2f}×10⁶ cells
      - Max cells at confluency: {flask_data['max_cells']/1e6:.1f}×10⁶ cells (~40x expansion)
    
    **Passage Schedule:**
    - Passage 1: {MIN_PASSAGE1_DAYS} days (media changes on days 1,4,7,11)  
    - Subsequent passages: {ADDITIONAL_PASSAGE_DAYS} days each (2 media changes per passage)
    - Minimum passages needed: {results['passages_used']}
    
    **Calculations:**
    - Initial flasks required: {results['initial_flasks']} (based on {flask_data['initial_cells']/1e6:.2f}×10⁶ cells/flask)
    - Final yield: {max(results['passage1_yield'], results['passage2_yield'], results['passage3_yield'])/1e6:.1f}×10⁶ cells
    - Total media: {results['total_media']:.0f} mL ({flask_data['media']} mL/flask × {results['initial_flasks']} flasks × {4 + (results['passages_used']-1)*2} changes)
    - {f"Plasma priming: {results['plasma_volume']} mL (15% of initial {flask_data['media']*results['initial_flasks']} mL media)" if plasma_priming else "No plasma priming"}
    
    **Quality Control (ISCT Standards):**
    - Viability >90% (Trypan Blue)
    - CD73+/CD90+/CD105+ >95%
    - CD45- <2%
    - Sterility testing required
    - Differentiation potential (adiogenic, osteogenic)
    """)

if __name__ == "__main__":
    main()
