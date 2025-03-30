import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math

# Constants
FLASK_TYPES = {
    'T-25': {'seeding': 0.7e6, 'confluent': 2.8e6, 'media': 5},
    'T-75': {'seeding': 2.1e6, 'confluent': 8.4e6, 'media': 15}, 
    'T-160': {'seeding': 4.6e6, 'confluent': 18.4e6, 'media': 30}
}

GVHD_RESPONSE = {
    'Grade I': {'min_dose': 0.5, 'max_dose': 1.0, 'response': [60, 80]},
    'Grade II': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]},
    'Grade III-IV': {'min_dose': 1.5, 'max_dose': 2.0, 'response': [40, 60]}
}

SEPARATOR_YIELD = {'Haemonetics': 1e6, 'Spectra Optia': 2e6}

SAFETY_FACTOR = 1.2
MAX_PASSAGES = 3
PASSAGE1_DAYS = 14  # Days for first passage
PASSAGE_DAYS = 7    # Days for subsequent passages
MEDIA_CHANGE_FREQ = 2  # Change media every 2 days
PLASMA_PERCENT = 0.15
MAX_PBSC = 50  # mL

def calculate_therapy(weight, dose, separator, flask_type, plasma_priming):
    # Validate inputs
    weight = max(8.0, min(weight, 120.0))
    dose = max(0.5, min(dose, 2.0))
    
    # PBSC calculation
    target_cells = dose * weight * 1e6
    pbsc_ml = min(target_cells / SEPARATOR_YIELD[separator], MAX_PBSC)
    
    # Passage calculations
    passages = []
    current_flasks = math.ceil(target_cells / FLASK_TYPES[flask_type]['confluent'])
    total_days = 0
    total_media = 0
    
    # Initial seeding (Passage 0)
    passages.append({
        'passage_num': 0,
        'flasks': current_flasks,
        'input': current_flasks * FLASK_TYPES[flask_type]['seeding'],
        'output': current_flasks * FLASK_TYPES[flask_type]['confluent'],
        'days': 0,
        'media_changes': 0
    })
    
    # Subsequent passages
    for passage_num in range(1, MAX_PASSAGES + 1):
        if passages[-1]['output'] >= target_cells:
            break
            
        # Calculate required flasks
        current_flasks = math.ceil(
            (passages[-1]['output'] / FLASK_TYPES[flask_type]['seeding']) * SAFETY_FACTOR
        )
        
        # Calculate passage duration and media changes
        days = PASSAGE1_DAYS if passage_num == 1 else PASSAGE_DAYS
        media_changes = (days // MEDIA_CHANGE_FREQ) + 1  # Include initial media
        
        passages.append({
            'passage_num': passage_num,
            'flasks': current_flasks,
            'input': current_flasks * FLASK_TYPES[flask_type]['seeding'],
            'output': current_flasks * FLASK_TYPES[flask_type]['confluent'],
            'days': days,
            'media_changes': media_changes
        })
        
        total_days += days
        total_media += current_flasks * FLASK_TYPES[flask_type]['media'] * media_changes
    
    # Plasma calculation
    plasma_vol = passages[0]['flasks'] * FLASK_TYPES[flask_type]['media'] * PLASMA_PERCENT if plasma_priming else 0
    
    return {
        'pbsc_ml': pbsc_ml,
        'passages': passages,
        'total_days': total_days,
        'total_media': total_media,
        'plasma_vol': plasma_vol,
        'target_cells': target_cells,
        'final_yield': passages[-1]['output']
    }

def plot_growth(passages, target_cells):
    fig, ax = plt.subplots(figsize=(10,6))
    
    # Build timeline
    days = [0]
    cells = [0]
    cumulative_days = 0
    
    for passage in passages[1:]:  # Skip seeding passage
        # Point after seeding
        days.append(cumulative_days)
        cells.append(passage['input'])
        
        # Point at confluency
        cumulative_days += passage['days']
        days.append(cumulative_days)
        cells.append(passage['output'])
    
    # Plot in millions
    ax.plot(days, [x/1e6 for x in cells], 'go-', markersize=8, linewidth=2)
    ax.axhline(target_cells/1e6, color='r', linestyle='--', 
              label=f'Target: {target_cells/1e6:.1f}×10⁶')
    
    ax.set_xlabel('Culture Days', fontsize=12)
    ax.set_ylabel('Cells (×10⁶)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, cumulative_days * 1.05)
    ax.set_ylim(0, max(max(cells), target_cells)*1.1/1e6)
    
    return fig

def plot_remission_probability(grade, dose):
    data = GVHD_RESPONSE[grade]
    x = np.linspace(0.5, 2.5, 100)
    opt = (data['min_dose'] + data['max_dose'])/2
    y = data['response'][0] + (data['response'][1] - data['response'][0]) * np.exp(-((x - opt)/0.3)**2
    
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(x, y, 'b-', linewidth=2)
    ax.axvline(dose, color='r', linestyle='--')
    ax.set_ylim(0, 100)
    ax.set_xlabel('Dose (×10⁶ cells/kg)', fontsize=12)
    ax.set_ylabel('Remission Probability (%)', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    return fig

def main():
    st.set_page_config(page_title="MSC Therapy Calculator", layout="wide")
    st.title("MSC Therapy Calculator")
    
    with st.sidebar:
        st.header("Patient Parameters")
        weight = st.slider("Weight (kg)", 8.0, 120.0, 70.0, 0.1, format="%.1f")
        grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        dose = st.slider("Dose (×10⁶/kg)", 0.5, 2.0, 1.0, 0.1)
        
        st.header("Lab Parameters")
        separator = st.selectbox("Cell Separator", list(SEPARATOR_YIELD.keys()))
        flask_type = st.selectbox("Flask Type", list(FLASK_TYPES.keys()))
        plasma_priming = st.checkbox("Plasma Priming (15% of initial media)")
    
    # Calculations
    results = calculate_therapy(weight, dose, separator, flask_type, plasma_priming)
    grade_data = GVHD_RESPONSE[grade]
    
    # Display results
    st.header("Therapy Parameters")
    cols = st.columns(4)
    cols[0].metric("PBSC Volume", f"{results['pbsc_ml']:.1f} mL")
    cols[1].metric("Initial Flasks", results['passages'][0]['flasks'])
    cols[2].metric("Total Media", f"{results['total_media']} mL")
    cols[3].metric("Culture Days", results['total_days'])
    
    if plasma_priming:
        st.info(f"**Plasma Volume Needed:** {results['plasma_vol']:.1f} mL")
    
    # Plots
    col1, col2 = st.columns(2)
    with col1:
        st.header("Cell Expansion")
        st.pyplot(plot_growth(results['passages'], results['target_cells']))
    
    with col2:
        st.header("GVHD Remission Probability")
        st.pyplot(plot_remission_probability(grade, dose))
    
    # Protocol details
    st.header("Protocol Summary")
    st.markdown(f"""
    **For {grade} GVHD:**
    - Recommended dose: {grade_data['min_dose']}-{grade_data['max_dose']}×10⁶/kg
    - Expected response: {grade_data['response'][0]}–{grade_data['response'][1]}%
    
    **Culture Parameters:**
    - Flask: {flask_type}
    - Seeding density: {FLASK_TYPES[flask_type]['seeding']/1e6:.1f}×10⁶ cells/flask
    - Confluent yield: {FLASK_TYPES[flask_type]['confluent']/1e6:.1f}×10⁶ cells/flask
    - Media changes: Every {MEDIA_CHANGE_FREQ} days
    
    **Passage Schedule:**
    - Passage 1: {PASSAGE1_DAYS} days
    - Subsequent passages: {PASSAGE_DAYS} days each
    - Total passages: {len(results['passages'])-1}
    
    **Final Yield:** {results['final_yield']/1e6:.1f}×10⁶ cells
    """)

if __name__ == "__main__":
    main()
