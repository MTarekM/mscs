import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math

# Flask data from reference table
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

# Constants
SAFETY_FACTOR = 1.2
MAX_PASSAGES = 3
PASSAGE1_DAYS = 14
PASSAGE_DAYS = 7
MEDIA_CHANGE_INTERVAL = 2
PLASMA_PERCENT = 0.15
MAX_PBSC = 50

def calculate_therapy(weight, dose, separator, flask_type, plasma_priming):
    target = dose * weight * 1e6
    pbsc = min(target / SEPARATOR_YIELD[separator], MAX_PBSC)
    
    # Passage calculations
    passages = []
    current_flasks = math.ceil(target / FLASK_TYPES[flask_type]['confluent'])
    total_days = 0
    total_media = 0
    
    # Initial seeding (Passage 0)
    passages.append({
        'flasks': current_flasks,
        'days': 0,
        'media_changes': 0,
        'output': current_flasks * FLASK_TYPES[flask_type]['confluent']
    })
    
    # Subsequent passages
    for passage_num in range(1, MAX_PASSAGES + 1):
        if passages[-1]['output'] >= target:
            break
            
        # Calculate required flasks
        current_flasks = math.ceil((passages[-1]['output'] / 
                                  FLASK_TYPES[flask_type]['seeding']) * 
                                  SAFETY_FACTOR)
        
        # Calculate passage duration and media changes
        days = PASSAGE1_DAYS if passage_num == 1 else PASSAGE_DAYS
        media_changes = (days // MEDIA_CHANGE_INTERVAL) + 1
        
        # Track cumulative totals
        total_days += days
        total_media += current_flasks * FLASK_TYPES[flask_type]['media'] * media_changes
        
        passages.append({
            'flasks': current_flasks,
            'days': days,
            'media_changes': media_changes,
            'output': current_flasks * FLASK_TYPES[flask_type]['confluent']
        })
    
    # Plasma calculation
    plasma_vol = passages[0]['flasks'] * FLASK_TYPES[flask_type]['media'] * PLASMA_PERCENT if plasma_priming else 0
    
    return {
        'pbsc_ml': pbsc,
        'total_days': total_days,
        'total_media': total_media,
        'plasma_vol': plasma_vol,
        'passages': passages,
        'target': target,
        'final_yield': passages[-1]['output']
    }

def plot_growth(results):
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Calculate timeline points
    days = [0]
    cells = [0]
    cumulative_days = 0
    
    for i, passage in enumerate(results['passages'][1:]):
        cumulative_days += passage['days']
        
        # Starting point (after previous passage)
        days.append(cumulative_days - passage['days'])
        cells.append(passage['output'] if i == 0 else cells[-1])
        
        # Ending point
        days.append(cumulative_days)
        cells.append(passage['output'])
    
    # Plot in millions
    ax.plot(days, [x/1e6 for x in cells], 'go-', markersize=8, linewidth=2)
    ax.axhline(results['target']/1e6, color='r', linestyle='--', 
              label=f'Target: {results["target"]/1e6:.1f}×10⁶')
    
    ax.set_xlabel('Culture Days', fontsize=12)
    ax.set_ylabel('Cells (×10⁶)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    return fig

def main():
    st.set_page_config(page_title="MSC Therapy Calculator", layout="wide")
    st.title("MSC Therapy Calculator")
    
    with st.sidebar:
        weight = st.slider("Weight (kg)", 8.0, 120.0, 8.0, 0.1)
        grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        dose = st.slider("Dose (×10⁶/kg)", 0.5, 2.0, 1.0, 0.1)
        separator = st.selectbox("Cell Separator", list(SEPARATOR_YIELD.keys()))
        flask_type = st.selectbox("Flask Type", list(FLASK_TYPES.keys()))
        plasma_priming = st.checkbox("Plasma Priming")
    
    results = calculate_therapy(weight, dose, separator, flask_type, plasma_priming)
    
    # Display results
    cols = st.columns(4)
    cols[0].metric("PBSC Volume", f"{results['pbsc_ml']:.1f} mL")
    cols[1].metric("Initial Flasks", results['passages'][0]['flasks'])
    cols[2].metric("Total Media", f"{results['total_media']} mL")
    cols[3].metric("Culture Days", results['total_days'])
    
    if plasma_priming:
        st.info(f"Plasma Volume Needed: {results['plasma_vol']:.1f} mL")
    
    st.pyplot(plot_growth(results))

if __name__ == "__main__":
    main()
