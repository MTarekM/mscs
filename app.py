import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import math

# Constants from reference tables
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

# Protocol constants
SAFETY_FACTOR = 1.2
MAX_PASSAGES = 3
PASSAGE0_DAYS = 2  # Initial seeding takes time too
PASSAGE1_DAYS = 14
PASSAGE_DAYS = 7
PLASMA_PERCENT = 0.15
MAX_PBSC = 50  # mL

# Fixed media change schedule for first passage (days 1, 4, 7, 11, 14)
PASSAGE1_MEDIA_CHANGES = [1, 4, 7, 11, 14]  # Day numbers when media is changed

def calculate_media_changes(days, media_freq, is_first_passage=False):
    """Calculate media changes based on days and frequency"""
    if is_first_passage:
        # For first passage, use the fixed schedule
        # Return the number of media changes (excluding day 0 which is initial media)
        return len([day for day in PASSAGE1_MEDIA_CHANGES if day <= days])
    else:
        # For other passages, change every media_freq days
        # Calculate how many times we'll change media during the culture period
        changes = [i for i in range(1, days+1) if i % media_freq == 0]
        return len(changes)

def calculate_therapy(weight, dose, separator, flask_type, plasma_priming, media_freq):
    # Validate inputs
    weight = max(8.0, min(weight, 120.0))
    dose = max(0.5, min(dose, 2.0))
    
    # Calculate target cells and PBSC volume
    target_cells = dose * weight * 1e6
    pbsc_ml = min(target_cells / SEPARATOR_YIELD[separator], MAX_PBSC)
    
    # Initialize passage tracking
    passages = []
    total_days = 0
    total_media = 0
    
    # Initial seeding (Passage 0)
    initial_flasks = math.ceil(target_cells / FLASK_TYPES[flask_type]['confluent'])
    
    # Initial passage also takes time and needs media
    p0_days = PASSAGE0_DAYS
    p0_media_changes = 1  # Initial media only for seeding
    
    passages.append({
        'passage_num': 0,
        'flasks': initial_flasks,
        'input': pbsc_ml * SEPARATOR_YIELD[separator],  # Starting cell count from PBSC
        'output': initial_flasks * FLASK_TYPES[flask_type]['seeding'],
        'days': p0_days,
        'media_changes': p0_media_changes,
        'media_schedule': [0]  # Day 0 (initial)
    })
    
    total_days += p0_days
    total_media += initial_flasks * FLASK_TYPES[flask_type]['media'] * p0_media_changes
    
    # Subsequent passages
    for passage_num in range(1, MAX_PASSAGES + 1):
        if passages[-1]['output'] >= target_cells:
            break
            
        current_flasks = math.ceil(
            (passages[-1]['output'] / FLASK_TYPES[flask_type]['seeding']) * SAFETY_FACTOR
        )
        
        days = PASSAGE1_DAYS if passage_num == 1 else PASSAGE_DAYS
        
        # Calculate media changes based on our schedule
        is_first = (passage_num == 1)
        media_changes = calculate_media_changes(days, media_freq, is_first)
        
        # Ensure at least one media change
        if media_changes == 0:
            media_changes = 1
            
        # Generate media change schedule
        if is_first:
            # First passage has fixed schedule
            media_schedule = [day for day in PASSAGE1_MEDIA_CHANGES if day <= days]
            media_schedule.insert(0, 0)  # Add initial media (day 0)
        else:
            # Other passages change every media_freq days
            media_schedule = [0]  # Initial media
            media_schedule.extend([i for i in range(1, days+1) if i % media_freq == 0])
        
        passages.append({
            'passage_num': passage_num,
            'flasks': current_flasks,
            'input': current_flasks * FLASK_TYPES[flask_type]['seeding'],
            'output': current_flasks * FLASK_TYPES[flask_type]['confluent'],
            'days': days,
            'media_changes': media_changes + 1,  # +1 for initial media
            'media_schedule': media_schedule
        })
        
        total_days += days
        total_media += current_flasks * FLASK_TYPES[flask_type]['media'] * (media_changes + 1)
    
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
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Build complete timeline including initial seeding
    days = [0]
    cells = [passages[0]['input']]  # Start with initial PBSC cells
    cumulative_days = 0
    
    # Add all passages to timeline
    for passage in passages:
        cumulative_days += passage['days']
        days.append(cumulative_days)
        cells.append(passage['output'])
    
    ax.plot(days, [x/1e6 for x in cells], 'go-', markersize=8, linewidth=2)
    ax.axhline(target_cells/1e6, color='r', linestyle='--', label=f'Target: {target_cells/1e6:.1f}×10⁶')
    
    # Add annotations for passages
    for i, passage in enumerate(passages):
        day_point = sum(p['days'] for p in passages[:i+1])
        ax.annotate(f"P{passage['passage_num']}", 
                   (day_point, passage['output']/1e6),
                   xytext=(5, 5), textcoords='offset points')
    
    # Add media change markers
    cumulative_days = 0
    for passage in passages:
        for change_day in passage['media_schedule'][1:]:  # Skip initial media (day 0)
            media_day = cumulative_days + change_day
            ax.axvline(media_day, color='gray', linestyle=':', alpha=0.5)
        cumulative_days += passage['days']
    
    ax.set_xlabel('Culture Days', fontsize=12)
    ax.set_ylabel('Cells (×10⁶)', fontsize=12)
    ax.set_title('MSC Expansion Timeline', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, cumulative_days * 1.05)
    ax.set_ylim(0, max(max(cells), target_cells)*1.1/1e6)
    
    return fig

def plot_remission_probability(grade, dose):
    data = GVHD_RESPONSE[grade]
    x = np.linspace(0.5, 2.5, 100)
    opt = (data['min_dose'] + data['max_dose'])/2
    y = data['response'][0] + (data['response'][1] - data['response'][0]) * np.exp(-((x - opt)/0.3)**2)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x, y, 'b-', linewidth=2, label='Remission Probability')
    ax.axvline(dose, color='r', linestyle='--', label=f'Selected Dose: {dose}×10⁶/kg')
    
    # Calculate actual probability at selected dose
    selected_prob = data['response'][0] + (data['response'][1] - data['response'][0]) * np.exp(-((dose - opt)/0.3)**2)
    ax.plot(dose, selected_prob, 'ro', markersize=8)
    ax.annotate(f"{selected_prob:.1f}%", (dose, selected_prob), xytext=(5, 5), textcoords='offset points')
    
    ax.set_ylim(0, 100)
    ax.set_xlabel('Dose (×10⁶ cells/kg)', fontsize=12)
    ax.set_ylabel('Probability (%)', fontsize=12)
    ax.set_title(f'GVHD {grade} Remission Probability', fontsize=14)
    ax.legend()
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
        media_freq = st.slider("Media Change Frequency (days)", 1, 7, 2, help="For passages after first. First passage uses fixed schedule.")
        
        st.info("First passage media changes on days 1, 4, 7, 11, and 14")
    
    results = calculate_therapy(weight, dose, separator, flask_type, plasma_priming, media_freq)
    
    st.header("Therapy Parameters")
    cols = st.columns(4)
    cols[0].metric("PBSC Volume", f"{results['pbsc_ml']:.1f} mL")
    cols[1].metric("Initial Flasks", results['passages'][0]['flasks'])
    cols[2].metric("Total Media", f"{results['total_media']} mL")
    cols[3].metric("Culture Days", results['total_days'])
    
    if plasma_priming:
        st.info(f"**Plasma Volume Needed:** {results['plasma_vol']:.1f} mL")
    
    # Display passages details
    st.subheader("Passage Details")
    passage_data = []
    for p in results['passages']:
        media_days = ", ".join([f"Day {d}" for d in p['media_schedule']])
        passage_data.append({
            "Passage": p['passage_num'],
            "Flasks": p['flasks'],
            "Input Cells": f"{p['input']/1e6:.1f}×10⁶",
            "Output Cells": f"{p['output']/1e6:.1f}×10⁶",
            "Days": p['days'],
            "Media Changes": p['media_changes'],
            "Media Schedule": media_days,
            "Media Volume": f"{p['flasks'] * FLASK_TYPES[flask_type]['media'] * p['media_changes']} mL"
        })
    st.table(passage_data)
    
    col1, col2 = st.columns(2)
    with col1:
        st.pyplot(plot_growth(results['passages'], results['target_cells']))
    
    with col2:
        st.pyplot(plot_remission_probability(grade, dose))
    
    st.subheader("Summary")
    st.markdown(f"""
    - **Target dose**: {dose}×10⁶ cells/kg for a {weight} kg patient ({results['target_cells']/1e6:.1f}×10⁶ total cells)
    - **Final yield**: {results['final_yield']/1e6:.1f}×10⁶ cells ({(results['final_yield']/results['target_cells'])*100:.1f}% of target)
    - **Total culture time**: {results['total_days']} days
    - **Total media required**: {results['total_media']} mL
    """)

if __name__ == "__main__":
    main()
