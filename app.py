import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

# Constants
GROWTH_RATE = 0.5  # doublings per day
CONFLUENCY_DENSITY = 15000  # cells/cm² at full confluency
HARVEST_CONFLUENCY = 0.8  # 80% confluency for harvest

# Flask specifications
FLASK_TYPES = {
    'T25': {
        'surface': 25,  # cm²
        'media_volume': (5, 10),  # mL range
        'seeding_density': 3000  # cells/cm²
    },
    'T75': {
        'surface': 75,
        'media_volume': (15, 20),
        'seeding_density': 3000
    },
    'T175': {
        'surface': 175,
        'media_volume': (30, 40),
        'seeding_density': 2500
    }
}

# GVHD response parameters
GVHD_RESPONSE = {
    'Grade I': {'min_dose': 0.5, 'max_dose': 1.0, 'response': [60, 80]},
    'Grade II': {'min_dose': 1.0, 'max_dose': 1.5, 'response': [50, 70]},
    'Grade III-IV': {'min_dose': 1.5, 'max_dose': 2.0, 'response': [40, 60]}
}

def calculate_msc_therapy(weight, desired_dose, flask_type, media_volume, media_freq):
    # Calculate total MSCs needed (×10⁶ cells)
    total_cells_needed = desired_dose * weight  # ×10⁶ cells
    total_cells = total_cells_needed * 1e6  # Convert to absolute cell count
    
    # Get flask parameters
    flask = FLASK_TYPES[flask_type]
    surface_area = flask['surface']
    seeding_density = flask['seeding_density']
    
    # Calculate initial cells per flask
    initial_cells = surface_area * seeding_density  # cells
    
    # Calculate maximum cells per flask at harvest confluency
    max_cells_per_flask = surface_area * CONFLUENCY_DENSITY * HARVEST_CONFLUENCY
    
    # Calculate passages needed
    passages = 0
    current_cells = initial_cells
    while current_cells < total_cells:
        passages += 1
        current_cells = min(current_cells * 5, max_cells_per_flask)  # 5-fold expansion
    
    # Calculate culture duration (days)
    days_per_passage = np.log(max_cells_per_flask/initial_cells) / np.log(2) / GROWTH_RATE
    culture_days = max(3, int(np.ceil(days_per_passage * (passages + 1))))
    
    # Calculate PBSC volume needed (based on initial cell yield)
    # Assuming 1×10⁶ MSCs per 50mL of PBSC
    pbsc_volume = max(0.05, (total_cells_needed / 1.0) * 0.05)  # in liters
    
    # Calculate flasks needed (always start with 1 flask)
    flasks = 1
    
    # Calculate total media needed (mL)
    total_media = media_volume * (culture_days // media_freq) * flasks * (passages + 1)
    
    return {
        'pbsc_volume': pbsc_volume,
        'flasks': flasks,
        'passages': passages,
        'culture_days': culture_days,
        'total_media': total_media,
        'initial_cells': initial_cells / 1e6,  # Convert back to ×10⁶ cells
        'target_cells': total_cells_needed,
        'max_cells_per_flask': max_cells_per_flask / 1e6
    }

def plot_growth_curve(initial_cells, target_cells, max_flask_cells, days):
    # Ensure valid parameters
    days = max(1, days)
    x = np.linspace(0, days, 100)
    
    # Calculate growth per passage
    passages = int(np.ceil(np.log(target_cells/initial_cells) / np.log(5)))
    days_per_passage = days / (passages + 1)
    
    # Create growth curve
    growth = []
    for day in x:
        passage = int(day // days_per_passage)
        growth.append(min(initial_cells * (5 ** passage) * np.exp(GROWTH_RATE * (day % days_per_passage)), 
                      max_flask_cells))
    
    fig, ax = plt.subplots()
    ax.plot(x, growth, 'g-', linewidth=2)
    ax.axhline(y=target_cells, color='r', linestyle='--', label='Target Dose')
    ax.set_xlabel('Culture Days')
    ax.set_ylabel('MSC Count (×10⁶)')
    ax.set_title('MSC Growth Curve')
    ax.legend()
    ax.grid(True, alpha=0.3)
    return fig

def plot_dose_response(grade, desired_dose):
    grade_data = GVHD_RESPONSE[grade]
    doses = np.linspace(grade_data['min_dose'], grade_data['max_dose'], 100)
    response = np.interp(doses, 
                        [grade_data['min_dose'], grade_data['max_dose']],
                        grade_data['response'])
    
    fig, ax = plt.subplots()
    ax.plot(doses, response, 'b-', linewidth=2)
    ax.axvline(x=desired_dose, color='r', linestyle='--', label='Desired Dose')
    ax.set_xlabel('Dose (×10⁶ MSCs/kg)')
    ax.set_ylabel('Response Probability (%)')
    ax.set_title(f'GVHD {grade} Response Curve')
    ax.legend()
    ax.grid(True, alpha=0.3)
    return fig

def main():
    st.set_page_config(page_title="MSC Therapy Calculator", layout="wide")
    st.title("Advanced MSC Therapy Calculator for GVHD")
    
    # Inputs
    with st.sidebar:
        st.header("Patient Parameters")
        weight = st.number_input("Weight (kg)", min_value=30, max_value=120, value=70)
        gvhd_grade = st.selectbox("GVHD Grade", list(GVHD_RESPONSE.keys()))
        desired_dose = st.slider("Desired Dose (×10⁶ MSCs/kg)", 
                               min_value=0.5, max_value=2.0, value=1.0, step=0.1)
        
        st.header("Culture Parameters")
        flask_type = st.selectbox("Flask Type", list(FLASK_TYPES.keys()))
        media_volume = st.slider(f"Media Volume per {flask_type} (mL)",
                               min_value=FLASK_TYPES[flask_type]['media_volume'][0],
                               max_value=FLASK_TYPES[flask_type]['media_volume'][1],
                               value=FLASK_TYPES[flask_type]['media_volume'][0])
        media_freq = st.slider("Media Change Frequency (days)", 
                              min_value=1, max_value=4, value=2)
    
    # Calculations
    results = calculate_msc_therapy(weight, desired_dose, flask_type, media_volume, media_freq)
    recommended_dose = GVHD_RESPONSE[gvhd_grade]
    
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
        st.metric("Total Media Needed", f"{results['total_media']:.0f} mL")
        st.metric("Recommended Dose Range", 
                 f"{recommended_dose['min_dose']}-{recommended_dose['max_dose']} ×10⁶/kg")
    
    # Plots
    st.header("Biological Projections")
    col1, col2 = st.columns(2)
    
    with col1:
        st.pyplot(plot_growth_curve(results['initial_cells'], 
                                  results['target_cells'],
                                  results['max_cells_per_flask'],
                                  results['culture_days']))
    
    with col2:
        st.pyplot(plot_dose_response(gvhd_grade, desired_dose))
    
    # Protocol Notes
    st.header("Protocol Details")
    st.markdown(f"""
    **Culture Specifications:**
    - Seeding density: {FLASK_TYPES[flask_type]['seeding_density']} cells/cm²
    - Harvest confluency: {HARVEST_CONFLUENCY*100:.0f}%
    - Maximum cells per {flask_type}: {results['max_cells_per_flask']:.1f}×10⁶
    - Media changes: Every {media_freq} days
    
    **Quality Control:**
    - Viability >90% (Trypan Blue)
    - CD73+/CD90+/CD105+ >95%
    - CD45- <2%
    - Sterility testing required
    """)

if __name__ == "__main__":
    main()
