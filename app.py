"""
Streamlit web interface for Monte Carlo radiation shielding simulator.

This app provides an interactive interface for simulating gamma photon
transport through shielding materials and analyzing the results.
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import json
from pathlib import Path

from simulator import run_simulation
from analysis import (
    compute_transmission,
    compute_hvl,
    extract_interaction_histogram,
    extract_energy_spectrum,
    summarize_results,
    compare_to_theory
)
from physics import CrossSections

# Load configuration
CONFIG_PATH = Path(__file__).parent / 'config.json'
with open(CONFIG_PATH, 'r') as f:
    CONFIG = json.load(f)

# Page configuration
st.set_page_config(
    page_title="Radiation Shielding Simulator",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better appearance
st.markdown("""
    <style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 1rem;
    }
    .metric-container {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
    </style>
""", unsafe_allow_html=True)

# Title
st.markdown('<h1 class="main-header">🔬 Monte Carlo Radiation Shielding Simulator</h1>', 
            unsafe_allow_html=True)
st.markdown("*Gamma photon transport through lead and concrete*")
st.markdown("---")

# Sidebar - Simulation Parameters
st.sidebar.header("⚙️ Simulation Parameters")

# Material selection
material = st.sidebar.selectbox(
    "Material",
    options=list(CONFIG['materials'].keys()),
    help="Select shielding material"
)

material_info = CONFIG['materials'][material]
st.sidebar.info(f"**Density:** {material_info['density_g_cm3']} g/cm³  \n{material_info['description']}")

# Source selection
st.sidebar.markdown("---")
st.sidebar.subheader("Source Configuration")

source_preset = st.sidebar.selectbox(
    "Source Preset",
    options=list(CONFIG['sources'].keys()),
    help="Common gamma sources"
)

if source_preset == "Custom":
    E0_keV = st.sidebar.number_input(
        "Custom Energy (keV)",
        min_value=50,
        max_value=3000,
        value=500,
        step=50,
        help="Photon energy in keV"
    )
else:
    E0_keV = CONFIG['sources'][source_preset]
    st.sidebar.info(f"**Energy:** {E0_keV} keV")

# Number of photons
st.sidebar.markdown("---")
n_photons = st.sidebar.number_input(
    "Number of Photons",
    min_value=100,
    max_value=CONFIG['simulation']['max_photons'],
    value=CONFIG['simulation']['default_n_photons'],
    step=1000,
    help="More photons = better statistics but slower"
)

# Random seed
seed = st.sidebar.number_input(
    "Random Seed",
    min_value=0,
    max_value=99999,
    value=CONFIG['simulation']['default_seed'],
    step=1,
    help="Set seed for reproducibility"
)

# Analysis mode selection
st.sidebar.markdown("---")
st.sidebar.header("📊 Analysis Mode")

mode = st.sidebar.radio(
    "Select Mode",
    options=["Single Thickness", "Thickness Sweep"],
    help="Single: Quick analysis at one thickness\nSweep: Full transmission curve and HVL"
)

# Main content area
if mode == "Single Thickness":
    st.header("📏 Single Thickness Analysis")
    
    # Thickness input
    col1, col2 = st.columns([2, 1])
    with col1:
        thickness_cm = st.slider(
            "Slab Thickness (cm)",
            min_value=0.1,
            max_value=10.0,
            value=2.0,
            step=0.1,
            help="Thickness of shielding material"
        )
    
    with col2:
        st.metric("Thickness", f"{thickness_cm} cm")
    
    # Run button
    if st.button("▶️ Run Simulation", type="primary", use_container_width=True):
        
        # Progress indicator
        with st.spinner(f"Simulating {n_photons} photons through {thickness_cm} cm of {material}..."):
            results = run_simulation(n_photons, E0_keV, thickness_cm, material, seed)
        
        st.success("✓ Simulation complete!")
        
        # Summary statistics
        summary = summarize_results(results)
        
        # Display key metrics
        st.subheader("📈 Results")
        
        col1, col2, col3, col4 = st.columns(4)
        
        T, T_err = compute_transmission(results)
        
        with col1:
            st.metric(
                "Transmission",
                f"{T:.4f}",
                delta=f"± {T_err:.4f}",
                help="Fraction of photons transmitted"
            )
        
        with col2:
            st.metric(
                "Transmitted",
                f"{summary['n_transmitted']:,}",
                delta=f"{summary['transmission_fraction']*100:.1f}%",
                help="Number of photons that passed through"
            )
        
        with col3:
            st.metric(
                "Absorbed",
                f"{summary['n_absorbed']:,}",
                help="Number of photons absorbed"
            )
        
        with col4:
            st.metric(
                "Backscattered",
                f"{summary['n_backscattered']:,}",
                help="Number of photons scattered back"
            )
        
        # Additional statistics
        st.markdown("---")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric(
                "Total Interactions",
                f"{summary['total_interactions']:,}",
                help="Total number of interactions (all photons)"
            )
        
        with col2:
            st.metric(
                "Compton Scattering",
                f"{summary['total_compton']:,}",
                help="Number of Compton scattering events"
            )
        
        with col3:
            st.metric(
                "Photoelectric",
                f"{summary['total_photoelectric']:,}",
                help="Number of photoelectric absorption events"
            )
        
        # Interaction depth histogram
        st.markdown("---")
        st.subheader("📊 Interaction Depth Distribution")
        
        bins, hist_compton, hist_pe = extract_interaction_histogram(results, thickness_cm)
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        bin_centers = (bins[:-1] + bins[1:]) / 2.0
        width = bins[1] - bins[0]
        
        ax.bar(bin_centers, hist_compton, width=width*0.95, 
               label='Compton Scattering', alpha=0.7, color='steelblue')
        ax.bar(bin_centers, hist_pe, width=width*0.95, bottom=hist_compton,
               label='Photoelectric Absorption', alpha=0.7, color='coral')
        
        ax.set_xlabel("Depth (cm)", fontsize=12, fontweight='bold')
        ax.set_ylabel("Number of Interactions", fontsize=12, fontweight='bold')
        ax.set_title(f"Interaction Distribution: {material.title()}, E₀ = {E0_keV} keV, d = {thickness_cm} cm",
                    fontsize=14, fontweight='bold')
        ax.legend(fontsize=11, framealpha=0.9)
        ax.grid(alpha=0.3, linestyle='--')
        ax.set_xlim(0, thickness_cm)
        
        st.pyplot(fig)
        plt.close()
        
        # Energy spectrum of transmitted photons
        if summary['n_transmitted'] > 0:
            st.markdown("---")
            st.subheader("⚡ Transmitted Photon Energy Spectrum")
            
            transmitted_energies = [
                r.final_energy for r in results if r.outcome == "transmitted"
            ]
            
            fig2, ax2 = plt.subplots(figsize=(12, 6))
            
            ax2.hist(transmitted_energies, bins=50, alpha=0.7, 
                    edgecolor='black', color='mediumseagreen', label='Transmitted photons')
            ax2.axvline(E0_keV, color='red', linestyle='--', linewidth=2,
                       label=f'Source Energy ({E0_keV} keV)')
            
            if transmitted_energies:
                mean_E = np.mean(transmitted_energies)
                ax2.axvline(mean_E, color='orange', linestyle=':', linewidth=2,
                           label=f'Mean ({mean_E:.1f} keV)')
            
            ax2.set_xlabel("Energy (keV)", fontsize=12, fontweight='bold')
            ax2.set_ylabel("Count", fontsize=12, fontweight='bold')
            ax2.set_title("Energy Distribution of Transmitted Photons (Compton Downscatter)",
                         fontsize=14, fontweight='bold')
            ax2.legend(fontsize=11, framealpha=0.9)
            ax2.grid(alpha=0.3, linestyle='--')
            
            st.pyplot(fig2)
            plt.close()
            
            # Energy statistics
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Mean Energy", f"{summary['mean_transmitted_energy_keV']:.1f} keV")
            with col2:
                st.metric("Std Dev", f"{summary['std_transmitted_energy_keV']:.1f} keV")
            with col3:
                energy_loss_percent = (E0_keV - summary['mean_transmitted_energy_keV'])/E0_keV * 100
                st.metric("Avg Energy Loss", f"{energy_loss_percent:.1f}%")
        
        else:
            st.warning("⚠️ No photons transmitted - slab too thick or energy too low")

else:  # Thickness Sweep mode
    st.header("📈 Thickness Sweep Analysis")
    
    st.markdown("Generate transmission curve and extract Half-Value Layer (HVL)")
    
    # Thickness range inputs
    col1, col2, col3 = st.columns(3)
    
    with col1:
        d_min = st.number_input(
            "Min Thickness (cm)",
            min_value=0.1,
            max_value=20.0,
            value=0.5,
            step=0.5,
            help="Starting thickness"
        )
    
    with col2:
        d_max = st.number_input(
            "Max Thickness (cm)",
            min_value=0.1,
            max_value=20.0,
            value=5.0,
            step=0.5,
            help="Ending thickness"
        )
    
    with col3:
        n_points = st.number_input(
            "Number of Points",
            min_value=3,
            max_value=20,
            value=10,
            step=1,
            help="Points in thickness sweep"
        )
    
    # Validate inputs
    if d_min >= d_max:
        st.error("❌ Min thickness must be less than max thickness")
    
    # Run sweep button
    if st.button("▶️ Run Thickness Sweep", type="primary", use_container_width=True):
        
        if d_min >= d_max:
            st.error("Please fix thickness range")
        else:
            # Generate thickness array
            thicknesses = np.linspace(d_min, d_max, n_points)
            transmissions = []
            transmission_errors = []
            
            # Progress tracking
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            # Run simulations
            for i, d in enumerate(thicknesses):
                status_text.text(f"Simulating thickness {d:.2f} cm ({i+1}/{n_points})...")
                
                # Run simulation with unique seed for each thickness
                results = run_simulation(n_photons, E0_keV, d, material, seed + i)
                
                T, T_err = compute_transmission(results)
                transmissions.append(T)
                transmission_errors.append(T_err)
                
                progress_bar.progress((i + 1) / n_points)
            
            status_text.text("✓ All simulations complete!")
            st.success("Simulation sweep finished!")
            
            # Convert to arrays
            transmissions = np.array(transmissions)
            transmission_errors = np.array(transmission_errors)
            
            # Compute HVL
            hvl_data = compute_hvl(thicknesses, transmissions, transmission_errors)
            
            # Get theoretical values
            xsec = CrossSections(material)
            mu_theory = xsec.mu_total(E0_keV)
            hvl_theory = np.log(2) / mu_theory
            
            # Display HVL results
            st.markdown("---")
            st.subheader("🎯 Half-Value Layer (HVL) Results")
            
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                if hvl_data['hvl_interpolation']:
                    st.metric(
                        "HVL (Interpolation)",
                        f"{hvl_data['hvl_interpolation']:.3f} cm",
                        help="HVL by direct interpolation at T=0.5"
                    )
                else:
                    st.metric("HVL (Interpolation)", "N/A")
            
            with col2:
                if hvl_data['hvl_fit']:
                    st.metric(
                        "HVL (Exp Fit)",
                        f"{hvl_data['hvl_fit']:.3f} cm",
                        help="HVL from exponential fit"
                    )
                else:
                    st.metric("HVL (Exp Fit)", "N/A")
            
            with col3:
                st.metric(
                    "HVL (Theory)",
                    f"{hvl_theory:.3f} cm",
                    help="Theoretical narrow-beam HVL"
                )
            
            with col4:
                if hvl_data['mu_effective']:
                    st.metric(
                        "μ_eff",
                        f"{hvl_data['mu_effective']:.3f} cm⁻¹",
                        help="Effective attenuation coefficient"
                    )
                else:
                    st.metric("μ_eff", "N/A")
            
            # Comparison to theory
            if hvl_data['hvl_fit']:
                comparison = compare_to_theory(hvl_data['hvl_fit'], E0_keV, material, mu_theory)
                
                st.info(f"""
                **Comparison to Theory:**
                - Simulated HVL: {comparison['hvl_simulated_cm']:.3f} cm
                - Theoretical HVL: {comparison['hvl_theory_cm']:.3f} cm
                - Difference: {comparison['absolute_difference_cm']:.3f} cm ({comparison['relative_difference_percent']:.1f}%)
                
                *Note: Simulated HVL is typically larger than narrow-beam theory due to build-up (scattered photons).*
                """)
            
            if hvl_data['fit_r_squared'] is not None:
                st.metric("Fit Quality (R²)", f"{hvl_data['fit_r_squared']:.4f}")
            
            # Transmission curve plot
            st.markdown("---")
            st.subheader("📉 Transmission vs. Thickness")
            
            fig, ax = plt.subplots(figsize=(12, 7))
            
            # Simulation data with error bars
            ax.errorbar(thicknesses, transmissions, yerr=transmission_errors,
                       fmt='o', capsize=5, markersize=8, linewidth=2,
                       label='Simulation Data', color='darkblue', ecolor='gray')
            
            # Exponential fit
            if hvl_data['mu_effective']:
                d_fit = np.linspace(d_min, d_max, 200)
                T_fit = np.exp(-hvl_data['mu_effective'] * d_fit)
                ax.plot(d_fit, T_fit, '--', linewidth=2.5, 
                       label=f"Exponential Fit: T = exp(-{hvl_data['mu_effective']:.3f}·d)", 
                       color='red')
            
            # Theoretical curve
            d_theory = np.linspace(d_min, d_max, 200)
            T_theory = np.exp(-mu_theory * d_theory)
            ax.plot(d_theory, T_theory, ':', linewidth=2.5,
                   label=f"Theory (Narrow Beam): T = exp(-{mu_theory:.3f}·d)",
                   color='green')
            
            # HVL reference line
            ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1.5,
                      label='T = 0.5 (HVL)')
            
            if hvl_data['hvl_interpolation']:
                ax.axvline(hvl_data['hvl_interpolation'], color='orange', 
                          linestyle='-.', alpha=0.5, linewidth=1.5,
                          label=f"HVL = {hvl_data['hvl_interpolation']:.2f} cm")
            
            ax.set_xlabel("Thickness (cm)", fontsize=13, fontweight='bold')
            ax.set_ylabel("Transmission Fraction", fontsize=13, fontweight='bold')
            ax.set_title(f"Transmission Curve: {material.title()}, E₀ = {E0_keV} keV, N = {n_photons:,}",
                        fontsize=15, fontweight='bold')
            ax.legend(fontsize=10, loc='best', framealpha=0.9)
            ax.grid(alpha=0.3, linestyle='--')
            ax.set_ylim(0, 1.05)
            ax.set_xlim(d_min, d_max)
            
            st.pyplot(fig)
            plt.close()
            
            # Log-scale plot
            st.markdown("---")
            st.subheader("📊 Transmission (Logarithmic Scale)")
            
            fig2, ax2 = plt.subplots(figsize=(12, 7))
            
            # Filter out zero transmissions for log plot
            mask = transmissions > 0
            
            ax2.semilogy(thicknesses[mask], transmissions[mask], 
                        'o', markersize=8, label='Simulation Data', color='darkblue')
            
            if hvl_data['mu_effective']:
                ax2.semilogy(d_fit, T_fit, '--', linewidth=2.5, 
                           label='Exponential Fit', color='red')
            
            ax2.semilogy(d_theory, T_theory, ':', linewidth=2.5,
                        label='Theory (Narrow Beam)', color='green')
            
            ax2.axhline(0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)
            
            ax2.set_xlabel("Thickness (cm)", fontsize=13, fontweight='bold')
            ax2.set_ylabel("Transmission (log scale)", fontsize=13, fontweight='bold')
            ax2.set_title("Transmission Curve (Logarithmic Scale)", 
                         fontsize=15, fontweight='bold')
            ax2.legend(fontsize=10, loc='best', framealpha=0.9)
            ax2.grid(alpha=0.3, which='both', linestyle='--')
            ax2.set_xlim(d_min, d_max)
            
            st.pyplot(fig2)
            plt.close()
            
            # Download results
            st.markdown("---")
            st.subheader("💾 Export Results")
            
            # Create CSV data
            import io
            csv_buffer = io.StringIO()
            csv_buffer.write("thickness_cm,transmission,transmission_error\n")
            for d, T, T_err in zip(thicknesses, transmissions, transmission_errors):
                csv_buffer.write(f"{d:.3f},{T:.6f},{T_err:.6f}\n")
            
            csv_data = csv_buffer.getvalue()
            
            st.download_button(
                label="📥 Download Data (CSV)",
                data=csv_data,
                file_name=f"transmission_data_{material}_{E0_keV}keV.csv",
                mime="text/csv"
            )

# Footer
st.sidebar.markdown("---")
st.sidebar.markdown("### 📚 About")
st.sidebar.info("""
**Monte Carlo Simulator**

Physics Model:
- Energy-dependent cross-sections
- Klein-Nishina Compton scattering
- Photoelectric absorption
- Free-path exponential sampling

**Developer Info:**
- Built with Python & Streamlit
- Monte Carlo transport engine
- Scientific visualization
""")

st.sidebar.markdown("---")
st.sidebar.markdown("### ⚠️ Limitations")
st.sidebar.warning("""
- Simplified cross-sections (use NIST data for accuracy)
- No pair production
- No secondary particles
- Planar geometry only
""")

# Tips expander
with st.sidebar.expander("💡 Tips & Guidelines"):
    st.markdown("""
    **For best results:**
    
    1. **Start small**: Test with 1,000 photons first
    2. **Increase gradually**: Use 10,000+ for publication
    3. **Thickness range**: 
       - Lead: 0.5-5 cm
       - Concrete: 5-30 cm
    4. **HVL extraction**: Need at least 3 points with 0.1 < T < 0.9
    
    **Common energies:**
    - Cs-137: 662 keV (medical)
    - Co-60: 1250 keV (industrial)
    - Ir-192: 380 keV (NDT)
    """)
    