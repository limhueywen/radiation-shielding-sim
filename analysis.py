"""
Analysis module for Monte Carlo simulation results.

Contains functions for:
- Computing transmission fractions
- Extracting half-value layers (HVL)
- Creating interaction depth histograms
- Statistical analysis and uncertainty quantification
"""

import numpy as np
import json
from typing import List, Tuple, Dict, Optional
from pathlib import Path
from simulator import SimulationResult

# Load configuration
CONFIG_PATH = Path(__file__).parent / 'config.json'
with open(CONFIG_PATH, 'r') as f:
    CONFIG = json.load(f)


def compute_transmission(results: List[SimulationResult]) -> Tuple[float, float]:
    """
    Compute transmission fraction and uncertainty.
    
    Transmission is defined as the fraction of photons that exit
    through the far boundary (z > thickness).
    
    Args:
        results: List of SimulationResult objects
        
    Returns:
        Tuple of (transmission_fraction, standard_error)
        
    Example:
        >>> results = run_simulation(1000, 662, 2.0, 'lead')
        >>> T, T_err = compute_transmission(results)
        >>> print(f"T = {T:.3f} ± {T_err:.3f}")
    """
    n_total = len(results)
    if n_total == 0:
        return 0.0, 0.0
    
    n_transmitted = sum(1 for r in results if r.outcome == "transmitted")
    
    # Transmission fraction
    T = n_transmitted / n_total
    
    # Binomial standard error: σ_T = sqrt(T(1-T)/N)
    T_err = np.sqrt(T * (1.0 - T) / n_total) if n_total > 0 else 0.0
    
    return T, T_err


def compute_energy_weighted_transmission(results: List[SimulationResult], 
                                        E0_keV: float) -> Tuple[float, float]:
    """
    Compute energy-weighted transmission (accounts for Compton downscatter).
    
    T_E = (Σ E_transmitted) / (N * E0)
    
    This metric is useful for dose calculations where energy matters.
    
    Args:
        results: List of SimulationResult objects
        E0_keV: Initial photon energy
        
    Returns:
        Tuple of (energy_weighted_transmission, standard_error)
    """
    n_total = len(results)
    if n_total == 0:
        return 0.0, 0.0
    
    transmitted_energies = [
        r.final_energy for r in results if r.outcome == "transmitted"
    ]
    
    total_energy_transmitted = sum(transmitted_energies)
    total_energy_incident = n_total * E0_keV
    
    T_E = total_energy_transmitted / total_energy_incident if total_energy_incident > 0 else 0.0
    
    # Approximate uncertainty (simplified)
    T_E_err = T_E * np.sqrt(1.0 / len(transmitted_energies)) if transmitted_energies else 0.0
    
    return T_E, T_E_err


def compute_hvl(thicknesses: np.ndarray, 
                transmissions: np.ndarray,
                transmission_errors: Optional[np.ndarray] = None) -> Dict[str, float]:
    """
    Extract half-value layer using multiple methods.
    
    Methods:
    1. Direct interpolation: Find thickness where T = 0.5
    2. Exponential fit: Fit T = exp(-μ_eff * d) and compute HVL = ln(2)/μ_eff
    
    Args:
        thicknesses: Array of slab thicknesses (cm)
        transmissions: Array of transmission fractions
        transmission_errors: Optional array of uncertainties
        
    Returns:
        Dictionary with keys:
            - 'hvl_interpolation': HVL from direct interpolation
            - 'hvl_fit': HVL from exponential fit
            - 'mu_effective': Effective attenuation coefficient
            - 'fit_r_squared': Quality of exponential fit
            - 'fit_intercept': Y-intercept of ln(T) vs d
    """
    results = {}
    
    # Method 1: Direct interpolation
    if np.any(transmissions <= 0.5) and np.any(transmissions >= 0.5):
        # Find where transmission crosses 0.5
        idx = np.where(transmissions <= 0.5)[0]
        if len(idx) > 0 and idx[0] > 0:
            # Linear interpolation between two points
            i = idx[0]
            t1, t2 = thicknesses[i-1], thicknesses[i]
            T1, T2 = transmissions[i-1], transmissions[i]
            
            # Interpolate: d_hvl = t1 + (0.5 - T1) * (t2 - t1) / (T2 - T1)
            hvl_interp = t1 + (0.5 - T1) * (t2 - t1) / (T2 - T1)
            results['hvl_interpolation'] = hvl_interp
        else:
            results['hvl_interpolation'] = None
    else:
        results['hvl_interpolation'] = None
    
    # Method 2: Exponential fit
    # Use only the "good" range (avoid very thin and very thick)
    min_T = CONFIG['analysis']['min_transmission_fit']
    max_T = CONFIG['analysis']['max_transmission_fit']
    
    mask = (transmissions > min_T) & (transmissions < max_T) & (transmissions > 0)
    
    if mask.sum() >= 3:  # Need at least 3 points
        log_T = np.log(transmissions[mask])
        d_fit = thicknesses[mask]
        
        # Linear regression: ln(T) = -μ_eff * d + b
        # Use weighted least squares if errors provided
        if transmission_errors is not None:
            weights = 1.0 / (transmission_errors[mask] / transmissions[mask])**2
            weights = weights / weights.sum()
            coeffs = np.polyfit(d_fit, log_T, 1, w=weights)
        else:
            coeffs = np.polyfit(d_fit, log_T, 1)
        
        slope, intercept = coeffs
        mu_eff = -slope  # Attenuation coefficient
        
        # HVL = ln(2) / μ_eff
        if mu_eff > 0:
            hvl_fit = np.log(2.0) / mu_eff
            results['hvl_fit'] = hvl_fit
            results['mu_effective'] = mu_eff
            results['fit_intercept'] = intercept
            
            # Compute R²
            log_T_pred = slope * d_fit + intercept
            ss_res = np.sum((log_T - log_T_pred)**2)
            ss_tot = np.sum((log_T - log_T.mean())**2)
            r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
            results['fit_r_squared'] = r_squared
        else:
            results['hvl_fit'] = None
            results['mu_effective'] = None
            results['fit_r_squared'] = None
            results['fit_intercept'] = None
    else:
        results['hvl_fit'] = None
        results['mu_effective'] = None
        results['fit_r_squared'] = None
        results['fit_intercept'] = None
    
    return results


def extract_interaction_histogram(results: List[SimulationResult],
                                  thickness_cm: float,
                                  n_bins: int = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Create histogram of interaction depths by type.
    
    Args:
        results: List of SimulationResult objects
        thickness_cm: Slab thickness in cm
        n_bins: Number of depth bins (default from config)
        
    Returns:
        Tuple of (bin_edges, compton_counts, photoelectric_counts)
        
    Example:
        >>> bins, hist_c, hist_pe = extract_interaction_histogram(results, 5.0)
        >>> plt.bar(bins[:-1], hist_c, label='Compton')
        >>> plt.bar(bins[:-1], hist_pe, bottom=hist_c, label='PE')
    """
    if n_bins is None:
        n_bins = CONFIG['analysis']['depth_bins']
    
    # Collect all interactions
    compton_depths = []
    pe_depths = []
    
    for result in results:
        for depth, itype in zip(result.interaction_depths, result.interaction_types):
            if 0 <= depth <= thickness_cm:  # Only count interactions inside slab
                if itype == "compton":
                    compton_depths.append(depth)
                elif itype == "photoelectric":
                    pe_depths.append(depth)
    
    # Create bins
    bins = np.linspace(0, thickness_cm, n_bins + 1)
    
    # Histogram each type
    hist_compton, _ = np.histogram(compton_depths, bins=bins)
    hist_pe, _ = np.histogram(pe_depths, bins=bins)
    
    return bins, hist_compton, hist_pe


def extract_energy_spectrum(results: List[SimulationResult],
                           outcome_filter: str = "transmitted",
                           n_bins: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract energy spectrum of photons with specific outcome.
    
    Args:
        results: List of SimulationResult objects
        outcome_filter: "transmitted", "absorbed", or "backscattered"
        n_bins: Number of energy bins
        
    Returns:
        Tuple of (bin_centers, counts)
    """
    energies = [
        r.final_energy for r in results 
        if r.outcome == outcome_filter and r.final_energy > 0
    ]
    
    if not energies:
        return np.array([]), np.array([])
    
    counts, bin_edges = np.histogram(energies, bins=n_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    
    return bin_centers, counts


def extract_angular_distribution(results: List[SimulationResult],
                                outcome_filter: str = "transmitted",
                                n_bins: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract angular distribution (cos θ) of final photon directions.
    
    θ is measured from the +z axis (forward direction).
    
    Args:
        results: List of SimulationResult objects
        outcome_filter: "transmitted", "absorbed", or "backscattered"
        n_bins: Number of angular bins
        
    Returns:
        Tuple of (cos_theta_centers, counts)
    """
    cos_theta_values = []
    
    for r in results:
        if r.outcome == outcome_filter:
            # Calculate cos(θ) from final direction
            # Assume direction was stored in final_position (should be in result)
            # For transmitted: use z-component of direction at exit
            # Simplified: just note this needs tracking in simulator
            pass
    
    # This requires storing final direction in SimulationResult
    # For now, return empty arrays
    return np.array([]), np.array([])


def compute_build_up_factor(thicknesses: np.ndarray,
                            transmissions: np.ndarray,
                            transmissions_energy_weighted: np.ndarray,
                            mu_theory: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute empirical build-up factor.
    
    Build-up factor accounts for scattered radiation contribution:
    B(d) = T_measured / T_narrow_beam
    
    Where T_narrow_beam = exp(-μ * d) assumes no scattering.
    
    Args:
        thicknesses: Array of thicknesses
        transmissions: Measured transmission (includes scatter)
        transmissions_energy_weighted: Energy-weighted transmission
        mu_theory: Theoretical attenuation coefficient (cm⁻¹)
        
    Returns:
        Tuple of (build_up_number, build_up_energy)
    """
    # Narrow beam (theoretical, no scatter)
    T_narrow = np.exp(-mu_theory * thicknesses)
    
    # Build-up factor for photon number
    B_number = transmissions / T_narrow
    
    # Build-up factor for energy
    B_energy = transmissions_energy_weighted / T_narrow
    
    return B_number, B_energy


def summarize_results(results: List[SimulationResult]) -> Dict[str, any]:
    """
    Generate summary statistics for a simulation run.
    
    Args:
        results: List of SimulationResult objects
        
    Returns:
        Dictionary with summary statistics
    """
    n_total = len(results)
    
    # Count outcomes
    n_transmitted = sum(1 for r in results if r.outcome == "transmitted")
    n_absorbed = sum(1 for r in results if r.outcome == "absorbed")
    n_backscattered = sum(1 for r in results if r.outcome == "backscattered")
    
    # Count interactions
    total_interactions = sum(len(r.interaction_depths) for r in results)
    total_compton = sum(r.interaction_types.count("compton") for r in results)
    total_pe = sum(r.interaction_types.count("photoelectric") for r in results)
    
    # Average interactions per photon
    avg_interactions = total_interactions / n_total if n_total > 0 else 0
    
    # Energy statistics for transmitted photons
    transmitted_energies = [
        r.final_energy for r in results if r.outcome == "transmitted"
    ]
    
    if transmitted_energies:
        mean_transmitted_energy = np.mean(transmitted_energies)
        std_transmitted_energy = np.std(transmitted_energies)
    else:
        mean_transmitted_energy = 0.0
        std_transmitted_energy = 0.0
    
    summary = {
        'n_photons': n_total,
        'n_transmitted': n_transmitted,
        'n_absorbed': n_absorbed,
        'n_backscattered': n_backscattered,
        'transmission_fraction': n_transmitted / n_total if n_total > 0 else 0.0,
        'total_interactions': total_interactions,
        'total_compton': total_compton,
        'total_photoelectric': total_pe,
        'avg_interactions_per_photon': avg_interactions,
        'mean_transmitted_energy_keV': mean_transmitted_energy,
        'std_transmitted_energy_keV': std_transmitted_energy
    }
    
    return summary


def compare_to_theory(simulated_hvl: float,
                     E0_keV: float,
                     material: str,
                     mu_theory: float) -> Dict[str, float]:
    """
    Compare simulated HVL to theoretical narrow-beam value.
    
    Args:
        simulated_hvl: HVL from simulation (cm)
        E0_keV: Source energy (keV)
        material: Material name
        mu_theory: Theoretical total attenuation coefficient (cm⁻¹)
        
    Returns:
        Dictionary with comparison metrics
    """
    # Theoretical HVL (narrow beam)
    hvl_theory = np.log(2.0) / mu_theory
    
    # Relative difference
    rel_diff = (simulated_hvl - hvl_theory) / hvl_theory * 100.0
    
    comparison = {
        'hvl_simulated_cm': simulated_hvl,
        'hvl_theory_cm': hvl_theory,
        'absolute_difference_cm': simulated_hvl - hvl_theory,
        'relative_difference_percent': rel_diff,
        'mu_theory_cm_inv': mu_theory,
        'material': material,
        'energy_keV': E0_keV
    }
    
    return comparison


# Utility function for plotting (optional)
def prepare_transmission_plot_data(thicknesses: np.ndarray,
                                  transmissions: np.ndarray,
                                  transmission_errors: np.ndarray,
                                  hvl_data: Dict[str, float],
                                  mu_theory: float) -> Dict[str, np.ndarray]:
    """
    Prepare data arrays for transmission curve plotting.
    
    Returns dictionary with x and y data for:
    - Simulation points with error bars
    - Exponential fit line
    - Theoretical curve
    """
    plot_data = {
        'thickness': thicknesses,
        'transmission': transmissions,
        'transmission_err': transmission_errors
    }
    
    # Generate smooth curves for fits
    d_smooth = np.linspace(thicknesses.min(), thicknesses.max(), 200)
    
    # Exponential fit
    if hvl_data['mu_effective']:
        T_fit = np.exp(-hvl_data['mu_effective'] * d_smooth)
        plot_data['thickness_fit'] = d_smooth
        plot_data['transmission_fit'] = T_fit
    
    # Theory
    T_theory = np.exp(-mu_theory * d_smooth)
    plot_data['thickness_theory'] = d_smooth
    plot_data['transmission_theory'] = T_theory
    
    return plot_data

