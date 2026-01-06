"""
Monte Carlo photon transport simulator.

This module implements the core Monte Carlo algorithm for simulating
gamma photon transport through shielding materials.
"""

import numpy as np
import json
from dataclasses import dataclass, field
from typing import List, Tuple
from pathlib import Path

from physics import (
    PhotonState, 
    CrossSections, 
    sample_compton_angle, 
    compton_energy,
    rotate_direction
)

# Load configuration
CONFIG_PATH = Path(__file__).parent / 'config.json'
with open(CONFIG_PATH, 'r') as f:
    CONFIG = json.load(f)

@dataclass
class SimulationResult:
    """Results from simulating a single photon."""
    outcome: str  # "transmitted", "absorbed", or "backscattered"
    interaction_depths: List[float]
    interaction_types: List[str]  # "compton" or "photoelectric"
    final_energy: float
    final_position: np.ndarray

def simulate_photon(
    E0_keV: float,
    thickness_cm: float,
    xsec: CrossSections,
    rng: np.random.Generator
) -> SimulationResult:
    """
    Simulate a single photon through the slab (Monte Carlo).
    
    This is the core Monte Carlo algorithm that:
    1. Samples free path from exponential distribution
    2. Advances photon position
    3. Checks boundary conditions
    4. Samples interaction type (photoelectric vs Compton)
    5. Updates photon state based on physics
    6. Repeats until photon is absorbed or escapes
    
    Args:
        E0_keV: Initial photon energy in keV
        thickness_cm: Slab thickness in cm
        xsec: Cross-section object for the material
        rng: NumPy random number generator
        
    Returns:
        SimulationResult containing outcome and interaction history
    """
    # Initialize photon state
    state = PhotonState(
        position=np.array([0.0, 0.0, 0.0]),
        direction=np.array([0.0, 0.0, 1.0]),
        energy=E0_keV
    )
    
    interaction_depths = []
    interaction_types = []
    
    max_interactions = CONFIG['simulation']['max_interactions_per_photon']
    energy_cutoff = CONFIG['simulation']['energy_cutoff_keV']

    # Monte Carlo loop 
    for interaction_count in range(max_interactions): 
        # Step 1: Sample free path (exponential distribution)
        mu_tot = xsec.mu_total(state.energy)
        if mu_tot <= 0: 
            break
        u = rng.random()
        s = -np.log(u) / mu_tot

        # Step 2: Advance position
        new_pos = state.position + s * state.direction

        # Step 3: Check boundaries 
        if new_pos[2] > thickness_cm: 
            # Transmitted through exit force
            return SimulationResult(
                outcome="transmitted", 
                interaction_depths=interaction_depths, 
                interaction_types=interaction_types, 
                final_energy=state.energy, 
                final_position=new_pos
            )
        
        if new_pos[2] < 0: 
            # Backscattered through entrance force
            return SimulationResult(
                outcome="backscattered", 
                interaction_depths=interaction_depths, 
                interaction_types=interaction_types,
                final_energy=state.energy, 
                final_position=new_pos
            )
        
        # Update position
        state.position = new_pos

        # Step 4: Select interaction type
        mu_pe = xsec.mu_photoelectric(state.energy)
        p_photoelectric = mu_pe / mu_tot

        if rng.random() < p_photoelectric: 
            # Photoelectric absorption
            interaction_depths.append(state.position[2])
            interaction_types.append("photoelectric")
            
            return SimulationResult(
                outcome="absorbed",
                interaction_depths=interaction_depths,
                interaction_types=interaction_types,
                final_energy=0.0,
                final_position=state.position
            )
        
        # Compton scattering
        interaction_depths.append(state.position[2])
        interaction_types.append("compton")
        
        # Step 5: Sample scattering angles
        cos_theta = sample_compton_angle(state.energy, rng)
        phi = rng.uniform(0, 2*np.pi)
        
        # Step 6: Update photon state
        new_energy = compton_energy(state.energy, cos_theta)
        new_direction = rotate_direction(state.direction, cos_theta, phi)
        
        state = PhotonState(
            position=state.position,
            direction=new_direction,
            energy=new_energy
        )
        
        # Energy cutoff (absorb very low energy photons)
        if state.energy < energy_cutoff:
            return SimulationResult(
                outcome="absorbed",
                interaction_depths=interaction_depths,
                interaction_types=interaction_types,
                final_energy=state.energy,
                final_position=state.position
            )
    
    # Max interactions reached (treat as absorbed)
    return SimulationResult(
        outcome="absorbed",
        interaction_depths=interaction_depths,
        interaction_types=interaction_types,
        final_energy=state.energy,
        final_position=state.position
    )

def run_simulation(
    n_photons: int,
    E0_keV: float,
    thickness_cm: float,
    material: str,
    seed: int = None
) -> List[SimulationResult]:
    """
    Run Monte Carlo simulation for N photons.
    
    Args:
        n_photons: Number of photons to simulate
        E0_keV: Source energy in keV
        thickness_cm: Slab thickness in cm
        material: Material name ("lead" or "concrete")
        seed: Random seed for reproducibility
        
    Returns:
        List of SimulationResult objects
    """
    if seed is None:
        seed = CONFIG['simulation']['default_seed']
    
    # Initialize random number generator
    rng = np.random.default_rng(seed)
    
    # Load cross-sections
    xsec = CrossSections(material)
    
    # Run simulations
    results = []
    for i in range(n_photons):
        result = simulate_photon(E0_keV, thickness_cm, xsec, rng)
        results.append(result)
    
    return results
