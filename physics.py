"""
Physics module for photon interactions.

Contains: 
- Cross-section calculations 
- Klein-Nishina Compton scattering
- Photon kinematics and rotations
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple

@dataclass
class PhotonState: 
    position: np.ndarray
    direction: np.ndarray
    energy: float

class CrossSections: 
    def __init__(self, material: str, density: float = None):
        self.material = material
        self.density = density if density else self._get_default_density()

    def _get_default_density(self) -> float: 
        densities = {'lead':11.34, 'concrete':2.3}
        return densities.get(self.material, 1.0)
    
    def mu_photoelectric(self, E_keV: float) -> float: 
        if self.material == 'lead':
            return 5.0 * (100.0 / E_keV)**3
        else: # Concrete
            return 0.15 * (100.0 / E_keV)**3   

    def mu_compton(self, E_keV: float) -> float:      
        if self.material == 'lead': 
            return 0.5
        else: # Concrete
            return 0.2
        
    def mu_total(self, E_keV: float) -> float: 
        return self.mu_photoelectric(E_keV) + self.mu_compton(E_keV)
    
    def sample_compton_angle(E_keV: float, rng: np.random.Generator) -> float: 
        alpha = E_keV / 511.0
        max_attempts = 10000
        for _ in range (max_attempts):
            mu = rng.uniform (-1.0, 1.0)

            ratio = 1.0 / (1.0 + alpha * (1.0 - mu))
            kn_value = ratio**2 * (ratio + 1.0/ratio - (1.0 - mu**2))

            max_kn = (1.0 + alpha)**2 / (1.0 + 2.0*alpha)
            envelope = 2.0 * max_kn

            if rng.uniform(0.0, envelope) < kn_value:
                return mu
            
            return rng.uniform(-1.0, 1.0)
        
    def compton_energy(E_keV: float, cos_theta: float) -> float:
        alpha = E_keV / 511.0  # E / m_e c²
        return E_keV / (1.0 + alpha * (1.0 - cos_theta))
    
    def rotate_direction(u: np.ndarray, cos_theta: float, phi: float) -> np.ndarray:
        sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta**2))

        # Build local coordinate system
        # v is perpendicular to u
        if abs(u[2]) < 0.99:
            v = np.array([-u[1], u[0], 0.0])
        else:
            v = np.array([1.0, 0.0, 0.0])
        v = v / np.linalg.norm(v)

        # w completes the right-handed system
        w = np.cross(u, v)

        # New direction in local frame
        u_new = (sin_theta * np.cos(phi) * v +
                 sin_theta * np.sin(phi) * w +
                 cos_theta * u)
    
    # Normalize to ensure unit vector
        return u_new / np.linalg.norm(u_new)