# Monte Carlo Radiation Shielding Simulator

Gamma photon transport through lead and concrete using Monte Carlo methods.

## Features
- Klein-Nishina Compton scattering with rejection sampling
- Energy-dependent photoelectric absorption
- Transmission curve analysis and HVL extraction
- Interactive Streamlit interface

## Installation
```bash
git clone https://github.com/YOUR_USERNAME/radiation-shielding-sim.git
cd radiation-shielding-sim
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Usage

### Run Streamlit App
```bash
streamlit run app.py
```

### Run Simulation Programmatically
```python
from simulator import run_simulation

results = run_simulation(
    n_photons=10000,
    E0_keV=662,  # Cs-137
    thickness_cm=2.0,
    material='lead',
    seed=42
)

# Analyze
n_transmitted = sum(1 for r in results if r.outcome == "transmitted")
print(f"Transmission: {n_transmitted/len(results):.3f}")
```

## Physics Model
- **Free path sampling**: Exponential distribution with μ_tot(E)
- **Compton scattering**: Klein-Nishina differential cross-section
- **Energy loss**: E' = E / (1 + (E/m_e c²)(1 - cos θ))
- **Angular distribution**: Rejection sampling from KN

## File Structure
- `physics.py` — Cross-sections, Klein-Nishina, kinematics
- `simulator.py` — Monte Carlo engine
- `analysis.py` — HVL extraction, histograms
- `app.py` — Streamlit interface

## References
- NIST XCOM database (cross-sections)
- Klein-Nishina formula (Compton scattering)