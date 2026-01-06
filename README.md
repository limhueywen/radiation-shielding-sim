# 🔬 Monte Carlo Radiation Shielding Simulator

A physics-based Monte Carlo simulator for gamma photon transport through shielding materials (lead and concrete). Features interactive visualization, half-value layer (HVL) extraction, and comprehensive statistical analysis.

## 📚 Course Information

- **Subject**: SIF2018 Radiation Physics  
- **Mini-Project**: Radiation Shielding Simulator
- **Author**: Lim Huey Wen  
- **Matric Number**: 23005873  
- **Institution**: Universiti Malaya  

---

## 🌟 Features

- **Physically Accurate Simulation**
  - Klein-Nishina Compton scattering with rejection sampling
  - Energy-dependent photoelectric absorption
  - Free-path sampling from exponential distributions
  
- **Dual Analysis Modes**
  - **Single Thickness**: Quick analysis with depth histograms and energy spectra
  - **Thickness Sweep**: Full transmission curves with HVL extraction

- **Interactive Visualizations**
  - Interaction depth distributions
  - Transmitted photon energy spectra
  - Transmission vs. thickness curves (linear and log scales)
  - Real-time progress tracking

- **Comprehensive Analysis**
  - Half-Value Layer (HVL) extraction via interpolation and exponential fitting
  - Comparison with theoretical narrow-beam attenuation
  - Statistical uncertainties (binomial errors)
  - Data export (CSV)

---

## 📋 Table of Contents

- [Quick Start](#-quick-start)
- [Installation](#-installation)
- [Usage](#-usage)
- [Physics Model](#-physics-model)
- [Project Structure](#-project-structure)
- [Examples](#-examples)
- [Advanced Usage](#-advanced-usage)
- [Contributing](#-contributing)

---

## 🚀 Quick Start
```bash
# Clone the repository
git clone https://github.com/limhueywen/radiation-shielding-sim/settings
cd radiation-shielding-sim

# Create virtual environment
python -m venv venv

# Activate virtual environment
source venv/bin/activate  # Mac/Linux
# OR
venv\Scripts\activate     # Windows

# Install dependencies
pip install -r requirements.txt

# Run the Streamlit app
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`

---

## 💾 Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager
- Git (for cloning)

### Step-by-Step Installation

#### Option 1: Using Virtual Environment (Recommended)
```bash
# 1. Clone repository
git clone https://github.com/limhueywen/radiation-shielding-sim/settings
cd radiation-shielding-sim

# 2. Create virtual environment
python -m venv venv

# 3. Activate virtual environment
# On macOS/Linux:
source venv/bin/activate
# On Windows:
venv\Scripts\activate

# 4. Install dependencies
pip install -r requirements.txt

# 5. Verify installation
python -c "import numpy, scipy, matplotlib, streamlit; print('✓ All dependencies installed')"
```

#### Option 2: Using Conda
```bash
# 1. Clone repository
git clone https://github.com/limhueywen/radiation-shielding-sim/settings
cd radiation-shielding-sim

# 2. Create conda environment
conda create -n radshield python=3.10
conda activate radshield

# 3. Install dependencies
pip install -r requirements.txt
```

### Dependencies

All dependencies are listed in `requirements.txt`:
```
numpy>=1.24.0       # Numerical computations
scipy>=1.10.0       # Scientific functions
matplotlib>=3.7.0   # Plotting
pandas>=2.0.0       # Data manipulation
streamlit>=1.28.0   # Web interface
```

---

## 🎯 Usage

### Running the Web Interface
```bash
streamlit run app.py
```

This launches an interactive web interface with:
- **Sidebar controls**: Material, energy, photon count, seed
- **Analysis modes**: Single thickness or thickness sweep
- **Real-time visualization**: Plots update as simulations complete

### Using the Simulator Programmatically

#### Basic Example
```python
from simulator import run_simulation
from analysis import compute_transmission, summarize_results

# Run simulation
results = run_simulation(
    n_photons=10000,
    E0_keV=662,           # Cs-137 source
    thickness_cm=2.0,
    material='lead',
    seed=42
)

# Analyze results
T, T_err = compute_transmission(results)
print(f"Transmission: {T:.4f} ± {T_err:.4f}")

summary = summarize_results(results)
print(f"Absorbed: {summary['n_absorbed']}")
print(f"Transmitted: {summary['n_transmitted']}")
```

#### Thickness Sweep Example
```python
import numpy as np
from simulator import run_simulation
from analysis import compute_transmission, compute_hvl

# Define thickness range
thicknesses = np.linspace(0.5, 5.0, 10)
transmissions = []
errors = []

# Run sweep
for d in thicknesses:
    results = run_simulation(10000, 662, d, 'lead', seed=42)
    T, T_err = compute_transmission(results)
    transmissions.append(T)
    errors.append(T_err)

# Extract HVL
hvl_data = compute_hvl(
    np.array(thicknesses),
    np.array(transmissions),
    np.array(errors)
)

print(f"HVL (interpolation): {hvl_data['hvl_interpolation']:.3f} cm")
print(f"HVL (exponential fit): {hvl_data['hvl_fit']:.3f} cm")
print(f"μ_effective: {hvl_data['mu_effective']:.3f} cm⁻¹")
```

#### Custom Configuration

Edit `config.json` to change defaults:
```json
{
  "simulation": {
    "default_n_photons": 10000,
    "max_photons": 100000,
    "energy_cutoff_keV": 10.0
  },
  "materials": {
    "lead": {"density_g_cm3": 11.34},
    "concrete": {"density_g_cm3": 2.3}
  },
  "sources": {
    "Cs-137": 662,
    "Co-60": 1250,
    "Custom": 500
  }
}
```

---

## 📐 Physics Model

### Geometry

- **Planar slab** of thickness *d* with normal along +z axis
- **Monoenergetic, collimated beam** incident at z = 0
- **Boundary conditions**:
  - z > d: transmitted
  - z < 0: backscattered
  - 0 ≤ z ≤ d: inside slab

### Interaction Physics

#### 1. Free Path Sampling

Path length between interactions sampled from exponential distribution:
```
s = -ln(u) / μ_tot(E)
```

where *u* ~ Uniform(0,1) and μ_tot(E) is the total attenuation coefficient.

#### 2. Interaction Selection

Interaction type selected via probability:
```
p_photoelectric = μ_pe(E) / μ_tot(E)
p_compton = μ_c(E) / μ_tot(E)
```

#### 3. Compton Scattering

**Energy Loss:**
```
E' = E / (1 + (E/m_e c²)(1 - cos θ))
```

**Angular Distribution:** Klein-Nishina differential cross-section:
```
dσ/dΩ ∝ (E'/E)² × (E/E' + E'/E - sin²θ)
```

Sampled via **rejection method** with uniform envelope.

#### 4. Photoelectric Absorption

Photon is completely absorbed. History terminates.

### Cross-Section Model

Current implementation uses **simplified power-law approximations**:

- **Lead**: μ_pe ∝ E⁻³, μ_c ≈ constant
- **Concrete**: Similar scaling with lower density

**For production use**: Replace with NIST XCOM database values 

---

## 📁 Project Structure
```
radiation-shielding-sim/
├── README.md                 # This file
├── requirements.txt          # Python dependencies
├── config.json              # Simulation parameters
├── .gitignore               # Git ignore rules
│
├── physics.py               # Physics calculations
│   ├── PhotonState          # Photon state dataclass
│   ├── CrossSections        # Material cross-sections
│   ├── sample_compton_angle # Klein-Nishina sampling
│   ├── compton_energy       # Compton kinematics
│   └── rotate_direction     # Vector rotation
│
├── simulator.py             # Monte Carlo engine
│   ├── SimulationResult     # Result dataclass
│   ├── simulate_photon      # Core MC loop (single photon)
│   └── run_simulation       # Ensemble simulation
│
├── analysis.py              # Post-processing
│   ├── compute_transmission # Transmission calculation
│   ├── compute_hvl          # HVL extraction
│   ├── extract_interaction_histogram
│   ├── extract_energy_spectrum
│   └── summarize_results
│
├── app.py                   # Streamlit web interface
│
├── data/                    # Cross-section data files
│   └── (optional NIST XCOM data)
│
└── figures/                 # Generated plots
    └── (saved visualizations)
```

### Module Descriptions

| Module | Purpose | Key Functions |
|--------|---------|---------------|
| `physics.py` | Fundamental physics calculations | Cross-sections, Klein-Nishina, rotations |
| `simulator.py` | Monte Carlo transport engine | Photon tracking, interaction sampling |
| `analysis.py` | Statistical analysis & HVL extraction | Transmission, HVL, histograms |
| `app.py` | Interactive web interface | Streamlit UI, visualization |
| `config.json` | Configuration parameters | Materials, sources, defaults |

---

## 📊 Examples

### Example 1: Lead Shielding for Cs-137
```python
from simulator import run_simulation
from analysis import compute_transmission

results = run_simulation(
    n_photons=50000,
    E0_keV=662,        # Cs-137
    thickness_cm=2.0,
    material='lead',
    seed=42
)

T, T_err = compute_transmission(results)
print(f"Transmission through 2 cm lead: {T:.4f} ± {T_err:.4f}")
# Expected output: ~0.52 ± 0.002 (near HVL)
```

### Example 2: HVL Determination
```python
import numpy as np
import matplotlib.pyplot as plt
from simulator import run_simulation
from analysis import compute_hvl, compute_transmission

# Sweep thickness
thicknesses = np.linspace(0.5, 5.0, 10)
T_values = []

for d in thicknesses:
    results = run_simulation(10000, 662, d, 'lead')
    T, _ = compute_transmission(results)
    T_values.append(T)

# Extract HVL
hvl_data = compute_hvl(thicknesses, np.array(T_values))
print(f"HVL: {hvl_data['hvl_fit']:.2f} cm")

# Plot
plt.figure(figsize=(10, 6))
plt.plot(thicknesses, T_values, 'o-', label='Simulation')
plt.axhline(0.5, color='r', linestyle='--', label='T = 0.5')
plt.axvline(hvl_data['hvl_fit'], color='g', linestyle='--', 
            label=f"HVL = {hvl_data['hvl_fit']:.2f} cm")
plt.xlabel('Thickness (cm)')
plt.ylabel('Transmission')
plt.legend()
plt.grid(True)
plt.show()
```

### Example 3: Energy Spectrum Analysis
```python
from simulator import run_simulation
import matplotlib.pyplot as plt

results = run_simulation(10000, 662, 3.0, 'lead')

# Extract transmitted energies
transmitted = [r.final_energy for r in results if r.outcome == "transmitted"]

plt.figure(figsize=(10, 6))
plt.hist(transmitted, bins=50, alpha=0.7, edgecolor='black')
plt.axvline(662, color='r', linestyle='--', label='Source (662 keV)')
plt.xlabel('Energy (keV)')
plt.ylabel('Count')
plt.title('Transmitted Photon Energy Spectrum')
plt.legend()
plt.grid(True)
plt.show()
```

---

## 🔬 Advanced Usage

### Custom Cross-Sections

To use NIST XCOM data instead of approximations:

1. Download data from [NIST XCOM](https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html)
2. Save as CSV in `data/` directory
3. Modify `CrossSections` class in `physics.py`:
```python
import pandas as pd
from scipy.interpolate import interp1d

class CrossSections:
    def __init__(self, material: str):
        # Load NIST data
        df = pd.read_csv(f'data/xcom_{material}.csv')
        self.E_grid = df['Energy_MeV'].values * 1000  # Convert to keV
        
        # Create interpolators (log-log for better behavior)
        self.mu_pe_interp = interp1d(
            np.log(self.E_grid), 
            np.log(df['mu_pe'].values),
            kind='cubic'
        )
        # ... similar for mu_compton
```

### Parallel Processing

For large-scale simulations:
```python
from multiprocessing import Pool
from functools import partial

def run_single(seed_offset, n_photons, E0, d, material):
    return run_simulation(n_photons, E0, d, material, seed=42+seed_offset)

# Parallel execution
with Pool(8) as pool:
    func = partial(run_single, n_photons=10000, E0=662, d=2.0, material='lead')
    results_list = pool.map(func, range(10))  # 10 independent runs
    
# Combine results
all_results = [r for results in results_list for r in results]
```

### Adding New Materials

Edit `config.json`:
```json
{
  "materials": {
    "lead": {"density_g_cm3": 11.34, "description": "Pure lead"},
    "concrete": {"density_g_cm3": 2.3, "description": "Standard concrete"},
    "tungsten": {"density_g_cm3": 19.3, "description": "Pure tungsten"}
  }
}
```

Then add cross-section data in `physics.py` for the new material.

---

## 🧪 Testing

### Run Unit Tests
```bash
# Test physics module
python -c "from physics import CrossSections; x = CrossSections('lead'); assert x.mu_total(662) > 0"

# Test simulator
python -c "from simulator import run_simulation; r = run_simulation(100, 662, 2.0, 'lead'); assert len(r) == 100"

# Test analysis
python -c "from analysis import compute_transmission; from simulator import run_simulation; r = run_simulation(100, 662, 2.0, 'lead'); T, err = compute_transmission(r); assert 0 <= T <= 1"
```

### Validation Tests

Compare against known results:
```python
# Test 1: Thin slab should have high transmission
results = run_simulation(10000, 662, 0.1, 'lead')
T, _ = compute_transmission(results)
assert T > 0.9, "Thin slab transmission too low"

# Test 2: Thick slab should block most photons
results = run_simulation(10000, 662, 10.0, 'lead')
T, _ = compute_transmission(results)
assert T < 0.01, "Thick slab transmission too high"

# Test 3: Energy conservation (transmitted E ≤ incident E)
results = run_simulation(10000, 662, 2.0, 'lead')
for r in results:
    if r.outcome == "transmitted":
        assert r.final_energy <= 662, "Energy increased!"
```

---

## 🤝 Contributing

Contributions are welcome! Please follow these guidelines:

1. **Fork** the repository
2. **Create a branch**: `git checkout -b feature/amazing-feature`
3. **Make changes** with clear commit messages
4. **Test** your changes thoroughly
5. **Submit a pull request**

### Code Style

- Follow PEP 8 style guidelines
- Add docstrings to all functions
- Include type hints where appropriate
- Write descriptive variable names

### Areas for Improvement

- [ ] Add pair production for E > 1.022 MeV
- [ ] Implement secondary electron tracking
- [ ] Add buildup factor calculations
- [ ] Support for layered geometries
- [ ] Variance reduction techniques
- [ ] GPU acceleration with CuPy
- [ ] Real-time 3D visualization

---

## ⚠️ Known Limitations

- **Simplified cross-sections**: Uses power-law approximations instead of NIST data
- **No pair production**: Only valid for E < 1 MeV (Cs-137, Co-60)
- **No secondary particles**: Electrons from photoelectric/Compton not tracked
- **Planar geometry only**: Cannot model spheres, cylinders, or complex shapes
- **Single material**: Layered shields not supported

For production use, consider using established codes like MCNP, GEANT4, or EGSnrc.

---

## 👤 Author

**[Lim Huey Wen]**
- GitHub: [@limhueywen](https://github.com/limhueywen)
- Email: hwen2004@gmail.com

---

## 🙏 Acknowledgments

- NIST for cross-section databases
- Streamlit for the excellent web framework
- NumPy/SciPy communities for scientific computing tools

---

**Current Version**: 1.0.0

**Roadmap**:
- ✅ Basic Monte Carlo engine
- ✅ Klein-Nishina sampling
- ✅ HVL extraction
- ✅ Streamlit interface
- ⬜ NIST data integration
- ⬜ Pair production
- ⬜ GPU acceleration

---

<div align="center">

**⚛️ Built with Python for Radiation Physics Education ⚛️**

*Last Updated: January 2026*

</div>