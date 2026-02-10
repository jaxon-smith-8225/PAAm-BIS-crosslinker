# BIS Crosslinker Simulation

A modular Python package for simulating crosslinked polymer networks using BIS crosslinkers and polyacrylamide (PAAm) chains

## Structure

```
biscrosslinker/
├── __init__.py          # Package initialization and exports
├── config.py            # Configuration parameters and constants
├── geometry.py          # Geometric utilities and collision detection
├── chemistry.py         # Reactive site identification and manipulation
├── network.py           # Network management classes
├── crosslinking.py      # Core crosslinking operations
└── main.py              # Main simulation script
```

## Module Overview

### `config.py`
Contains all configuration parameters:
- Output settings (file names)
- Rotation parameters (angles, random vs systematic)
- Collision detection settings (chain radius)
- Simulation limits (max chains, min distances)

### `geometry.py`
Geometric utilities for molecular manipulation:
- `ChainCylinder`: Bounding volume class for collision detection
- `segment_segment_distance()`: 3D line segment distance calculation
- `cylinders_collide()`: Collision detection between cylinders
- `angle_between_vectors()`: Vector angle calculation
- `rotate_by_angle()`: Rotate molecules around an axis
- `align_molecule()`: Align molecules using vector alignment
- `apply_random_rotation()`: Apply random rotations

### `chemistry.py`
Chemistry-specific functions:
- `scan_chain()`: Identify reactive sites (amide groups) in polymer chains
- `unpack_linking_site()`: Extract atom references from site dictionaries
- `get_reactive_group()`: Get atoms to remove during crosslinking
- `calculate_reactive_direction()`: Compute reaction direction vectors
- `create_chain_cylinder()`: Create bounding volumes for chains
- `is_site_too_close()`: Check crosslink spacing constraints

### `network.py`
Classes for managing the growing network:
- `PAAmTemplate`: Reusable template with pre-computed reactive sites
  - Caches reactive site information
  - Provides fresh Universe instances
  - Optimizes repeated chain additions
- `CrosslinkNetwork`: Manages the entire crosslinked structure
  - Tracks all chains and cylinders
  - Handles collision detection
  - Finds valid crosslinking sites
  - Maintains crosslink positions and metadata

### `crosslinking.py`
Core crosslinking operations:
- `replace_with_bis()`: Replace reactive site with BIS crosslinker
- `connect_PAAm()`: Connect new PAAm chain to existing crosslinker
  - Handles site selection
  - Performs rotation attempts
  - Checks for collisions
  - Returns merged structure

### `main.py`
Main execution script:
- Initializes template and network
- Runs main simulation loop
- Handles errors and retries
- Outputs final structure

## Usage

### Running the simulation

```bash
python -m biscrosslinker.main
```

Or if you're in the same directory:

```bash
python main.py
```

### Using as a library

```python
from biscrosslinker import PAAmTemplate, CrosslinkNetwork
from biscrosslinker import replace_with_bis, connect_PAAm
import MDAnalysis as mda

# Load template
template = PAAmTemplate("PAAm25mer.pdb")

# Initialize network
initial_chain = mda.Universe("PAAm25mer.pdb", context='default')
network = CrosslinkNetwork(initial_chain)

# Add chains programmatically
# ... (see main.py for full example)
```

### Customizing parameters

Edit `config.py` to adjust:
- `MAX_NUM_CHAINS`: Maximum number of chains to add
- `CHAIN_RADIUS`: Collision detection radius
- `MIN_RAND_ANGLE`, `MAX_RAND_ANGLE`: Rotation angle range
- `RAND_CHAIN_ROTATION`: Random vs systematic rotation
- `MIN_CROSSLINK_DISTANCE`: Minimum spacing between crosslinks

## Requirements

- MDAnalysis
- NumPy

## Input Files Required

- `PAAm25mer.pdb`: Template PAAm chain structure
- `BIS.pdb`: BIS crosslinker molecule structure

## Output

- `output.pdb`: Final crosslinked network structure (configurable in `config.py`)

## Benefits of Modular Structure

1. **Easier testing**: Each module can be tested independently
2. **Better organization**: Related functions grouped logically
3. **Reusability**: Import specific functions in other projects
4. **Maintainability**: Changes isolated to relevant modules
5. **Collaboration**: Multiple people can work on different modules
6. **Documentation**: Clearer purpose for each component