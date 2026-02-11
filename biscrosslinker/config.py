"""
Configuration parameters for BIS crosslinker simulation
"""

# Output settings
OUTPUT_STRUCTURE = "output.pdb"

# Rotation parameters
MIN_RAND_ANGLE = 30
MAX_RAND_ANGLE = 200
RAND_CHAIN_ROTATION = True
MAX_ROTATION_ATTEMPTS = 25

# Collision detection
CHAIN_RADIUS = 4.2

# Simulation limits
MAX_NUM_CHAINS = 20

# Crosslink spacing
MIN_CROSSLINK_DISTANCE = 20.0

# Loop formation
MAX_BIS_ORIENTATION_ANGLE = 180.0  # degrees
ENABLE_LOOPING = True