"""
BIS Crosslinker - A package for simulating crosslinked polymer networks

Main components:
- network: PAAmTemplate and CrosslinkNetwork for managing structure
- crosslinking: Core operations for connecting chains
- chemistry: Reactive site identification and manipulation
- geometry: Collision detection and molecular alignment
- config: Configuration parameters
"""

from .network import PAAmTemplate, CrosslinkNetwork
from .crosslinking import replace_with_bis, connect_PAAm
from .chemistry import scan_chain, create_chain_cylinder
from .geometry import ChainCylinder
from .config import *

__all__ = [
    'PAAmTemplate',
    'CrosslinkNetwork',
    'replace_with_bis',
    'connect_PAAm',
    'scan_chain',
    'create_chain_cylinder',
    'ChainCylinder',
]