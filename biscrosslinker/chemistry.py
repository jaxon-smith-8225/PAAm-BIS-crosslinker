"""
Chemistry utilities for identifying and manipulating reactive sites
"""

import numpy as np
import MDAnalysis as mda
from .config import CHAIN_RADIUS
from .geometry import ChainCylinder


def is_site_too_close(site, crosslink_positions, min_distance):
    """
    Check if a reactive site is too close to existing crosslinks
    
    Args:
        site: Reactive site dictionary
        crosslink_positions: List of crosslink position vectors
        min_distance: Minimum allowed distance
        
    Returns:
        True if site is too close to any crosslink
    """
    site_pos = site['backbone_carbon'][0].position
    for crosslink_pos in crosslink_positions:
        distance = np.linalg.norm(site_pos - crosslink_pos)
        if distance < min_distance:
            return True
    return False


def scan_chain(universe):
    """
    Scan a polymer chain for reactive sites (amide groups)
    
    Identifies carbon atoms with the pattern: C-O, C-N(H2), where C is also
    bonded to a backbone carbon.
    
    Args:
        universe: MDAnalysis Universe containing the polymer chain
        
    Returns:
        List of dictionaries, each containing:
            - central_carbon: The reactive carbon atom
            - backbone_carbon: List with backbone carbon atom
            - oxygen: List with oxygen atom
            - nitrogen: List with nitrogen atom
            - hydrogens: List with two hydrogen atoms
    """
    # Select all carbon atoms
    carbons = universe.select_atoms('element C')
    reactive_groups = []
    
    for c_atom in carbons:
        # Get bonded atoms
        bonded = c_atom.bonded_atoms
        
        # Look for oxygen and nitrogen among bonded atoms
        oxygen = bonded.select_atoms('element O')
        nitrogen = bonded.select_atoms('element N')
        
        # Must have exactly 1 oxygen and 1 nitrogen
        if len(oxygen) != 1 or len(nitrogen) != 1:
            continue
        
        # Find the backbone carbon (another carbon bonded to this one)
        bonded_carbons = bonded.select_atoms('element C')
        if len(bonded_carbons) != 1:
            continue
        backbone_carbon = bonded_carbons
        
        # Get hydrogens bonded to the nitrogen
        n_atom = nitrogen[0]
        hydrogens = n_atom.bonded_atoms.select_atoms('element H')
        
        # Must have exactly 2 hydrogens on nitrogen
        if len(hydrogens) != 2:
            continue
        
        reactive_groups.append({
            'central_carbon': c_atom,
            'backbone_carbon': backbone_carbon,
            'oxygen': oxygen,
            'nitrogen': nitrogen,
            'hydrogens': hydrogens
        })
    
    return reactive_groups


def unpack_linking_site(linking_site):
    """
    Extract atoms from linking site dictionary into a more convenient format
    
    Args:
        linking_site: Dictionary from scan_chain()
        
    Returns:
        Dictionary with individual atom references
    """
    return {
        'reactive_carbon': linking_site['central_carbon'],
        'backbone_carbon': linking_site['backbone_carbon'][0],
        'oxygen': linking_site['oxygen'][0],
        'nitrogen': linking_site['nitrogen'][0],
        'hydrogen1': linking_site['hydrogens'][0],
        'hydrogen2': linking_site['hydrogens'][1]
    }


def get_reactive_group(unpacked_site):
    """
    Create AtomGroup of atoms to remove during crosslinking
    
    Args:
        unpacked_site: Dictionary from unpack_linking_site()
        
    Returns:
        MDAnalysis AtomGroup containing atoms to be removed
    """
    return mda.AtomGroup([
        unpacked_site['reactive_carbon'],
        unpacked_site['oxygen'],
        unpacked_site['nitrogen'],
        unpacked_site['hydrogen1'],
        unpacked_site['hydrogen2']
    ])


def calculate_reactive_direction(unpacked_site):
    """
    Calculate direction vector from backbone to reactive carbon
    
    Args:
        unpacked_site: Dictionary from unpack_linking_site()
        
    Returns:
        Direction vector (numpy array)
    """
    return (unpacked_site['reactive_carbon'].position - 
            unpacked_site['backbone_carbon'].position)


def create_chain_cylinder(chain_atom_group, radius=CHAIN_RADIUS,
                          start_index=None, end_index=None):
    """
    Create cylindrical bounding volume for a chain.

    If start_index and end_index are both provided, they are used to look up
    the terminal atom positions directly — useful when you know the exact
    indices of the end carbons. Positions are looked up fresh from the atom
    group each call, so they are always current after translations/rotations.
    Otherwise, the function falls back to the first and last carbon atoms in
    the atom group.

    Args:
        chain_atom_group: MDAnalysis AtomGroup representing the chain
        radius: Radius of the cylinder
        start_index: Optional integer index of the start terminal atom
        end_index: Optional integer index of the end terminal atom

    Returns:
        ChainCylinder object
    """
    if start_index is not None and end_index is not None:
        start_point = chain_atom_group.atoms[start_index].position.copy()
        end_point   = chain_atom_group.atoms[end_index].position.copy()
        return ChainCylinder(start_point, end_point, radius)

    # Fallback: find terminal carbons by their unique name CTS
    terminal_carbons = chain_atom_group.select_atoms('name CTS')

    if len(terminal_carbons) == 2:
        start_point = terminal_carbons[0].position.copy()
        end_point   = terminal_carbons[1].position.copy()
    elif len(terminal_carbons) < 2:
        # Last resort if CTS atoms aren't present (e.g. BIS or partial structures)
        start_point = chain_atom_group[0].position.copy()
        end_point   = chain_atom_group[-1].position.copy()
    else:
        raise ValueError(f"Expected exactly 2 CTS atoms, found {len(terminal_carbons)}. "
                         f"Check your PDB template.")

    return ChainCylinder(start_point, end_point, radius)