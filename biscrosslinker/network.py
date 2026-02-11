"""
Network management classes for building crosslinked polymer structures
"""

import MDAnalysis as mda
from .chemistry import scan_chain, is_site_too_close, create_chain_cylinder
from .geometry import cylinders_collide


class PAAmTemplate:
    """Reusable PAAm template with pre-computed reactive sites"""
    
    def __init__(self, pdb_file="PAAm25mer.pdb"):
        """
        Initialize template from PDB file
        
        Args:
            pdb_file: Path to PAAm template PDB file
        """
        # Load the template 
        self.universe = mda.Universe(pdb_file, context='default', 
                                     to_guess=['elements', 'bonds'])
        self.original_positions = self.universe.atoms.positions.copy()
        
        # Scan for reactive sites
        self.reactive_sites = scan_chain(self.universe)
        
        # Pre-compute which sites are valid (not first or last)
        self.valid_site_indices = list(range(1, len(self.reactive_sites) - 1))
        
        print(f"Template loaded: {len(self.reactive_sites)} reactive sites found, "
              f"{len(self.valid_site_indices)} are valid for crosslinking")
    
    def get_fresh_universe(self):
        """
        Return a new universe with positions reset to original
        
        Returns:
            New MDAnalysis Universe instance
        """
        # Create new universe instance
        new_univ = mda.Universe(self.universe.filename, context='default', 
                                to_guess=['elements', 'bonds'])
        # Reset to original positions
        new_univ.atoms.positions = self.original_positions.copy()
        return new_univ
    
    def get_site_info(self, site_index):
        """
        Get reactive site info without re-scanning
        
        Args:
            site_index: Index of the reactive site
            
        Returns:
            Site dictionary
        """
        return self.reactive_sites[site_index]


class CrosslinkNetwork:
    """Manages the growing crosslinked structure and its metadata"""
    
    def __init__(self, initial_chain_universe):
        """
        Initialize network with first chain
        
        Args:
            initial_chain_universe: MDAnalysis Universe of the initial chain
        """
        self.structure = initial_chain_universe  # Current mda.Universe
        self.cylinders = []  # Collision detection cylinders
        self.crosslink_positions = []  # Positions of crosslinks
        self.chain_info = []  # Metadata about each chain
        
        # Add the initial chain
        initial_cylinder = create_chain_cylinder(initial_chain_universe.atoms)
        self.cylinders.append(initial_cylinder)
        self.chain_info.append({
            'chain_id': 0,
            'cylinder': initial_cylinder,
            'added_at_step': 0
        })
    
    def add_chain(self, new_structure, new_cylinder, crosslink_position):
        """
        Add a new chain to the network
        
        Args:
            new_structure: Updated MDAnalysis universe with new chain
            new_cylinder: ChainCylinder for the new chain
            crosslink_position: Position where crosslink was made
        """
        self.structure = new_structure
        self.cylinders.append(new_cylinder)
        self.crosslink_positions.append(crosslink_position)
        
        chain_id = len(self.chain_info)
        self.chain_info.append({
            'chain_id': chain_id,
            'cylinder': new_cylinder,
            'added_at_step': chain_id
        })
    
    def num_chains(self):
        """
        Return total number of chains
        
        Returns:
            Number of chains in the network
        """
        return len(self.chain_info)
    
    def check_collision(self, new_cylinder, skip_cylinder=None):
        """
        Check if new cylinder collides with any existing cylinders

        Args:
            new_cylinder: ChainCylinder to check
            skip_cylinder: ChainCylinder to skip (e.g. the parent cylinder the
                           new chain is directly attached to, which will always
                           share an endpoint and should not count as a collision)

        Returns:
            True if collision detected, False otherwise
        """
        for cyl in self.cylinders:
            if skip_cylinder is not None and cyl is skip_cylinder:
                continue
            if cylinders_collide(new_cylinder, cyl):
                return True
        return False
    
    def find_valid_sites(self, min_distance=10.0):
        """
        Find all reactive sites that are far enough from existing crosslinks
        
        Args:
            min_distance: Minimum distance from existing crosslinks
            
        Returns:
            List of valid site dictionaries
        """
        linking_groups = scan_chain(self.structure)
        
        # Filter sites: exclude first/last and too close to crosslinks
        valid_sites = []
        for site in linking_groups[1:-1]:
            if not is_site_too_close(site, self.crosslink_positions, min_distance):
                valid_sites.append(site)
        
        return valid_sites